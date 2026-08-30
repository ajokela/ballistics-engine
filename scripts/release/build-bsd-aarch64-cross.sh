#!/usr/bin/env bash
# Cross-build the three aarch64 BSD release binaries on the x86_64 build server.
#
# This is the BUILD half of a two-step lane. It produces binaries and NEVER
# executes them - the build host is x86_64 and physically cannot. The binaries
# are not releasable until `validate-bsd-aarch64.sh` has run each one inside the
# matching aarch64 KVM guest. Running only this script and shipping the output
# would reintroduce exactly the "it compiled, therefore it works" gap that the
# old in-guest builds closed by accident.
#
#   ./build-bsd-aarch64-cross.sh 0.35.1 [OUTDIR]
#   ./validate-bsd-aarch64.sh    0.35.1  OUTDIR      # <- mandatory second step
#
# Why cross-build at all: the in-guest lane (build-k3s-bsds.sh) spends ~45 min
# per OS compiling inside an emulated-storage KVM guest on ONE ARM node, and
# that node has been the release critical path - and its single point of
# failure - repeatedly. Cross-compiling takes ~3 min per OS on 32 x86_64 cores.
#
# TOOLCHAIN CONSTRAINTS - these are load-bearing, not preferences:
#
#   * clang-22 + lld-22, no older. clang 14 AND clang 19 both SEGFAULT (exit
#     139, frontend crash) compiling ring's crypto/curve25519/curve25519.c for
#     aarch64-unknown-openbsd. Debian bookworm-security carries clang-22
#     (1:22.1.8-1~deb12u1), so plain `apt-get install clang-22 lld-22` on
#     rust:1-bookworm is enough - no apt.llvm.org repo needed.
#
#   * OpenBSD ships ONLY versioned shared libraries (libc.so.102.0) with no
#     plain libc.so symlink. lld does not implement OpenBSD's versioned-library
#     search. Given -lc it finds no libc.so, SILENTLY falls back to the static
#     libc.a, and then dies with "duplicate symbol: atexit" against crtbeginS.o.
#     The fix is to create an unversioned symlink for every lib*.so.* in the
#     sysroot's usr/lib. Without it the link fails; with it there are zero
#     errors. This is done by the sysroot prep below - do not remove it.
#
#   * FreeBSD is Tier 2 with prebuilt std (`rustup target add`). OpenBSD and
#     NetBSD are Tier 3 with NO prebuilt std, so they need
#     `cargo +nightly build -Z build-std=std,panic_abort` plus a nightly
#     toolchain carrying rust-src.
#
#   * The FreeBSD sysroot must match the GUEST release (14.3), not whatever is
#     newest. A binary linked against an older base runs on a newer system; the
#     reverse does not.
#
# Everything is cached and idempotent: the builder image is content-addressed by
# its Dockerfile, the three sysroots are downloaded once and stamped, and the
# cargo registry plus target dir persist between runs.
#
# Env overrides:
#   CROSS_HOST      ssh destination of the x86_64 build server (default alex@10.1.1.27)
#   CROSS_WORKDIR   cache/scratch dir on that host          (default ~/bsd-cross)
#   CROSS_ONLY_OS   build just one of freebsd|openbsd|netbsd (default all three)
#   CROSS_REBUILD_IMAGE=1   force a builder-image rebuild
set -Eeuo pipefail

die() { echo "error: $*" >&2; exit 1; }

V="${1:?usage: build-bsd-aarch64-cross.sh VERSION [OUTDIR]}"
OUT="${2:-$HOME/release-$V}"
[[ "$V" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-][A-Za-z0-9.]+)?$ ]] || die "bad VERSION: $V"

CROSS_HOST="${CROSS_HOST:-alex@10.1.1.27}"
CROSS_WORKDIR="${CROSS_WORKDIR:-bsd-cross}"
CROSS_ONLY_OS="${CROSS_ONLY_OS:-}"
CROSS_REBUILD_IMAGE="${CROSS_REBUILD_IMAGE:-0}"
if [[ -n "$CROSS_ONLY_OS" ]]; then
  case "$CROSS_ONLY_OS" in freebsd|openbsd|netbsd) ;; *) die "CROSS_ONLY_OS must be freebsd, openbsd or netbsd" ;; esac
  OS_LIST="$CROSS_ONLY_OS"
else
  OS_LIST="freebsd openbsd netbsd"
fi

REPO="$(git rev-parse --show-toplevel)" || die "not inside a git repo"
TAG="v$V"
git -C "$REPO" rev-parse -q --verify "refs/tags/$TAG" >/dev/null ||
  die "tag $TAG does not exist locally (git fetch --tags?)"
SOURCE_SHA="$(git -C "$REPO" rev-parse "${TAG}^{commit}")"
SOURCE_DATE_EPOCH="$(git -C "$REPO" log -1 --format=%ct "${TAG}^{commit}")"

# The version in the TAGGED tree is what ships. Reading Cargo.toml from the
# working tree instead would happily build a mismatched binary whenever the
# checkout has moved on, which is how stale binaries reach a release.
TAG_VERSION="$(git -C "$REPO" show "${TAG}:Cargo.toml" | awk '
  /^\[package\][[:space:]]*$/ { in_package = 1; next }
  /^\[/ && !/^\[package\]/    { in_package = 0 }
  in_package && /^version[[:space:]]*=/ {
    v = $0; sub(/^[^=]*=[[:space:]]*"/, "", v); sub(/".*/, "", v); print v; exit
  }')"
[[ "$TAG_VERSION" == "$V" ]] ||
  die "Cargo.toml at $TAG says version $TAG_VERSION, not $V"

RUNID="$(date +%s)-$(printf '%04x' $((RANDOM % 65536)))"
STAGE="$(mktemp -d "${TMPDIR:-/tmp}/bsd-cross.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT HUP INT TERM

echo "==> staging source from $TAG ($SOURCE_SHA)"
git -C "$REPO" archive --format=tar.gz -o "$STAGE/src.tar.gz" "$TAG"
if command -v sha256sum >/dev/null 2>&1; then
  SRC_SHA256="$(sha256sum "$STAGE/src.tar.gz" | awk '{print $1}')"
  SUMCHECK() { sha256sum -c "$@"; }
else
  SRC_SHA256="$(shasum -a 256 "$STAGE/src.tar.gz" | awk '{print $1}')"
  SUMCHECK() { shasum -a 256 -c "$@"; }
fi
# A release build must be reproducible from a committed lockfile, not from
# whatever cargo resolves today.
tar -tzf "$STAGE/src.tar.gz" | grep -qx 'Cargo.lock' ||
  die "$TAG does not contain Cargo.lock; a release cross-build requires one"

echo "==> host $CROSS_HOST  runid $RUNID  targets: $OS_LIST"
ssh -o ConnectTimeout=15 "$CROSS_HOST" "mkdir -p '$CROSS_WORKDIR/staging/$RUNID'" ||
  die "cannot reach build host $CROSS_HOST"
scp -q "$STAGE/src.tar.gz" "$CROSS_HOST:$CROSS_WORKDIR/staging/$RUNID/src.tar.gz" ||
  die "failed to stage source on $CROSS_HOST"

# ---------------------------------------------------------------------------
# Remote driver. Quoted heredoc: nothing here is expanded locally; every value
# arrives as a positional argument.
# ---------------------------------------------------------------------------
cat > "$STAGE/remote.sh" <<'REMOTE_DRIVER'
set -Eeuo pipefail
die() { echo "error: $*" >&2; exit 1; }

V="$1"; RUNID="$2"; SRC_SHA256="$3"; SOURCE_SHA="$4"; SOURCE_DATE_EPOCH="$5"
TAG="$6"; OS_LIST="$7"; WORKDIR="$8"; REBUILD_IMAGE="$9"

case "$WORKDIR" in /*) ;; *) WORKDIR="$HOME/$WORKDIR" ;; esac
STAGE="$WORKDIR/staging/$RUNID"
OUTD="$WORKDIR/out/$RUNID"
SYSROOTS="$WORKDIR/sysroots"
mkdir -p "$SYSROOTS" "$WORKDIR/cargo" "$WORKDIR/target" "$OUTD"

[[ "$(sha256sum "$STAGE/src.tar.gz" | awk '{print $1}')" == "$SRC_SHA256" ]] ||
  die "staged source digest does not match what was sent"

command -v docker >/dev/null 2>&1 || die "docker not on PATH on this host"
DOCKER=docker
docker info >/dev/null 2>&1 || DOCKER="sudo docker"
$DOCKER info >/dev/null 2>&1 || die "cannot talk to the docker daemon (tried docker and sudo docker)"

# --- builder image, content-addressed by its Dockerfile ---------------------
mkdir -p "$WORKDIR/image"
cat > "$WORKDIR/image/Dockerfile" <<'DOCKERFILE'
FROM rust:1-bookworm
# clang-22/lld-22 come from bookworm-security. clang 14 and 19 both segfault
# compiling ring's curve25519.c for aarch64-unknown-openbsd; 22 does not.
RUN apt-get update && apt-get install -y --no-install-recommends \
      clang-22 lld-22 curl xz-utils ca-certificates jq file \
 && rm -rf /var/lib/apt/lists/*
# All three targets are built with NIGHTLY, for two different reasons:
#   * aarch64-unknown-openbsd / -netbsd are Tier 3: no std is distributed at
#     all, so cargo has to compile it from rust-src via -Z build-std, which is
#     a nightly-only flag.
#   * aarch64-unknown-freebsd DOES have a prebuilt std, but as of rust 1.97.1 it
#     is published on the nightly channel only - `rustup target add
#     aarch64-unknown-freebsd` on stable fails with "has no prebuilt artifacts
#     available for target". So FreeBSD uses nightly's prebuilt std (fast, no
#     build-std). Re-check on future stables; if it lands there, FreeBSD can
#     move back to stable without touching anything else in this script.
RUN rustup toolchain install nightly --profile minimal --component rust-src \
 && rustup +nightly target add aarch64-unknown-freebsd
DOCKERFILE
IMG_TAG="$(sha1sum "$WORKDIR/image/Dockerfile" | cut -c1-12)"
IMG="ballistics-bsd-cross:$IMG_TAG"
if [[ "$REBUILD_IMAGE" == 1 ]] || ! $DOCKER image inspect "$IMG" >/dev/null 2>&1; then
  echo "==> building builder image $IMG"
  $DOCKER build -t "$IMG" "$WORKDIR/image" > "$WORKDIR/image/build.log" 2>&1 ||
    { tail -40 "$WORKDIR/image/build.log"; die "builder image build failed (full log: $WORKDIR/image/build.log)"; }
else
  echo "==> reusing builder image $IMG"
fi
IMG_ID="$($DOCKER image inspect -f '{{.Id}}' "$IMG")"

# --- sysroots: fetched once, stamped, reused --------------------------------
FREEBSD_REL=14.3
OPENBSD_REL=7.8
NETBSD_REL=10.1
# Bump when the extraction/symlink recipe below changes, so already-cached
# sysroots are rebuilt instead of silently keeping the old, broken layout.
SYSROOT_RECIPE=2
sysroot_dir() {
  case "$1" in
    freebsd) echo "$SYSROOTS/freebsd-$FREEBSD_REL" ;;
    openbsd) echo "$SYSROOTS/openbsd-$OPENBSD_REL" ;;
    netbsd)  echo "$SYSROOTS/netbsd-$NETBSD_REL" ;;
  esac
}

cat > "$WORKDIR/image/prep-sysroot.sh" <<'PREPSYSROOT'
#!/usr/bin/env bash
# Runs inside the builder container. $1 = os, $2 = release, $3 = sysroot dir,
# $4 = stamp string to write once the sysroot is known good.
set -Eeuo pipefail
die() { echo "error: $*" >&2; exit 1; }
OS="$1"; REL="$2"; SR="$3"; STAMP="$4"
DL="$(dirname "$SR")/.downloads"
mkdir -p "$SR" "$DL"

fetch() { # url dest
  [[ -s "$2" ]] && return 0
  echo "    fetching $(basename "$2")"
  curl -fSL --retry 3 --retry-delay 5 -o "$2.part" "$1" || die "download failed: $1"
  mv "$2.part" "$2"
}

need() { [[ -e "$SR/$1" ]] || die "$OS sysroot is missing $1 - extraction did not produce a usable sysroot"; }

case "$OS" in
freebsd)
  fetch "https://download.freebsd.org/releases/arm64/aarch64/${REL}-RELEASE/base.txz" "$DL/freebsd-$REL-base.txz"
  # Match the GUEST release: an older base links binaries that run on newer
  # systems, never the other way round.
  tar -xf "$DL/freebsd-$REL-base.txz" -C "$SR" ./lib ./usr/lib ./usr/include
  need lib/libc.so.7
  need usr/lib/libc.so       # an ld script GROUPing libc.so.7 + libc_nonshared.a
  need usr/lib/crt1.o
  need usr/include/stdio.h
  ;;
openbsd)
  fetch "https://cdn.openbsd.org/pub/OpenBSD/${REL}/arm64/base${REL//./}.tgz" "$DL/openbsd-$REL-base.tgz"
  fetch "https://cdn.openbsd.org/pub/OpenBSD/${REL}/arm64/comp${REL//./}.tgz" "$DL/openbsd-$REL-comp.tgz"
  tar -xzf "$DL/openbsd-$REL-base.tgz" -C "$SR" ./usr/lib
  tar -xzf "$DL/openbsd-$REL-comp.tgz" -C "$SR" ./usr/lib ./usr/include
  # THE OPENBSD SYMLINK STEP. OpenBSD ships only versioned shared objects
  # (libc.so.102.0) and no plain libc.so. lld has no OpenBSD versioned-library
  # search, so `-lc` finds nothing shared, silently falls back to libc.a, and
  # the link dies with "duplicate symbol: atexit" against crtbeginS.o. Give lld
  # the unversioned names it actually looks for.
  n=0
  for so in "$SR"/usr/lib/lib*.so.*; do
    [[ -e "$so" ]] || continue
    base="$(basename "$so")"
    link="${base%%.so.*}.so"
    [[ -e "$SR/usr/lib/$link" ]] || ln -s "$base" "$SR/usr/lib/$link"
    n=$((n + 1))
  done
  echo "    created unversioned .so symlinks for $n versioned libraries"
  (( n > 0 )) || die "openbsd sysroot has no lib*.so.* - the symlink fix found nothing to link"
  need usr/lib/libc.so
  need usr/lib/crtbeginS.o
  need usr/include/stdio.h
  ;;
netbsd)
  B="https://cdn.netbsd.org/pub/NetBSD/NetBSD-${REL}/evbarm-aarch64/binary/sets"
  fetch "$B/base.tar.xz" "$DL/netbsd-$REL-base.tar.xz"
  fetch "$B/comp.tar.xz" "$DL/netbsd-$REL-comp.tar.xz"
  # NetBSD already ships plain libc.so - it needs no OpenBSD-style .so symlink
  # farm. It does, however, need a DIFFERENT symlink (see below).
  tar -xf "$DL/netbsd-$REL-base.tar.xz" -C "$SR" ./lib ./usr/lib
  tar -xf "$DL/netbsd-$REL-comp.tar.xz" -C "$SR" ./usr/lib ./usr/include
  # THE NETBSD MACHINE-INCLUDE SYMLINK. The comp set ships the arch headers as
  # usr/include/aarch64/, but /usr/include/machine -> aarch64 is created by the
  # installer, not carried in the tarball. Without it <sys/cdefs.h> fails on
  # `#include <machine/cdefs.h>` and every C dependency (ring's curve25519.c
  # first) dies with "'machine/cdefs.h' file not found".
  [[ -e "$SR/usr/include/machine" ]] || ln -s aarch64 "$SR/usr/include/machine"
  need usr/lib/libc.so
  need usr/lib/crt0.o
  need usr/include/stdio.h
  need usr/include/machine/cdefs.h
  ;;
*) die "unknown os $OS" ;;
esac
echo "$STAMP" > "$SR/.stamp"
[[ -n "${HOST_UID:-}" ]] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" "$SR" "$DL"
exit 0
PREPSYSROOT
chmod +x "$WORKDIR/image/prep-sysroot.sh"

for os in $OS_LIST; do
  SR="$(sysroot_dir "$os")"
  case "$os" in freebsd) REL=$FREEBSD_REL ;; openbsd) REL=$OPENBSD_REL ;; netbsd) REL=$NETBSD_REL ;; esac
  WANT_STAMP="$REL/recipe$SYSROOT_RECIPE"
  if [[ -f "$SR/.stamp" && "$(cat "$SR/.stamp")" == "$WANT_STAMP" ]]; then
    echo "==> $os sysroot $WANT_STAMP cached"
    continue
  fi
  echo "==> preparing $os $WANT_STAMP sysroot"
  rm -rf "$SR"
  $DOCKER run --rm \
    -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
    -v "$SYSROOTS:/sysroots" \
    -v "$WORKDIR/image/prep-sysroot.sh:/prep.sh:ro" \
    "$IMG" bash /prep.sh "$os" "$REL" "/sysroots/$(basename "$SR")" "$WANT_STAMP" ||
    die "$os sysroot preparation failed"
done

# --- the cross build itself -------------------------------------------------
cat > "$WORKDIR/image/cross-build.sh" <<'CROSSBUILD'
#!/usr/bin/env bash
# Runs inside the builder container, as root. $1 = version, $2 = os.
set -Eeuo pipefail
die() { echo "error: $*" >&2; exit 1; }
V="$1"; OS="$2"

case "$OS" in
  # BUILD_STD=1 means "no std is distributed for this triple, compile it from
  # rust-src". FreeBSD has a prebuilt std on nightly, so it does not need it.
  freebsd) T=aarch64-unknown-freebsd; SR=$(echo /sysroots/freebsd-*); BUILD_STD_NEEDED=0 ;;
  openbsd) T=aarch64-unknown-openbsd; SR=$(echo /sysroots/openbsd-*); BUILD_STD_NEEDED=1 ;;
  netbsd)  T=aarch64-unknown-netbsd;  SR=$(echo /sysroots/netbsd-*);  BUILD_STD_NEEDED=1 ;;
  *) die "unknown os $OS" ;;
esac
[[ -d "$SR" ]] || die "no sysroot for $OS"

TU="$(echo "$T" | tr 'a-z-' 'A-Z_')"   # aarch64-unknown-openbsd -> AARCH64_UNKNOWN_OPENBSD
TL="$(echo "$T" | tr '-' '_')"         # aarch64-unknown-openbsd -> aarch64_unknown_openbsd

# Point BOTH rustc's linker driver and the cc-rs C compilations (ring, etc.) at
# clang-22 with the right target triple and sysroot. lld-22 is selected
# explicitly; the host's default GNU ld cannot link these targets.
export "CARGO_TARGET_${TU}_LINKER=clang-22"
export "CARGO_TARGET_${TU}_RUSTFLAGS=-C link-arg=--target=$T -C link-arg=--sysroot=$SR -C link-arg=-fuse-ld=lld-22"
export "CC_${TL}=clang-22"
export "CFLAGS_${TL}=--target=$T --sysroot=$SR"
export CARGO_HOME=/cargo
export CARGO_TARGET_DIR=/target
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-0}"

mkdir -p /work/src
tar -xzf /stage/src.tar.gz -C /work/src
cd /work/src

TOOLCHAIN="+nightly"
if [[ "$BUILD_STD_NEEDED" == 1 ]]; then
  BUILD_STD=(-Z build-std=std,panic_abort)
else
  BUILD_STD=()
fi

echo "--- toolchain for $OS"
cargo "$TOOLCHAIN" -V
rustc "$TOOLCHAIN" -V
clang-22 --version | head -1

# Deliberately NOT `set -o pipefail` here: we want tee's exit status to keep the
# pipeline alive so the specific cargo status below can be reported.
cargo "$TOOLCHAIN" build --locked --release --target "$T" "${BUILD_STD[@]}" 2>&1 | tee "/out/$OS-build.log"
rc=${PIPESTATUS[0]}
(( rc == 0 )) || die "cargo build failed for $OS (rc=$rc); see $OS-build.log"

BIN="/target/$T/release/ballistics"
[[ -x "$BIN" ]] || die "no binary at $BIN after a successful build for $OS"

# Prove we produced a foreign aarch64 object for the RIGHT OS, not a host
# binary that happens to sit in the target dir.
FILEDESC="$(file -b "$BIN")"
echo "$FILEDESC" | grep -q 'ARM aarch64' || die "$OS: not an aarch64 binary: $FILEDESC"
case "$OS" in
  freebsd) echo "$FILEDESC" | grep -qi 'freebsd' || die "$OS: binary is not FreeBSD-flavoured: $FILEDESC" ;;
  openbsd) echo "$FILEDESC" | grep -qi 'openbsd' || die "$OS: binary is not OpenBSD-flavoured: $FILEDESC" ;;
  netbsd)  echo "$FILEDESC" | grep -qi 'netbsd'  || die "$OS: binary is not NetBSD-flavoured: $FILEDESC" ;;
esac

mkdir -p "/out/$OS"
cp "$BIN" "/out/$OS/ballistics-$V-$OS-aarch64"
cp Cargo.lock "/out/$OS/Cargo.lock.used"
rustc "$TOOLCHAIN" -V > "/out/$OS/rustc.txt"
cargo "$TOOLCHAIN" -V > "/out/$OS/cargo.txt"
clang-22 --version | head -1 > "/out/$OS/compiler.txt"
printf '%s\n' "$FILEDESC" > "/out/$OS/file.txt"
printf '%s\n' "$T" > "/out/$OS/target.txt"
printf '%s\n' "$(basename "$SR")" > "/out/$OS/sysroot.txt"
# The container runs as root so rustup/cargo can write freely; hand the caches
# and outputs back to the invoking user or the next (non-root) step cannot even
# write the .sha256 next to the binary.
[[ -n "${HOST_UID:-}" ]] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" /out /cargo /target
echo "OK $OS: $FILEDESC"
CROSSBUILD
chmod +x "$WORKDIR/image/cross-build.sh"

HOSTNAME_REMOTE="$(hostname)"
for os in $OS_LIST; do
  echo "==> cross-building $os"
  start=$(date +%s)
  $DOCKER run --rm \
    -e SOURCE_DATE_EPOCH="$SOURCE_DATE_EPOCH" \
    -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
    -v "$SYSROOTS:/sysroots:ro" \
    -v "$STAGE:/stage:ro" \
    -v "$WORKDIR/cargo:/cargo" \
    -v "$WORKDIR/target:/target" \
    -v "$OUTD:/out" \
    -v "$WORKDIR/image/cross-build.sh:/cross-build.sh:ro" \
    "$IMG" bash /cross-build.sh "$V" "$os" ||
    die "cross build failed for $os"
  echo "    $os built in $(( $(date +%s) - start ))s"
done

# --- checksums + provenance -------------------------------------------------
BUILT_AT_UTC="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
CLANG_VER="$($DOCKER run --rm "$IMG" clang-22 --version | head -1)"
for os in $OS_LIST; do
  d="$OUTD/$os"
  name="ballistics-$V-$os-aarch64"
  [[ -f "$d/$name" ]] || die "missing built artifact for $os"
  digest="$(sha256sum "$d/$name" | awk '{print $1}')"
  size_bytes="$(stat -c '%s' "$d/$name")"
  # Two-space "<hash>  <file>" form. `shasum -c` rejects a bare digest, which
  # is exactly what the in-guest builds emit and build-k3s-bsds.sh has to
  # normalize afterwards. Emit the right thing at the source instead.
  printf '%s  %s\n' "$digest" "$name" > "$d/$name.sha256"
  lock_sha="$(sha256sum "$d/Cargo.lock.used" | awk '{print $1}')"

  jq -n \
    --arg package_version "$V" \
    --arg release_tag "$TAG" \
    --arg source_ref "$TAG" \
    --arg source_sha "$SOURCE_SHA" \
    --arg source_archive_sha256 "$SRC_SHA256" \
    --arg cargo_lock_sha256 "$lock_sha" \
    --arg os "$os" \
    --arg target "$(cat "$d/target.txt")" \
    --arg sysroot "$(cat "$d/sysroot.txt")" \
    --arg rustc "$(cat "$d/rustc.txt")" \
    --arg cargo "$(cat "$d/cargo.txt")" \
    --arg cc "$(cat "$d/compiler.txt")" \
    --arg image "$IMG@$IMG_ID" \
    --arg node "$HOSTNAME_REMOTE" \
    --arg run_id "$RUNID" \
    --arg built_at_utc "$BUILT_AT_UTC" \
    --arg name "$name" \
    --arg sha256 "$digest" \
    --arg filedesc "$(cat "$d/file.txt")" \
    --argjson size_bytes "$size_bytes" \
    --argjson source_date_epoch "$SOURCE_DATE_EPOCH" \
    '{
      schema_version: 1,
      project: "ballistics-engine",
      package_version: $package_version,
      release_tag: $release_tag,
      source_ref: $source_ref,
      source_sha: $source_sha,
      source_date_epoch: $source_date_epoch,
      source_archive_sha256: $source_archive_sha256,
      cargo_lock_sha256: $cargo_lock_sha256,
      lock_origin: "source",
      os: $os,
      arch: "aarch64",
      target: $target,
      profile: "release",
      locked: true,
      mode: "release",
      toolchain: {rustc: $rustc, cargo: $cargo, cc: $cc},
      builder: {
        kind: "docker-cross-x86_64",
        image: $image,
        node: $node,
        guest_os: null,
        run_id: $run_id,
        sysroot: $sysroot,
        cross_note: "Built on an x86_64 host; the build host cannot execute this binary. Runtime behaviour is established by validate-bsd-aarch64.sh inside the matching aarch64 KVM guest."
      },
      built_at_utc: $built_at_utc,
      tests: [
        {name: "locked-release-cross-build", status: "passed",
         command: "cargo build --locked --release --target \($target)"},
        {name: "target-object-format-check", status: "passed", command: ("file(1): " + $filedesc)},
        {name: "release-library-tests", status: "skipped",
         command: "not runnable on an x86_64 cross host"},
        {name: "cli-version-and-trajectory-smoke", status: "pending",
         command: "run validate-bsd-aarch64.sh to execute this binary on aarch64"}
      ],
      validation: null,
      artifacts: [{name: $name, sha256: $sha256, size_bytes: $size_bytes}]
    }' > "$d/$name.provenance.json"
  rm -f "$d/Cargo.lock.used" "$d/rustc.txt" "$d/cargo.txt" "$d/compiler.txt" \
        "$d/file.txt" "$d/target.txt" "$d/sysroot.txt"
done

rm -rf "$STAGE"
echo "REMOTE_OUT=$OUTD"
REMOTE_DRIVER

echo "==> running the cross build on $CROSS_HOST"
REMOTE_LOG="$STAGE/remote.log"
set +e
ssh -o ConnectTimeout=15 "$CROSS_HOST" 'bash -s' -- \
  "$V" "$RUNID" "$SRC_SHA256" "$SOURCE_SHA" "$SOURCE_DATE_EPOCH" "$TAG" \
  "$OS_LIST" "$CROSS_WORKDIR" "$CROSS_REBUILD_IMAGE" \
  < "$STAGE/remote.sh" 2>&1 | tee "$REMOTE_LOG"
rc=${PIPESTATUS[0]}
set -e
(( rc == 0 )) || die "remote cross build failed (rc=$rc)"

REMOTE_OUT="$(awk -F= '/^REMOTE_OUT=/{print $2}' "$REMOTE_LOG" | tail -1)"
[[ -n "$REMOTE_OUT" ]] || die "remote driver did not report an output directory"

echo "==> collecting artifacts into $OUT"
mkdir -p "$OUT"
for os in $OS_LIST; do
  name="ballistics-$V-$os-aarch64"
  for ext in "" .sha256 .provenance.json; do
    scp -q "$CROSS_HOST:$REMOTE_OUT/$os/$name$ext" "$OUT/" ||
      die "failed to fetch $name$ext from $CROSS_HOST"
  done
done

echo "==> verifying checksums locally"
( cd "$OUT" && for os in $OS_LIST; do SUMCHECK "ballistics-$V-$os-aarch64.sha256"; done ) ||
  die "checksum verification failed in $OUT"

n=0; for os in $OS_LIST; do n=$((n + 1)); done
cat <<EOF

OK: $n cross-built aarch64 BSD binaries + provenance in $OUT

NOT YET RELEASABLE. These binaries have never been executed. Run:

    scripts/release/validate-bsd-aarch64.sh $V $OUT

which boots each aarch64 BSD guest on the KVM node, gates --version == $V, and
runs a real trajectory solve. Only after it exits 0 are these release assets.
EOF
