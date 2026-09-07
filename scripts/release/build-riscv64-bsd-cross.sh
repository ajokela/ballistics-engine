#!/usr/bin/env bash
# freebsd/netbsd/openbsd on riscv64, cross-compiled on the x86_64 fleet runner.
#
# The aarch64 BSDs boot a golden guest under KVM and build natively inside it.
# That is not available here for two reasons: riscv64 guests get no KVM on an
# amd64 host, so an in-guest build is TCG emulation; and OpenBSD/riscv64 has no
# rustc in ports at all, so for that target a native build is impossible rather
# than merely slow. All three are therefore cross-compiled, ~3 min each.
#
# These are tier-3 Rust targets: rustup ships no std for them, so std is built
# from source with -Z build-std against a sysroot taken from the matching
# verification host. Harvest the sysroots with k8s-cluster's
# scripts/riscv64-sysroots.sh; they change only when a host OS is upgraded.
#
# Building is not shipping: validate-riscv64.sh must run each binary on the real
# OS before these become release assets.
set -euo pipefail
V="${1:?usage: build-riscv64-bsd-cross.sh VERSION [OUTDIR]}"
OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"

SYSROOTS="${RISCV64_SYSROOTS:-/home/alex/riscv-poc/sysroots}"
# Pinned by DIGEST, not by the floating :nightly-bookworm tag. -Z build-std is an
# unstable interface that has broken across nightlies before, so "whatever nightly
# the tag points at today" is not a build input a release can rest on. Today this
# lane only works because the runner happens to have this image cached; a
# `docker system prune` would silently change the compiler under it. Pinning also
# makes the shipped binary attributable: this digest is rustc 1.100.0-nightly
# (a69a63265 2026-09-03). Bump deliberately, and re-run the lane when you do.
IMAGE="${RISCV64_IMAGE:-rustlang/rust@sha256:429e94b9fc4c29a16ac6546821a749401620652c14b8dd5033064267b37ab9aa}"
W="${RISCV64_WORKDIR:-/tmp/riscv64-bsd-$V}"
PUB="${PUB:-ajokela/ballistics-engine}"

for os in freebsd netbsd openbsd; do
  [ -d "$SYSROOTS/$os" ] || { echo "missing sysroot $SYSROOTS/$os (run k8s-cluster scripts/riscv64-sysroots.sh)" >&2; exit 1; }
done

# REF defaults to the release tag. Overridable so the lane can be exercised on a
# branch before that tag exists -- otherwise this script is untestable until the
# moment it has to work.
REF="${RISCV64_REF:-v$V}"
# sudo, because the container runs as root (it apt-gets a toolchain) and leaves
# root-owned files under $W/target. A plain rm then fails on the SECOND run --
# the first release works and every rebuild after it dies in cleanup.
sudo rm -rf "$W"
git clone -q --depth 1 --branch "$REF" "https://github.com/$PUB" "$W"
cp "$W/scripts/release/cross-riscv64.sh" "$W/cross-riscv64.sh"

command -v jq >/dev/null 2>&1 || { echo "jq not on PATH (needed for provenance)" >&2; exit 1; }

# Reproducibility inputs, recorded the way the aarch64 lane records them.
SOURCE_SHA="$(git -C "$W" rev-parse HEAD)"
SOURCE_DATE_EPOCH="$(git -C "$W" log -1 --pretty=%ct)"
LOCK_SHA="$(sha256sum "$W/Cargo.lock" | awk '{print $1}')"
# The aarch64 lane pins by TAG and appends the resolved image id to make the
# provenance reference exact. This lane already pins by digest, so appending
# would emit "rustlang/rust@sha256:429e...@sha256:db89..." -- two digests and not
# a resolvable reference. Append only when the pin is not already a digest.
IMG_ID="$(sudo docker image inspect --format '{{.Id}}' "$IMAGE" 2>/dev/null || echo unknown)"
case "$IMAGE" in
  *@sha256:*) IMAGE_REF="$IMAGE" ;;
  *)          IMAGE_REF="$IMAGE@$IMG_ID" ;;
esac
RUNID="${GITHUB_RUN_ID:-local-$SOURCE_DATE_EPOCH}"
NODE="$(hostname)"
BUILT_AT_UTC="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"

for os in freebsd netbsd openbsd; do
  echo "==> $os riscv64"
  sudo docker run --rm \
    -v "$W":/src -v "$SYSROOTS":/sysroots:ro -w /src \
    -e SYSROOTS=/sysroots -e SRC=/src \
    "$IMAGE" bash -c "
      set -e
      apt-get update -qq >/dev/null
      apt-get install -y -qq clang lld llvm >/dev/null
      rustup component add rust-src >/dev/null 2>&1 || true
      bash /src/cross-riscv64.sh $os
      # Record the exact toolchain that produced THIS binary, from inside the
      # container that produced it -- asking the host afterwards would describe a
      # different machine.
      rustc -V  > /src/rustc.txt
      cargo -V  > /src/cargo.txt
      clang --version | head -1 > /src/compiler.txt" >/dev/null

  T="riscv64gc-unknown-$os"
  NAME="ballistics-$V-$os-riscv64"
  install -m 0755 "$W/target/$T/release/ballistics" "$OUT/$NAME"
  DIGEST="$(sha256sum "$OUT/$NAME" | awk '{print $1}')"
  SIZE="$(stat -c '%s' "$OUT/$NAME")"
  printf '%s  %s\n' "$DIGEST" "$NAME" > "$OUT/$NAME.sha256"

  # Same schema as build-bsd-aarch64-cross.sh, so both fleet lanes are readable
  # by one consumer. The runtime test stays "pending" and validation stays null
  # until validate-riscv64.sh actually runs the binary on the target OS -- a
  # cross-build on x86_64 cannot make any claim about runtime behaviour.
  jq -n \
    --arg package_version "$V" --arg release_tag "v$V" \
    --arg source_ref "$REF" --arg source_sha "$SOURCE_SHA" \
    --arg cargo_lock_sha256 "$LOCK_SHA" \
    --arg os "$os" --arg target "$T" \
    --arg sysroot "$SYSROOTS/$os" \
    --arg rustc "$(cat "$W/rustc.txt")" --arg cargo "$(cat "$W/cargo.txt")" \
    --arg cc "$(cat "$W/compiler.txt")" \
    --arg image "$IMAGE_REF" --arg node "$NODE" --arg run_id "$RUNID" \
    --arg built_at_utc "$BUILT_AT_UTC" --arg name "$NAME" --arg sha256 "$DIGEST" \
    --arg filedesc "$(file -b "$OUT/$NAME")" \
    --argjson size_bytes "$SIZE" --argjson source_date_epoch "$SOURCE_DATE_EPOCH" \
    '{
      schema_version: 1,
      project: "ballistics-engine",
      package_version: $package_version,
      release_tag: $release_tag,
      source_ref: $source_ref,
      source_sha: $source_sha,
      source_date_epoch: $source_date_epoch,
      source_archive_sha256: null,
      cargo_lock_sha256: $cargo_lock_sha256,
      lock_origin: "source",
      os: $os,
      arch: "riscv64",
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
        cross_note: "Built on an x86_64 host, which cannot execute this binary. These are tier-3 Rust targets with no prebuilt std, so std was compiled from source with -Z build-std against a sysroot taken from the matching verification host. Runtime behaviour is established by validate-riscv64.sh on the real OS."
      },
      built_at_utc: $built_at_utc,
      tests: [
        {name: "locked-release-cross-build", status: "passed",
         command: ("cargo build -Z build-std=std,panic_abort --locked --release --target " + $target)},
        {name: "target-object-format-check", status: "passed", command: ("file(1): " + $filedesc)},
        {name: "release-library-tests", status: "skipped",
         command: "not runnable on an x86_64 cross host"},
        {name: "cli-version-and-trajectory-smoke", status: "pending",
         command: "run validate-riscv64.sh to execute this binary on riscv64"}
      ],
      validation: null,
      artifacts: [{name: $name, sha256: $sha256, size_bytes: $size_bytes}]
    }' > "$OUT/$NAME.provenance.json"
done
rm -f "$W/rustc.txt" "$W/cargo.txt" "$W/compiler.txt"

echo "==> built:"
ls -l "$OUT"/ballistics-"$V"-*-riscv64
