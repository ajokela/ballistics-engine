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
IMAGE="${RISCV64_IMAGE:-rustlang/rust:nightly-bookworm}"
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
      bash /src/cross-riscv64.sh $os" >/dev/null

  T="riscv64gc-unknown-$os"
  install -m 0755 "$W/target/$T/release/ballistics" "$OUT/ballistics-$V-$os-riscv64"
  ( cd "$OUT" && sha256sum "ballistics-$V-$os-riscv64" > "ballistics-$V-$os-riscv64.sha256" )
done

echo "==> built:"
ls -l "$OUT"/ballistics-"$V"-*-riscv64
