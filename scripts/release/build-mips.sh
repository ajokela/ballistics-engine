#!/usr/bin/env bash
# linux-mips64el, cross-compiled on the Orange Pi (alex@10.1.1.10) per MIPS_BUILD_NOTES.md.
# Build under /opt/mips-build (NVMe) - the home eMMC runs near-full and has failed releases.
# nightly + -Z build-std because rustup ships no rust-std for this target.
# Smoke test via qemu-user-static on the Pi, NOT the flaky MIPS VM.
set -euo pipefail
V="${1:?usage: build-mips.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
H=alex@10.1.1.10
ssh "$H" "set -e; source ~/.cargo/env
  export CARGO_TARGET_MIPS64EL_UNKNOWN_LINUX_GNUABI64_LINKER=mips64el-linux-gnuabi64-gcc
  # Exact path, NOT a glob: /opt/mips-build accumulates old per-version checkouts
  # (ballistics-engine-0.24.1 etc), and a multi-match glob makes cd fail with
  # 'too many arguments' (bit the 0.33.1 release).
  [ -d /opt/mips-build/ballistics-engine/.git ] || git clone https://github.com/ajokela/ballistics-engine /opt/mips-build/ballistics-engine
  cd /opt/mips-build/ballistics-engine
  git fetch --tags && git checkout v$V
  cargo +nightly build -Z build-std=std,panic_abort --target mips64el-unknown-linux-gnuabi64 --locked --release"
GOT=$(ssh "$H" "cd /opt/mips-build/ballistics-engine && qemu-mips64el -L /usr/mips64el-linux-gnuabi64 target/mips64el-unknown-linux-gnuabi64/release/ballistics --version")
[ "$GOT" = "ballistics $V" ] || { echo "VERSION GATE FAILED (qemu): '$GOT'"; exit 1; }
scp "$H":/opt/mips-build/ballistics-engine/target/mips64el-unknown-linux-gnuabi64/release/ballistics "$OUT/ballistics-$V-linux-mips64el"
( cd "$OUT" && shasum -a 256 "ballistics-$V-linux-mips64el" > "ballistics-$V-linux-mips64el.sha256" )
echo "OK: $OUT/ballistics-$V-linux-mips64el ($GOT under qemu)"
