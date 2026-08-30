#!/usr/bin/env bash
# linux-mips64el, per MIPS_BUILD_NOTES.md.
#
# This platform has ALWAYS been cross-compiled with a qemu-user smoke test - the
# notes' step 1 is literally "cross-compile on a fast Linux host with the
# MIPS64EL GNU cross toolchain". Only the host has moved: it was the Orange Pi
# (alex@10.1.1.10), and defaults to the x86_64 build-server since 2026-08-30,
# which is both faster (32 cores) and not a single point of failure with the
# BSD KVM host. Set MIPS_HOST to move it back.
#
# nightly + -Z build-std because rustup ships no rust-std for this target.
set -euo pipefail
V="${1:?usage: build-mips.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
H="${MIPS_HOST:-alex@10.1.1.27}"
W="${MIPS_WORKDIR:-/tmp/mips-build-$V}"
TARGET=mips64el-unknown-linux-gnuabi64
BIN="target/$TARGET/release/ballistics"

ssh "$H" "set -e
  rm -rf '$W' && git clone -q https://github.com/ajokela/ballistics-engine '$W'
  cd '$W' && git fetch --tags -q && git checkout -q 'v$V'
  sudo docker run --rm -v '$W':/src -w /src rust:1-bookworm bash -c '
    set -e
    apt-get update -qq >/dev/null && apt-get install -y -qq gcc-mips64el-linux-gnuabi64 >/dev/null
    rustup toolchain install nightly --profile minimal --component rust-src >/dev/null 2>&1
    export CARGO_TARGET_MIPS64EL_UNKNOWN_LINUX_GNUABI64_LINKER=mips64el-linux-gnuabi64-gcc
    cargo +nightly build -Z build-std=std,panic_abort --target $TARGET --locked --release'"

# Smoke-test under qemu-user, NOT the flaky MIPS VM.
GOT=$(ssh "$H" "sudo docker run --rm -v '$W':/src -w /src debian:bookworm-slim bash -c '
    apt-get update -qq >/dev/null && apt-get install -y -qq qemu-user-static libc6-mips64el-cross libgcc-s1-mips64el-cross >/dev/null 2>&1
    qemu-mips64el-static -L /usr/mips64el-linux-gnuabi64 $BIN --version'")
[ "$GOT" = "ballistics $V" ] || { echo "VERSION GATE FAILED (qemu): '$GOT'"; exit 1; }

scp -q "$H:$W/$BIN" "$OUT/ballistics-$V-linux-mips64el"
( cd "$OUT" && shasum -a 256 "ballistics-$V-linux-mips64el" > "ballistics-$V-linux-mips64el.sha256" )
echo "OK: $OUT/ballistics-$V-linux-mips64el ($GOT under qemu, host=$H)"
