#!/usr/bin/env bash
# linux-riscv64.
#
# TWO MODES. Default is `cross` while the real board is unavailable.
#   RISCV_MODE=native  - build ON real hardware (~35 min). The original path; use
#                        it whenever a board is up, because it is the only one
#                        that exercises real RISC-V silicon. Set RISCV_HOST to
#                        the board: user@10.1.1.25 (Milk-V, Debian, 4 cores) as
#                        of 2026-09-14; root@10.1.1.26 has been down since the
#                        2026-08-30 outage.
#   RISCV_MODE=cross   - cross-compile in Docker on a fast x86_64 host and gate
#                        the binary under qemu-user. Added 2026-08-30 when the
#                        board was down after a power outage.
#
# The --version gate is NOT optional in either mode: this platform once served a
# stale binary off a corrupted filesystem, which is the whole reason it exists.
# In cross mode the gate runs under qemu-riscv64, so it proves the binary
# executes and self-reports correctly, NOT that it runs on real silicon.
set -euo pipefail
V="${1:?usage: build-riscv.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
MODE="${RISCV_MODE:-cross}"
TARGET=riscv64gc-unknown-linux-gnu
BIN="target/$TARGET/release/ballistics"

if [ "$MODE" = native ]; then
  H="${RISCV_HOST:-user@10.1.1.25}"
  # Do not rely on the board's PATH. A non-interactive ssh runs no login shell,
  # and Debian's stock ~/.bashrc returns before any rustup line it might carry,
  # so `cargo` is simply absent unless we put it there ourselves. Measured on
  # user@10.1.1.25: `ssh <host> cargo --version` -> "command not found" while an
  # interactive login finds 1.98.1.
  ssh "$H" "set -e
    export PATH=\"\$HOME/.cargo/bin:\$PATH\"
    [ -f \"\$HOME/.cargo/env\" ] && . \"\$HOME/.cargo/env\" || true
    command -v cargo >/dev/null || { echo 'no cargo on $H' >&2; exit 1; }
    cd ~/ballistics-engine 2>/dev/null || git clone https://github.com/ajokela/ballistics-engine ~/ballistics-engine
    cd ~/ballistics-engine
    git fetch --tags --force && git checkout -q 'v$V' && cargo build --release --locked"
  GOT=$(ssh "$H" "~/ballistics-engine/target/release/ballistics --version")
  SRC="$H:~/ballistics-engine/target/release/ballistics"
else
  H="${CROSS_HOST:-alex@10.1.1.27}"
  W="${CROSS_WORKDIR:-/tmp/riscv-build-$V}"
  ssh "$H" "set -e
    rm -rf '$W' && git clone -q https://github.com/ajokela/ballistics-engine '$W'
    cd '$W' && git fetch --tags -q && git checkout -q 'v$V'
    sudo docker run --rm -v '$W':/src -w /src rust:1-bookworm bash -c '
      set -e
      apt-get update -qq >/dev/null && apt-get install -y -qq gcc-riscv64-linux-gnu >/dev/null
      rustup target add $TARGET >/dev/null 2>&1
      export CARGO_TARGET_RISCV64GC_UNKNOWN_LINUX_GNU_LINKER=riscv64-linux-gnu-gcc
      cargo build --release --locked --target $TARGET'"
  GOT=$(ssh "$H" "sudo docker run --rm -v '$W':/src -w /src debian:bookworm-slim bash -c '
      apt-get update -qq >/dev/null && apt-get install -y -qq qemu-user-static libc6-riscv64-cross libgcc-s1-riscv64-cross >/dev/null 2>&1
      qemu-riscv64-static -L /usr/riscv64-linux-gnu $BIN --version'")
  SRC="$H:$W/$BIN"
fi

[ "$GOT" = "ballistics $V" ] || { echo "VERSION GATE FAILED ($MODE): '$GOT'"; exit 1; }
scp -q "$SRC" "$OUT/ballistics-$V-linux-riscv64"
( cd "$OUT" && shasum -a 256 "ballistics-$V-linux-riscv64" > "ballistics-$V-linux-riscv64.sha256" )
echo "OK: $OUT/ballistics-$V-linux-riscv64 ($GOT, mode=$MODE)"
