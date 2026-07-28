#!/usr/bin/env bash
# linux-riscv64 on REAL hardware root@10.1.1.26 (native, ~35 min; cargo+git present).
# NOT the emulated QEMU VM. This host once served a stale binary off a corrupted
# filesystem, which is why the --version gate here is the whole point.
set -euo pipefail
V="${1:?usage: build-riscv.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
H=root@10.1.1.26
ssh "$H" "set -e; cd ~/ballistics-engine || git clone https://github.com/ajokela/ballistics-engine ~/ballistics-engine && cd ~/ballistics-engine; git fetch --tags && git checkout v$V && cargo build --release --locked"
GOT=$(ssh "$H" "~/ballistics-engine/target/release/ballistics --version")
[ "$GOT" = "ballistics $V" ] || { echo "VERSION GATE FAILED: '$GOT'"; exit 1; }
scp "$H:~/ballistics-engine/target/release/ballistics" "$OUT/ballistics-$V-linux-riscv64"
( cd "$OUT" && shasum -a 256 "ballistics-$V-linux-riscv64" > "ballistics-$V-linux-riscv64.sha256" )
echo "OK: $OUT/ballistics-$V-linux-riscv64 ($GOT on-device)"
