#!/usr/bin/env bash
# Run each riscv64 release binary on the OS it targets, before it ships.
#
# A cross-build never executes what it produced: the runner is x86_64 and
# physically cannot. So every riscv64 asset is gated here, the same way
# validate-bsd-aarch64.sh gates the ARM ones. Two checks per target:
#
#   1. --version, hard-gated to equal VERSION. This platform once served a stale
#      binary off a corrupted filesystem, which is why the gate exists at all.
#   2. A full trajectory solve. A wrong float ABI or a soft-float mismatch links
#      cleanly and passes --version; it only shows up in the arithmetic. All four
#      targets must agree to the digit.
#
# Hosts (all riscv64):
#   freebsd/openbsd/linux  QEMU guests on the emulator host, ports 2222/2223/2224
#   netbsd                 real silicon, a Milk-V Mars (StarFive JH7110)
set -euo pipefail
V="${1:?usage: validate-riscv64.sh VERSION [OUTDIR] [os...]}"
OUT="${2:-$HOME/release-$V}"
shift 2 2>/dev/null || shift $# 
# Default to all four; a caller may gate just the targets its own job built, so
# the linux job and the BSD job can run independently.
TARGETS=("$@"); [ ${#TARGETS[@]} -gt 0 ] || TARGETS=(linux freebsd netbsd openbsd)

VMHOST="${RISCV64_VMHOST:-alex@10.1.1.12}"
NETBSD="${RISCV64_NETBSD:-root@10.1.1.30}"
SOLVE="trajectory --velocity 2700 --bc 0.243 --drag-model g7 --mass 175 \
--diameter 0.308 --auto-zero 100 --max-range 300"
EXPECT_FPS="${RISCV64_EXPECT_FPS:-2162.56}"

fail=0
for os in "${TARGETS[@]}"; do
  BIN="$OUT/ballistics-$V-$os-riscv64"
  [ -f "$BIN" ] || { echo "::error::missing $BIN"; fail=1; continue; }

  case "$os" in
    netbsd)  RUN() { ssh -o BatchMode=yes "$NETBSD" "cat > ~/ballistics; chmod +x ~/ballistics; ~/ballistics $1"; } ;;
    linux)   RUN() { ssh -o BatchMode=yes "$VMHOST" "/home/alex/vms/linux-riscv64/ssh.sh 'cat > /root/ballistics; chmod +x /root/ballistics; /root/ballistics $1'"; } ;;
    freebsd) RUN() { ssh -o BatchMode=yes "$VMHOST" "/home/alex/vms/freebsd-riscv64/ssh.sh 'cat > /tmp/ballistics; chmod +x /tmp/ballistics; /tmp/ballistics $1'"; } ;;
    openbsd) RUN() { ssh -o BatchMode=yes -o HostKeyAlias=openbsd-riscv64-qemu -J "$VMHOST" -p 2223 alex@127.0.0.1 \
                       "cat > ~/ballistics; chmod +x ~/ballistics; ~/ballistics $1"; } ;;
  esac

  GOT=$(RUN --version < "$BIN" | tr -d '\r')
  if [ "$GOT" != "ballistics $V" ]; then
    echo "::error::$os riscv64 version gate: got '$GOT', want 'ballistics $V'"; fail=1; continue
  fi
  FPS=$(RUN "$SOLVE" < "$BIN" | grep -oE 'Impact Velocity: *[0-9.]+' | grep -oE '[0-9.]+$')
  if [ "$FPS" != "$EXPECT_FPS" ]; then
    echo "::error::$os riscv64 solve gate: impact velocity $FPS, want $EXPECT_FPS"; fail=1; continue
  fi
  echo "  ok  $os riscv64: ballistics $V, $FPS fps"
done

[ "$fail" -eq 0 ] || { echo "::error::riscv64 validation failed"; exit 1; }
echo "==> validated on real systems: ${TARGETS[*]}"
