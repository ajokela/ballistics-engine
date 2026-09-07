#!/usr/bin/env bash
# Run each riscv64 release binary on the OS it targets, before it ships.
#
# A cross-build never executes what it produced: the runner is x86_64 and
# physically cannot. So every riscv64 asset is gated here, the same way
# validate-bsd-aarch64.sh gates the ARM ones.
#
# Per target: upload to a per-version path, verify the uploaded bytes IN the
# guest against the build's own .sha256, then run --version and a full solve.
# The solve matters because a wrong float ABI links cleanly and self-reports the
# right version -- it only shows up in the arithmetic.
#
# Hosts (all riscv64):
#   freebsd/openbsd/linux  QEMU guests on the emulator host, each reached through
#                          its own ~/vms/<os>-riscv64/ssh.sh wrapper, which carries
#                          the host-guest key and known_hosts
#   netbsd                 real silicon, a Milk-V Mars (StarFive JH7110)
set -euo pipefail
V="${1:?usage: validate-riscv64.sh VERSION [OUTDIR] [os...]}"
OUT="${2:-$HOME/release-$V}"
shift 2 2>/dev/null || shift $#
TARGETS=("$@"); [ ${#TARGETS[@]} -gt 0 ] || TARGETS=(linux freebsd netbsd openbsd)

VMHOST="${RISCV64_VMHOST:-alex@10.1.1.12}"
NETBSD="${RISCV64_NETBSD:-root@10.1.1.30}"

# The board's HPN-patched sshd advertises sntrup761x25519-sha512@openssh.com but
# stalls mid-exchange on it: the connection hangs at SSH2_MSG_KEX_ECDH_REPLY and
# only dies after its LoginGraceTime (600s), which reads as a dead host. OpenSSH
# 9.6+ picks that algorithm by default. Pin the classic exchange HERE, in the
# code, rather than in a runner's ~/.ssh/config -- a correctness requirement of a
# release gate must travel with the gate, or a rebuilt runner silently loses it.
NETBSD_SSH=(ssh -o BatchMode=yes -o ConnectTimeout=20 -o KexAlgorithms=curve25519-sha256)
VM_SSH=(ssh -o BatchMode=yes -o ConnectTimeout=20)

# Per-version, so a re-dispatch of the same tag cannot find a previous run's file.
REMOTE="ballistics-validate-$V"
SOLVE="trajectory --velocity 2700 --bc 0.243 --drag-model g7 --mass 175 \
--diameter 0.308 --auto-zero 100 --max-range 300"
MUZZLE=2700   # the solve's muzzle velocity; impact must be below it

# Run a shell snippet on the target OS. $1 is the snippet; stdin is forwarded.
remote() {
  local os="$1" snip="$2"
  case "$os" in
    netbsd)  "${NETBSD_SSH[@]}" "$NETBSD" "$snip" ;;
    linux)   "${VM_SSH[@]}" "$VMHOST" "/home/alex/vms/linux-riscv64/ssh.sh   '$snip'" ;;
    freebsd) "${VM_SSH[@]}" "$VMHOST" "/home/alex/vms/freebsd-riscv64/ssh.sh '$snip'" ;;
    openbsd) "${VM_SSH[@]}" "$VMHOST" "/home/alex/vms/openbsd-riscv64/ssh.sh '$snip'" ;;
    *) return 2 ;;
  esac
}

# Capture without letting errexit abort before the ::error:: is printed. A bare
# `x=$(cmd | grep ...)` under `set -euo pipefail` exits the script the moment the
# command fails or the grep does not match, so every fail=1/continue below it is
# unreachable and the job ends with a bare non-zero status and no annotation.
capture() { set +e; CAP_OUT=$("$@" 2>&1); CAP_RC=$?; set -e; }

fail=0
declare -a SEEN_OS=() SEEN_FPS=()

for os in "${TARGETS[@]}"; do
  BIN="$OUT/ballistics-$V-$os-riscv64"
  SUM="$BIN.sha256"
  if [ ! -f "$BIN" ]; then echo "::error::missing $BIN"; fail=1; continue; fi

  WANT_SHA=$(awk '{print $1}' "$SUM" 2>/dev/null || true)
  [ -n "$WANT_SHA" ] || WANT_SHA=$(sha256sum "$BIN" | awk '{print $1}')

  # Upload and checksum in ONE chain. `&&` throughout, so a failed write cannot
  # fall through to executing whatever was already at that path -- with `;` and a
  # fixed filename, a re-dispatched release could validate the PREVIOUS run's
  # binary and pass both gates while uploading a binary nothing ever ran.
  capture remote "$os" "rm -f ~/$REMOTE && cat > ~/$REMOTE && chmod +x ~/$REMOTE && \
    (sha256sum ~/$REMOTE 2>/dev/null || sha256 -q ~/$REMOTE 2>/dev/null || cksum -a sha256 ~/$REMOTE 2>/dev/null)" < "$BIN"
  if [ "$CAP_RC" -ne 0 ]; then
    echo "::error::$os riscv64: upload failed (rc=$CAP_RC): $CAP_OUT"; fail=1; continue
  fi
  GOT_SHA=$(printf '%s' "$CAP_OUT" | tr -d '\r' | grep -oiE '[0-9a-f]{64}' | head -1 || true)
  if [ "$GOT_SHA" != "$WANT_SHA" ]; then
    echo "::error::$os riscv64: binary changed in transit (want $WANT_SHA, got '${GOT_SHA:-none}')"; fail=1; continue
  fi

  capture remote "$os" "~/$REMOTE --version"
  GOT=$(printf '%s' "$CAP_OUT" | tr -d '\r' | tail -1)
  if [ "$CAP_RC" -ne 0 ] || [ "$GOT" != "ballistics $V" ]; then
    echo "::error::$os riscv64 version gate: got '$GOT' (rc=$CAP_RC), want 'ballistics $V'"; fail=1; continue
  fi

  capture remote "$os" "~/$REMOTE $SOLVE"
  if [ "$CAP_RC" -ne 0 ]; then
    echo "::error::$os riscv64: solve did not run (rc=$CAP_RC): $(printf '%s' "$CAP_OUT" | tail -3)"; fail=1; continue
  fi
  FPS=$(printf '%s' "$CAP_OUT" | grep -oE 'Impact Velocity: *[0-9.]+' | grep -oE '[0-9.]+$' | head -1 || true)
  if [ -z "$FPS" ]; then
    echo "::error::$os riscv64: could not parse Impact Velocity from the solve output"; fail=1; continue
  fi
  # Structural, not a pinned magic number: the engine's numerics change between
  # releases by design, and an absolute constant here would block a legitimate
  # release on all four platforms with no fix short of re-cutting the tag. A
  # wrong float ABI still cannot pass -- it would land outside these bounds or
  # disagree with the other targets below.
  if ! awk -v v="$FPS" -v m="$MUZZLE" 'BEGIN{exit !(v>0 && v<m)}'; then
    echo "::error::$os riscv64: impact velocity $FPS implausible (expected 0 < v < $MUZZLE)"; fail=1; continue
  fi
  # Optional hard pin, for a release that wants to assert an exact figure.
  if [ -n "${RISCV64_EXPECT_FPS:-}" ] && [ "$FPS" != "$RISCV64_EXPECT_FPS" ]; then
    echo "::error::$os riscv64: impact velocity $FPS != pinned $RISCV64_EXPECT_FPS"; fail=1; continue
  fi

  SEEN_OS+=("$os"); SEEN_FPS+=("$FPS")
  echo "  ok  $os riscv64: ballistics $V, $FPS fps"

  # Remove it. These four guests are long-lived and mutable -- unlike the aarch64
  # lane, which boots a disposable overlay of a read-only golden image per run --
  # so anything left behind persists into the next release. That persistence is
  # what made a stale-binary false pass reachable in the first place; the
  # versioned path and the in-guest checksum above close the hole, and this keeps
  # the guests from accumulating release binaries besides.
  capture remote "$os" "rm -f ~/$REMOTE"
  [ "$CAP_RC" -eq 0 ] || echo "  note: could not remove ~/$REMOTE on $os (rc=$CAP_RC)"
done

# Cross-target agreement. This is the property that actually catches a bad float
# ABI on one platform: every target ran identical inputs through identical code,
# so any disagreement is a defect regardless of what the absolute number is.
if [ "${#SEEN_FPS[@]}" -gt 1 ]; then
  for i in "${!SEEN_FPS[@]}"; do
    if [ "${SEEN_FPS[$i]}" != "${SEEN_FPS[0]}" ]; then
      echo "::error::riscv64 targets disagree: ${SEEN_OS[0]}=${SEEN_FPS[0]} vs ${SEEN_OS[$i]}=${SEEN_FPS[$i]}"
      fail=1
    fi
  done
  [ "$fail" -ne 0 ] || echo "  ok  all ${#SEEN_FPS[@]} targets agree at ${SEEN_FPS[0]} fps"
fi

[ "$fail" -eq 0 ] || { echo "::error::riscv64 validation failed"; exit 1; }
echo "==> validated on real systems: ${TARGETS[*]}"
