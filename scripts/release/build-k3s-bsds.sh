#!/usr/bin/env bash
# aarch64 BSDs (freebsd/openbsd/netbsd) on the K3S cluster. MAKE TARGETS ONLY -
# never arm64-build.sh directly. Jobs run KVM guests pinned to BSD_NODE; do
# not co-schedule other host builds on that node.
#
# BSD_NODE defaults to nanopct6 (RK3588: 8 cores, GICv3, KVM) since 2026-08-30,
# when orangepi5-max went down in a power outage. Any replacement MUST be
# big.LITTLE-compatible with the builder's `taskset -c 4-7` pin and expose
# gic-version=3 under KVM - the Raspberry Pi 5s have only 4 cores, so that pin
# fails outright there. The pin is a CORRECTNESS requirement, not a tuning knob:
# without it the guest dies in early firmware with a Synchronous Exception,
# because a vCPU migrating between A76 and A55 cores sees different CPU features.
# Golden images live at /opt/vms/base/<os>-aarch64.qcow2 ON THAT NODE.
#
# BSD_IMAGE_DIGEST is resolved here rather than by arm64-build.sh, whose own
# resolution shells into $ORANGEPI (`docker pull`) and therefore fails whenever
# that host is down - even when the build itself would have run fine elsewhere. Each in-guest build hard-gates
# --version == VERSION and runs a trajectory smoke; provenance.json ships with
# each binary (0.29.0+ convention - keep shipping them).
# Known flake: netbsd guest networking; re-run just that OS:
#   make bsd-build PROJECT=ballistics-engine OS=netbsd MODE=release REF=v$V RELEASE_TAG=v$V VERSION=$V FETCH=0
set -euo pipefail
V="${1:?usage: build-k3s-bsds.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
K8S="${K8S_CLUSTER_DIR:-$HOME/projects/k8s-cluster}"
export BSD_NODE="${BSD_NODE:-nanopct6}"
BSD_NODE_SSH="${BSD_NODE_SSH:-pi@10.1.1.31}"
# crictl, not ctr: only crictl reads the CRI registry auth, so `ctr images pull`
# fails with "no basic auth credentials" against registry.localnet.
# Match the REPO digest (repo@sha256:...) - inspecti also prints the config
# digest and image id, and a bare sha256 grep picks the wrong one.
if [ -z "${BSD_IMAGE_DIGEST:-}" ]; then
  BSD_IMAGE_DIGEST=$(ssh "$BSD_NODE_SSH" \
    "sudo k3s crictl pull registry.localnet/bsd-builder:latest >/dev/null 2>&1
     sudo k3s crictl inspecti registry.localnet/bsd-builder:latest 2>/dev/null" \
    | tr ',' '\n' | grep -o 'bsd-builder@sha256:[0-9a-f]\{64\}' | head -1 | sed 's/.*@//')
fi
[[ "$BSD_IMAGE_DIGEST" =~ ^sha256:[0-9a-f]{64}$ ]] ||
  { echo "could not resolve bsd-builder digest from $BSD_NODE_SSH"; exit 1; }
export BSD_IMAGE_DIGEST
echo "==> BSD_NODE=$BSD_NODE  digest=${BSD_IMAGE_DIGEST:0:19}..."
( cd "$K8S" && make bsd-release-build PROJECT=ballistics-engine VERSION="$V" RELEASE_TAG="v$V" )
# arm64-build.sh files artifacts per project and per tag, NOT at the top of
# build-artifacts/ (which holds only ballistics-engine/ and lattice/).
ART="$K8S/build-artifacts/ballistics-engine/v$V"
for os in freebsd openbsd netbsd; do
  for ext in "" .sha256 .provenance.json; do
    cp "$ART/ballistics-$V-$os-aarch64$ext" "$OUT/"
  done
done
# The in-guest builds emit a bare 64-hex digest with no filename, which
# `shasum -c` rejects outright. The published convention - what 0.29.0 and
# 0.30.0 shipped, and what the hosted platforms emit - is "<hash>  <file>".
# Normalize here so verification and the uploaded asset agree.
( cd "$OUT" && for os in freebsd openbsd netbsd; do
    f="ballistics-$V-$os-aarch64"
    printf '%s  %s\n' "$(tr -d ' \n' < "$f.sha256")" "$f" > "$f.sha256.tmp"
    mv "$f.sha256.tmp" "$f.sha256"
  done
  shasum -c ballistics-$V-*bsd*-aarch64.sha256 )
echo "OK: 3 aarch64 BSD binaries + provenance in $OUT"
