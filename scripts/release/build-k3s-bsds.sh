#!/usr/bin/env bash
# aarch64 BSDs (freebsd/openbsd/netbsd) on the K3S cluster. MAKE TARGETS ONLY -
# never arm64-build.sh directly. Jobs run KVM guests pinned to orangepi5-max; do
# not co-schedule other host builds there. Each in-guest build hard-gates
# --version == VERSION and runs a trajectory smoke; provenance.json ships with
# each binary (0.29.0+ convention - keep shipping them).
# Known flake: netbsd guest networking; re-run just that OS:
#   make bsd-build PROJECT=ballistics-engine OS=netbsd MODE=release REF=v$V RELEASE_TAG=v$V VERSION=$V FETCH=0
set -euo pipefail
V="${1:?usage: build-k3s-bsds.sh VERSION}"; OUT="${2:-$HOME/release-$V}"; mkdir -p "$OUT"
K8S="${K8S_CLUSTER_DIR:-$HOME/projects/k8s-cluster}"
( cd "$K8S" && make bsd-release-build PROJECT=ballistics-engine VERSION="$V" RELEASE_TAG="v$V" )
for os in freebsd openbsd netbsd; do
  for ext in "" .sha256 .provenance.json; do
    cp "$K8S/build-artifacts/ballistics-$V-$os-aarch64$ext" "$OUT/"
  done
done
( cd "$OUT" && shasum -c ballistics-$V-*bsd*-aarch64.sha256 )
echo "OK: 3 aarch64 BSD binaries + provenance in $OUT"
