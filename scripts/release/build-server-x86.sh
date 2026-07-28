#!/usr/bin/env bash
# FALLBACK path: linux x86_64/aarch64, windows, and the three x86_64 BSDs via the
# build-server on 10.1.1.27. The PRIMARY path for these eight platforms is GitHub
# Actions (build-and-release.yml); use this when Actions is down or for local
# verification. Run the orchestrator ON the box (its remote mode is broken; and
# `python` is not on the non-interactive PATH there - use python3).
# KNOWN WART: docker builds have a hard 1800 s timeout in build_orchestrator.py
# (~line 558) while the QEMU-emulated linux-aarch64 build takes ~37 min - the
# orchestrator reports FAILURE while the container finishes and deposits a good,
# reproducible binary. Check the artifact before believing the failure.
# It clones MAIN HEAD (no tag pin): only run when main == the release tag.
set -euo pipefail
V="${1:?usage: build-server-x86.sh VERSION}"
ssh alex@10.1.1.27 'cd ~/projects/build-server && python3 scripts/agent_interface.py build --repo https://github.com/ajokela/ballistics-engine --platforms linux-x86_64,linux-aarch64,windows-x86_64,freebsd-x86_64,openbsd-x86_64,netbsd-x86_64'
echo "Artifacts land under alex@10.1.1.27:~/projects/build-server/builds/build_*/ - scp, rename to ballistics-$V-<platform>, verify --version where executable."
