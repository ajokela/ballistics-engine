#!/usr/bin/env bash
# Execute the three cross-built aarch64 BSD binaries on REAL ARM hardware.
#
# This is the VALIDATION half of the lane whose build half is
# build-bsd-aarch64-cross.sh. The split exists because cross-compiling on the
# x86_64 build server takes ~3 min per OS instead of ~45 min in-guest, and it
# takes the single ARM node off the build critical path. But the old in-guest
# build proved something the cross-build cannot: that the binary actually RUNS.
# This script is where that proof now happens, and it is not optional. A
# cross-built binary that has never executed is not a release asset.
#
#   ./validate-bsd-aarch64.sh 0.35.1 ~/release-0.35.1
#
# For each OS it boots the matching golden guest under KVM on the ARM node and,
# inside that guest, runs:
#   1. `--version`, hard-gated to equal VERSION. A mismatch fails the run.
#   2. a real trajectory solve, checked for a plausible results table.
# Both must pass for every OS or this script exits non-zero.
#
# The solve is deliberately NOT asserted against exact numbers. Physics output
# legitimately shifts between engine versions, and a numeric equality gate would
# turn every intentional change into a false release blocker. What is asserted
# is that the solver ran on aarch64 and produced a self-consistent table -
# which is what "the binary works on this platform" actually means.
#
# HARDWARE CONSTRAINT - `taskset -c 4-7` is a CORRECTNESS requirement:
#   RK3588 is big.LITTLE (4x A76 + 4x A55). A vCPU that migrates between an A76
#   and an A55 sees a different CPU feature set mid-execution, and the guest
#   dies in early firmware with a Synchronous Exception. Pinning QEMU to the
#   A76 cluster (cores 4-7) is what makes the guest boot at all. Any
#   replacement node must have >=8 cores in that layout - a 4-core Raspberry Pi
#   5 cannot satisfy this pin.
#
# The pods carry `app: bsd-build` with the same required podAntiAffinity as the
# real build Jobs, so a validation guest can never time-slice cores 4-7 against
# a build guest, or against another validation guest.
#
# Golden images at /opt/vms/base/<os>-aarch64.qcow2 are READ-ONLY. Every run
# boots a disposable qcow2 overlay; nothing here can modify a golden image.
#
# Env overrides:
#   VALIDATE_NODE       k8s node with /dev/kvm + golden images (default nanopct6)
#   VALIDATE_NAMESPACE  namespace for the transient pods       (default default)
#   VALIDATE_IMAGE      builder image                (default registry.localnet/bsd-builder:latest)
#   VALIDATE_ONLY_OS    validate just one of freebsd|openbsd|netbsd
#   SSH_WAIT_SECS       guest boot budget, seconds              (default 600)
#   KUBECONFIG          as usual (default ~/.kube/config)
set -Eeuo pipefail

die() { echo "error: $*" >&2; exit 1; }

V="${1:?usage: validate-bsd-aarch64.sh VERSION ARTDIR}"
ART="${2:?usage: validate-bsd-aarch64.sh VERSION ARTDIR}"
[[ "$V" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-][A-Za-z0-9.]+)?$ ]] || die "bad VERSION: $V"
[[ -d "$ART" ]] || die "artifact directory does not exist: $ART"

NODE="${VALIDATE_NODE:-nanopct6}"
NS="${VALIDATE_NAMESPACE:-default}"
IMAGE="${VALIDATE_IMAGE:-registry.localnet/bsd-builder:latest}"
SSH_WAIT_SECS="${SSH_WAIT_SECS:-600}"
# How many times a guest may fail to BOOT before the OS is called failed. Boot
# flakes only; a failed gate is never retried.
BOOT_ATTEMPTS="${BOOT_ATTEMPTS:-3}"

if [[ -n "${VALIDATE_ONLY_OS:-}" ]]; then
  case "$VALIDATE_ONLY_OS" in freebsd|openbsd|netbsd) ;; *) die "VALIDATE_ONLY_OS must be freebsd, openbsd or netbsd" ;; esac
  OS_LIST="$VALIDATE_ONLY_OS"
else
  OS_LIST="freebsd openbsd netbsd"
fi

command -v kubectl >/dev/null 2>&1 || die "kubectl not on PATH"
command -v jq >/dev/null 2>&1 || die "jq not on PATH (needed to stamp validation into provenance)"
if command -v sha256sum >/dev/null 2>&1; then
  sha256_of() { sha256sum "$1" | awk '{print $1}'; }
else
  sha256_of() { shasum -a 256 "$1" | awk '{print $1}'; }
fi

kubectl get node "$NODE" >/dev/null 2>&1 || die "node $NODE not found (KUBECONFIG=${KUBECONFIG:-$HOME/.kube/config})"
node_ready="$(kubectl get node "$NODE" -o jsonpath='{range .status.conditions[?(@.type=="Ready")]}{.status}{end}')"
[[ "$node_ready" == "True" ]] || die "node $NODE is not Ready (status=$node_ready)"
kubectl -n "$NS" get secret bsd-build-ssh >/dev/null 2>&1 ||
  die "secret bsd-build-ssh missing in namespace $NS - the guests are unreachable without it"

RUNID="$(date +%s)-$(printf '%04x' $((RANDOM % 65536)))"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/bsd-validate.XXXXXX")"
POD=""
cleanup() {
  local rc=$?
  trap - EXIT HUP INT TERM
  [[ -n "$POD" ]] && kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=false >/dev/null 2>&1 || true
  rm -rf "$TMP"
  exit "$rc"
}
trap cleanup EXIT HUP INT TERM

# Per-OS state lives in files under $TMP rather than associative arrays: this
# script has to run under the bash 3.2 that ships as /bin/bash on macOS, where
# `declare -A` does not exist.

# Refuse before booting anything if an artifact is missing or corrupt: an hour
# of guest boots to discover a truncated scp is the wrong order of operations.
echo "==> pre-flight: artifacts in $ART"
for os in $OS_LIST; do
  name="ballistics-$V-$os-aarch64"
  [[ -f "$ART/$name" ]] || die "missing binary: $ART/$name"
  [[ -f "$ART/$name.sha256" ]] || die "missing checksum: $ART/$name.sha256"
  want="$(awk '{print $1}' "$ART/$name.sha256")"
  [[ "$want" =~ ^[0-9a-f]{64}$ ]] || die "$name.sha256 does not hold a sha256 digest"
  got="$(sha256_of "$ART/$name")"
  [[ "$want" == "$got" ]] || die "$name does not match its .sha256 (want $want, got $got)"
  printf '%s\n' "$got" > "$TMP/$os.expectsha"
  echo "    $name  $(echo "$got" | cut -c1-16)...  OK"
done

# ---------------------------------------------------------------------------
# In-pod runner. Boots one guest, ships the binary in, runs the two gates.
# ---------------------------------------------------------------------------
cat > "$TMP/run-in-guest.sh" <<'RUNNER'
#!/usr/bin/env bash
set -Eeuo pipefail
die() { echo "GUEST-ERROR: $*" >&2; exit 1; }

OS="${OS:?}"; VERSION="${VERSION:?}"; EXPECT_SHA="${EXPECT_SHA:?}"
SSH_WAIT_SECS="${SSH_WAIT_SECS:-600}"
BIN=/work/ballistics-bin
WORK=/work
SSH_PORT=2222

[[ -f "$BIN" ]] || die "binary was never delivered to the pod"
got="$(sha256sum "$BIN" | awk '{print $1}')"
[[ "$got" == "$EXPECT_SHA" ]] ||
  die "binary corrupted in transit to the pod (want $EXPECT_SHA, got $got)"
chmod 0755 "$BIN"

BASE="/opt/vms/base/${OS}-aarch64.qcow2"
[[ -f "$BASE" ]] || BASE="/opt/vms/base/${OS}-aarch64.img"
[[ -f "$BASE" ]] || die "no golden image for $OS on this node"
BASE_FMT="$(qemu-img info "$BASE" | awk '/file format/{print $3}')"

# Disposable overlay. The golden image is never opened for writing.
echo "==> overlay from $BASE ($BASE_FMT)"
qemu-img create -f qcow2 -b "$BASE" -F "$BASE_FMT" "$WORK/overlay.qcow2" >/dev/null
cp /usr/share/AAVMF/AAVMF_VARS.fd "$WORK/efivars.fd"

DISK2_ARGS=()
BASE2="/opt/vms/base/${OS}-aarch64-disk2.qcow2"
if [[ -f "$BASE2" ]]; then
  qemu-img create -f qcow2 -b "$BASE2" -F qcow2 "$WORK/overlay2.qcow2" >/dev/null
  DISK2_ARGS=(-drive "if=none,file=$WORK/overlay2.qcow2,format=qcow2,id=hd1"
              -device virtio-blk-pci,drive=hd1,vectors=0,romfile=)
fi

# Cloned NetBSD guests have no reusable on-disk entropy seed; without a virtio
# RNG sshd stalls before it can even emit a banner.
RNG_ARGS=()
[[ "$OS" == netbsd ]] && RNG_ARGS=(-object rng-random,id=rng0,filename=/dev/urandom
                                   -device virtio-rng-pci,rng=rng0,vectors=0)

# LogLevel=ERROR is required, not cosmetic: without it ssh writes
# "Warning: Permanently added ... to the list of known hosts." to stderr on
# every single connection, and any 2>&1 capture of the guest's output then has
# that warning as its first line. That is how an earlier revision read the
# --version gate as reporting "hosts.".
SSH=(ssh -i /ssh/id_build -p "$SSH_PORT" -o StrictHostKeyChecking=no
     -o UserKnownHostsFile=/dev/null -o LogLevel=ERROR -o ConnectTimeout=10)
SCP=(scp -i /ssh/id_build -P "$SSH_PORT" -o StrictHostKeyChecking=no
     -o UserKnownHostsFile=/dev/null -o LogLevel=ERROR)
SSH_READY=0

shutdown_guest() {
  if (( SSH_READY )); then
    "${SSH[@]}" root@127.0.0.1 \
      'shutdown -p now 2>/dev/null || halt -p 2>/dev/null || poweroff 2>/dev/null || true' \
      >/dev/null 2>&1 || true
  fi
  [[ -s "$WORK/qemu.pid" ]] && kill "$(cat "$WORK/qemu.pid")" >/dev/null 2>&1 || true
}
trap 'rc=$?; trap - EXIT; shutdown_guest; exit $rc' EXIT HUP INT TERM

# taskset -c 4-7 pins QEMU to the RK3588 A76 cluster. NOT a tuning knob: a vCPU
# that migrates onto an A55 sees different CPU features and the guest dies in
# early firmware with a Synchronous Exception.
echo "==> booting $OS guest (pinned to A76 cores 4-7)"
taskset -c 4-7 qemu-system-aarch64 \
  -machine virt,gic-version=3,highmem=off -accel kvm -cpu host \
  -m 2G -smp 4 \
  -drive if=pflash,format=raw,unit=0,readonly=on,file=/usr/share/AAVMF/AAVMF_CODE.fd \
  -drive "if=pflash,format=raw,unit=1,file=$WORK/efivars.fd" \
  -drive "if=none,file=$WORK/overlay.qcow2,format=qcow2,id=hd0" \
  -device virtio-blk-pci,drive=hd0,bootindex=0,vectors=0,romfile= \
  "${DISK2_ARGS[@]}" "${RNG_ARGS[@]}" \
  -netdev "user,id=n0,hostfwd=tcp::${SSH_PORT}-:22" \
  -device virtio-net-pci,netdev=n0,vectors=0,romfile= \
  -display none -serial "file:$WORK/console.log" \
  -pidfile "$WORK/qemu.pid" -daemonize

# Bounded by WALL CLOCK, not attempt count. QEMU user-mode networking accepts
# on the forwarded port whether or not anything inside the guest is listening,
# so a hung guest burns ConnectTimeout on every probe and an "8 attempts x 5s"
# loop silently becomes many minutes.
echo "==> waiting for guest ssh (up to ${SSH_WAIT_SECS}s)"
start=$SECONDS
deadline=$(( start + SSH_WAIT_SECS ))
while (( SECONDS < deadline )); do
  if "${SSH[@]}" root@127.0.0.1 true 2>/dev/null; then
    SSH_READY=1; echo "    guest up after $(( SECONDS - start ))s"; break
  fi
  sleep 5
done
if (( ! SSH_READY )); then
  echo "==> guest console (last 40 lines):"; tail -40 "$WORK/console.log" || true
  echo "==> one ssh probe with stderr shown:"
  "${SSH[@]}" -v -o ConnectTimeout=15 root@127.0.0.1 true 2>&1 | tail -25 || true
  die "$OS guest never came up within ${SSH_WAIT_SECS}s"
fi

"${SSH[@]}" root@127.0.0.1 'uname -a' > "$WORK/guest-uname.txt"
echo "==> guest: $(cat "$WORK/guest-uname.txt")"

# / is small on the OpenBSD and NetBSD golden images; use the same roomy
# filesystems the in-guest builder uses.
case "$OS" in
  freebsd) GDIR=/root ;;
  openbsd) GDIR=/usr/obj ;;
  netbsd)  GDIR=/usr/pkg ;;
esac
GBIN="$GDIR/ballistics-validate-$VERSION"

echo "==> shipping binary to guest $GBIN"
"${SSH[@]}" root@127.0.0.1 "mkdir -p '$GDIR' && rm -f '$GBIN'"
"${SCP[@]}" "$BIN" "root@127.0.0.1:$GBIN"
"${SSH[@]}" root@127.0.0.1 "chmod 0755 '$GBIN'"
# Prove the bytes that will run are the bytes that were built.
guest_sha="$("${SSH[@]}" root@127.0.0.1 \
  "sha256 -q '$GBIN' 2>/dev/null || sha256sum '$GBIN' 2>/dev/null | awk '{print \$1}' || cksum -a sha256 -n '$GBIN' 2>/dev/null | awk '{print \$1}'")"
guest_sha="$(echo "$guest_sha" | tr -d '[:space:]')"
if [[ "$guest_sha" =~ ^[0-9a-f]{64}$ ]]; then
  [[ "$guest_sha" == "$EXPECT_SHA" ]] ||
    die "binary changed in transit into the $OS guest (want $EXPECT_SHA, got $guest_sha)"
  echo "    in-guest digest matches"
else
  echo "    WARN: guest has no usable sha256 tool; skipping in-guest digest check"
fi

# ---- gate 1: --version, hard equality -------------------------------------
echo "==> gate 1: --version"
set +e
ver_out="$("${SSH[@]}" root@127.0.0.1 "'$GBIN' --version" 2>&1)"
ver_rc=$?
set -e
echo "    output: $ver_out"
(( ver_rc == 0 )) || die "$OS: --version exited $ver_rc: $ver_out"
# Anchor on the binary's own output line rather than "first non-empty line":
# anything the transport prints must not be able to satisfy the version gate.
ver_got="$(echo "$ver_out" | awk '/^ballistics[[:space:]]/{print $2; exit}')"
[[ -n "$ver_got" ]] ||
  die "$OS: --version produced no 'ballistics <version>' line: $ver_out"
[[ "$ver_got" == "$VERSION" ]] ||
  die "$OS VERSION GATE FAILED: binary reports '$ver_got', expected '$VERSION'"
echo "    VERSION GATE PASSED: $ver_got"

# ---- gate 2: a real trajectory solve --------------------------------------
echo "==> gate 2: trajectory solve"
set +e
solve_out="$("${SSH[@]}" root@127.0.0.1 \
  "'$GBIN' trajectory --bc 0.475 --velocity 2700 --mass 168 --diameter 0.308 --bore-height 1.5 --angle 0" 2>&1)"
solve_rc=$?
set -e
echo "$solve_out" | sed 's/^/    | /'
(( solve_rc == 0 )) || die "$OS: trajectory solve exited $solve_rc"

# Structural, not numeric. Engine versions legitimately move the numbers; what
# must hold is that the solver ran and emitted a coherent table.
for needle in 'TRAJECTORY RESULTS' 'Max Range' 'Impact Velocity' 'Impact Energy'; do
  echo "$solve_out" | grep -q "$needle" ||
    die "$OS: trajectory output has no '$needle' - not a results table"
done
field() { echo "$solve_out" | sed -n "s/.*$1:[[:space:]]*\([0-9.]\{1,\}\).*/\1/p" | head -1; }
range="$(field 'Max Range')"
vimp="$(field 'Impact Velocity')"
eimp="$(field 'Impact Energy')"
[[ -n "$range" && -n "$vimp" && -n "$eimp" ]] ||
  die "$OS: could not parse Max Range / Impact Velocity / Impact Energy from the table"
awk -v r="$range" -v v="$vimp" -v e="$eimp" '
  BEGIN {
    if (r <= 0)          { print "Max Range " r " is not positive"; exit 1 }
    if (v <= 0)          { print "Impact Velocity " v " is not positive"; exit 1 }
    if (v >= 2700)       { print "Impact Velocity " v " is not below the 2700 fps muzzle velocity"; exit 1 }
    if (e <= 0)          { print "Impact Energy " e " is not positive"; exit 1 }
  }' || die "$OS: trajectory result is not physically plausible"
echo "    SOLVE GATE PASSED: range=$range yd  impact=$vimp fps  energy=$eimp ft-lb"

printf '%s\n' "$ver_got" > "$WORK/result-version.txt"
printf '%s|%s|%s\n' "$range" "$vimp" "$eimp" > "$WORK/result-solve.txt"
tr -d '\n' < "$WORK/guest-uname.txt" > "$WORK/result-uname.txt"
echo "VALIDATION-OK $OS $ver_got"
RUNNER

overall_rc=0

for os in $OS_LIST; do
  name="ballistics-$V-$os-aarch64"
  POD="bsd-validate-$os-$RUNID"
  echo
  echo "================================================================"
  echo "  validating $name on $NODE"
  echo "================================================================"
  started=$(date +%s)

  # The FreeBSD guest boots only about half the time: when it does not, AAVMF dies
  # immediately with "Synchronous Exception at 0x00000000BFAF1A14" and the console
  # holds nothing but the firmware banner. Measured directly - the SAME binary was
  # validated twice minutes apart on an idle node, failing at 614s then passing at
  # 29s - so it is the guest, not the artifact and not the release.
  #
  # Retried here rather than by re-running the whole job, which re-does a clean
  # cross-build to get back to a coin flip.
  #
  # ONLY a guest that never came up is retried. A version gate or solve gate that
  # fails is a real defect in the binary and must stay fatal on the first failure -
  # retrying those would turn this script from a gate into a slot machine.
  boot_attempt=1
  while :; do

  # Idempotent: a pod left behind by an interrupted run is removed first.
  kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=true >/dev/null 2>&1 || true

  # ...and so is any validation pod from an EARLIER run. Pod names carry $RUNID,
  # so a stale pod never collides by name -- but every validation pod carries
  # `app: bsd-build` under a REQUIRED podAntiAffinity on hostname, so one left
  # behind still occupies the node's only slot and the new pod cannot schedule.
  # In v0.36.3 a run that died in ImagePullBackOff left its pods behind, and the
  # next run's freebsd -- validated first -- could not schedule, while openbsd
  # and netbsd passed once the stale pod had been reaped. The runbook already
  # claimed this cleanup happened; only the same-name case actually did.
  #
  # Matched by NAME PREFIX, not by the label: `app: bsd-build` is shared with the
  # real build Jobs, and deleting by label alone would kill an in-progress build.
  stale=$(kubectl -n "$NS" get pods -o name 2>/dev/null \
            | sed 's|^pod/||' \
            | grep '^bsd-validate-' \
            | grep -v -- "-$RUNID\$" || true)
  if [[ -n "$stale" ]]; then
    echo "==> reaping stale validation pod(s) holding the anti-affinity slot:"
    echo "$stale" | sed 's/^/      /'
    echo "$stale" | xargs -n1 kubectl -n "$NS" delete pod --ignore-not-found --wait=true >/dev/null 2>&1 || true
  fi

  cat > "$TMP/pod-$os.yaml" <<POD_YAML
apiVersion: v1
kind: Pod
metadata:
  name: $POD
  labels:
    app: bsd-build
    role: validate
    os: $os
spec:
  restartPolicy: Never
  automountServiceAccountToken: false
  terminationGracePeriodSeconds: 30
  nodeSelector:
    kubernetes.io/hostname: $NODE
  # Same required anti-affinity as the real build Jobs: only ONE guest may hold
  # cores 4-7 on this node at a time.
  affinity:
    podAntiAffinity:
      requiredDuringSchedulingIgnoredDuringExecution:
      - labelSelector:
          matchExpressions:
          - { key: app, operator: In, values: [bsd-build] }
        topologyKey: kubernetes.io/hostname
  containers:
  - name: validator
    image: $IMAGE
    # IfNotPresent, not Always: validation must not depend on registry
    # uptime. In v0.36.3 the cluster registry (zot) was down -- its pod
    # stuck Unknown on a half-dead node, its Longhorn volume unable to
    # move because Longhorn's own CSI pods were stuck on the same node --
    # and an always-pull turned that into a release blocker even though
    # the correct builder image was ALREADY cached on the node.
    # The image is a stable test harness, not the artifact under test:
    # what is validated is the freshly staged binary. Refresh the image
    # explicitly when it actually changes.
    #
    # NB: this comment lives inside an unquoted heredoc, where backticks
    # and dollar signs are still evaluated -- keep both out of it.
    imagePullPolicy: IfNotPresent
    # Override the builder image ENTRYPOINT, which would otherwise start a full
    # in-guest build. We only want its qemu/ssh/AAVMF toolbox.
    command: ["sleep", "infinity"]
    securityContext:
      privileged: true          # required to open /dev/kvm
      readOnlyRootFilesystem: true
    env:
    - { name: OS,             value: "$os" }
    - { name: VERSION,        value: "$V" }
    - { name: EXPECT_SHA,     value: "$(cat "$TMP/$os.expectsha")" }
    - { name: SSH_WAIT_SECS,  value: "$SSH_WAIT_SECS" }
    resources:
      requests: { memory: "4Gi", cpu: "2", ephemeral-storage: "4Gi" }
      limits:   { memory: "4Gi", cpu: "4", ephemeral-storage: "8Gi" }
    volumeMounts:
    - { name: vms,  mountPath: /opt/vms, readOnly: true }
    - { name: ssh,  mountPath: /ssh,     readOnly: true }
    - { name: kvm,  mountPath: /dev/kvm }
    - { name: work, mountPath: /work }
    - { name: tmp,  mountPath: /tmp }
  volumes:
  - name: vms
    hostPath: { path: /opt/vms, type: Directory }   # golden images, read-only
  - name: ssh
    secret: { secretName: bsd-build-ssh, defaultMode: 0400 }
  - name: kvm
    hostPath: { path: /dev/kvm, type: CharDevice }
  - name: work
    emptyDir: { sizeLimit: 8Gi }
  - name: tmp
    emptyDir: { sizeLimit: 1Gi }
POD_YAML

  kubectl -n "$NS" apply -f "$TMP/pod-$os.yaml" >/dev/null ||
    { echo "!! $os: could not create pod $POD"; overall_rc=1; break; }

  if ! kubectl -n "$NS" wait --for=condition=Ready "pod/$POD" --timeout=300s >/dev/null 2>&1; then
    echo "!! $os: pod $POD never became Ready"
    kubectl -n "$NS" describe pod "$POD" 2>&1 | tail -25
    kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=true >/dev/null 2>&1 || true
    # A pod that never goes Ready has not run the binary at all, so this says
    # nothing about the artifact: it is the node, the image, or the scheduler
    # (anti-affinity contention, a slow pull, transient node pressure). Same
    # class as a guest that never boots, so it gets the same budget instead of
    # being fatal on the first try. The GATES stay fatal on first failure --
    # retrying a version or solve failure would turn this script from a gate
    # into a slot machine.
    if (( boot_attempt < BOOT_ATTEMPTS )); then
      echo "== $os pod never became Ready (attempt $boot_attempt/$BOOT_ATTEMPTS); retrying"
      boot_attempt=$(( boot_attempt + 1 ))
      started=$(date +%s)
      continue
    fi
    echo "!! $os: pod never became Ready after $BOOT_ATTEMPTS attempts"
    overall_rc=1; break
  fi

  # Deliver via tar on stdin - the same mechanism `kubectl cp` uses, but with
  # the failure visible here instead of swallowed. The runner re-checks the
  # digest inside the pod, so a truncated stream cannot pass unnoticed.
  cp "$ART/$name" "$TMP/ballistics-bin"
  cp "$TMP/run-in-guest.sh" "$TMP/run-in-guest.sh.payload"
  # COPYFILE_DISABLE stops macOS tar from injecting ._AppleDouble members.
  if ! COPYFILE_DISABLE=1 tar -C "$TMP" -cf - ballistics-bin run-in-guest.sh.payload |
       kubectl -n "$NS" exec -i "$POD" -- tar -xf - -C /work; then
    echo "!! $os: failed to deliver the binary into $POD"
    kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=false >/dev/null 2>&1 || true
    overall_rc=1; break
  fi

  set +e
  kubectl -n "$NS" exec -i "$POD" -- bash /work/run-in-guest.sh.payload 2>&1 |
    tee "$TMP/$os.log"
  exec_rc=${PIPESTATUS[0]}
  set -e
  elapsed=$(( $(date +%s) - started ))
  printf '%s\n' "$elapsed" > "$TMP/$os.secs"

  if (( exec_rc != 0 )) || ! grep -q "^VALIDATION-OK $os " "$TMP/$os.log"; then
    if grep -q 'never came up within' "$TMP/$os.log" && (( boot_attempt < BOOT_ATTEMPTS )); then
      echo "== $os guest never booted (attempt $boot_attempt/$BOOT_ATTEMPTS); retrying"
      kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=true >/dev/null 2>&1 || true
      boot_attempt=$(( boot_attempt + 1 ))
      started=$(date +%s)
      continue
    fi
    echo "!! $os VALIDATION FAILED after ${elapsed}s (exec rc=$exec_rc)"
    overall_rc=1
  else
    grep "^VALIDATION-OK $os " "$TMP/$os.log" | tail -1 | awk '{print $NF}' > "$TMP/$os.ver"
    sed -n 's/^==> guest: //p' "$TMP/$os.log" | tail -1 > "$TMP/$os.uname"
    echo "== $os validated in ${elapsed}s"
  fi

  break
  done

  kubectl -n "$NS" delete pod "$POD" --ignore-not-found --wait=true >/dev/null 2>&1 || true
  POD=""
done

# ---------------------------------------------------------------------------
# Stamp the outcome into provenance. A cross-build alone leaves the runtime
# tests "pending"; only a passing run on aarch64 may mark them passed, and it
# records WHERE that happened.
#
# Stamped PER OS, not gated on the overall result: an OS that genuinely ran and
# passed on aarch64 should say so even when a sibling OS failed, so a re-run can
# target only what is still outstanding. An OS that failed keeps "pending" and
# is therefore still visibly unvalidated.
# ---------------------------------------------------------------------------
if ls "$TMP"/*.ver >/dev/null 2>&1; then
  echo
  echo "==> recording validation in provenance"
  stamped_at="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  for os in $OS_LIST; do
    [[ -s "$TMP/$os.ver" ]] || continue
    p="$ART/ballistics-$V-$os-aarch64.provenance.json"
    [[ -f "$p" ]] || { echo "    WARN: no provenance for $os, skipping"; continue; }
    jq --arg node "$NODE" --arg guest "$(cat "$TMP/$os.uname" 2>/dev/null || true)" \
       --arg at "$stamped_at" --arg ver "$(cat "$TMP/$os.ver")" \
       --argjson secs "$(cat "$TMP/$os.secs" 2>/dev/null || echo 0)" '
      .validation = {
        method: "qemu-kvm-aarch64-guest",
        node: $node,
        guest_os: $guest,
        validated_at_utc: $at,
        version_reported: $ver,
        elapsed_seconds: $secs
      }
      | .builder.guest_os = $guest
      | .tests = ([ .tests[] | if .name == "cli-version-and-trajectory-smoke"
            then .status = "passed"
               | .command = "executed on aarch64 under KVM: --version gate + trajectory solve"
            else . end ])
    ' "$p" > "$p.tmp" && mv "$p.tmp" "$p"
    echo "    $os: provenance stamped"
  done
fi

echo
echo "================================================================"
for os in $OS_LIST; do
  if [[ -s "$TMP/$os.ver" ]]; then
    printf '  %-8s PASS  --version=%-10s %ss  (%s)\n' "$os" "$(cat "$TMP/$os.ver")" \
      "$(cat "$TMP/$os.secs" 2>/dev/null || echo '?')" \
      "$(cat "$TMP/$os.uname" 2>/dev/null || echo 'unknown guest')"
  else
    printf '  %-8s FAIL\n' "$os"
  fi
done
echo "================================================================"

if (( overall_rc != 0 )); then
  echo
  echo "VALIDATION FAILED. These binaries must NOT be released."
  exit 1
fi
echo
echo "All aarch64 BSD binaries executed on real ARM hardware and passed the"
echo "--version gate ($V) and a live trajectory solve. Release assets are good."
