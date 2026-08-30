# Release runbook

The releasing human's (or agent's) job is three decisions; everything else is a script
or a CI job. Codified from the 0.29.0/0.30.0 releases, which were driven by hand.

## The three human steps

1. **Decide**: `scripts/release/bump.sh X.Y.Z` -> commit -> push -> **wait for CI green**.
2. **Tag the green commit**: `git tag vX.Y.Z <sha> && git push origin vX.Y.Z`.
   The tag triggers `build-and-release.yml` (8 hosted platforms + GH release skeleton)
   and is the ref every fleet build checks out.
3. **Approve/triage**: watch the pipeline; re-run red jobs (netbsd guest-net flake
   clears on retry); run `verify-channels.sh X.Y.Z` last. Red verify = not released.

Asset count alone is not proof: 29 can hide a duplicate plus a gap. Check each of
the 13 platforms is present exactly once.

## Platform matrix (13)

| Where | Platforms | How |
|---|---|---|
| GitHub-hosted CI | macos x2, linux x86_64/aarch64, windows, 3x BSD x86_64 | `build-and-release.yml` on tag |
| **build** 10.1.1.27 cross (`CROSS_HOST`) + **validate** ARM KVM node (`VALIDATE_NODE`, default nanopct6) | 3x BSD aarch64 (+provenance) | `build-bsd-aarch64-cross.sh` **then** `validate-bsd-aarch64.sh` |
| K3S cluster (`BSD_NODE`, default nanopct6) | 3x BSD aarch64 — *fallback* | `build-k3s-bsds.sh` / fleet runner |
| 10.1.1.27 cross (or 10.1.1.26 native via `RISCV_MODE=native`) | linux-riscv64 | `build-riscv.sh` |
| 10.1.1.27 cross (`MIPS_HOST`) | linux-mips64el | `build-mips.sh` |
| 10.1.1.27 build-server | fallback for the 8 hosted ones | `build-server-x86.sh` |

Every binary must pass a `--version == X.Y.Z` gate ON its target (or under qemu for
MIPS). This is non-negotiable: stale binaries have shipped-adjacent twice.

**Degraded-fleet defaults (2026-08-30 power outage).** `orangepi5-max` (10.1.1.10)
and the RISC-V board (10.1.1.26) are down, so three lanes were re-pointed:
BSD aarch64 -> `BSD_NODE=nanopct6`, and RISC-V + MIPS -> cross-compiled in Docker
on 10.1.1.27 and gated under qemu-user. MIPS was ALWAYS cross-compiled, so only
its host moved; RISC-V previously built natively on real silicon, so its gate now
proves the binary self-reports correctly under emulation, NOT that it runs on real
hardware. Restore with `RISCV_MODE=native` and `BSD_NODE=orangepi5-max` once the
hosts are back.

## The aarch64 BSD lane: build on x86_64, validate on ARM

Two scripts, run in this order. Neither is optional.

```bash
scripts/release/build-bsd-aarch64-cross.sh X.Y.Z ~/release-X.Y.Z   # ~3 min/OS on 10.1.1.27
scripts/release/validate-bsd-aarch64.sh    X.Y.Z ~/release-X.Y.Z   # boots each guest on ARM
```

### Why it is split

The old lane (`build-k3s-bsds.sh`, still present as the fallback) compiles
*inside* each aarch64 BSD guest. That takes ~45 min per OS, and all three
guests run on one ARM KVM node, so the whole lane is ~2.5 h of critical path
gated on a single machine — a machine that has taken the release down more than
once (orangepi5-max went NotReady in the 2026-08-30 power outage; before that it
was CPU exhaustion during the Longhorn migration, and before that host-build
contention broke netbsd in 0.27.0).

Cross-compiling on the 32-core x86_64 build server takes ~3 min per OS. But
cross-compiling on its own would *lose* something the in-guest build gave us for
free: proof that the binary actually runs. A compiler that emits a valid ELF for
a target it cannot execute has told you nothing about whether the program works
there. So the runtime proof is not dropped, it is moved into its own step —
`validate-bsd-aarch64.sh` boots the matching golden guest under KVM on the ARM
node and runs the binary there:

1. `--version`, hard-gated to equal `X.Y.Z`. Mismatch fails the run.
2. a real trajectory solve, checked for a coherent results table.

The solve is deliberately **not** asserted against exact numbers — physics
output legitimately moves between engine versions, and a numeric-equality gate
would turn every intentional change into a false release blocker. The
`--version` gate is the strict one.

The build half never touches ARM hardware, so the validate half **must**. If you
find yourself tempted to ship cross-built binaries without running
`validate-bsd-aarch64.sh`, you have re-created the exact gap this split was
designed to make impossible.

### Toolchain constraints (all load-bearing)

- **clang-22 + lld-22, not older.** clang 14 *and* clang 19 both **segfault**
  (exit 139, frontend crash) compiling `ring`'s
  `crypto/curve25519/curve25519.c` for `aarch64-unknown-openbsd`. clang-22
  compiles it. Debian bookworm-security carries clang-22
  (`1:22.1.8-1~deb12u1`), so plain `apt-get install -y clang-22 lld-22` on
  `rust:1-bookworm` is enough — no apt.llvm.org repo needed.

- **All three targets build on the nightly toolchain**, for two different
  reasons. `aarch64-unknown-openbsd` and `aarch64-unknown-netbsd` are Tier 3:
  no std is distributed at all, so cargo compiles it from `rust-src` via
  `-Z build-std=std,panic_abort`, a nightly-only flag.
  `aarch64-unknown-freebsd` *does* have a prebuilt std, but as of rust 1.97.1 it
  is published on the **nightly channel only** — `rustup target add
  aarch64-unknown-freebsd` on stable fails with *"has no prebuilt artifacts
  available for target"*. Re-check on future stables; if it lands there,
  FreeBSD can move back to stable without touching anything else.

- **OpenBSD needs an unversioned-`.so` symlink farm in the sysroot.** OpenBSD
  ships *only* versioned shared objects (`libc.so.102.0`) and no plain
  `libc.so`. lld does not implement OpenBSD's versioned-library search, so given
  `-lc` it finds nothing shared, **silently falls back to the static `libc.a`**,
  and the link then dies with `duplicate symbol: atexit` against `crtbeginS.o` —
  an error that points nowhere near the actual cause. The sysroot prep creates
  an unversioned symlink for every `lib*.so.*` in `usr/lib` (42 of them at 7.8).
  With it, zero link errors. Do not remove this step to "clean up".

- **NetBSD needs `usr/include/machine -> aarch64` in the sysroot.** The comp set
  ships the arch headers as `usr/include/aarch64/`, but the `machine` symlink is
  created by the installer and is not carried in the tarball. Without it
  `<sys/cdefs.h>` fails on `#include <machine/cdefs.h>` and every C dependency
  dies with *"'machine/cdefs.h' file not found"*. NetBSD *does* already ship a
  plain `libc.so`, so it needs no OpenBSD-style symlink farm — the two BSDs need
  two different fixes, and neither implies the other.

- **Sysroots match the GUEST release, not the newest available**: FreeBSD 14.3,
  OpenBSD 7.8, NetBSD 10.1. A binary linked against an older base runs on a
  newer system; the reverse does not. If a golden image is upgraded, bump the
  matching `*_REL` in the build script.

- Changing the sysroot extraction recipe requires bumping `SYSROOT_RECIPE` in
  the build script, or already-cached sysroots keep their old, broken layout.

### Validation-host constraint

`taskset -c 4-7` in the QEMU invocation is a **correctness** requirement, not a
tuning knob. RK3588 is big.LITTLE (4x A76 + 4x A55); a vCPU that migrates from
an A76 onto an A55 sees a different CPU feature set mid-execution and the guest
dies in early firmware with a Synchronous Exception. Any replacement node must
have ≥8 cores in that layout — a 4-core Raspberry Pi 5 cannot satisfy the pin.
The validation pods carry `app: bsd-build` with the same required
podAntiAffinity as the build Jobs, so only one guest ever holds cores 4-7.

Golden images at `/opt/vms/base/<os>-aarch64.qcow2` are mounted read-only and
every run boots a disposable qcow2 overlay; nothing in this lane can modify one.

### Caching, idempotence, overrides

Both scripts are re-runnable. The builder image is content-addressed by its
Dockerfile, the three sysroots are downloaded once and stamped, and the cargo
registry plus target dir persist under `~/bsd-cross` on the build host. A repeat
run of one OS is ~30 s. Validation deletes any pod left over from an interrupted
run before creating its own.

| Env | Default | Effect |
|---|---|---|
| `CROSS_HOST` | `alex@10.1.1.27` | ssh destination of the x86_64 build server |
| `CROSS_WORKDIR` | `~/bsd-cross` | sysroot/cargo/target cache on that host |
| `CROSS_ONLY_OS` | all three | build one of `freebsd`/`openbsd`/`netbsd` |
| `CROSS_REBUILD_IMAGE` | `0` | force a builder-image rebuild |
| `VALIDATE_NODE` | `nanopct6` | k8s node with `/dev/kvm` + golden images |
| `VALIDATE_ONLY_OS` | all three | validate one OS |
| `SSH_WAIT_SECS` | `600` | guest boot budget |

### Falling back to the in-guest build

`build-k3s-bsds.sh X.Y.Z` still works and still produces the same artifact set.
Use it when the cross lane cannot: the build server is down, a new dependency
needs a sysroot piece that is not extracted, or a cross-built binary fails
validation in a way that suggests the cross toolchain rather than the code. It
is slower (~45 min/OS) and it builds *and* tests inside the guest, so it needs
no separate validation step — but it puts the ARM node back on the critical
path, which is the thing this lane exists to avoid. Note that it emits a **bare
64-hex digest** in its `.sha256` files and normalizes them afterwards; the cross
script emits the correct `<hash>  <file>` two-space form directly.

## Channels, in order

1. `prepublish-check.sh` FIRST, then `cargo publish --locked` (idempotence guard: skip
   if the sparse index already has it).

   The check refuses a non-pristine working tree, a HEAD that is not at a tag, a
   tag/Cargo.toml version disagreement, and any drift between the real packaged file
   set and `packaged-files.txt`. It exists because 0.33.1 shipped seven Applied
   Ballistics geometry files that were sitting UNTRACKED in the working directory:
   cargo refused the publish exactly as designed, and the refusal was overridden with
   `--allow-dirty`. **Never pass `--allow-dirty` to a publish.** If the tree is dirty,
   the fix is to commit, remove, or ignore the files — not to override the guard.

   When the shipped file set changes on purpose, regenerate the manifest in the same
   commit as the change:

   ```bash
   cargo package --list --locked | sort > scripts/release/packaged-files.txt
   ```
2. Binaries -> `assemble.sh` (refuses on <13 binaries or any checksum mismatch;
   GH release + gs://ballistics-releases/X.Y.Z/).
3. `deploy-wasm.sh` (4 copies + both firebase sites + badge sed + live byte-verify).
4. Docs: `cd ~/projects/ballistics/ballistics-docs-site && ./update-docs.sh`.
5. Bindings: bump + tag ballistics-engine-py (its CI publishes to PyPI on tag);
   bump ballistics-engine-rb (3 spots: Cargo.toml version, Cargo.toml dep, gemspec —
   plus lib/ballistics_engine.rb VERSION), `gem build` + container-validate + `gem push`.
6. npm: passkey-gated, human-only until an automation token / OIDC trusted publishing
   is set up.
7. **ballistics.tools** (the download hub) — the channel everyone forgets: it drifted from
   0.22.0 to 0.31.0 unnoticed because it was not listed here. Source is
   `~/projects/ballistics-tools-site` (NOT the stale duplicate under `ballistics.rs/`),
   deployed with `cd ~/projects/ballistics && firebase deploy --only hosting:tools`.
   Two real build artifacts live there, both of which must be REBUILT, not just relinked:
   - **Node WASM package** (`downloads/ballistics-wasm-node/` + its `.zip`):
     `scripts/build-wasm.sh --target nodejs --out-dir /tmp/wasm-node-X.Y.Z`,
     copy the five files over, re-zip. Smoke-test it with the site's OWN documented example
     before publishing — the API is a consuming builder, so every setter must be chained.
   - **Three BSD Python wheels** (FreeBSD/OpenBSD/NetBSD x86_64): PyPI cannot host these, so
     the hub is the only way BSD users get the Python binding. Built with maturin on the BSD
     hosts from `ballistics-engine-py` at the new version, so they are blocked on
     crates.io + the binding bump (step 5) and land LAST.
   Do NOT touch `downloads/bc5d*/manifest.json` — those carry the table/Flask versioning
   (0.34.x), not the engine's.
8. `verify-channels.sh X.Y.Z` — the release is done when this exits 0, not before.

## What the first automated release (0.30.1) taught

- **The three fleet jobs SERIALIZE, and mostly SHOULD.** One runner, one job slot:
  k3s-bsds + riscv + mips run back to back (~3.5 h). Do NOT "fix" this by adding
  runner slots: `orangepi5-max` is 10.1.1.10, i.e. the BSD KVM host and the MIPS
  cross-compile host are the SAME machine, and co-scheduling a host build with its
  VM jobs is what broke netbsd in 0.27.0. Only `riscv` (10.1.1.26, separate
  hardware) is safe to overlap, and it is ~19 min of a ~3.5 h critical path — so
  the serialization costs little and buys the invariant. Plan wall clock, not
  parallelism.
- **A partial fleet run does NOT need a rebuild.** Artifacts persist on a cluster
  volume reachable through the long-lived `artifacts-reader-*` pod:
  `kubectl exec`/`kubectl cp` from `/workspace/artifacts/ballistics-engine/v<tag>/`.
  This is how freebsd+openbsd survived a cancelled run intact.
- **Re-running one BSD: delete the stuck k8s Job first.** Two guests of the same OS
  must never co-schedule on orangepi5-max. Then
  `make bsd-build PROJECT=ballistics-engine OS=<os> MODE=release REF=vX.Y.Z RELEASE_TAG=vX.Y.Z VERSION=X.Y.Z`.
- **The cluster BSD builds emit a BARE 64-hex digest, no filename**, which
  `shasum -c`/`sha256sum -c` reject outright. Normalize to `<hash>  <file>` before
  verifying or uploading. Everything else (hosted, riscv, mips) already emits the
  two-space form; the release must be uniform.
- **npm lives in the wasm-pack output dir** that `deploy-wasm.sh` writes:
  `cd /tmp/wasm-X.Y.Z && npm publish` (after `npm login`; it is passkey-gated).
  Publish older versions first so `latest` ends up on the newest.

## Standing gotchas (hard-won; do not relearn)

- fuzz/Cargo.lock holds TWO ballistics-engine pins (local + the 0.21.5 reference
  engine); bump.sh moves only the local one.
- build-server: run ON .27, `python3` not `python`, 1800 s docker timeout falsely
  fails the ~37 min emulated aarch64 build (artifact is still good + reproducible).
- K3S: make targets only; never co-schedule host builds on orangepi5-max.
- Build WASM through `scripts/build-wasm.sh`, never `wasm-pack` directly. It defaults to
  the full terminal and verifies the emitted .wasm carries every gateable command. Built
  by hand you must pass `--no-default-features` AND `-- --features wasm-terminal`, and
  omitting the second ships a terminal with only `trajectory` and `version` — it builds,
  deploys, and passes the badge check, then fails on the first `zero`. Firebase project
  is `ballistics-rs` (the engine CLAUDE.md's ballistics-buddy reference is stale).
- Registry front pages cache-lag: verify via index.crates.io, PyPI /simple/,
  and the RubyGems JSON API only.
- Hermetic tests: no $HOME, no network. A network-dependent test failed CI ON the
  0.30.0 bump commit after four green runs, looking exactly like a bump regression.

## Phase 3 (private orchestration)

The private repo `ajokela/ballistics-release` owns the one-button pipeline: fleet
jobs on a self-hosted runner (10.1.1.27, labels `fleet`), channel publishing with
secrets that never enter this public repo, and the final verify job. See that
repo's README for runner care and required secrets.
