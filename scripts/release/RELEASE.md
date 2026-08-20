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
| K3S cluster | 3x BSD aarch64 (+provenance) | `build-k3s-bsds.sh` / fleet runner |
| 10.1.1.26 (real RISC-V) | linux-riscv64 | `build-riscv.sh` |
| 10.1.1.10 (Orange Pi, cross) | linux-mips64el | `build-mips.sh` |
| 10.1.1.27 build-server | fallback for the 8 hosted ones | `build-server-x86.sh` |

Every binary must pass a `--version == X.Y.Z` gate ON its target (or under qemu for
MIPS). This is non-negotiable: stale binaries have shipped-adjacent twice.

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
