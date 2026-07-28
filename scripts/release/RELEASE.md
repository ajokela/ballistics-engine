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

1. `cargo publish --locked` (idempotence guard: skip if the sparse index already has it).
2. Binaries -> `assemble.sh` (refuses on <13 binaries or any checksum mismatch;
   GH release + gs://ballistics-releases/X.Y.Z/).
3. `deploy-wasm.sh` (4 copies + both firebase sites + badge sed + live byte-verify).
4. Docs: `cd ~/projects/ballistics/ballistics-docs-site && ./update-docs.sh`.
5. Bindings: bump + tag ballistics-engine-py (its CI publishes to PyPI on tag);
   bump ballistics-engine-rb (3 spots: Cargo.toml version, Cargo.toml dep, gemspec —
   plus lib/ballistics_engine.rb VERSION), `gem build` + container-validate + `gem push`.
6. npm: passkey-gated, human-only until an automation token / OIDC trusted publishing
   is set up.
7. `verify-channels.sh X.Y.Z` — the release is done when this exits 0, not before.

## Standing gotchas (hard-won; do not relearn)

- fuzz/Cargo.lock holds TWO ballistics-engine pins (local + the 0.21.5 reference
  engine); bump.sh moves only the local one.
- build-server: run ON .27, `python3` not `python`, 1800 s docker timeout falsely
  fails the ~37 min emulated aarch64 build (artifact is still good + reproducible).
- K3S: make targets only; never co-schedule host builds on orangepi5-max.
- wasm-pack needs `--no-default-features`; firebase project is `ballistics-rs`
  (the engine CLAUDE.md's ballistics-buddy reference is stale).
- Registry front pages cache-lag: verify via index.crates.io, PyPI /simple/,
  and the RubyGems JSON API only.
- Hermetic tests: no $HOME, no network. A network-dependent test failed CI ON the
  0.30.0 bump commit after four green runs, looking exactly like a bump regression.

## Phase 3 (private orchestration)

The private repo `ajokela/ballistics-release` owns the one-button pipeline: fleet
jobs on a self-hosted runner (10.1.1.27, labels `fleet`), channel publishing with
secrets that never enter this public repo, and the final verify job. See that
repo's README for runner care and required secrets.
