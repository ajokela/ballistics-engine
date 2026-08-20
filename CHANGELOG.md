# Changelog

All notable changes to the ballistics-engine project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **The WASM module's terminal commands are individually removable**, so an app that embeds the
  engine to solve trajectories no longer ships the twelve commands it never calls. Each
  non-trajectory command of the browser terminal (`WasmBallistics.runCommand`) now sits behind
  its own cargo feature — `wasm-zero`, `wasm-lead`, `wasm-monte-carlo`, `wasm-truing`,
  `wasm-estimate-bc`, `wasm-bc-convert`, `wasm-drag-curve`, `wasm-reticle`, `wasm-powder`,
  `wasm-recoil`, `wasm-power-factor` — with `wasm-terminal` turning on the whole set. Each one
  takes its handler, its exclusive helpers, and its help text (including its `Examples:` lines)
  out of the binary with it.

  `trajectory`, `version`, and the entire `Calculator` builder API are NOT gated: `Calculator`
  composes a `trajectory` command line internally, so every setter it exposes (BC, drag model,
  wind and wind segments, spin drift, Coriolis, atmosphere, sight height, zero range) keeps
  working with all eleven features off. `setZeroRange` in particular is unaffected by
  `wasm-zero` — it emits `--auto-zero` on `trajectory`, a different code path from the `zero`
  command.

  Measured on 0.33.2 (`--target web`, default release profile), against 917,924 bytes raw /
  344,592 gzip for the full set: `wasm-monte-carlo` 115 KB raw / 46 KB gzip (it takes the
  `--wez` sweep with it, which is a flag on `monte-carlo` rather than a command of its own),
  `wasm-truing` 92/31, `wasm-bc-convert` 65/22, `wasm-reticle` 55/21, `wasm-lead` 21/7,
  `wasm-powder` 17/6, `wasm-estimate-bc` 17/7, `wasm-zero` 15/4, `wasm-recoil` 12/3,
  `wasm-power-factor` 11/4, `wasm-drag-curve` 7/3. The commands share almost nothing, so the
  rows are close to additive — they sum to 427 KB raw / 153 KB gzip against a measured
  all-removed saving of 434 KB / 153 KB. A trajectory-only module is 483,496 raw / 191,421
  gzip / 154,797 brotli, 44% off the wire.

  Removing a command changes nothing the remaining ones compute. The full-terminal build is
  byte-identical to the ungated build across every command (the only difference is plain
  `monte-carlo`, which is unseeded and already nondeterministic against itself; the seeded
  `--wez` sweep matches byte-for-byte), and `Calculator` output is byte-identical between the
  full and trajectory-only builds across a matrix of drag-model, wind, wind-segment,
  spin-drift, Coriolis, atmosphere and zero configurations. A command compiled out answers
  `Unknown command`, and `help` advertises only what is present.

  These features are deliberately **not** in `default`. Every wasm32 build passes
  `--no-default-features` (the default `pdf`/`online` features don't compile for the target),
  so a default-on gate would be stripped by the one flag every affected build already uses. The
  shipped builds therefore pass `--features wasm-terminal` explicitly, and
  `scripts/build-npm.sh`, `scripts/release/deploy-wasm.sh` and `scripts/release/RELEASE.md`
  were updated accordingly. Note the bare `--` in those invocations: `wasm-pack` forwards only
  post-`--` arguments to cargo, so `--features` placed before it is consumed as an invalid
  `wasm-pack` flag and the build fails rather than silently mis-configuring.

  Splitting the help text into per-command chunks costs the full build ~3 KB raw / ~1 KB
  gzipped (35 `push_str` calls where there was one literal) — the price of the trimming table.

  Not affected, because they were never in the WASM build: `explain`, `error-budget`,
  `tolerance`, `dial-plan` and `adaptive-card` (0.33.x decision-support) are native-CLI-only.
  They have no WASM dispatch entry, and dead-code elimination already keeps them out of the
  module.

  **Build entry point.** All WASM builds now go through `scripts/build-wasm.sh`, which exists
  because the gate has one dangerous failure mode: a build that forgets the feature flag
  compiles, deploys, and passes a version check while answering `Unknown command` to everything
  but `trajectory`. The script closes it from both ends — its DEFAULT preset is `full`, so a
  missing argument yields the complete terminal rather than a stripped one, and every build is
  verified against the command set its preset promised by inspecting the emitted `.wasm`.
  Crucially the expectation for `--preset full` is "all twelve commands" derived from the preset
  itself, not from the computed feature string: deriving it from the feature list would mean
  losing that list also lowers the bar, and the verifier would agree that a stripped module is
  what `full` asked for. `deploy-wasm.sh`, `build-npm.sh`, `RELEASE.md` and CI all route through
  it; the CI step builds and verifies both presets, so a dropped flag fails there rather than on
  ballistics.rs.

  CI now type-checks both ends of the range (`wasm-terminal` and the bare build) plus each of
  the eleven features in isolation — a helper gated one feature too tightly compiles fine in
  the full set and fails only standalone. `wasm_tests.rs`'s 170-test terminal-parity suite
  requires `wasm-terminal` (it drives the whole command surface); a new always-compiled
  `minimal_surface_tests` module covers the ungated surface — the `Calculator` fluent API, each
  physics setter measurably moving the answer, `getFullTrajectory`, and a gated-out command
  reporting itself unknown — so a trimmed build is not shipped untested.
- **Bridge command `card.pdf`**: the engine's printable field dope card, returned as
  `{pdf_base64, byte_length, page_count, row_count, kind, source,
  unprintable_title_chars}` (`byte_length` describes
  the decoded document, not the base64 text). The request is the SAME `CardRequestV1` the
  on-screen card commands take, plus an optional presentation-only `pdf` block (`title`,
  `location`, `powder`, `bullet`, `target_speed` for the Lead column, `font_scale` **or**
  `font_preset`, `bold_data`) — so an app stores one request per saved card and replays it
  against `card.range_table` for the screen and `card.pdf` for the printout. The Range column
  follows the request's own distance unit (yards/metres), unlike `trajectory -o pdf`, whose
  card is always yards; the header/footer block is imperial on all three surfaces. Gated on
  the `pdf` feature and listed in `meta.capabilities` only when compiled in (a pdf-less build
  reports it as an unknown command) — but the `pdf` request block itself is accepted in every
  build, so a stored request still drives the on-screen card on an engine that cannot print
  it.
- **`card.pdf` can PRINT STORED ROWS instead of solving**, via one `card.pdf`-only request
  key: `stored_card: {card: <a stored card.range_table result, verbatim>, engine_version,
  bc5d_table_version}`. With it, nothing is solved, no trajectory runs, and
  `bc5d_table_path` is never opened — so a saved card reprints identically after an engine
  bump, after the correction table at that path is overwritten in place by a table-set
  refresh, and even after that file is deleted (pinned by
  `stored_rows_print_without_a_solve_or_a_correction_table`, which exports successfully
  against a table path that does not exist and shows the same request failing without the
  block). The footer's `BC:` is the stored card's own `bc_for_solve`; `source` in the
  response is `stored_rows` so a caller can verify it got a reprint. Omitting the key (or
  sending `null`) keeps the previous solve-then-print behaviour unchanged.
  The stored response, its `units` block and its rows are read back through mirrors that
  IGNORE fields this build has not heard of — deliberately not `deny_unknown_fields`, unlike
  `CardRequestV1`, where the caller really is the author. What arrives there is the engine's
  own output round-tripped through an app's storage, so strictness bought no validation and
  cost a hard cross-version break: the first field ever added to `CardResponseV1` would have
  made every card saved by the platform that had already bumped its engine pin unexportable
  on the platform that had not, with a serde message listing internal struct fields. The
  stored-side "this is not a range table" guard also now refuses `wind_angles_deg`, which
  the request-side guard already refused.
- **`card.pdf` reports a title its font cannot draw**, in `unprintable_title_chars` (the
  distinct characters, in order of first use; empty when the whole title printed). The card
  is drawn in Liberation Sans, which covers Latin, Latin-1, Cyrillic, Greek and the usual
  punctuation and nothing else — and printpdf emits NOTHING for an unmapped codepoint, so a
  card titled `射撃カード 308` printed a header reading `308`, an Arabic or Thai title
  printed a blank header, and an emoji left a gap, all with `ok: true` and a normal byte
  length. Every caller-supplied header/footer string now prints
  `pdf_dope_card::UNPRINTABLE_SUBSTITUTE` (`?`) where such a character stood, so the gap is
  visible, and `pdf_dope_card::unprintable_chars` is public so a caller can ask before or
  after printing. Not an error: the document and every row are unaffected.
- **Provenance on the paper.** The dope card footer now prints `Engine:<version>` and, for a
  reprint, `Table:<version>` beside the timestamp (`DopeCardConfig::engine_version` /
  `table_version`; the CLI's cards state the build too). An empty string prints nothing
  rather than a placeholder, so a card that used no correction table says nothing about one.
  A reprint states the versions the rows came from, which is what makes a printed card and a
  screen reconcilable afterwards.
- **`card.range_table` / `card.come_ups` / `card.wind` responses gain `bc_for_solve`**: the
  scalar BC the rows were actually computed with (the published BC unless a BC5D table
  applied its muzzle correction). Additive; it is what a saved card stores so a reprint can
  state the BC its numbers used.
- `pdf_dope_card::dope_card_rows_per_page` / `dope_card_page_count` expose the dope card's
  pagination, and `generate_dope_card_pdf` now paginates by them, so a caller reporting a
  page count reads the generator's own arithmetic instead of a copy that could drift.
- **BC5D offline correction on the mobile bridge.** Solve-json v1 gains an optional,
  additive `corrections` block with `bc5d_table_path`: a local path to a caliber-specific
  BC5D correction table (`bc5d_<caliber>.bin`, the exact dual-CRC format the CLI's
  `--bc-table-dir` consumes; mobile apps download tables themselves — the mobile engine
  build has no `online` feature by design). The service CRC-verifies the table, generates
  the same velocity-keyed BC segment ladder as the CLI, folds the muzzle correction into
  the scalar fallback BC, and applies the schedule to BOTH the zero solve and the flight
  (the 0.22.11 auto-zero rule). The block is echoed on `resolved_request` and carried by
  the resolved→request round-trip; segments always regenerate from the published BC, so
  re-solving never compounds the correction. Non-G1/G7 drag models are typed as G1 with a
  new stable `bc5d_drag_model_coerced` warning; on wasm32 the field is a structured
  `invalid_value` (no filesystem), never a silent no-op. Documented in
  `docs/SOLVE_JSON_V1.md`.
- **Card requests learn velocity-keyed BC.** `CardRequestV1` (bridge `card.come_ups`,
  `card.range_table`, `card.wind`) gains optional `bc_segments`
  (`[{velocity_min, velocity_max, bc}]`, velocities in the request's units system like the
  CLI's `--bc-segment`) and `bc5d_table_path` (same semantics as solve's corrections
  block). Explicit `bc_segments` win over the table, exactly like the CLI's precedence;
  the resolved schedule feeds the zero solve and every sampled run on all three card
  surfaces consistently.
- **Bridge command `bc5d.info`**: `{"path": ...}` → header metadata (`valid`, `crc_ok`,
  format version, caliber, generator API version, timestamp, per-axis bin counts, cell
  count, weight/velocity coverage) for a downloaded table, using the exact same
  load-with-verification path the solve/card fields use; corrupt or missing files are a
  clean `command_failed` envelope. Listed in `meta.capabilities` only on builds with
  filesystem access (not wasm32).
- Process-wide BC5D table cache (`bc_table_5d::path_cache`, non-wasm32): tables loaded by
  path are parsed and CRC-verified once, cached behind an `RwLock` keyed by
  (canonical path, file size, mtime) with a small bounded capacity, so repeated card/solve
  calls do not re-parse multi-MB tables; replacing a downloaded file invalidates naturally.
- **Bridge card commands** `card.come_ups`, `card.range_table`, and `card.wind`, backed by
  the new transport-free `card_service` module (`CardRequestV1`/`CardResponseV1`). These
  replicate the CLI card commands row-for-row — same zero solve, sampled trajectory,
  nearest-sample selection, and adjustment/bias/CF/click ordering — and the
  `card_bridge_golden` test suite executes the actual CLI binary against the bridge to
  prove it. Requests are denominated in the declared `units` system (imperial/metric),
  mirroring the CLI flags they wrap.
- The adjustment display layer (`AdjustmentUnit`, `drop_to_adjustment`,
  `adjustment_display`, `windage_adjustment_display`, `adjustment_unit_label`) moved from
  the CLI into the library's `adjustment` module so every surface — CLI, WASM terminal,
  bridge — shares one conversion/bias/CF/quantization boundary. CLI behavior unchanged.
- **Bridge profile commands** `profile.validate`, `profile.normalize`, and (behind the
  `profile-import` feature, and listed in `meta.capabilities` only when compiled in)
  `profile.import_a7p`. `validate` runs exactly the cheap invariants the CLI applies when
  loading/saving a profile (units string, tracking-CF band, turret/optic assembly +
  click-value parse) and returns `{valid, warnings, normalized}`; `normalize` returns the
  document re-serialized by this engine — the supported way for an app to round-trip a
  stored blob through a newer engine; `import_a7p` takes `{a7p_base64, zero_click?,
  strict?}` (decoded payload capped at 1 MiB, strict RFC 4648 decoding) and returns the
  mapped `ProfileData` plus the full import report (`warnings`/`mapped`/`unmapped`/
  `unknown_fields`) — the same mapping, warnings, and `--zero-click`/`--strict` semantics
  as the CLI's `profile import`, byte-for-byte (asserted by `tests/profile_bridge.rs`).
- `ProfileData` (with `ProfileZeroSet`/`ProfileBcSegment`/`ProfileDragPoint`, the
  turret-pair validators, and `optic_profile` assembly) moved verbatim from the CLI binary
  into the new library module `ballistics_engine::profile`, and the `.a7p` →
  `ProfileData` mapping (`map_a7p_to_profile` + `ImportReport`) into
  `ballistics_engine::profile_import`, so the CLI and the bridge share one definition.
  The serde wire shape of `~/.ballistics/profiles/*.json` is unchanged — same field
  names, defaults, null-vs-absent key behavior, and unknown-key tolerance, locked in by a
  round-trip fixture test in `src/profile.rs`. File persistence and profile unit
  conversion stay in the CLI; the new module is fs-free and compiles for wasm32.

### Changed
- **`card.pdf` states which card it prints, and refuses one it would misrepresent.** The
  response carries `kind: "range_table"`, and a request carrying a wind card's
  `wind_speeds`/`wind_angles_deg`, or a `stored_card` whose `kind` is not `range_table`, is
  now refused. Previously a stored `card.wind` request exported `ok=true` with a Wind column
  of `0.0` on every row while the screen showed drift to −0.42 MIL. A stored card whose
  `units` labels disagree with the request's own axes is refused for the same reason: the two
  are not the same card.
- **A card too big to print is refused from its row and page count**, before any document is
  built (`card_service::MAX_PDF_ROWS` = 5000, `MAX_PDF_PAGES` = 60 → `resource_limit`); the
  4 MiB byte cap remains the backstop for a card made huge by its labels. The refusal text
  now states what is true of the card — how many rows, how many pages — instead of advising
  "coarsen step, shorten the range domain", controls a saved card's snapshot does not have.
- The dope card's one-decimal-place rule for MIL/MOA/SMOA/IPHY cells is now documented as the
  contract for every surface that shows these rows (a turret's resolution is 0.1; two
  decimals on screen against one on paper put `2.45` against `2.4`, half a click apart), with
  the rounding vectors pinned — near-ties included, since rounding a double's shortest
  decimal spelling disagrees with rounding its binary value at values like `7.35`.

### Fixed
- **A BC5D table is now refused for a caliber it was not built for.** Nothing bound a
  table's CONTENT to the shot it corrected: the CRC only proves the bytes are intact, the
  CLI picks the file by NAME (`bc5d_<key>.bin`), and the bridge/solve-json/WASM surfaces
  take a caller-supplied path or raw bytes. Because every lookup clamps to the edge bins,
  a foreign table never failed — it returned a full, plausible-looking ladder that
  silently biased every row (the published `bc5d_224.bin` on a 175 gr G1 0.505 .308 at
  2600 fps gives segment BCs 0.4710..0.5114 where the .308 table gives 0.4989..0.5072,
  ~6.7 % off in the low bands; a .308 table on a .243 shot measures -17.9 % on effective
  BC). Every surface that has a shot in hand — bridge cards, bridge `solve`, `solve-json`,
  the CLI's `--bc-table-dir`, and the WASM `loadBc5dTable` path — now compares the table's
  header caliber against the bullet diameter and refuses a mismatch with an error naming
  both values ("table is for 0.224, shot is 0.308"), rather than correcting with foreign
  data or quietly downgrading to an uncorrected solve the caller cannot distinguish from a
  corrected one. Matching rule: equality of the 3-digit BC5D caliber key
  (`round(caliber * 1000)`) — the same key that selects `bc5d_<key>.bin`, so the CLI and
  every path/bytes consumer agree by construction; in practice a half-thousandth bucket
  (a `308` table accepts `[0.3075, 0.3085)`), which also absorbs the `f32` header's
  representation error. A MATCHING table is bit-identical to before (asserted against
  pre-change golden ladders and card/solve rows). The published tables all carry the right
  header caliber today, so this closes a latent identity gap (a rotated manifest, a
  hand-copied `.bin`, a future generator bug) rather than a live mis-serve.
  `bc5d.info` additionally reports `caliber_key`, the exact integer the guard compares, so
  an app can pre-check a downloaded table and show its own message. Library API:
  `Bc5dError` gains a `CaliberMismatch` variant, `Bc5dTable` gains `caliber_key()` and
  `ensure_caliber_matches()`, `bc_table_5d::caliber_to_key` is now public, and
  `path_cache::load_verified_for_caliber` is the guarded loader every shot-bearing
  path consumer uses (plain `load_verified` remains for `bc5d.info`, which describes a
  file rather than judging a shot).

## [0.33.2] - 2026-08-18

### Added
- **JSON command bridge** for embedded (mobile/FFI) consumers: one C entry-point set
  (`ballistics_bridge_call` / `_call_n` / `_free`) dispatching a versioned JSON envelope
  to transport-free library services — commands `meta.capabilities`, `meta.version`, and
  `solve` (delegating verbatim to the solve-json v1 service; the test suite asserts
  byte-identical output through both paths). All failures are in-band error envelopes:
  the calls never return NULL and never panic across the ABI. New default-on feature
  `bridge`; shipped header `include/ballistics_bridge.h`; mobile artifact builds via
  `scripts/build-mobile-ios.sh` / `scripts/build-mobile-android.sh`
  (`--no-default-features --features bridge,ffi,pdf,profile-import`, 16 KiB
  `max-page-size` on Android).

### Fixed
- serde_json now uses the `float_roundtrip` feature crate-wide: the default float parser
  can be one ULP off re-reading its own shortest output, which any JSON hop (bridge,
  solve-json round-trips) would silently inflict on solver inputs.
- Release tooling: `build-mips.sh` no longer trips over old per-version checkouts on the
  cross-build host; `extract-notes.sh` resolves the repo from its own location instead of
  the caller's working directory.

### Removed
- The Applied Ballistics geometry datasets and associated scripts that were inadvertently
  swept into the 0.33.1 release (and its crates.io package, which is yanked). This
  material is not part of the engine and has been withdrawn; `data/validation/` again
  contains only golden-physics case files.

## [0.33.1] - 2026-08-17

### Fixed
- Tip-off yaw drag is now additive per McCoy — `CD = CD0 + CD_delta2 * delta^2` — replacing
  the multiplicative `Cd * (1 + delta^2)` form whose implicit quadratic yaw-drag coefficient
  equaled `CD0` itself (~0.3 per rad^2), an order of magnitude below literature values
  (~4-20 per rad^2 for spitzer rifle bullets; 7.62 M80 = 9.6). New `cd_delta2` input
  (default 7.5 per rad^2, documented mid-range constant) controls the coefficient; the term
  is denominator-correct on both the G-model and custom-drag-table paths and is exactly
  inert when `tipoff_yaw` is zero, leaving all baseline trajectories bit-identical.
  Tip-off yaw remains modeled only in the RK4/RK45 generic integrator kernel
  (`derivatives::compute_derivatives`); the `TrajectorySolver` Euler/RK4 paths do not
  consume it (pre-existing scope). (MBA-1227)

## [Unreleased]

### Added
- **G1/G7 ballistic-coefficient drag-family conversion** (MBA-1375): the new `bc-convert`
  command converts a published scalar BC at an explicit Mach or velocity using the live
  reference curves (`BC_target = BC_source × Cd_target / Cd_source`). Repeated
  `--bc-segment VMIN:VMAX:BC` inputs convert an existing velocity-banded ladder and report
  fixed one-BC G1 and G7 least-squares fits plus the better-fitting family. Native and browser
  terminal surfaces share the same deterministic library calculation and table/CSV/JSON
  rendering; phase one intentionally accepts only G1 and G7.

### Fixed
- **WEZ variance attribution supports every built-in reference drag model** (MBA-1442):
  solve-json v1 and MCP now accept and advertise the full `G1`/`G2`/`G5`/`G6`/`G7`/`G8`/
  `GI`/`GS`/`RA4` set, restoring the `monte-carlo --wez` dominant-source and variance-share
  columns for the five models that previously reported `n/a`. Attribution remains unavailable
  for custom drag decks, ranges the nominal trajectory cannot reach, and structural sensitivity
  refusals; `p_hit` remains available in those cases. The JSON contract change is additive, but
  Rust consumers with an exhaustive match over the public `DragModelV1` enum must add arms for
  the five new variants.

## [0.33.0] - 2026-08-05

### Added
- **Decision-support analysis train**: three CLI subcommands built on a shared perturbation
  kernel over solve-json v1 requests:
  - `explain --a A.json --b B.json --ranges R1,R2,...` (MBA-1345) attributes the difference
    between two resolved firing solutions to a symmetric counterfactual swap of each of seven
    input groups (drag, muzzle velocity, zero/sight geometry, atmosphere, wind, shot geometry,
    effects), with any nonlinear interaction between groups reported once per range as an
    explicit interaction remainder rather than distributed across groups.
  - `error-budget --request FILE --ranges R1,... --sigma AXIS=VALUE [--target rect:WxH|circle:D]`
    (MBA-1347) propagates each declared per-input one-sigma uncertainty to impact covariance via
    central differences, ranks the sources by their own share of impact variance (never
    collapsed into an "other" bucket), and, with `--target`, reports the hit probability and
    each source's hit-probability gain if it alone were measured perfectly -- the
    value-of-information a shooter with time to improve exactly one input needs.
  - `tolerance --request FILE --range R --target ... --axis AXIS --domain AXIS=LO:HI`
    (MBA-1350) bisects one input at a time outward from its nominal value until the impact
    leaves the target, reporting a one-variable, no-probability-attached bound (or an explicit
    "no measurable effect" / "unbounded in domain" / "unavailable" when a bound does not apply).

  See [CLI_USAGE.md](CLI_USAGE.md#solution-diff-attribution-explain--mba-1345) for all three,
  worked examples, and each report's own stated limits.

  Underneath the CLI this is a new library layer: a shared input taxonomy of 32 perturbable axes
  (`ballistics_engine::perturbation`) with typed axis access and a central-difference/bisection
  kernel that re-solves the real trajectory through the same `solve_v1` path every other
  consumer of a solve-json v1 request uses; `ballistics_engine::explain`,
  `ballistics_engine::error_budget`, and `ballistics_engine::tolerance` publish versioned
  (`schema_version`) JSON reports on top of it. The kernel depends on a resolved solve-json v1
  request now being round-trippable back into a solvable request
  (`From<&ResolvedSolveRequestV1> for SolveRequestV1`) -- previously the resolved form was
  output-only.
- **Turret/hold decision-support train** (0.33.0 decision-support Plan B): saved profiles can
  now store a shooter's actual turret mechanics and reticle hold extent, and two new CLI
  subcommands turn that model plus a computed correction into ranked, executable plans or a
  compact, error-bounded range card:
  - **Turret-mechanics profile fields** (MBA-1348): `profile save` gains twelve turret/hold
    fields -- elevation/windage click graduation, clicks-per-revolution, zero-stop, elevation/
    windage travel remaining from the current zero, the turret's current dialed state, and the
    reticle's four-direction usable hold extent -- carried forward unchanged on re-save like
    every other optional profile field, with `--clear-turret`/`--clear-click` for intentional
    removal.
  - **`dial-plan --elevation VAL [--windage VAL] --range R (--profile NAME |
    --elevation-click ...)`** (MBA-1348) ranks whole-click dial, reticle-hold, and hybrid
    execution plans for a TRUE angular correction, on top of a new residual-carrying planner
    core, `ballistics_engine::optic::plan_corrections`: each plan reports its exact dial-space/
    true-angular split and its linear miss at range, and is never silently clamped -- an
    optic whose declared travel or hold bounds can't fully execute a plan gets an explicit
    `feasible: false` with the limiting mechanism named instead.
  - **`adaptive-card <load/profile args> --start R1 --end R2 [--anchor R]...
    [--elevation-budget VAL] [--windage-budget VAL] [--max-rows N] [--adjustment mil|moa]
    [-o table|csv|json|pdf]`** (MBA-1351) builds a range card that PROVABLY
    reconstructs the trajectory within a stated elevation/windage error budget: greedy
    worst-point insertion over the solved trajectory's own sample grid, then a separate dense
    pass MEASURES the finished card's true worst-case error and reports it -- together with
    the grid it was verified against -- in a footer (table/csv/pdf) or as report fields
    (`-o json`, `AdaptiveCardReportV1` pretty-printed verbatim). A budget tighter than the
    optic's own half-click quantization floor is honestly reported `budget_met: false`, never
    silently relaxed, and every requested `--anchor` always appears in the card regardless of
    what the measured error says. **Not a claim of fewer rows than a well-chosen fixed step**:
    measured against the most generous uniform-step baseline on a smooth trajectory, adaptive
    placement lost 10, tied 5, and won 5 across five domains and four budgets -- the value here
    is the measured error bound, the always-present anchors, and no step size to guess, never a
    shorter card.

  See [CLI_USAGE.md](CLI_USAGE.md#constrained-dial--hold-planning-dial-plan--mba-1348) for
  both commands, worked examples (including the `come-ups` -> `dial-plan` pairing), and each
  report's own stated limits.

  Underneath the CLI, four modules are new or newly public in the library:
  `ballistics_engine::optic` (the turret/hold profile model and `plan_corrections`),
  `ballistics_engine::hold_curve` (`HoldCurve`, the solved/sampled drop-vs-range curve
  `adaptive-card` and the existing reticle/BDC inverse solvers now share),
  `ballistics_engine::card` (the unified `CardRow` display-row type, see Changed below, plus
  the `adaptive_card` engine), and `ballistics_engine::pdf_dope_card` (the PDF dope card, now
  `CardRow`-based -- also see Changed below for what moved).

  None of the card, optic, or PDF changes in this bullet touch solve-json v1: no field, no
  accepted request shape, and no resolved-request behavior changed by anything here (see the
  Changed entries below for what solve-json v1 *did* gain this release).
- **Confidence-controlled Monte Carlo** (0.33.0 decision-support Plan C, MBA-1352): every
  *sampled* `monte-carlo` hit-probability estimate -- the fixed-count run and the new opt-in
  `--adaptive` mode, not the separate `--wez` sweep -- now states its sample count, method,
  confidence level, and interval:
  - The existing fixed-`n` `monte-carlo --target-distance` run gains an additive Wilson
    interval alongside the unchanged point estimate: a `Hit probability NN% CI: [lo, hi]
    (Wilson, n=N)` line under the existing `Hit Probability:` line, and an additive
    `hit_probability_ci` key (`method`, `confidence_percent`, `low`, `high`, `samples`) next to
    the unchanged `hit_probability` on `-o json`/`-o full`.
  - **`monte-carlo --adaptive [--confidence 90|95|99] [--target-ci-half-width P]
    [--min-samples N] [--max-samples N] [--mc-batch-size N] [--seed S]`** samples in batches
    until an anytime-valid confidence sequence's half-width reaches `--target-ci-half-width`
    (or `--max-samples` is exhausted) instead of a fixed `--num-sims` guess, and reports the
    point estimate, its interval, method, confidence level, stop reason, and three sample
    cardinalities (`attempts` drawn, `samples` that produced an outcome, `arrivals` that
    reached the target plane) as a versioned report, `AdaptiveMcReportV1` (`-o json`/`-o full`
    verbatim, `full` accepted as an alias for `json`), or a plain-text summary block (table
    output). Ignores `--num-sims` and `--wind-direction-std`; incompatible with `--wez`
    (`--adaptive --wez` is a usage error), which is a different report entirely and states no
    interval of its own, per spec. New `--seed` flag (reused by the fixed-count path above too)
    makes any run reproducible. `--confidence`/`--seed` are no-ops under `--wez` (its per-range
    `p_hit` carries no interval and it seeds itself internally) -- each prints a one-line stderr
    note there instead of doing nothing silently, and the sweep still runs and exits
    unchanged.

  See [CLI_USAGE.md](CLI_USAGE.md#confidence-controlled-sampling---adaptive--mba-1352) for the
  full flag table, captured examples, and the reported assumptions (sampling error only, never
  model error; anytime-valid stopping is specifically what makes stopping early on the
  confidence sequence legitimate, unlike repeated peeks at a fixed-`n` Wilson interval).

  Underneath the CLI this is a new public library layer: `ballistics_engine::mc_stats`
  (`Welford` streaming moments, `ConfidenceLevel`, `wilson_interval`, and the anytime-valid
  `BernoulliConfidenceSequence`) and `ballistics_engine::special::{ln_gamma, ln_beta}` (the
  log-gamma/log-beta functions the confidence sequence's mixture martingale is built on), plus
  `run_monte_carlo_adaptive_seeded`, `AdaptiveMcReportV1`, `McConvergence`, `McStopReason`,
  `ConfidenceLevel`, `wilson_interval`, `MC_ADAPTIVE_SCHEMA_VERSION_V1`,
  `MC_ADAPTIVE_METHOD_V1`, and `MC_ADAPTIVE_ASSUMPTIONS_V1` re-exported at the crate root, and
  two new methods on the existing `MonteCarloResults`: `hit_probability_wilson` and
  `position_is_hit`.

### Changed
- **WEZ (`monte-carlo --wez`) attribution now runs on the shared central-difference kernel**
  instead of its own one-sided finite differences. Published share values (`wind_call_share`,
  `mv_sd_share`, `other_share`) move slightly as a result -- about 5e-4 absolute on the
  project's own characterization fixture -- but bucket names and JSON field names are
  unchanged.
- **WEZ attribution now reports `n/a` (`attribution_unavailable: true`) for three
  configurations**: a loaded custom drag table or a drag model outside solve-json v1's
  `G1`/`G6`/`G7`/`G8` set (`G2`, `G5`, `GI`, `GS`, `RA4` -- neither has a solve-json v1
  representation at all), and now also whenever the shared attribution kernel hits a structural
  refusal on one of the shot's error sources for that particular range (an input combination the
  underlying sensitivity solve can't evaluate). `p_hit` is unaffected in all three cases; only
  the per-source variance-share attribution is unavailable.
- **solve-json v1 requests may now supply `shot.zero_distance_m` and `shot.muzzle_angle_rad`
  together**; this previously failed as `conflicting_fields`. The explicit angle now wins for
  elevation, the zero search is skipped, and a `zero_distance_elevation_not_resolved` warning
  notice is emitted -- `zero_distance_m` is still retained on the resolved request and still
  affects the windage-convergence bias, `summary.equivalent_horizontal_range_m`, and downrange
  wind-coverage validation. This is an additive widening of the accepted request shape scoped to
  these two fields only: no request that was previously valid because it did NOT combine
  `zero_distance_m` and `muzzle_angle_rad` changes behavior.
- **`resolved_request` now echoes seven more fields when the request explicitly supplies
  them**: `sight_offset_lateral_m`, `zero_poi_up_m`, `zero_poi_right_m`, `drops_reference`,
  `pressure_reference`, `wind_reference`, and `reticle` (from a `reticle hold` request) each gain
  a resolved counterpart -- needed so a resolved request is a complete, round-trippable
  description of the solve, which the new decision-support kernel above depends on. Every field
  skip-serializes when absent, so a request that never supplies any of the seven keeps a
  byte-identical response; a request that DOES supply one gains a key it did not have before --
  explicitly supplying a value equal to its own behavioral default (e.g. `wind_reference:
  "shooter"`, `sight_offset_lateral_m: 0.0`) is now distinguished from never having supplied the
  field at all, a distinction the response did not previously make.
- **Five hand-rolled bracket/interpolation searches (hold curves, wind-scenario corridors,
  reticle holds, the CLI's regular-interval trajectory sampler, and the inclined-shot
  equivalent-horizontal-range solver) now share one search-and-lerp core.** The regular-interval
  sampler's interpolated values -- `TrajectoryResult::sampled_points`, and therefore
  `come-ups`/range-table/wind-card/compare stdout, `wind_scenarios`' inputs, and the
  Python/Ruby/WASM bindings -- may shift by about 1 ULP, and up to 2 ULP in a pathological
  worst case (the re-association can move the delta term by roughly a ULP on each side): the
  shared core divides to get an interpolation fraction and then multiplies, where this one site
  previously multiplied then divided. Invisible at the display precision every existing fixture
  pins; the other four sites' arithmetic is unchanged.
- **`come-ups`/`range-table`/`wind-card`/`compare`'s four output surfaces share one internal
  row type**, `ballistics_engine::card::CardRow`, replacing four separate function-local
  structs that each said the same handful of things (range, drop, wind, velocity, energy,
  time) a different way. Internal only: **the four surfaces' output is byte-identical; only
  the internal row types were unified** -- pinned by a 12-case golden suite
  (`tests/card_golden_cli.rs`) covering table/CSV/JSON for all four, so binding maintainers do
  not need to re-verify any of their wire formats -- separately from the ULP-scale sampler
  shift noted above, which is invisible at these four surfaces' pinned display precision.
- **The PDF dope card moved from a `ballistics`-binary-private module into the library**,
  `ballistics_engine::pdf_dope_card` (behind the existing `pdf` feature), and rewritten onto
  `CardRow` above: `generate_dope_card_pdf` now takes `&[CardRow]` plus an explicit
  `RangeUnit` (`Yards`/`Meters`) parameter instead of the module's own binary-only
  `DopeCardRow { range_yd: u32, .. }`. The existing call site (`trajectory -o pdf`) is
  unaffected -- its rows are still rounded to whole yards before rendering and its PDF output
  is unchanged -- but the relocated function is now reachable from any library consumer
  (language bindings; `adaptive-card -o pdf` above) that enables the `pdf` feature.
- `pdf` is a default feature, so the relocated PDF dope card now compiles into the library for
  default-feature consumers (Python/Ruby/FFI builds), adding printpdf and ~825 KB of embedded
  fonts to those artifacts; `--no-default-features` builds are unaffected.
- **Fixed-count `monte-carlo` output gains an additive confidence-interval line**: the unseeded
  path is byte-identical to before, and the new `--seed` option's numeric estimates are pinned
  bit-for-bit reproducible for a given seed (see Added above for the new line/key itself). Note
  the two Monte Carlo paths' `std_*` fields are different estimators over different
  populations, not merely different numbers: the legacy fixed-count box/JSON `std_range`/
  `std_impact_velocity` are population standard deviations (divide by `n`) over all trials,
  while `AdaptiveMcReportV1`'s `std_*` fields are sample standard deviations (divide by `n-1`,
  `assumptions[3]`) over arrivals only.

### Fixed
- **`$.atmosphere.pressure_reference` was rejected as `unknown_field`** by
  `decode_solve_request_v1`'s hand-maintained shape validator ever since the field shipped
  (MBA-1397): its allowlist was never updated to match, even though direct `SolveRequestV1`
  construction and `resolve_atmosphere` always handled the field correctly. The whole
  QNH-pressure feature was unreachable through the JSON wire path as a result. Now accepted and
  performs the QNH-to-station-pressure reduction as documented.

## [0.32.0] - 2026-07-31

### Added
- **Import third-party Ventum reticle specs** (MBA-1440). `reticle import <ventum.json>` converts
  a Ventum reticle definition into the engine's `ReticleDescription` and emits the canonical
  schema that `reticle hold --reticle-json` and `profile save --reticle-json` already consume
  (table output previews the marks; `-o json` produces the interchange form). The converted
  description is validated before emit, so a decoration-only spec fails as `NoMarks` rather than
  producing an empty reticle.
- **Online reverse-solver subcommands** for the `ballistics` CLI, behind the default-on
  `online` feature. `login` saves a CLI access token from your ballisticsinsight.com account
  (the `BALLISTICS_API_TOKEN` environment variable overrides the on-disk
  `~/.ballistics/credentials.toml`, which is written `0600`) and `logout` clears it; the saved
  token is attached as a bearer credential when calling the hosted service:
  - `recommend-powder` — powder and charge recommendations for a cartridge, bullet weight, and
    desired muzzle velocity (`--barrel-length`, `--velocity-tolerance` optional).
  - `recommend-twist` — an optimal barrel twist rate from bullet caliber, weight, overall
    length, and nose length.
  - `recommend-col` — a recommended cartridge overall length for a cartridge (and optional
    bullet weight/type).
  - `calibrate-bc` — a bullet's ballistic coefficient estimated from its diameter, length, and
    weight.
  Each targets `https://api.ballistics.7.62x51mm.sh` by default (override with `--api-url`) and
  prints an actionable hint when the service asks the caller to sign in. These subcommands call
  an online service; every local trajectory and analysis subcommand is unchanged and needs no
  token or network access.

## [0.31.0] - 2026-07-29

### Added
- **Robust hold corridors across named segmented-wind scenarios** (MBA-1349). Field users
  usually have several concrete plausible wind calls, not a distribution they believe; a
  nominal trajectory hides that ambiguity and Monte Carlo demands a probability model they
  do not have. New `hold-corridor` subcommand and a new `wind_scenarios` module: give it a
  bounded set of NAMED scenarios (`WindScenarioSetV1`, whose segments are the same
  `SPEED:ANGLE:UNTIL[:VERTICAL]` tokens `--wind-segment` takes) plus the ranges you care
  about, and each scenario is solved ONCE through the existing segmented-wind machinery.
  At every range it reports every scenario's hold, the min/max corridor they span, the
  minimax hold, the worst-case miss from it (naming the scenario responsible), what
  holding the designated `nominal` would have cost instead, and whether one hold keeps
  every scenario on the target.
  **Two documented metrics**, because the shapes genuinely differ: a rectangle
  (`--target rect:WxH`, or no target at all) treats the axes independently, so the minimax
  hold is the per-axis midpoint of the extremes; a circle (`--target circle:D`) minimizes
  the largest Euclidean distance, so the hold is the center of the minimum enclosing
  circle — computed exhaustively over the candidate circles, which is exact and
  deterministic where the usual randomized construction is not. **Boundary contact counts
  as a fit.**
  **The zero is solved once, in calm air, and reused by every scenario** — a rifle has one
  zero, and re-zeroing per scenario would collapse the very corridor this exists to show.
  **Caps are enforced before any solving**: ≤ 8 scenarios and ≤ 64 ranges, alongside
  malformed segments, unknown `nominal` names, duplicate names and ranges, and unsupported
  versions — all structured, typed errors raised before a single trajectory runs.
  **Reordering scenarios cannot change the answer**: they are sorted by name internally
  before anything is solved.
  **No probabilities anywhere** (an explicit non-goal): nothing is weighted, nothing is
  interpolated between scenarios, and the finite set is never folded into a standard
  deviation. A three-scenario corridor is not a confidence interval and the output says
  so; statistical dispersion remains `monte-carlo --wez`'s job. `-o json` emits the
  versioned `RobustHoldReportV1`. Native-only this train, alongside MBA-1362's three
  solvers — the core and its shared formatter are already structured for the WASM
  follow-up.
- **Three reticle/BDC inverse solvers: `mark-to-range`, `bdc-match`, `optimal-zero`**
  (MBA-1362). Read-only solvers over an existing load — no new physics, no feature flag —
  sharing **one** drop-vs-range helper (a single solved, finely sampled angular-drop curve
  with a forward lookup and a monotone bisection inverse), so none of them can disagree
  with each other, or with `reticle hold --range`, about what "the drop at 500 yards" is.
  - **`mark-to-range`** (Nightforce / Nikon Spot On / Swarovski / TRACT) maps each
    subtension to the range where it lands, converting nominal marks to TRUE angular for
    the focal plane and magnification in use first. Marks that do not map to a range are
    **reported** — `not_a_holdover` for a mark at or above center, `beyond_max_range` with
    the deepest hold the load actually reaches — never silently dropped, so the table
    always carries one row per mark. Takes `--mark <MIL>` flags or a whole
    `--reticle-json` description.
  - **`bdc-match`** (Zeiss Rapid-Z) fits the magnification that makes an SFP BDC ladder
    match the load, from `--mark-range MIL:RANGE` pairs. Substituting `u = ref_mag ÷ mag`
    makes the residual linear, so the fit is the closed-form least-squares slope through
    the origin — exact and deterministic, with no starting guess. FFP is rejected with an
    explanation (its subtensions do not depend on magnification, so nothing is solvable)
    rather than returning a meaningless number, and a residual above `--residual-warn`
    (default 0.2 mil) prints a warning that the reticle does not fit this load at ANY
    magnification.
  - **`optimal-zero`** (GeoBallistics HDZ) min-max searches the one zero that minimizes
    the LARGEST hold a whole target list needs (`--target RANGE[:HEIGHT]`, 2 to 16, with
    `--vital` as the default height), then reports each target's hold, the signed miss a
    dead-center hold would produce, and whether every target stays inside its vital zone.
    Golden-section over a bracketed zero range with hard-coded constants, so the answer is
    reproducible bit for bit; targets are sorted internally so entry order cannot change
    it.
  All three are **native CLI only this train** — WASM parity for them is a tracked
  follow-up, matching the precedent that `mpbr`/`come-ups`/`range-table` are CLI-only.
- **Reticle hold points: schema, parametric generators, and a hold-point API**
  (MBA-1361; the Strelok Pro / Nikon Spot On / Hawke ChairGun / AB Quantum class, and the
  biggest UX gap we had against them for shooters who HOLD rather than dial). No reticle
  concept existed anywhere in the stack before this. New `reticle` subcommand with two
  verbs: `reticle generate {mil-grid|tree|bdc}` builds a description, and `reticle hold`
  places a firing solution in one and names the nearest mark. Take the solution directly
  (`--drop-mil` / `--wind-mil`) or let `--range` solve the trajectory for it.
  **FFP/SFP aware, by published optics-manual math:** a first-focal-plane mark subtends
  the same angle at every magnification, while a second-focal-plane mark's TRUE subtension
  is `nominal × reference-mag ÷ magnification` — a 2 mil mark on a reticle calibrated at
  10x covers 4 mil of target at 5x. The hold point is a property of the trajectory, so it
  stays true angular; only the MARKS are rescaled, and the nearest-mark distance is
  measured after that rescaling. An `off_reticle` flag fires when the hold runs outside the
  marks' bounding box grown by 20 % of its span per axis. Non-physical magnifications are
  rejected on both planes, with typed errors throughout.
  **One schema, every surface:** `ReticleDescription` (`serde`, FFP/SFP, reference
  magnification, marks with angular offsets and kinds) is shared by the CLI, an optional
  saved-profile `reticle` field (`profile save --reticle-json` / `--clear-reticle`,
  carried forward by an unrelated re-save like the DSF table), the browser terminal, and a
  new optional solve-json v1 request block (`reticle: {range_m, magnification,
  description}`) whose response adds `reticle_hold` — and only then, so every existing
  response stays byte-identical. `reticle generate -o json` emits exactly what
  `--reticle-json` consumes.
  **FFI is append-only:** a new `FFIReticleHold` struct and
  `ballistics_hold_point_in_reticle(...)`, with no change to any existing `repr(C)` layout,
  so no consumer needs a recompile. `marks_len` is validated against a stated bound before
  a single element is read (the MBA-1407 drag-table lesson applied).
  **Hard IP exclusions, permanent unless licensed:** no Horus/TREMOR-family grid layouts
  (`mil-grid` builds a plain mil-hash CROSS, not a filled 2-D grid), no Time-of-Flight
  wind-dot calibration (wind enters only as an angular deflection the caller already
  solved), and no vendor reticle catalog. The `bdc` generator takes ALREADY-SOLVED drops
  and runs no solve of its own, so a ladder's provenance stays the caller's.
  The browser terminal has the same two verbs through the *same* formatter (identical bytes
  for identical inputs), with two documented differences: `--reticle-json` there takes the
  description as inline JSON text rather than a path (no filesystem), and
  `reticle hold --range` / `--profile` are native-only.
- **Wind-call truing: back-solve the effective crosswind from an observed horizontal
  miss** (MBA-1392; Vortex Ace class, and the only competitor that had it). New
  `true-wind` subcommand on the native CLI and the WASM terminal: give it where your
  groups actually landed left/right of the aim point
  (`--miss RANGE:RIGHT[:SIGMA]`, repeatable) and it reports the constant crosswind that
  reproduces that miss through the real forward trajectory model, plus a wind-call
  correction factor (`solved ÷ called`) when you supply `--called-wind`. Each
  observation is fitted independently by a bracketed root find over
  ±100 mph — deflection is monotone and near-linear in crosswind speed (Didion) — and
  the combined answer is the plain mean, or an inverse-variance weighted mean with its
  own 1σ when every `--miss` carries a sigma (mixing the two is a hard error, not a
  silent blend). A miss no wind in that band can produce is rejected with a diagnostic
  naming the band rather than clamped into a plausible-looking number.
  **A horizontal miss is not purely wind, so the command separates it:** `--twist-rate`
  is REQUIRED and gyroscopic spin drift is always modelled and subtracted (without it a
  1:11" .308's ~3.5 in of right drift at 700 yd would be read as wind); `--latitude` and
  `--shot-direction` are optional but must be supplied TOGETHER, and add Coriolis.
  Whatever the model had no data for stays absorbed in the solved wind and is named in
  the report's "NOT subtracted" line, so a contaminated number is never presented as pure
  wind. Signs, in one documented block: `--miss` positive = impact RIGHT of aim; solved
  wind positive = wind FROM the shooter's LEFT (9 o'clock) pushing impacts right — the
  wind-FROM convention (0 = headwind) established by the 0.19.0 sign fix — and a
  right-hand twist drifts right. `RANGE` follows `--units`; the linear miss and its sigma
  are INCHES in both unit systems (the `--drop-unit in` precedent); wind speeds follow
  `--units` (mph/m·s⁻¹). Because `--miss` is a linear tape measurement and not a dial
  reading, the MBA-1358 scope tracking correction factors do NOT apply and the command
  exposes none. Fully offline (no API path exists), table/JSON/CSV, and native and WASM
  render through ONE shared formatter so the two surfaces cannot drift.
- **Truing forward model extended to wind, spin drift, and Coriolis** (MBA-1392,
  enabling the above). The two duplicated truing solver assemblies now take an optional
  `TruingEnvironment` (wind + `TruingTwist` + `TruingEarthFrame`), and the forward model
  samples the LATERAL axis (`position.z`) alongside drop — the old return value named
  `z` actually held the DOWNRANGE distance, which is now named `downrange_m` so the two
  cannot be confused again. Every effect is opt-in and only enabled when the data that
  makes it meaningful is present; `TruingEnvironment::default()` is the historical calm,
  effects-off model, so `true-velocity` (single-observation, multi-observation joint
  MV+BC, and uncertainty-aware), `plan-truing`, and `dsf` are byte-identical.
- **Earth-fixed compass wind bearings** (MBA-1368; Vortex Wind Bearing Capture /
  Lapua Ballistics class). New `--wind-ref {shooter|compass}` on `trajectory` and
  `monte-carlo` (native CLI and the WASM terminal's trajectory command): `shooter`
  (the default) is today's behavior byte-identically; `compass` treats EVERY wind
  direction the run consumes — `--wind-direction`, a location CSV's WIND_DIR, and
  each `--wind-segment` angle — as an absolute bearing (0 = north), re-referenced
  once at the input boundary as `relative = bearing - shot azimuth` (wind-FROM both
  sides, normalized to [0, 360)) so downstream physics never sees a bearing. Compass
  mode requires `--shot-direction` (hard error naming the flag — never a silent
  treat-as-shooter-relative; `monte-carlo` gains the flag for exactly this) and
  rejects clock positions (shooter-relative by definition). Monte Carlo converts the
  base direction BEFORE any dispersion sampling (the direction sigma disperses
  around the converted value). solve-json v1 gains the optional
  `wind.wind_reference` field (`"shooter"`/`"compass"`; compass requires an explicit
  `shot.shot_azimuth_rad` and folds the conversion into the resolved wind echo, the
  QNH precedent), mirrored in the MCP schema. The WASM params builder gains additive
  `setWindReference(mode)`/`setShotDirection(deg)` methods (no existing signature
  changed). FFI deliberately unchanged: bindings stay shooter-relative numeric
  degrees, documented on `FFIWindConditions`. Regression-pinned: wind FROM north on
  a shot due north is a pure headwind (the 0.19.0 sign convention).
- **Clock-position wind direction entry** (MBA-1367; Kestrel/ATrag/AB-class field
  convention). Every `--wind-direction` flag (trajectory, monte-carlo, come-ups, lead,
  range-table, compare, `profile save`) and both WASM terminal wind-direction args now
  accept marked clock positions alongside plain degrees: `3oc` = 90°, `10h30` = 315°
  (`(H % 12) × 30 + MM × 0.5`, wind-FROM, 12 o'clock = headwind = 0° per the
  post-0.19.0 convention), and — on standalone flags only — the colon form `10:30`.
  Inside `--wind-segment` the ANGLE field takes the colon-free marked forms
  (`10:3oc:400`); `10:30:400` keeps its historical SPEED:ANGLE:DIST meaning. Bare
  numbers remain degrees everywhere, so every existing invocation is unchanged.
  One shared `wind::parse_wind_direction` helper (plus `parse_wind_direction_standalone`
  and `parse_wind_segment_str_detailed`) backs the clap parsers, the segment grammar,
  and the WASM arg loops; malformed clock tokens (hour outside 1-12, minutes outside
  0-59) fail loudly with the rule in the message.

- **Multiple named zeroes and per-load POI offsets in saved profiles** (MBA-1360;
  Lapua Sight-In POI / ATrag zero zones / Strelok multi-zero class). Profiles gain an
  optional `zero_sets` list — each set is a name plus an optional `zero_distance`
  (display units, rescaled by unit conversion like the profile's own zero) and
  optional `poi_up_mil`/`poi_right_mil` per-load DIAL corrections (constant angular
  mils; positive = dial up/right more, so a load impacting high stores a negative
  value). Managed with `profile zero-set add|remove|list`; selected per run with
  `--zero-set NAME` on `trajectory`, `come-ups`, `wind-card`, `range-table`, `dsf`,
  and `plan-truing`. A selected set's zero_distance feeds the auto-zero exactly as
  the profile's zero would (explicit CLI zero flags always win); its dial corrections
  are added to total-correction dial outputs (come-ups, wind-card drift, range-table
  Drop/Wind, the PDF dope card's Drop/Wind) at the shared conversion boundary, BEFORE
  the MBA-1358 tracking-CF division (`dial = (true + bias) / CF`); component holds
  (mover Ring/lead) deliberately stay bias-free. Unknown set names are a hard error
  listing the available sets. Data sources wired in: the `--profile` CSV row's
  previously allowlisted-but-unused `V_OFFSET_MIL`/`H_OFFSET_MIL` columns now form a
  selectable ephemeral set named after the row, and `.a7p` import with `--zero-click`
  additionally records the file's zero click state as a set named `a7p-zero`
  (dial-correction convention; the MBA-1359 engine-field mapping is unchanged).
  Forward-compat follows the documented `bc_segments` pattern: old readers load and
  solve the master zero identically (and reject `--zero-set` loudly as an unknown
  flag); only an old reader's re-save drops the stored sets. Defaults (no `--zero-set`)
  are byte-identical everywhere.
- **Scope tracking correction factors from a tall-target test** (MBA-1358, Litz).
  Per-axis `--elevation-cf`/`--windage-cf` store the published tall-target ratio
  `actual measured travel / dialed travel` (0.95 = under-tracks 5%; strictly between
  0.5 and 1.5; CLI flag over the new validated-on-load `elevation_cf`/`windage_cf`
  saved-profile fields) and DIVIDE every dial-unit output by it exactly once at the
  shared conversion boundary — an under-tracking scope needs more dial. Surfaces:
  come-ups, range-table/compare Drop and Wind columns, wind-card drift, mover lead and
  the Ring column/fields (a dialed quantity; native and WASM, including the WASM `lead`
  command), the PDF dope-card rows, and the zero command's MOA/mrad dial outputs on
  both the native CLI and the WASM terminal's banners — while raw drop/drift inches
  (and degree bore-angle echoes) are never scaled. New `tall-target` subcommand
  computes the CF with no trajectory solve. `true-velocity --elevation-cf` MULTIPLIES
  dialed (mil/moa) observations by the CF before the fit (scope-dial units -> true
  angular), so scope tracking error is not baked into trued MV/BC, and shows dial-unit
  report values back in scope units (÷CF); linear `in` drops never scale. Defaults
  (no flag/fields) are byte-identical. FFI is deliberately excluded (no dial outputs).
- **Equivalent horizontal range (BDC shoot-to) for inclined shots** (MBA-1395). When a
  look angle is set and a zero was solved, `trajectory` prints one summary line —
  `Equivalent horizontal range: X yd (shoot-to for BDC)` — on the native CLI and the
  WASM terminal, and solve-json v1 gains `summary.equivalent_horizontal_range_m`
  (absent when undefined). Definition: the flat-fire range whose ANGULAR elevation
  correction against the same zero matches the inclined solution's (SIG AMR / Leica
  EHR / Gunwerks style; McDonald/Litz angular-match inversion over one extra flat
  re-solve), deliberately not the rifleman's-rule cosine approximation. Omitted for
  flat shots (byte-identical), targets at/inside the zero range, and non-positive
  corrections. New public API `TrajectorySolver::equivalent_horizontal_range`.
- **Drops-reference toggle: LOS vs target plane** (MBA-1403). New
  `--drops-reference {los|target}` on `trajectory` (native CLI and the WASM terminal)
  plus the solve-json v1 request field `shot.drops_reference`. `los` (the default) keeps
  the historical LOS-perpendicular sampled drop byte-identically; `target` reports drop
  as vertical in the target plane — the LOS-perpendicular drop divided by
  `cos(shooting angle)` (JBM's "target plane" reference) — relabels the sampled
  table/CSV drop column (`Drop (target)` / `drop_target_in`), and lets a supplied
  target height slope the sampler's LOS datum. Output-mode toggle only: the solved
  trajectory, drift, and zeroing semantics are unchanged; rejected at shooting angles
  of 90° or beyond; FFI is deliberately not plumbed (raw-kinematics surface). The three
  per-integrator trajectory-sampling blocks in `cli_api.rs` were consolidated into one
  shared helper (byte-identical refactor) so the sampler's LOS datum is chosen in
  exactly one place.
- **Lateral sight offset input for offset-mounted optics** (MBA-1396). New
  `BallisticInputs.sight_offset_lateral_m` (positive = sight RIGHT of bore) models the
  windage half of sight geometry the way `sight_height` models the vertical half: the
  bullet starts that far LEFT of the sight line (`initial_position` lateral term,
  composing with the cant displacement), and a zero solve applies the windage-zero
  convergence (`offset / zero_distance`) so the trajectory crosses the sight line
  laterally at the zero range — not a constant parallel offset. Without a zero solve
  (explicit muzzle angle) only the physical displacement applies. Distinct from and
  additive with MBA-1359's zero POI offset (angular zero state vs physical mount
  geometry). Surfaces: `--sight-offset` (inches imperial / mm metric, the
  `--sight-height` units) on `trajectory`, `mpbr`, `come-ups`, `wind-card`,
  `range-table`, `compare`, and `monte-carlo` (mirroring `--cant`; MC solves no zero, so
  only the physical displacement applies there) plus the WASM terminal's `trajectory`;
  solve-json v1 request field `rifle.sight_offset_lateral_m`; a third appended
  `c_double` on `FFIBallisticInputs`; `compute_wez` gained a trailing
  `sight_offset_lateral_m` parameter (cant precedent); back-compatible
  `sight_offset_lateral_m` saved-profile field (always meters). Omitting every new
  input is byte-identical.
- **Deliberate zero POI offset inputs** (MBA-1359, Kestrel ZH/ZO semantics). New
  `BallisticInputs` fields `zero_poi_vertical_m`/`zero_poi_horizontal_m` record that the
  rifle is deliberately zeroed off — e.g. 0.1 in high / 0.2 in left at the zero range —
  and shift the whole solution by the equivalent angular bias (offset / zero distance),
  applied post-solve to the zero elevation and azimuth (the zero search itself still
  solves perfect convergence). Sign convention: positive = impacts HIGH / impacts RIGHT
  at the zero distance. Surfaces: `--zero-poi-up`/`--zero-poi-right` (inches imperial /
  cm metric) on `trajectory`, `mpbr`, `come-ups`, `wind-card`, `range-table`, and
  `compare` plus the WASM terminal's `trajectory`; solve-json v1 request fields
  `shot.zero_poi_up_m`/`shot.zero_poi_right_m` (SI, additive, no response-shape change);
  two appended `c_double`s on `FFIBallisticInputs` (ABI append-only convention);
  back-compatible `zero_poi_up_m`/`zero_poi_right_m` saved-profile fields (always
  meters); and `.a7p` import can now convert the file's `zero_x`/`zero_y` device click
  counts via the new `profile import --zero-click <e.g. 0.1mil>` flag (the `.a7p` format
  does not store the device's click size, so without the flag the counts remain reported
  as unmapped, unchanged). Omitting every new input is byte-identical to 0.30.x.

### Compatibility
- **Two intentional default-output deltas** (everything else in this release is
  byte-identical without its new flag). (1) The MBA-1367 sentinel fix below changes
  output only for the explicit-`--wind-direction 0`-plus-location-CSV combination.
  (2) `summary.equivalent_horizontal_range_m` and a matching CLI/WASM line (MBA-1395)
  appear on ANY inclined shot (`shooting_angle != 0`) with no opt-in flag, because a
  shoot-to range is meaningful exactly when a look angle is set; flat-fire output is
  unchanged. A consumer that must reject unknown response fields will see this on
  inclined-shot responses.
- **solve-json v1 grows only by addition** (see `docs/SOLVE_JSON_V1.md`): new request
  fields are optional and default to prior behavior; new response fields are omitted when
  their feature is inactive (the one exception is the EHR field above). Removing,
  renaming, re-signing, or re-unitting a field still requires v2.

### Fixed
- **`trajectory --wind-direction 0` no longer loses to a location CSV** (MBA-1367
  sentinel fix). The dispatch treated `wind_direction != 0.0` as "user set", so an
  explicit `--wind-direction 0` (and now `12oc`, which maps to 0°) was silently
  replaced by a `--location` CSV's WIND_DIR column. Explicit flag presence now decides:
  an explicit zero wins over the CSV; omitting the flag still inherits the CSV value
  exactly as before (the only behavioral delta is the explicit-zero-plus-CSV
  combination).

## [0.30.1] - 2026-07-28

### Fixed
- **A supplied `--location`/`--site` (or `--profile`/`--profile-row`) now loads or stops the
  run** (MBA-1425). Three silences existed: `--location` without `--site` was ignored entirely
  (no message at all — the file was never read); an unreadable path warned and continued; a
  missing site row warned and continued. All three produced a clean trajectory computed against
  air or a gun the user never chose. The FILE flags now require their row selectors at the
  parser level (deliberately asymmetric: a standalone `--site`/`--profile-row` keeps its
  pre-existing PDF-labeling job), and both load failures are hard errors naming the path/row.
- **`zero --from-angle` reports real apexes it used to withhold** (MBA-1419 item 3, resolved).
  The max ordinate is now gated on whether the diagnostics trajectory actually TURNED OVER
  inside its envelope, replacing the far-crossing proxy — which was over-conservative (a
  45° G1 launch turns over at ~2040 yd with a genuine ~1432 yd apex the proxy refused to
  print). When the solve is still climbing at its envelope, the table says so, and JSON/CSV
  now OMIT `max_ordinate` instead of emitting the truncation height as if it were an apex.
- **The browser terminal's `zero` now has an atmosphere at all** (MBA-1426 items 3-4). Both of
  its solve branches previously ran against the default atmosphere and REJECTED
  `--temperature`/`--pressure`/`--altitude` as unknown flags, while native zero has carried
  them since 0.29.0. It now accepts the full set including `--pressure-type` and
  `--density-altitude` (with the supersede notice riding the output text, since this surface
  has no stderr), feeding both the target-distance and `--from-angle` solves.
- **Terminal parity sweep** (MBA-1426 items 3-5): `--pressure-type` on the terminal's
  `true-velocity`, `estimate-bc`, and `lead`; `--zero-density-altitude` on the terminal's
  trajectory (naming itself, not `--density-altitude`, in the refusal notice); and the
  from-angle output now names the primary crossing, matching native's `primary_crossing`.

### Added
- **`zero --drag-model`** (MBA-1419). The zero command lacked drag-model selection while every
  sibling — including the browser terminal's own zero — had it, so a G7 BC silently ran
  against the G1 reference table. Default `g1` preserves existing behavior byte-for-byte.

## [0.30.0] - 2026-07-26

The theme of this release is reaching every surface: the numbers the engine
already computes now arrive where consumers actually are — bindings that call
the fast path, scripts that parse the calculator subcommands, saved profiles
that carry their own weather, and above all the browser terminal, which gains
the release's two headline features on the same bytes as native.


### Fixed
- **The browser build now emits `zero_angle_degrees` in JSON** (MBA-1402 parity). 0.29.0 added
  the solved auto-zero angle to native's table, JSON and CSV, but the WASM build shipped it in
  the text banner only — so a consumer parsing the browser JSON could not read back the angle the
  same run had just printed. Reported against a real 0.29.0 WASM build. Present only when
  auto-zero ran, absent (not null) otherwise, matching native. Note the WASM CSV remains
  structurally different from native's (it has no summary-row form at all), which is tracked
  separately rather than bolted on here.
- **Browser-terminal parity gaps from 0.29.0** (MBA-1418). Three places where the hand-written
  WASM layer quietly disagreed with native. The zero-day `--zero-pressure-type` refusal is now
  disclosed in the terminal (it forced `Absolute` silently; native prints a warning, and the
  browser has no stderr, so it rides the table-only notice block its siblings already use). The
  `recoil` CSV header emitted `velocity_m/s` in metric where native emits `velocity_mps` — the
  imperial spelling already agreed, which is why it survived. And `zero --from-angle` is now
  listed in the terminal help and covered by tests, having worked but been undiscoverable.
- **Pinned what `--pressure-type` means for a `--location` CSV pressure** (MBA-1421). No
  behaviour change: a mode typed on the command line already applied to a CSV-supplied pressure,
  which is the defensible reading — unlike a profile-stored mode (which describes the value the
  profile stored alongside it), a typed mode has exactly one pressure in the run to describe.
  That is now covered by tests and stated in CLI_USAGE instead of being an accident of
  implementation. Note the CSV row selector is `--site`, not `--location-name`.
- **Docs and comments that described behaviour the code does not have** (MBA-1420). The QNH
  section of CLI_USAGE said a saved profile's `pressure_reference` "round-trips with the
  profile", which oversold it: the field round-trips through the profile file but is never
  applied to a solve. `--plot`'s help and its code comment both described two panels; it renders
  four (drop, drift, velocity, energy). A `recoil.rs` comment showed arithmetic that does not
  reproduce its own cited figure (752122.5 + 63877.5 is 816000, not 816001, and the shown steps
  give 30.161 ft-lb, not the 30.22 the SAAMI document states) — the pinned test value was always
  correct. The `zero --from-angle` handoff is now documented with its preconditions: a recovered
  range reproduces the original zero only under the same wind, bore height, and zero-day
  atmosphere. And a comment now records which of `atmosphere` and `inputs.pressure` is
  authoritative after the QNH reduction, since the two differ and only one reaches the solve.
- **`--density-altitude` out-of-range errors now name the flag** (MBA-1420). An extreme value
  drove the pressure inverse negative and surfaced as `atmosphere.pressure must be finite`,
  pointing at a field the user never typed; NaN surfaced as a temperature error. The flag now
  carries its own wide bound (-5000 to 120000) and fails naming itself.
- **`zero --from-angle` now says which crossing `zero_range` is** (MBA-1419). A bore angle
  generally crosses the line of sight twice, and the single `zero_range` value is
  far-when-present-else-near — a consumer reading it alone could not tell which of two valid
  answers it received. JSON gains a `primary_crossing` key (`"far"` or `"near"`) and CSV a
  matching row. The diagnostics solve envelope is also floored at the distance the crossing
  search itself covers, so a near-only fallback no longer sizes the follow-up solve off a very
  short near zero. Nothing was silently wrong before this; the rule was documented and both
  crossings already printed.
- **solve-json v1 rejects projectile values that cannot describe a projectile** (MBA-1413).
  `projectile.mass_kg`, `diameter_m`, `length_m` and `ballistic_coefficient` now carry physical
  bounds and return a typed `invalid_value` error with the exact JSONPath. Found by fuzzing,
  which got the engine to accept and solve a bullet weighing 1.1e-66 kg. The visible symptom was
  narrower — at that magnitude the request did not survive a JSON round trip, breaking the
  protocol's stability invariant — but the real defect was accepting the number at all. The
  bounds are enormous (1 mg to 100 kg, 0.1 mm to 1 m), so no real load is affected.
  `rifle.muzzle_velocity_mps` and `shot.max_range_m` are deliberately left unbounded: the MCP
  server documents a split where a structurally valid but unsolvable request returns a tool
  error rather than a protocol error, and bounding those would reclassify it.
- **Input conditioning now reaches every entry point into integration** (MBA-1415). The BC
  reference-standard conversion, the SI-derived imperial fields, and the powder-resolved muzzle
  velocity were applied only in `TrajectorySolver::new`. `fast_integrate` and
  `fast_integrate_with_segments` are public and are called directly by the Python binding, so
  those consumers silently bypassed all three: set a conditioned field, get no error and no
  effect. All three now live in one idempotent `BallisticInputs::normalize_for_solve`, called by
  the constructor and both fast entry points. Latent rather than live — no shipped binding
  exposed the affected field yet — but this is the third instance of the shape (MBA-1296,
  MBA-1356), so it is fixed structurally and pinned by a cross-entry-point parity test.

### Changed
- **printpdf 0.10 → 0.12.4.** PDF dope cards are now self-contained: 0.10 emitted no font
  resources, no `/Encoding` and no `/ToUnicode`, leaving text to whatever the viewer defaulted
  to; 0.12 embeds subset fonts with a ToUnicode map, so a card renders identically everywhere
  and stays searchable and copyable. **The cost is file size** — a sample dope card goes from
  ~6 KB to ~836 KB, which is the embedded font data. `default-features = false` is preserved, so
  the Linux-only `mmapio`/`rust-fontconfig` chain that breaks the BSD targets stays out.

### Added
- **Density altitude reaches the zero surfaces** (MBA-1422): `zero --density-altitude`, and
  `trajectory --zero-density-altitude` for declaring a zero-day density altitude alongside the
  existing `--zero-pressure`/`--zero-pressure-type`. MBA-1366 scoped the flag to `trajectory`
  only, so a shooter carrying conditions as a density altitude could not solve a zero the same
  way they solve a trajectory. Precedence matches the shot-day rule: DA supersedes altitude and
  pressure, an explicit temperature still wins, and the pressure mode is refused with a notice
  because DA leaves it no value of its own to describe.
- **`trajectory --saved-profile` now loads a profile's atmosphere** (MBA-1417). `profile save` has
  always stored temperature, pressure, humidity and altitude — and since 0.29.0 `pressure_reference`
  and `density_altitude` — but the trajectory path loaded none of them, so a profile carried its
  ballistics and silently dropped its conditions. They now load as a set, with per-value precedence
  CLI flag > `--location` CSV > profile > standard. **If you have profiles that stored a
  non-standard atmosphere, runs against them will change** — that is the fix, but re-baseline
  before relying on it. Profiles that never set an atmosphere are byte-identical, pinned by a test.
  The stored pressure mode applies only when the profile's own pressure value is in use, and a
  stored density altitude only when the run supplies no pressure or altitude of its own.
- **`--pressure-type` now works on every subcommand that accepts `--pressure`** (MBA-1416):
  `estimate-bc`, `true-velocity`, `plan-truing`, `mpbr`, `come-ups`, `lead`, `wind-card`,
  `stability`, `range-table` and `compare`. These ten previously accepted a pressure and silently
  treated it as absolute station pressure, so a weather-report barometer value entered at
  elevation was wrong by the ISA reduction on every command except `trajectory` and `zero`. The
  mode reduces against each command's own `--altitude`. Omitting it is exactly `absolute`, so
  existing invocations are byte-identical — pinned by a test per command, alongside one asserting
  the flag actually reaches the physics rather than merely parsing.
- **`ballistics drag-curve`** (MBA-1424) prints a built-in reference drag function as
  `(Mach, Cd)` data in table, CSV, or JSON form, for all nine models (G1/G2/G5/G6/G7/G8/GI/GS/
  RA4). Points come out verbatim from the table — no resampling — so a consumer charts exactly
  what the solver interpolates rather than a re-vendored copy that drifts. Note the tables do not
  share a Mach domain: most run to Mach 5, while GS and RA4 stop at Mach 4. Library callers can
  use `drag::reference_drag_table`. These are the public-domain Aberdeen/BRL functions as
  tabulated in McCoy, plus the British RA 1929 function. Available in the browser terminal too
  (MBA-1426), where native and WASM emit **the same shared formatter's bytes** — the sharing is
  deliberate, so the two surfaces cannot drift the way the recoil CSV header once did.
- **`trajectory --with-drag-coefficient`** (MBA-1423) adds a `drag_coefficient` key to each
  point of `--full` JSON. This is the projectile's *own* Cd — the reference-table value scaled
  by the form factor `SD / BC` — not the G-model curve, which would be identical for every
  bullet sharing a drag model. A segmented BC therefore shows its band steps, and a custom
  drag table passes through unscaled. Covers only the speeds actually flown, and pairs Cd with
  the station speed of sound the document already reports Mach against. Off by default: JSON
  without the flag is byte-identical to 0.29.0. Mirrored in the browser terminal's `-o json`
  (MBA-1427) — the build the requester of this feature actually runs — reading the same
  solver-annotated value rather than re-deriving it.

## [0.29.0] - 2026-07-26

The theme of this release is inputs that mean what the shooter means: several
values you can already type were being silently reinterpreted, and several
flags did nothing without saying so.

### Changed
- **`trajectory --auto-zero` output now includes the solved zero angle**
  (MBA-1402) on the table, JSON, CSV and browser-terminal surfaces. This is a
  deliberate, user-visible change to a pre-existing flag: if you parse
  `--auto-zero` output, expect the new line/key/column.
- `--adjustment-unit` now warns when it cannot affect anything (MBA-1414). On
  `trajectory` it reaches only the mover Ring column (`--target-speed`) and the
  PDF dope card (`-o pdf`); supplying it with neither used to be silent.

### Added
- Declare which reference atmosphere a BC is corrected to (MBA-1365):
  `--bc-reference icao|army-standard-metro`. Many Sierra, Hornady and Barnes
  BCs are Army Standard Metro referenced, which reads about 1.8% too little
  drag if fed in unconverted. ICAO remains the default and is unchanged.
- Enter a sea-level (QNH / weather-report "barometer") pressure and have it
  reduced to station pressure (MBA-1397): `--pressure-type absolute|qnh`, plus
  `--zero-pressure-type` for the zero day. Typing a QNH value at elevation and
  having it treated as station pressure was a classic source of vertical error.
- Carry conditions as a single density altitude (MBA-1366):
  `--density-altitude`, the entry mode Shooter, Nosler, AB Analytics, Ballistic
  AE and TRASOL all support. Implemented by back-solving an ISA-equivalent
  atmosphere, so Mach, lapse rate and segmented-atmosphere behaviour are all
  preserved. It supersedes `--altitude`/`--pressure`; an explicit
  `--temperature` still wins, with the pressure re-solved so the same density
  altitude is reproduced. The `DA`/`DENSITY_ALTITUDE` CSV columns now work.
- Correct a downrange chronograph reading back to true muzzle velocity
  (MBA-1377): `--chrono-distance`. Most chronographs read 10-15 ft out.
- Recover the zero range a stored bore angle implies (MBA-1402):
  `zero --from-angle`. A bore angle generally implies TWO zeros — the same
  relationship that makes a 25 yd and a ~300 yd zero the same angle — so both
  the near and far crossing are reported.
- `recoil` and `power-factor` calculators (MBA-1372), using SAAMI's own recoil
  formulae and published USPSA/IDPA/SASS thresholds.
- Velocity and energy panels on `trajectory --plot` (MBA-1394), alongside drop
  and drift.
- Adjustment units on `wind-card`, `lead`, `range-table` and `compare`, plus
  independent per-axis selection via `--windage-unit` (MBA-1410).

### Fixed
- `--windage-unit moa` with `--adjustment-unit clicks` printed whole clicks
  under an MOA header on `range-table` and `compare` (MBA-1410).
- `come-ups --windage-click-value` was inert and silent; it now warns
  (MBA-1410). The PDF dope card no longer prints clicks with a trailing `.0`.

## [0.28.1] - 2026-07-25

### Fixed
- Inclined zero solves work again (MBA-1412): since 0.25.0, any nonzero
  `shooting_angle` made `calculate_zero_angle_with_conditions` (and every
  surface behind it — bindings, browser terminal, FFI, CLI) fail with "Zero
  angle did not converge". Zeroing now solves a level rifle (`SightLine`
  contract, matching the cant doctrine); solve-json v1's documented
  world-vertical zero contract is unchanged (`WorldVertical`).
- Browser-terminal `lead` honors `--cd-scale` on loaded custom drag decks
  (previously applied the deck but silently ignored the scale) (MBA-1411).
- Browser-terminal BC5D coercion warnings survive the solve error path
  (previously lost on failures).
- `generate-bc-segments` warns (once) when a non-G1/G7 family is coerced to
  G1, matching the other auxiliary lookups.
- BC-fit reporting labels the actual drag family (G2/G5/GI/GS/RA4) instead of
  a `G?` wildcard, and `estimate-bc --drag-model` accepts (and documents) any
  single family.
- The `--cd-scale` pairing error only suggests `--bc-adjustment` on commands
  that have that flag.

### Added
- Truing validity windows (MBA-1405): the truing report prints the
  MV-calibration window (90-100% of the Mach 1.2 crossing distance), warns
  per observation outside it, prints cause-aware notes when no window exists
  (supersonic throughout vs never supersonic), and cross-references
  `plan-truing`; the `dsf` verb prints the DSF window (at or beyond 90% of
  the Mach 0.9 crossing). Browser terminal in byte-parity lockstep; truing
  JSON gains additive `mv_window_start_m`/`mv_window_end_m` (null when
  absent). `TrajectoryResult` gains additive labeled crossing fields
  (`mach_1_2_distance_m`/`mach_1_0_distance_m`/`mach_0_9_distance_m`).
- Browser-terminal DSF parity (MBA-1411): repeatable `--dsf-point MACH:DSF`
  on `trajectory` passes a per-call table (no profile storage in the
  terminal), validated by the engine, applied drop-only with the native note.

### Performance
- Velocity-band BC transition boundaries are computed once per segment table
  instead of per integration step (no numeric change).

## [0.28.0] - 2026-07-24

### Breaking
- `wez::compute_wez` gained a positional `cd_scale` parameter (MBA-1356).
- `TrajectoryParams` gained a required `cd_scale` field — external constructors
  must set it (`1.0` = previous behavior) (MBA-1356).
- `BallisticInputs` gained a public `cd_scale` field, breaking exhaustive
  struct literals (soften with `..Default::default()`) (MBA-1356).

### Added
- Turret adjustment units (MBA-1355): `--adjustment-unit smoa|iphy|clicks` on
  trajectory/come-ups (native + browser terminal), suffixed click graduations
  (`0.25moa`, `0.1mil` — bare numbers are an error), whole-click display
  rounding (ties away from zero), unit-suffixed column headers, new
  `adjustment` library module. Follow-up for remaining commands and per-axis
  units: MBA-1410.
- `.drg` drag-file support (MBA-1409): `--drag-table` on trajectory,
  monte-carlo, and zero accepts QuickTARGET-style `.drg` decks; `(cd, mach)`
  vs `(mach, cd)` column order is auto-detected, with a hard "ambiguous
  columns" error (and CSV workaround hint) when detection cannot decide.
- Real G2, G5, GI, and GS reference drag tables plus the new RA4 model
  (MBA-1386): the native CLI now accepts the full nine-family `--drag-model`
  set; FFI adds additive `bc_type` slot 8 for RA4; six new stderr warnings
  cover the remaining G1/G7-only auxiliary lookups. This retires 0.27.1's
  G1-fallback notes ("real tables remain future work" no longer applies).
- `--cd-scale` whole-curve drag multiplier for custom decks (MBA-1356) on
  trajectory, monte-carlo, and zero (native + browser terminal), with a
  [0.5, 2.0] plausibility warning and a pairing error without `--drag-table`;
  FFI adds additive `_scaled` export variants.
- Mach-keyed DSF truing table (MBA-1357): new `dsf` verb derives drop-scale
  factors from observed drops (staging gates steer supersonic observations to
  muzzle-velocity truing first), stores up to six points per saved profile
  (supersede within 0.05 Mach, `profile save --clear-dsf` to reset), and
  trajectory/come-ups auto-apply the table as a strictly drop-only
  post-processing correction (velocity/energy/time outputs unchanged;
  table-format-only "DSF table active" note). New `truing_dsf` library
  module; profile field is additive (older profiles load unchanged). Browser
  terminal parity is tracked as MBA-1411.

### Changed
- Numeric shift: `-d g2|g5|gi|gs` now solve against the real reference tables
  (previously a G1 approximation) — affects library, browser terminal, and
  FFI consumers using those families (MBA-1386).
- Numeric shift: profile/CSV `drag_model` strings `G6`/`G8`, previously
  warn-coerced to G1 on the native CLI, now run real G6/G8 physics
  (MBA-1386).
- Numeric shift: velocity-band BC transitions are smoothstep-blended at band
  boundaries (margin = min(50 fps, 25% of the narrower adjacent band)), so
  banded-BC trajectories are continuous across band edges; mid-band values
  are unchanged (MBA-1404).
- Come-ups lead column emits a `lead_smoa` JSON key and the Ring header
  carries its unit (`Ring(mil)`) under the default and `--target-speed`
  paths (MBA-1355).

## [0.27.1] - 2026-07-19

### Added
- Truing experiment design (MBA-1346): new native `plan-truing` command and
  `truing_plan` library API choose an exact-size, minimum-separation-compliant set of
  target ranges from an explicit candidate list/grid. The deterministic
  exhaustive/greedy-exchange optimizer reuses batched production-solver Jacobians,
  reports sensitivity, conditioning, singular values, weak-axis uncertainty,
  station information gain, and rejected/unreachable candidates, and returns an
  honest MV-only recommendation when the available facility cannot identify BC.
- Uncertainty-aware MV/BC truing (MBA-1353): opt-in known-1σ observation weighting,
  per-reading sigma overrides, explicit independent MV/BC priors, a constrained
  weighted joint MAP, local-Gaussian covariance/95% parameter intervals/correlation,
  effective degrees of freedom, and propagated model/future-observation drop bands.
  Approximation failures and weak/prior-dominated identification are structured and
  visible; the legacy point-estimate path/schema remains unchanged when uncertainty
  is not requested.

### Fixed
- Restored the MBA-1346/MBA-1353 truing analysis features (plan-truing command,
  uncertainty-aware true-velocity, truing_plan/truing_uncertainty library
  modules) whose implementation had not been committed; recovered and
  validated against their surviving test suites (884-test baseline).
- FFI: `ballistics_calculate_trajectory_with_drag_table` (and the zero-angle
  variant) now reject custom drag tables longer than
  `MAX_FFI_DRAG_TABLE_LEN` (4096 rows) instead of attempting unbounded
  allocations from a caller-supplied length (MBA-1407).
- Truing: the remaining drifted grain-to-kg literals (including the recovered predict_many_in_unit site) now use the canonical
  `constants::GRAINS_TO_KG` (~1.5e-7 relative correction, invisible at output
  precision; completes MBA-1327/MBA-1333) (MBA-1408).
- WASM terminal: `--auto-zero` no longer inherits the shot's Coriolis
  conditions (latitude/azimuth), matching the native CLI's zero solve
  (fix-half of MBA-1384; zero-day opt-in flags remain future work).
- Drag models: requesting a family without a dedicated table now warns
  instead of silently using G1 — native CLI profile/CSV strings other than
  G1/G7 warn on stderr; the browser terminal prints a note for G2/G5/GI/GS
  (fix-half of MBA-1386; real tables remain future work).

## [0.27.0] - 2026-07-18

### Added
- `powder` command (MBA-737): resolve the powder-temperature-adjusted muzzle velocity
  without running a trajectory — linear fps-per-degree model or a measured
  `--powder-temp-curve` (clamped interpolation, overrides linear), `--sweep
  START:END:STEP` velocity ladders, optional muzzle energy via `-m/--mass`, and
  table/json/csv output. Same command in the WASM terminal. The physics is one shared
  implementation (`resolve_powder_adjusted_velocity`, now public with
  `interpolate_powder_temp_curve`): the solver's inline powder resolution was extracted
  to it, so the command prints exactly what the trajectory/lead solvers fly. NOTE:
  `powder` always applies the linear model; `trajectory`/`lead` still gate theirs
  behind `--use-powder-sensitivity` (a curve applies there unconditionally).
- WASM terminal: `trajectory --plot` (MBA-1337) — the native terminal charts (braille
  or ascii, drop + lateral drift vs range) render in the browser terminal;
  `terminal_plot` moved into the library unchanged.
- WASM terminal: full `true-velocity` command (MBA-1343) — single-observation truing
  AND the multi-observation `--observed RANGE:DROP` joint MV+BC calibration, with
  table/json/csv output character-identical to the native CLI (multi-obs JSON floats
  may differ in the last ~3 of 17 significant digits across ISAs). The truing core
  moved to a new `ballistics_engine::truing` module.
- WASM terminal: `monte-carlo --wez` (MBA-1343) — the full Weapon Employment Zone
  sweep (rect/circular targets, wind-call error in quadrature, seeded per-step MC,
  variance attribution) with summary/statistics/full output byte-identical to native.
  The WEZ compute moved to a new `ballistics_engine::wez` module; the per-step seed
  (reproducibility contract) lives in the library.
- Truing quality/robustness (MBA-1337): quality bands now judge a mil-equivalent RMS
  so the verdict no longer changes with `--drop-unit`; identifiability diagnostics
  differentiate in mil space; the MV-only convergence flag is reported instead of a
  hardcoded `true`; exactly-determined fits (observations == fitted parameters) say so
  instead of reading "excellent".
- MCP server hardening (MBA-1337): fractional JSON-RPC ids rejected (integral floats
  like `2.0` still accepted); `engine_info` enforces its advertised
  `additionalProperties: false`; requests before `initialize` are rejected with
  `-32002` (ping exempt).
- `solve-json` v1: named cross-interface conformance fixtures (MBA-1309).
- Docs: units-by-output-surface reference table in CLI_USAGE (native world-frame JSON
  vs table/CSV deflection vs solve-json vs the WASM terminal's `drop_inches`/`drop_cm`
  keys), plus the default-bore-height (60 in vs 1500 mm) and mover-ring unit notes.
- `compare` command (MBA-735): side-by-side multi-load comparison. Repeatable
  `--load "NAME:DRAG:BC:MASS:VELOCITY[:DIAMETER]"` specs (2-8 loads, mixable with
  `--profile`), each load zeroed independently at the shared `--zero-distance`, then run
  through identical wind/atmosphere. Table output shows per-load drop/drift/velocity per
  range; JSON adds energy, time of flight, zero angles, and per-row deltas vs load #1;
  CSV emits one sanitized column group per load. Saved-profile loads carry their
  velocity-BC segments and custom Cd(Mach) drag curves into both the zeroing and the
  trajectory runs (tagged in the table legend and flagged in JSON) — the range-table
  scalar-BC caveat does not apply to `compare`.

### Changed
- WASM terminal: the new `true-velocity` and `monte-carlo` handlers reject bare
  (non-flag) tokens and value-less flags with clear errors, matching native clap —
  previously such input was silently ignored, which in multi-observation truing could
  fit fewer observations than the user supplied (MBA-1343 review).
- Bin-only dependencies (clap, clap_complete, csv, dirs, strsim) moved behind a
  default-on `cli` feature (MBA-1331); the `ballistics` binary declares
  `required-features = ["cli"]`. Plain builds and `cargo install` are unchanged;
  lib-only consumers can drop the five with `default-features = false`. NOTE any
  `--no-default-features` build that needs the BINARY must add `--features cli` or the
  bin target is silently skipped.
- WEZ dominant-source labels are now the snake_case spellings (`wind_call`, `mv_sd`,
  `other`) on every surface (MBA-1337); the summary table and statistics CSV
  previously used display spellings ("wind call", "MV SD", "other/group"). The
  `-o full` JSON contract is unchanged.

### Fixed
- WASM terminal: `--auto-zero` with a nonzero `--shooting-angle` no longer fails
  ("Zero angle did not converge … not bracketed") or bakes the incline into the zero
  (MBA-1344) — the zero solve now runs on a level range like the native CLI, so the
  rifle's zero is independent of the shot's slope and cant.
- The drag-table CSV fallback loader no longer silently loses the first data row of a
  headerless file (a side effect of the old csv-crate headers mode), and tolerates
  quoted fields, stray non-UTF-8 bytes, and bare-CR line endings (MBA-1331).

## [0.26.0] - 2026-07-17

### Changed
- **BREAKING (MBA-1339):** the native CLI `--bore-height` flag now takes **inches** in
  imperial and **millimetres** in metric (previously feet / metres), unifying it with
  `--sight-height` and the WASM CLI's `--muzzle-height`. Defaults are unchanged
  (60 in / 1500 mm = 5 ft / 1.5 m), so any run that does not set `--bore-height` is
  byte-identical — but a run passing e.g. `--bore-height 5` now means 5 inches, not 5 feet.
  Both surfaces now accept both flag names: native `--bore-height` gains a `--muzzle-height`
  alias and the WASM `--muzzle-height` gains a `--bore-height` alias, with identical units
  and defaults on each.
- Removed the unused `rayon` dependency (MBA-1332). No code referenced it; downstream
  builds shed the transitive dependency tree.
- Internal: all 57 duplicated kg-scale grain-conversion literals now reference
  `constants::GRAINS_TO_KG`, enforced by the constants-guard test (MBA-1333). The value
  was exact and internally consistent everywhere, so trajectory output is unchanged.

### Added
- Fallible segmented-wind construction and raw integration APIs (MBA-1338):
  `wind::WindSegmentError` — a typed validation error carrying the caller's segment
  index, the offending field, and the violated rule — plus `WindSock::try_new`,
  `TryFrom<Vec<WindSegment>> for WindSock`, `try_integrate_trajectory`, and
  `try_solve_trajectory_rust`. Checked APIs reject non-finite or negative `speed_kmh`,
  non-finite `angle_deg`, non-finite or non-positive `until_m`, and non-finite
  `vertical_mps` **before** sorting, precomputation, or producing any trajectory points.
  `wind::validate_wind_segments` is now public. The existing infallible signatures are
  unchanged and documented as legacy/unchecked entry points; their error text is
  byte-identical to previous releases.
- WASM tests for `--print-bc-segments` table-gating and `lead -o` invalid-value rejection.

### Fixed
- WASM (MBA-1294): the auto-zero "Rifle zeroed at …" banner is now table-only; `-o json`
  and `-o csv` output stay pure machine payloads. The native "non-finite state" integration
  error now names the likely cause (extreme bore/muzzle height → air density) and points to
  `--altitude`.

## [0.25.2] - 2026-07-16

### Fixed
- Build: `printpdf` is now depended on with `default-features = false`, dropping its
  default `html` feature. That feature was the sole path to `azul-layout` →
  `rust-fontconfig` → `mmapio`, whose `MAP_LOCKED` usage is Linux-only and failed to
  compile on every BSD target (`cannot find value MAP_LOCKED in crate libc`). This had
  broken FreeBSD/OpenBSD/NetBSD binary builds since the printpdf 0.7 → 0.10 bump. The
  `pdf` feature still generates dope-card PDFs via core printpdf + `allsorts-azul` font
  embedding — unchanged output. As a side effect the `azul-*`, `rust-fontconfig`, and
  `mmapio` transitive subtree is removed from **all** targets, slimming build times and
  binary size. No physics, API, or numerical-output change; PDF export is byte-compatible.

## [0.25.1] - 2026-07-16

### Added
- WASM `clearDragTable()`: unload a custom drag table previously installed via
  `loadDragTable()`, reverting every `trajectory`/`zero`/`lead`/`monte-carlo` run to the
  standard G1/G7 + BC drag. The exact inverse of `loadDragTable` — it lets a single engine
  instance compare a G7-BC solve against a measured Cd curve (CDM) without a second instance
  (`load → solve CDM → clear → solve G7`), the same load/clear pattern `clearWindSegments`
  provides for segmented wind. Returns whether a table was loaded; idempotent.

### Fixed
- Build: gate the `criterion` benchmark dependency to non-wasm targets. `criterion` 0.8
  pulls in Rayon, which hard-`compile_error!`s on `wasm32`, which had broken the entire
  `wasm-pack test` harness (every WASM unit and integration test was un-runnable).
  Benchmarks only ever run natively.

## [0.25.0] - 2026-07-15

### Added
- `.a7p` profile import (`profile import file.a7p` — ArcherBC2 format): zero-dependency
  MD5 + proto3 wire reader, full honest mapping report, multi-BC velocity bands into
  profile BC segments, CUSTOM Mach:Cd drag curves end-to-end (profile schema v2, backward
  compatible).
- Mover ring: `trajectory --target-speed` per-point ring column (speed x ToF) in
  table/CSV/JSON (`mover_ring_m`/`mover_ring_mil`), `--adjustment-unit moa` honored;
  target speed bounded [0,300] mph on every surface.
- `lead`: powder-temperature flags (native + WASM); WASM `lead` gains the environmental
  flags (`--temperature`/`--pressure`/`--humidity`/`--altitude`/wind) it never had.
- WASM `loadDragTable(bytes)`: custom Mach:Cd drag curves in the browser, same CSV parser
  as the native CLI, auto-applied to trajectory/zero/lead/monte-carlo.
- Truing v2: `true-velocity --observed RANGE:DROP` (repeatable) — joint MV+BC fit against
  the real forward model with identifiability gating and per-observation residuals.
- WEZ mode: `monte-carlo --wez` — P(hit)-vs-range sweep for WxH targets with
  `--wind-call-error` and per-range variance attribution.
- MCP server: `ballistics mcp` — zero-dependency Model Context Protocol server over stdio
  (tools: solve on the solve-json v1 schema, engine_info).
- solve-json v1: versioned single-solve JSON protocol (library service + `solve-json`
  stdin transport + compatibility lock).
- `trajectory --plot [braille|ascii]`: inline terminal charts of drop and drift.
- Trajectory JSON self-description: top-level `legend` block documenting units and axis
  conventions (x = lateral, z = downrange, yards) — existing keys unchanged.
- Golden physics validation harness (`cargo test --features validation`): versioned
  reference cases (analytic + pinned cross-implementation) with per-source tolerances.
- npm packaging: `scripts/build-npm.sh` produces the publish-ready `ballistics-engine`
  package (resumed from 0.13.4).

### Changed
- Integration divergence guard: `solve()` now returns `Err` for stiff inputs whose
  accepted steps exceed the physical speed budget (previously could return Ok with
  non-physical results, e.g. negative terminal range); result postcondition also rejects
  negative range/time/velocity/energy (MBA-1293).
- Grain<->gram conversions unified on exact SI constants (~1e-7 relative shift in a few
  display paths); a guard test bans drifted literals.
- PDF `--target-speed` now respects `--units metric` (was always mph).

### Fixed
- wasm_tests.rs brace nesting (tests were lexically nested inside one function).
- WEZ baseline solve errors propagate instead of panicking.
- no_std feasibility spike report committed under docs/ (qualified GO).

## [0.22.14] - 2026-07-08

### Added
- `--print-bc-segments` (trajectory command): prints the BC5D-generated
  velocity:BC segment ladder as ready-to-paste `--bc-segment` arguments, in the
  active `--units`. Lets devices that cannot hold the 5D tables (e.g. the WASM
  CLI) run BC5D-equivalent corrections from a transcribed ladder.

### Fixed
- BC5D info line labeled the bullet length "(est)" even when supplied via
  `--bullet-length`/CSV; it now reports the source ("user"/"est") and notes
  that length is not a v2 table lookup axis (axes: drag type, weight, BC,
  muzzle velocity, current velocity).

## [0.22.13] - 2026-07-08

### Fixed

- **CIPM air density ignored humidity.** `calculate_air_density_cipm` formed the water-vapor
  mole fraction by dividing the vapor pressure (hPa) by the total pressure in Pa, making it
  100x too small and returning essentially dry-air density regardless of relative humidity
  (e.g. 15 C / 1013.25 hPa / 50% RH gave ~1.2254 kg/m³ instead of the CIPM-2007 moist value
  ~1.2211 — moist air is lighter than dry air). The vapor pressure is now converted to Pa
  before the division. Added a regression test pinning the result to the reference Python
  implementation (0.1% tolerance) at three temperature/pressure/humidity combinations.

- **Custom drag tables used the wrong retardation denominator.** All three solver paths
  (`derivatives.rs`, `fast_trajectory.rs`, `cli_api.rs`) divided a custom drag table's Cd by
  `bc_value`. A custom curve supplies the projectile's ACTUAL drag coefficient, so the
  physically correct denominator is the SECTIONAL DENSITY in lb/in²
  (`weight_grains / 7000 / diameter_in²`): `Cd_own / SD == Cd_ref / BC`. Trajectories with a
  custom table therefore wrongly scaled with whatever `bc_value` happened to be set. The
  denominator is now the sectional density derived from the bullet's mass/diameter (imperial
  mirror fields, with SI fallback); if both are unavailable it falls back to `bc_value` with a
  one-time warning instead of panicking. New `BallisticInputs::sectional_density_lb_in2` /
  `custom_drag_denominator` helpers; regression test asserts custom-table trajectories are
  invariant to `bc_value` and still differ from the G-model run.

## [0.22.12] - 2026-07-07

### Fixed

- **WASM: `--units metric` (and `-u`) now works.** `run_command` pre-scanned `--units` to
  select the unit system but left the flag in the argument list, and the `trajectory` / `zero`
  handlers' unknown-flag rejection then threw `Unknown flag: --units`. Metric commands were
  therefore unusable from WASM (including metric `--bc-segment`), despite `--help` advertising
  the flag. The handlers now skip `--units`/`-u` (its value is already consumed globally).
  Native CLI was unaffected. WASM-only change.

## [0.22.11] - 2026-07-07

### Added

- **`--bc-segment VMIN:VMAX:BC` — manual velocity-keyed BC segments** (CLI + WASM, repeatable).
  Supply your own velocity-dependent BC ladder directly instead of relying on the auto-estimator
  or a BC5D table: the given BC applies while the bullet's current speed is in `[VMIN, VMAX)`.
  Requested by an external user doing downrange BC work.
  - VMIN/VMAX follow `--units` (fps imperial, m/s metric); converted to fps internally.
  - Keyed to **velocity**, orthogonal to distance-keyed `--wind-segment` — combine both and each
    applies on its own axis. A velocity outside every segment falls back to the base `--bc`.
  - Passing any `--bc-segment` implies `--use-bc-segments` and **overrides** `--bc-table-dir`.
  - WASM: `--bc-segment` accepted by `runCommand`; validation errors thrown as JS exceptions.

### Fixed

- **`--auto-zero` now solves the launch angle with the active velocity-keyed BC** (manual
  `--bc-segment`, `--use-bc-segments`, or a `--bc-table-dir` table) instead of the base `--bc`.
  Previously the native CLI zero-angle solver ignored BC segments, so a segment that changed
  early-flight drag mis-zeroed the shot (it grounded short of the requested zero distance and
  reported wrong drop/range). Now matches the WASM build. (Pre-existing bug, surfaced by
  `--bc-segment`.)
- **`--bc-segment` + `--bc-table-dir` together**: the out-of-segment fallback base BC is no
  longer scaled by the BC5D table's muzzle correction — manual segments are a clean override.

## [0.22.10] - 2026-07-06

### Added

- **BC5D correction tables now usable from WASM (Node/browser).** The 5-dimensional BC
  correction tables (hosted at https://ballistics.tools/downloads/bc5d/) could previously
  only be applied by the native CLI via `--bc-table-dir`, because loading them required
  filesystem access. WASM has no filesystem, so a JS/Node host now fetches the `.bin` and
  hands the raw bytes to the engine:
  - New WASM binding `loadBc5dTable(bytes: Uint8Array)` parses a `bc5d_<caliber>.bin` in
    memory and returns a short summary. `hasBc5dTable()` reports whether one is loaded.
  - Once a table is loaded, any `trajectory` run with `--use-bc-segments` synthesizes
    velocity-dependent BC segments from the table and applies them — the same offline
    parity with the online solver (ClusterBCDegradation + BC segments + weather) that the
    native `--bc-table-dir` path produces. Load a table matching the bullet's caliber.
  - Public library API: `Bc5dTable::from_bytes(&[u8])` (a filesystem-free counterpart to
    `Bc5dTable::load`) and `Bc5dTable::generate_segments(base_bc, drag_type, weight_grains,
    muzzle_velocity_fps)` (the table→`BCSegmentData` synthesis, previously CLI-only).

## [0.22.9] - 2026-07-06

### Fixed / Added

- **`estimate-bc` now fits real dope-card data correctly.** Two defects made it return
  wrong BCs from dope cards (reported by an external user):
  - It hard-coded a fixed ICAO-standard atmosphere with no way to change it, so a card
    measured at a different air density fit to the wrong BC. Added `--temperature`,
    `--pressure`, `--humidity`, `--altitude` (same units as the trajectory command) so the
    fit runs at the conditions the data actually came from.
  - It fit **bore-referenced** (flat-fire) drop, but dope cards are **zeroed** (drop below
    line of sight). A zeroed card has a ~0-drop point at the zero range, which is impossible
    for a flat-fire model and drove the fit to nonsense (a sparse 3-point card pinned both
    G1 and G7 at the sweep ceiling). Added `--zero-range` (+ `--sight-height`): given it,
    drop is fit as drop-below-LOS from a rifle zeroed there — exactly what a dope card
    prints. Omitted, drop stays bore-referenced (unchanged, backward compatible).
  - `estimate-bc` now **warns** when drop data looks zeroed (a point near 0 drop) but
    `--zero-range` wasn't supplied, and **flags a fit as UNRELIABLE** when it runs to the BC
    search limit instead of silently returning a boundary value. G7 fits are capped at a
    physical 0.70 (they were sharing G1's 1.2 ceiling).

  Native CLI (table/json/csv; JSON/CSV gain a `reliable` field) and WASM. Public API:
  `estimate_bc_fit` gains `atmosphere` / `zero_range` / `sight_height` params, and
  `BcEstimate` gains an `at_bound` flag.

## [0.22.8] - 2026-07-06

### Added

- **`estimate-bc` now estimates G1 *and* G7 BCs, and can fit velocity data as well as drop
  data.** Previously the command fit a single G1 BC from a drop curve only. New options:
  - `--drag-model g1|g7|both` (default `both`) — estimate the BC referenced to either
    standard drag model, or both at once. (A bullet's G7 BC is a different number from its G1
    BC — for a boat-tail, roughly half; both are now reported side by side.)
  - `--velocity-data "dist,vel;..."` — fit against a downrange velocity-retention curve
    instead of, or in addition to, drop. A velocity fit doesn't depend on the zero, sight
    height, or launch angle, so it is often the more reliable basis.
  - `--data "dist,drop;..."` on the native CLI too (n-point drop series), alongside the
    existing `--distance1/--drop1/--distance2/--drop2` input (now optional).

  The command prints one row per (drag model × supplied data basis), each with a fit-quality
  RMS. Supplying both a drop series and a velocity series with `--drag-model both` yields all
  four variants (G1/drop, G1/velocity, G7/drop, G7/velocity). Native CLI (table/json/csv) and
  WASM. New public API `estimate_bc_fit` / `BcFitMode` / `BcEstimate`; the existing
  `estimate_bc_from_trajectory` remains as a G1/drop convenience wrapper. Requested by an
  external WASM user.

## [0.22.7] - 2026-07-05

### Changed

- **Dependency updates.** Bumped core dependencies to their current major
  versions. No engine API or physics changes; the native binary, the WASM
  build, and the full test suite are all unchanged and green.
  - `nalgebra` 0.34 → 0.35
  - `rand` 0.9 → 0.10
  - `rand_distr` 0.5 → 0.6
  - `getrandom` 0.3 → 0.4 (the wasm32 `wasm_js` backend)

- **Monte Carlo RNG stream changed (seeded-output note).** `rand` 0.10 revised
  its internal generator, so a given seed now produces a different sequence of
  draws than 0.22.6 and earlier. Monte Carlo *statistical* results — means,
  spreads, confidence intervals — are unchanged, but if you depend on a fixed
  seed reproducing an exact set of simulated shots, those specific values will
  differ from prior versions. This is a one-time shift from the dependency
  update, not a change in the physics.

  (`printpdf` 0.9 was evaluated and deliberately **not** taken: it is a full
  API rewrite that would require rewriting the dope-card PDF layer for no
  user-facing gain, so `printpdf` stays pinned at 0.7.)

## [0.22.6] - 2026-07-05

### Fixed

- **Zero-day atmosphere appeared to have no effect at a short auto-zero.** With
  `--auto-zero 100`, a base command and one adding extreme zero-day atmosphere
  (`--zero-temperature 20 --zero-pressure 20 --zero-humidity 90 --zero-altitude 12000`)
  produced byte-identical output. The atmosphere *was* reaching the zero solve, but the
  zero-angle search converged as soon as the height error at the zero distance was under
  1 mm (~0.037 MOA at 100 yd) — coarser than the small zero-day-density effect there — so
  two very different atmospheres rounded to the same zero angle. Tightened the convergence
  to 0.1 mm (angle tolerance 1e-7 rad), so the effect now resolves even at short zeros. It
  remains physically small (a 100 yd zero is nearly density-independent: ~0.04" at 100 yd,
  ~0.2" at 490 yd for those extreme conditions) but is no longer quantized to zero. This
  also makes the zero angle itself more precise. Reported by an external user.

## [0.22.5] - 2026-07-05

### Added

- **Decoupled powder temperature for the curve, symmetric on shot and zero days.** The
  `--powder-temp-curve` maps *powder* temperature to velocity, but was interpolated at the
  ambient *air* temperature. Powder isn't always at air temperature (warmed in a chamber,
  cooled in a pocket), so velocity and air density now use separate temperatures:
  - `--powder-temp` (now accepts no value / is optional) sets the powder temperature the
    curve is interpolated at on the shot side; it defaults to `--temperature` when omitted.
    (With the linear `--powder-temp-sensitivity` model instead, `--powder-temp` keeps its
    existing "reference temperature" meaning, default 70°F/21°C — the models are mutually
    exclusive.)
  - `--zero-powder-temp` (new) sets the zero-day powder temperature for the curve on the
    `--auto-zero` side, defaulting to `--zero-temperature`. An explicit `--zero-velocity`
    still takes precedence.
  - Air temperature (`--temperature` / `--zero-temperature`) continues to exclusively drive
    air density. Fully backward compatible: with neither flag, the curve is looked up at the
    air temperature exactly as before. Native CLI + WASM; integration tests added.

### Fixed

- **WASM: an explicit `--zero-velocity` was overridden by a shot-day `--powder-temp-curve`
  during the zero solve.** The WASM zero solve cloned the shot inputs (carrying the curve),
  which re-interpolated over the supplied zero velocity. `--zero-velocity` now disables the
  curve for the zero solve, matching the native CLI.

## [0.22.4] - 2026-07-05

### Fixed

- **The WASM command surface now rejects unknown flags instead of silently ignoring
  them.** The hand-rolled WASM argument parser used a catch-all that dropped any
  unrecognized token, so a typo or a flag that isn't wired into the WASM surface looked
  like a successful no-op — making "no error" an unreliable signal that a flag is active.
  Unrecognized `--flags` now return an `Unknown flag: <name>` error, matching the native
  CLI's clap behavior. Applied to all WASM subcommands (trajectory, zero, monte-carlo,
  estimate-bc).

## [0.22.3] - 2026-07-04

### Added

- **MOA support for the PDF dope card (`--adjustment-unit mil|moa`).** The printable
  dope card was MIL-only; its Drop/Wind/Lead columns now render in MIL (default) or MOA,
  with the column sub-headers labeled accordingly. The MIL and MOA paths share one
  converter (`drop_to_adjustment`), and the moving-target lead is computed in the same
  unit. Backward compatible — the default is MIL. (MBA-724)

### Changed

- Unified the dope-card angular conversion on `drop_to_adjustment`, removing the
  now-redundant MIL-only `yards_to_mil` / `calculate_lead_mil` helpers.

## [0.22.2] - 2026-07-04

### Added

- **`--powder-temp-curve`: a measured powder-temperature → muzzle-velocity table.** The
  linear `--powder-temp-sensitivity` model assumes a constant fps/°F, but real powders are
  non-linear (temperature-stable powders flatten; others steepen when hot). This flag takes
  comma-separated `TEMP:VELOCITY` points (e.g. `"40:2620,70:2700,100:2760"`, in °F/fps or
  °C/m·s⁻¹ per `--units`) and interpolates the muzzle velocity at the ambient `--temperature`,
  clamped at the endpoints (no extrapolation beyond measured data). When supplied it
  supersedes the linear sensitivity model. It is data-driven — the shooter enters points they
  actually chronographed — rather than a guessed curve shape, and mirrors the existing
  `bc_segments` interpolation design.
- The curve composes with `--auto-zero`: with `--zero-temperature`, the zero angle is solved
  using the curve's velocity at the zeroing-day temperature (an explicit `--zero-velocity`
  still takes precedence). Wired into both the native CLI and the WASM command surface, with
  integration tests for exact points, interpolation, endpoint clamping, and backward
  compatibility.

## [0.22.1] - 2026-07-04

### Added

- **Zero-day condition overrides for `--auto-zero`.** New optional flags let the zero
  **angle** be solved under the conditions the rifle was actually zeroed in, while the
  trajectory runs under the current shot-day conditions — correctly modeling the point-of-
  impact shift when you sight in on one day/velocity and shoot on another (e.g. a cold vs.
  warm powder temperature drawn from a powder-temp/velocity table):
  - `--zero-velocity` — zero-day muzzle velocity (fps / m·s⁻¹)
  - `--zero-temperature` — zero-day air temperature (°F / °C)
  - `--zero-pressure` — zero-day barometric pressure (inHg / hPa)
  - `--zero-humidity` — zero-day relative humidity (percent)
  - `--zero-altitude` — zero-day altitude (feet / meters)

  Each flag is independently optional; any omitted flag falls back to the corresponding
  shot-day value, so supplying none of them reproduces the previous single-condition
  behavior exactly (verified: zero output diff). Wired into both the native CLI and the
  WASM command surface.

## [0.22.0] - 2026-07-01

A large correctness pass from an adversarial logic/math audit of the whole crate (44
confirmed findings across drag, atmosphere, ground/zero geometry, wind, spin/rotational
physics, Monte Carlo, and cross-surface consistency). This release changes numeric output
for most trajectories — re-validate any downstream reference tables. Highest-impact items
marked ⚠️.

### Fixed

**Drag & atmosphere**
- **⚠️ Ground plane no longer double-counts bore height.** `ground_threshold` was set to
  `-bore_height` while the solver already starts the bullet at a ground-referenced
  `muzzle_height`, so flat-fire shots flew until `2×bore_height` below the muzzle instead
  of 1 — ground-impact range and time-of-flight were overstated (~37%/~45% in a typical
  flat-fire case). Ground is now `y = 0` in all CLI/FFI/WASM/Monte-Carlo call sites.
- **⚠️ Removed a transonic drag-rise double-count.** The G1/G7 table Cd already embeds the
  full McCoy transonic rise; an additional multiplier (up to 2.5×, applied unconditionally
  near Mach 1) was stacking a second rise on top of it, inflating drop and reducing retained
  velocity for any shot touching Mach 0.85–1.2 — and for the entire flight of subsonic loads
  (.22LR, subsonic .300BLK, etc).
- **⚠️ Corrected `CD_TO_RETARD` to the ICAO-referenced value** (`2.08551e-4`, was the Army
  Standard Metro constant `2.049e-4`), matching density normalized to ICAO sea level
  (1.225 kg/m³). Drag was uniformly ~1.75% low for any ICAO-referenced BC.
- **Mach number and air density now use the same resolved atmosphere.** Speed-of-sound was
  computed from the raw input temperature while density used the ICAO-lapse-resolved
  station temperature at `--altitude`, so drag was evaluated at the wrong Mach for
  altitude shots with default atmosphere (a few percent error near transonic).
- **Miller Sg density correction, and the fast/Monte-Carlo Mach lookup, now use the same
  resolved station atmosphere** instead of the raw (often default sea-level) input —
  affects stability classification, spin drift, and Magnus at altitude.
- Moist-air speed-of-sound and CIPM-2007 enhancement/compressibility coefficients corrected
  (were using a specific-humidity coefficient with a mole fraction, and hPa where the CIPM
  constants expect Pa); both errors were small (<1.5%) but real.
- Magnus moment coefficient is now continuous at Mach 1.2 (was a 33% discontinuous jump
  between the transonic and supersonic branches).

**Zeroing & geometry**
- **`zero` subcommand now zeroes to the line of sight**, not the bore line — was ignoring
  `--sight-height` entirely (~1.9 MOA error at typical sight heights).
- **`--auto-zero` now solves with the user's actual drag model and atmosphere**, not a
  hardcoded G1/ISA-sea-level zero — a G7 auto-zero previously struck several inches off
  the intended zero distance.
- **`sight_adjustment_moa` sign fixed** (was subtracting the sight-height correction instead
  of adding it — told users to dial the wrong direction).
- **`point_blank_range` fixed** — was reading the lateral (always-zero) axis with a
  bore-referenced threshold; now downrange, referenced to the line of sight.
- **`true-velocity` drop formula no longer double-counts barrel elevation** and no longer
  flips the sight-height sign (was up to +11% effective-velocity error).
- **`--shooting-angle` (inclined fire) is now applied** — previously accepted by every
  surface (CLI/FFI/WASM) and by the low-level fast/library-binding integration path, but
  read by no solver; gravity is now rotated into the shot-aligned frame everywhere,
  including `fast_integrate_with_segments`.
- Trajectory `bullet_length` now honors `--bullet-length` instead of a hardcoded 4.5-caliber
  estimate, so stability/spin-drift/aerodynamic-jump match the `stability` subcommand for
  the same bullet.

**Monte Carlo**
- **CEP is now the median target-plane radial miss** (was `range_std_dev × 1.1774`, off by
  up to ~100×).
- **Fell-short samples are no longer misclassified as hits** — a downrange shortfall was
  being compared against the target-plane hit radius.
- **Vertical angular dispersion is no longer sampled twice** — muzzle-angle jitter and an
  independent "pointing error" of the same magnitude were both applied, inflating vertical
  spread ~41%.
- **`max_range` is now sized to the target**, so simulations no longer silently truncate at
  the 1000 m integrator default on longer shots.
- **FFI and WASM Monte Carlo now use the same ground/bore convention as the CLI** (previously
  only the CLI had a real ground-impact range; FFI/WASM reported the integrator cap with
  zero standard deviation).
- **FFI Monte Carlo now honors its `atmosphere` argument** instead of discarding it and
  running every sample at ICAO sea level.

**Wind**
- Wind beyond the last `--wind-segment` now correctly goes to zero on the shear-enabled fast
  path (was falling back to the first segment's wind).
- `fast_integrate_with_segments` (the library/binding entry point) now actually steps
  through wind segments by downrange distance instead of applying only the first segment
  for the whole flight.

**Cross-surface consistency (CLI vs FFI vs WASM)**
- `estimate-bc` drop values are now inches on both CLI and WASM (CLI was parsing them as
  yards — a 36× unit mismatch).
- `--twist-rate` is mm/turn in metric mode on WASM, matching the CLI (MBA-970).
- Powder-temperature modeling now matches between the native solver and WASM (WASM
  previously used a hardcoded 70°F/21°C reference and ignored `--temperature` entirely).
- `--powder-temp`'s default no longer gets reinterpreted as 70°C in metric mode.
- WASM now sets `ground_threshold` on `trajectory` (was the -100 m library default, a ~3.5×
  disagreement with the CLI on the same command).
- WASM's Drop column (table/JSON/CSV) is now relative to the line of sight instead of the
  bore line.
- WASM's default wind direction is now 0° (headwind), matching the CLI, instead of 90°
  (full crosswind).
- FFI trajectory samples now report real `mach` and `spin_rate_rps` instead of hardcoded 0.0.
- FFI `ballistics_quick_trajectory` drop is now referenced to the sight line, not the bore.
- Metric-mode `--bc-table`/`--bc-table-dir` lookups now convert to the tables' native
  inches/fps units instead of feeding them mm/m-s.
- `--bc-adjustment` now actually affects the trajectory when a `--bc-table` is also
  supplied (previously affected only the displayed BC).
- Metric-mode PDF dope cards now convert to the card's imperial fields instead of writing
  raw metric values under imperial labels.
- Density-altitude on the PDF dope card no longer applies the 120 ft/°C rule to a Fahrenheit
  delta (was overstating the temperature correction ~1.8×).
- `monte-carlo` help text now documents per-`--units` values instead of claiming SI
  regardless of the active unit system.

**Numerical**
- `fast_integrate`'s `FastSolution::sol()` no longer returns NaN when queried at its own
  reported event time (was pushing a duplicate final (time, state) sample).
- The sampled-trajectory `Apex` flag is no longer placed on the first sample for flat or
  descending shots that never rise above the muzzle.

### Performance
- Speed-of-sound/temperature/pressure are now resolved once per trajectory solve and
  threaded through drag and Magnus calculations, instead of being re-resolved (including a
  vapor-pressure polynomial) on every RK4/RK45 stage.

## [0.21.5] - 2026-06-27

Two monte-carlo correctness fixes (MBA-967, MBA-971).

### Fixed
- **monte-carlo "Mean Range" is now meaningful and target-independent (MBA-967).**
  Simulations used `BallisticInputs::default()` (muzzle_height 0, ground_threshold -100 m)
  so each sim flew to the integrator's range cap instead of a real ground impact — "Mean
  Range" reported the cap (~1000 m), not the trajectory's ground-impact range. And a
  per-sample `max_range < target_distance` skip tied the range statistic to
  `--target-distance`. MC now uses the same bore-height/ground convention as `trajectory`
  (1.5 m) and computes range/velocity over all samples, so Mean Range matches `trajectory`
  for the same physics and no longer changes with the target distance.
- **monte-carlo `hit_probability` unified across the CLI and FFI (MBA-971).** The CLI
  counted ground-impact ranges within 1 m of the target (a range-precision notion that read
  0% for any target short of the impact range) while the FFI used a position notion with a
  redundant clause. Both now use a single `MonteCarloResults::hit_probability(radius)`
  method — the fraction of samples landing within the hit radius of the point of aim at the
  target plane. Samples that fall short are recorded as definite misses, so an unreachable
  target correctly yields ~0% (was a spurious ~100%).

### Added
- **`monte-carlo --target-radius` flag** (default 0.3 m, unit-aware) to define the hit zone
  used by hit probability, instead of a hardcoded radius.

## [0.21.4] - 2026-06-26

Ten CLI/correctness fixes surfaced by an adversarial edge-case sweep (MBA-960..966,
968, 969, 970). Two are behavioral changes existing consumers should note:

### Fixed
- **⚠️ Metric mode no longer runs in a near-vacuum (MBA-960/961).** The CLI atmosphere
  defaults (`29.92` / `59.0`) were imperial literals applied before `--units` was parsed,
  so in metric mode they were read as 29.92 hPa (~24 km) and 59 °C — drag ~30× too low,
  impact velocity ~42% high. `--pressure`/`--temperature` are now resolved and validated
  *after* `--units` (per-unit defaults + ranges). **Metric trajectories without explicit
  atmosphere flags now produce different (correct) output than before.** Also fixes the
  `--pressure` validator (was an hPa range applied to imperial inHg: rejected valid `14`,
  accepted an hPa typo like `1013`).
- **⚠️ JSON output now honors `--units` (MBA-962).** The JSON formatter always emitted raw
  SI regardless of `--units` (table/csv were already correct). It now converts to the
  requested units and adds an explicit `"units"` field. **JSON consumers using `--units
  imperial` previously received metric numbers; they now receive imperial.**
- **`--use-powder-sensitivity` now affects the trajectory (MBA-963).** It was a no-op on
  the CLI solve path (the muzzle-velocity adjustment was computed but never applied); the
  per-degree sensitivity delta was also mis-converted. Now wired into `TrajectorySolver`.
- **`spin_drift` summary is consistent with the applied path (MBA-964).** It reported
  `null` while the path silently applied drift from a default 1:12 twist.
- **`--wind-shear-model` is now honored (MBA-965).** The selector was discarded (everything
  ran as PowerLaw) and invalid values were accepted; models are now distinct and validated.
- **Pitch-damping / precession diagnostics are serialized (MBA-966).** `--enable-pitch-damping`
  / `--enable-precession` computed diagnostics that the CLI never emitted; now in JSON.
- **Default RK45 reaches `--max-range` exactly (MBA-968).** It stopped ~2% short (no
  boundary interpolation), disagreeing with Euler/RK4-fixed; the final state is now
  interpolated to the boundary.
- **Time-cap runs are labeled "no impact" (MBA-969).** An upward shot with
  `--ignore-ground-impact` reported a bogus impact at the hidden 100 s integration cap.
- **`--twist-rate` is mm/turn in metric for `trajectory` (MBA-970).** It was read as inches
  regardless of `--units`, so a metric 1:10 twist gave zero stability; now matches the
  `stability` subcommand.

## [0.21.3] - 2026-06-25

### Fixed
- **Direct-atmosphere mode no longer rejected as degenerate (0.21.2 regression).** The
  0.21.2 degenerate-atmosphere guard rejected the legitimate direct-atmosphere tuple
  `(air_density, speed_of_sound, 0.0, 0.0)` because slots 2 and 3 are `0.0` sentinels and
  the guard treated "pressure ≤ 0" as non-physical. The guard now recognizes direct mode
  (real density < 2.0 kg/m³ and speed of sound > 200 m/s with zero pressure/ratio slots)
  and only rejects genuinely non-physical input. Test: `fast_path_rejects_degenerate_atmosphere`
  now also asserts direct mode succeeds.
- **Transonic Mach transitions are now flagged in sampled output.** The trajectory sampler's
  `transonic_distances` was a `Vec::new()` placeholder (a TODO), so no sampled point ever
  received a `MachTransition` flag regardless of trajectory. The Euler / RK4 / RK45 solvers
  now record the downrange distances where the projectile crosses Mach 1.2 (transonic) and
  Mach 1.0 (subsonic), and `add_trajectory_flags` marks the nearest sample at each crossing.
  Test: `transonic_crossing_flags_a_sampled_point`.
- **MBA-717: fast/MC path no longer hardcodes bullet geometry.** `build_inputs` used fixed
  placeholders (0.308" diameter / 1.24" length / 10" twist) because `TrajectoryParams`
  didn't carry the geometry, so spin-drift / Magnus / stability on that path ignored the
  real bullet. `TrajectoryParams` now carries `bullet_diameter` / `bullet_length` /
  `twist_rate`, plumbed from the inputs. (Behavioral impact is small today — the fast-path
  yaw-of-repose has no angular state, so its Magnus/spin-drift is ~0 — but the data is now
  correct for callers that do use it and for future work on that path.)

### Changed
- **MBA-722: humidity-scale convention documented + centralized.** `BallisticInputs.humidity`
  is a 0–1 fraction while `AtmosphericConditions.humidity` is a 0–100 percent (the scale the
  atmosphere density helpers expect) — an easy 100× footgun. Added a `humidity_percent()`
  helper that does the clamped 0–1 → 0–100 conversion, used it at the Monte-Carlo boundary,
  and cross-documented both fields. No numerical change (existing paths already converted).
- **Documented the dual-mode `atmo_params` tuple.** The atmosphere 4-tuple has two modes
  (standard `(alt, temp_c, pressure_hPa, density_ratio)` vs direct `(density, sound, 0, 0)`)
  and slot 3 is a density RATIO that rides in the `humidity` field — an easy footgun. Added
  doc comments at `TrajectoryParams.atmos_params`, `FastIntegrationParams.atmo_params`, and
  the `atmo_is_physical` guard spelling out both modes. Documentation only.

## [0.21.2] - 2026-06-25

### Changed
- **Coriolis is now independent of the advanced-effects umbrella on the fast/MC path.**
  `fast_integrate_with_segments` previously gated its Earth-rotation ω — and spin-drift /
  Magnus — together on `enable_advanced_effects`, so a caller could not request
  Coriolis-only. The ω is now gated on `enable_coriolis` (+ a latitude), and Magnus on
  `enable_magnus`, matching the RK4 solver. A caller can enable Coriolis with
  `enable_advanced_effects = false`. Tests: `fast_path_coriolis_independent_of_advanced_effects`.

### Fixed
- **Degenerate atmosphere no longer reports success.** A non-physical `atmo_params`
  (pressure ≤ 0 or non-finite — typically a unit mix-up, e.g. inHg passed where hPa is
  expected) yielded zero/NaN air density and silently truncated the fast path to a
  single-point stub trajectory marked `success = true`. The fast entrypoints now validate
  the atmosphere up front and return `success = false` for degenerate input. Test:
  `fast_path_rejects_degenerate_atmosphere`.

## [0.21.1] - 2026-06-24

### Fixed
- **Coriolis shot direction on the fast / Monte-Carlo path.** 0.21.0 fixed the RK4
  solver but missed a third Earth-rotation construction site: `fast_integrate_with_segments`
  (the path the Python binding uses) still built ω from `azimuth_angle` instead of
  `shot_azimuth`, so east and west shots were identical there. Now uses `shot_azimuth`,
  matching the RK4 and `fast_integrate` paths — directional Coriolis (Eötvös) is applied
  consistently across all solver paths. Regression test
  `fast_trajectory::tests::fast_path_coriolis_uses_shot_direction`. (This corrects the
  "fast/MC path is still a no-op" caveat noted in 0.21.0.)

### Changed
- **DOWNSTREAM:** the fast/MC path now changes output for Coriolis shots with a non-North
  bearing (north / unset unchanged). The Python binding accepts a `shot_direction` key
  (degrees, 0=N) to drive it; pass it alongside `latitude` + advanced effects.

## [0.21.0] - 2026-06-24

### Fixed
- **Coriolis now responds to shot direction (Eötvös effect).** `--shot-direction` was
  stored but never reached the trajectory solver, so the Coriolis correction was always
  computed as a due-North shot — east- and west-fired shots produced identical results.
  The solver now carries the firing compass bearing in a dedicated `shot_azimuth` input
  (distinct from `azimuth_angle`, the small aiming offset that rotates the launch
  velocity), and both the RK4 solver and the fast/Monte-Carlo path use it for the
  Earth-rotation vector. An east shot now gets the upward Eötvös term (`+2Ω·cosφ·v_east`,
  shoots higher) and a west shot the downward one (shoots lower); lateral drift is
  unchanged. Regression test `coriolis_direction_tests::eotvos_east_higher_than_west`.

### Changed
- **DOWNSTREAM:** numerical output changes for Coriolis shots with a non-North
  `--shot-direction` (previously those were silently treated as North; North / unset is
  unchanged). FFI `FFIBallisticInputs` gains an appended `shot_azimuth` field (radians,
  0=N) — existing field offsets are preserved, but FFI consumers that construct the full
  struct must add it. WASM `runCommand` now accepts `--shot-direction`. The fast/MC path
  still does not plumb latitude/bearing on its own (`build_inputs` hardwires
  `latitude: None`), so directional Coriolis there remains out of scope for this release.

## [0.20.0] - 2026-06-21

Dependency modernization — a coordinated upgrade of the major dependencies. No new
features, CLI flags, or physics; numerical trajectory output is unchanged. The one
observable difference is in Monte Carlo: under the new RNG the per-run sample
*sequence* differs, but the statistical results (mean, std-dev, CEP, confidence
ellipse) are unchanged.

### Changed
- **`ndarray` 0.15 → 0.17** and **`ndarray-npy` 0.8 → 0.10** (upgraded together; the
  drag-table `.npy` loader is unaffected).
- **`rand` 0.8 → 0.9** and **`rand_distr` 0.4 → 0.5**. Monte Carlo uses an unseeded
  RNG, so this only changes the sample sequence, not the statistical behavior.
- **`ureq` 2.12 → 3.1** (online feature). Migrated to the 3.x API; the public error
  mapping (`ServerError(code, body)`, timeout vs. network) and the unbounded
  BC-table read are preserved.
- **`getrandom` 0.2 → 0.3** with the `wasm_js` backend feature (WASM builds).
- CI: `actions/upload-artifact` 6 → 7, `actions/download-artifact` 7 → 8,
  `cross-platform-actions/action` 0.30 → 0.32.

## [0.19.0] - 2026-06-20

Adds downrange-segmented wind to the CLI and WASM, and **corrects the wind-direction
sign on the `cli_api`/`TrajectorySolver` path** (a breaking output change — see below).

### Added
- **Downrange-segmented wind** — wind can now vary along the trajectory (e.g. a muzzle
  reading plus downrange sensor stations).
  - CLI: repeatable `--wind-segment SPEED:ANGLE:UNTIL_DISTANCE` on the `trajectory`
    command. SPEED and UNTIL_DISTANCE follow `--units` (mph & yards imperial, m/s &
    meters metric); ANGLE is degrees in the wind-FROM convention (same as
    `--wind-direction`). Each segment applies from the previous boundary out to its
    UNTIL_DISTANCE; wind is a step function with **zero wind beyond the last segment**
    (a coverage warning is printed when the segments don't reach `--max-range`).
    Segments **override** scalar `--wind-speed`/`--wind-direction`, and are **not
    compatible with `--enable-wind-shear`** (rejected with an error).
  - WASM: `runCommand` accepts the same repeatable `--wind-segment` flag, and the
    `Calculator` gains `addWindSegment(speedMph, directionDeg, untilYards)` and
    `clearWindSegments()`.
  - API: `TrajectorySolver::set_wind_segments(Vec<WindSegment>)` performs a per-step
    downrange lookup via `WindSock`, preserving all solver physics (Coriolis, spin
    drift, aerodynamic jump, zero-finding, full `TrajectoryResult` output). A solve
    with no segments is numerically identical to before. New public
    `wind::parse_wind_segment_str(&str, imperial)`.

### Fixed
- **CLI `--twist-right` can now select left-hand twist.** It was a presence-only flag
  (always true; `--twist-right false` errored), so left-hand twist — and the
  aerodynamic-jump direction reversal it produces — was unreachable from the CLI. It now
  takes an optional value: `--twist-right` / `--twist-right true` = right-hand (default),
  `--twist-right false` = left-hand (the bare flag still works, so existing usage is
  unaffected).

### Changed (BREAKING)
- **Wind direction sign corrected on the `cli_api`/`TrajectorySolver` path.** The scalar
  wind vector was built with the opposite sign to the standard wind-clock convention
  (and to the engine's own `WindSock`, the Monte-Carlo path, and the Python/Ruby
  bindings): `--wind-direction 0` produced a *tailwind* and `90°` drifted the bullet the
  wrong way. It is now correct: **0° = headwind, 90° = from the right (drifts left),
  180° = tailwind, 270° = from the left** — matching `WindSock`/the bindings.
  **Trajectories computed with a non-zero wind direction will change** (head/tail effect
  and crosswind drift direction flip); no change when wind speed is 0. This affects
  **every consumer of the `TrajectorySolver`/cli_api path**: the CLI `trajectory`,
  `wind-card`, `range-table`, and `come-ups` subcommands; the WASM `runCommand` and
  `Calculator`; the C FFI (`FFIWindConditions.direction`, for iOS/Android callers); and
  any direct `TrajectorySolver` user. The Monte-Carlo path and the Python/Ruby bindings
  already used the correct convention and are unchanged. The aerodynamic-jump crosswind sign
  on this path was aligned to match (right twist + wind-from-right → impact up, drift
  left), bringing `cli_api` into agreement with the fast-integrate path. The wind-shear
  path was already sign-aligned with the scalar path and flips with it.
- The `--wind-direction` help text was corrected (previously the inaccurate
  "0=North, 90=East").

## [0.18.3] - 2026-06-19

Extends aerodynamic jump (MBA-959) to the fast-integrate path. Default-off; no change to existing trajectories.

### Added
- **Aerodynamic jump in `fast_trajectory::fast_integrate`** — the fast fixed-step kernel (used by the `ballistics-engine-py` `fast_integrate` binding and the Monte-Carlo path) now applies aerodynamic jump when `enable_aerodynamic_jump` is set, by perturbing the prebuilt launch velocity vertically. New public `fast_trajectory::aerodynamic_jump_launch_offset_rad(inputs, atmo_params)` returns the Litz vertical offset (radians) from the engine's Miller Sg; crosswind is taken from `wind_speed`/`wind_angle` (BallisticInputs convention: 0 = headwind, +90° = from the right). Previously only the `cli_api::TrajectorySolver` (Euler/RK4/RK45) path applied AJ; the fast path ignored the flag.

## [0.18.2] - 2026-06-19

Adds an opt-in aerodynamic-jump effect (MBA-959). Default-off and additive — no change to existing trajectories.

### Added
- **Aerodynamic jump** (`enable_aerodynamic_jump`, default off): the vertical point-of-impact shift a crosswind imparts to a spin-stabilized bullet, applied as a muzzle launch-angle perturbation in the solver. Uses Bryan Litz's estimator `Y = 0.01·Sg − 0.0024·L + 0.032` (MOA per mph of crosswind), fed by the engine's own Miller stability factor (Sg) and bullet length in calibers. Exposed on `TrajectoryResult.aerodynamic_jump`, via the CLI `--enable-aerodynamic-jump` flag, and as the public `aerodynamic_jump::litz_crosswind_jump_moa()`. The jump is purely vertical; direction follows Litz (right twist: crosswind from the right → impact up). The zero is found on the bare bore, so the jump appears as an additive POI shift rather than being absorbed by the zero. Most accurate near Sg ≈ 1.75 (Litz regression validity).

### Notes
- The legacy heuristic `aerodynamic_jump::calculate_aerodynamic_jump` is retained for backward compatibility but is superseded by the Litz estimator and is no longer used by the solver.

## [0.18.1] - 2026-06-18

Two follow-up fixes after 0.18.0. Neither changes trajectory output.

### Fixed
- **BCCR drag-correction loader** now verifies the stored CRC32 (a standard IEEE/zlib CRC over the data section) instead of reading and discarding it — corrupt-but-in-range files are rejected with `ChecksumMismatch` instead of loading silently. Gated on a non-zero stored checksum so older checksum-less files still load (MBA-953).
- **Precession/nutation angular diagnostics** — corrected four dimensional errors: the precession and nutation frequencies are now rad/s via the standard epicyclic form `φ = (Iₓp/2Iy)[1 ± √(1 − 1/Sg)]` (with the dimensionless Miller stability factor), and the yaw no longer random-walks. Diagnostic-only — these feed only the opt-in `max_yaw_angle` / `max_precession_angle` outputs (`--enable-precession-nutation`); no trajectory change (MBA-941).

## [0.18.0] - 2026-06-18

A cross-solver consistency and physics-correction pass (JIRA MBA-938..958). 17 commits since v0.17.0; build and full test suite (210+) green. Several changes correct real bugs that **shift numerical results** on the Coriolis, transonic/subsonic drag, spin-drift-at-altitude, and WASM-zero paths — see Changed and validate against a reference before deploying downstream.

### Fixed
- **Coriolis deflection direction** — corrected to the physical `-2 Ω×v` with the right lateral ω sign; a Northern-hemisphere shot now drifts right (was left), validated vs py_ballisticcalc. The fast/Monte-Carlo path (used by the Python/Ruby bindings) applied no Coriolis at all and now matches the default solver to ~1%.
- **Transonic drag resolution** — the engine silently used a coarse 21-point G1/G7 fallback (a runtime loader looked for a directory that never exists at runtime); the high-resolution 79/84-point tables are now baked in via `include_str!`. Transonic drop error vs py_ballisticcalc fell from +43" to +3" at 1300 yd.
- **Custom drag tables** (`custom_drag_table`) are now honored by all three solvers (were plumbed onto the inputs but never read — the G-model was always used).
- **Named bullet shapes** (boat/round/flat in `bullet_model`) now drive the transonic shape on every solver (the integrate_trajectory path ignored the name and used a caliber/weight heuristic).
- **`use_form_factor`** is honored consistently across all three solvers (was derivatives-only; no-op by default).
- **WASM standalone zero** targets the line of sight (sight height), not the bore line — it was off by the sight-height angle (~2 MOA at 100 yd / 2" sight).
- **PDF output** is gated behind the `pdf` feature so `--no-default-features` builds the binary.
- **Sampled `--full` CSV** reports drop/drift in inches (Imperial) / meters (Metric), matching the table and API (was yards).
- Removed a dead WASM `--wind-dir-std` setter; documented the wind-direction-spread approximation.

### Changed
> **FLAGGED for downstream:** these correct real bugs but **shift numerical output**. Re-validate dope and any cached/expected outputs.
- **Coriolis** — direction fix plus the previously-absent Monte-Carlo/binding application; affects any Coriolis-enabled trajectory.
- **Transonic/subsonic drag** — the high-resolution G1/G7 tables shift drag below ~Mach 1.3; the low-velocity **Reynolds correction was removed** (it ran on only one of three solvers, py_ballisticcalc does not model it, and the default solver already matches pbc subsonically), shifting subsonic (<200 m/s) shots on the binding path.
- **Spin drift at altitude** — the Miller stability density correction is now the canonical linear `(T/T0)(P0/P)` (was `sqrt`), raising spin drift at altitude (sea level unchanged). The same correction now also feeds the Magnus / yaw-of-repose Sg (off by default).
- **WASM zero** shifts by the sight-height angle (~2 MOA at 100 yd) — now a proper sight-line zero.

### Added
- `transonic_drag::resolve_projectile_shape(...)` — shared bullet-shape resolution (name first, then heuristic) used by all three solvers.
- `TrajectoryParams.ground_threshold` — the integrate_trajectory ground plane is configurable (was hardcoded `-1000.0`; default preserved).

### Performance (output-identical)
- Velocity-BC segments are estimated **once** at solver setup on the integrate_trajectory path (the Python-binding hot path) instead of being rebuilt — allocating a string + a vector — on every derivative evaluation. Byte-identical output.

## [0.17.0] - 2026-06-17

An autonomous adversarial hardening pass plus a follow-up altitude/temperature fix. 82 fixes across 78 commits since v0.16.0; build and full test suite green. Most changes are correctness/robustness fixes that preserve output, but several correct real bugs that **shift numerical results** — see Changed below and validate against a reference before deploying downstream.

### Fixed
- **RK45 time-of-flight desync** — the default solver advanced `time` by the *next* step's dt, corrupting time_of_flight and per-point/sample times.
- **Arden-Buck vapor pressure** — the linear factor sat outside `exp()`, over-estimating saturation pressure (~7x) and corrupting humidity-dependent air density.
- **Gyroscopic spin-drift twist sign** — applied twice (canceling), so left-twist barrels drifted the wrong way (right-twist/default unchanged).
- **API twist-rate unit** — sent in meters where inches/turn was expected (~33x off); now sent verbatim as inches/turn across API/FFI/WASM.
- **Wind-layer extrapolation** — above the top custom layer the profile extrapolated garbage instead of clamping; WindSock cursor no longer skips segments on large query jumps.
- **Transonic drag** applied exactly once everywhere (was double-applied in `derivatives.rs`, missing in `fast_trajectory.rs`, shape-misclassified in `cli_api.rs`); **Reynolds drag** applied once and capped at 5x (was double-applied / uncapped ~16x spike).
- **Euler solver** now uses the shared acceleration kernel — previously ignored the ballistic coefficient (diverging up to ~2.3x) and omitted Magnus/Coriolis.
- **Projectile shape (default solver)** received caliber in meters instead of inches, under-applying transonic drag for G1/G2/G5/G6/G8/GI/GS.
- **Monte Carlo** — fixed long-range truncation and a 100x humidity mis-scale (fraction vs percent); BC-correction table no longer silently zeros the BC when queried with the other drag model.
- **Drag model selection** — G2/G5/GI/GS now explicitly alias to G1 (documented).
- **Output** — corrected dope-card Drop MIL sign, `--full` table downrange/drift column swap, dope-card footer weekday (was off by 4 days), stability min-twist units label, and estimate-bc verification range.
- **WASM** — auto-zero solves to line-of-sight (not bore); wind direction consumed as radians; Monte Carlo honors `--drag-model`; full-trajectory JSON keys fixed; wasm32 target now compiles; help text synced to the parser and dead zero flags dropped.
- **Edge cases / divide-by-zero / NaN guards** — sectional-density and BC zero divisors, aero-jump zero/negative twist and `mass_kg<=0`, moist-air speed of sound at zero pressure, air density at non-positive pressure, zero-angle solver reporting non-convergence as success.
- **Robustness** — BC5D/BCCR table loaders use checked-arithmetic dimension bounds and NaN-safe lookups; FFI stops leaking a `CString` per version call and caps oversized/negative Monte Carlo sim counts; fixed panics on PDF UTF-8 header truncation, `percentile(p>1)` OOB, and Euler zero-relative-velocity.
- **Input validation** — CLI rejects non-finite numeric arguments; `target_distance` must be finite and positive (was `> 0`, letting `+inf` through); only `-inf` maps to the ignore-ground sentinel.

### Changed
> **FLAGGED for downstream:** these correct real bugs but **shift numerical output** — altitude, air-density, gravity, and transonic results all move. Re-validate dope and any cached/expected outputs.
- **Air density unified** across all three integrators — `cli_api` no longer double-counts altitude (`exp(-alt/8000)` atop station pressure) and now applies humidity; Monte Carlo no longer reprojects shooter-altitude density to sea level. Shifts *every* trajectory (~0.32% at sea level, larger at altitude).
- **Altitude now drives both station pressure and temperature** — with a default (15 °C / sea-level) input, `--altitude` derives ICAO station pressure *and* applies the −6.5 °C/km lapse rate; an explicit pressure or temperature stays authoritative. Previously altitude-only ran too warm/thin. Validated against py_ballisticcalc to ~0.04% (0–3000 m).
- **Gravity** is the SI-canonical `9.80665` in the default RK4/RK45 kernel (was `9.81`).
- **Transonic single-application** and **Reynolds once + capped** change near-Mach-1 and Monte Carlo drag.
- **Spin twist sign** and **aero-jump Sg** (kg→grains, was ~1000x off) affect opt-in spin/jump diagnostics; left-twist spin-drift direction changes.
- **Default `bullet_length` 4.0→4.5x caliber** (now unified across library/FFI/WASM) — shifts opt-in spin/Magnus/Miller diagnostics only, not the base trajectory.
- **Zero solver returns `Err` on non-convergence** (was a best-effort `Ok(angle)`) — downstream callers must handle `Err`.
- **WASM default drag model G7→G1** to match the G1-scale default BC.
- **Monte Carlo range/velocity means exclude short-falling draws** when no explicit target is given (keeps CEP measurable and FFI arrays equal-length); pass an explicit target for unbiased per-target stats.

### Added
- `atmosphere::resolve_station_pressure(pressure_hpa, altitude_m) -> Option<f64>` — an explicit station pressure stays authoritative; a default pressure with real altitude derives the ICAO station pressure.
- `atmosphere::resolve_station_temperature(temperature_c, altitude_m) -> Option<f64>` — mirrors the above for temperature (returns `None` at the 15 °C default so the ICAO lapse rate applies).

### Performance (output-identical)
- Build `BallisticInputs` once per integration instead of per derivative eval (eliminates per-step String alloc and Vec/DragTable clones on the Monte Carlo hot path).
- Binary-search the drag Mach axis; precompute WindSock segment vectors; hoist RK45 air-density/wind and per-step pitch-damping/drag-coefficient allocations; skip BC-segment re-sort when already sorted; lowercase the bullet model once.

## [0.13.30] - 2026-01-26

### Added
- **Terms of Service Acceptance** - When using `--online` for the first time, users must accept the Terms of Service from https://ballistics.rs/terms.txt
  - TOS is fetched from the server and displayed in the terminal
  - User must type 'y' or 'yes' to accept
  - Acceptance is stored in `~/.ballistics/tos.json`
  - Subsequent runs skip the prompt unless TOS version changes

### Dependencies
- Added `dirs` crate for cross-platform home directory detection

## [0.13.29] - 2026-01-26

### Added
- **Expanded Test Coverage** - Added 36 new unit tests to modules with lighter coverage:
  - `trajectory_integration`: RK4/RK45 consistency, ground impact detection, target distance, wind handling (+6 tests)
  - `fast_trajectory`: Solution interpolation edge cases, BC segment boundaries, event arrays (+6 tests)
  - `stability_advanced`: Bullet type parameters, atmospheric corrections, dynamic stability (+11 tests)
  - `spin_drift_advanced`: Drift direction, edge cases, transonic correction, yaw of repose (+13 tests)

### Changed
- Total test count increased from 156 to 192 tests

## [0.13.28] - 2026-01-26

### Fixed
- **Wind Direction CSV Override** - Fixed default check comparing against 90.0 instead of 0.0 (the actual CLI default), which prevented location CSV wind direction overrides from working

## [0.13.27] - 2026-01-26

### Fixed
- **Location CSV Overrides** - Fixed bug where `final_humidity` and `final_wind_direction` were set but not used, causing location CSV overrides for humidity and wind direction to be ignored
- **Compiler Warnings** - Eliminated all unused variable warnings in main.rs

### Changed
- Location CSV now properly overrides humidity (`HUMIDITY` column) and wind direction (`WIND_DIR` column)

**Related Tickets:**
- MBA-600: Fix unused variable warnings in ballistics-engine CLI

## [0.13.26] - 2026-01-26

### Added
- **Longitude Parameter** - `--longitude` for weather zone features (degrees, -180 to 180)
- **Shot Direction Parameter** - `--shot-direction` for azimuth (0=North, 90=East)

### Notes
- Negative values require equals format: `--longitude=-115.2`
- These parameters are required when using `--enable-weather-zones` or `--enable-3d-weather`

**Related Tickets:**
- MBA-601: Add weather integration control parameter for --online CLI mode

## [0.13.25] - 2026-01-26

### Added
- **Weather Zone Control** - `--enable-weather-zones` flag for automatic weather zone generation
- **3D Weather Control** - `--enable-3d-weather` flag for altitude-dependent atmospheric corrections
- **Wind Shear Model Selection** - `--wind-shear-model` parameter (none, logarithmic, power_law, ekman_spiral)
- **Weather Zone Interpolation** - `--weather-zone-interpolation` parameter (linear, cubic, step)

### Notes
- All weather parameters are passed through to Flask API when using `--online` mode
- Weather features require latitude, longitude, and shot-direction to be specified

**Related Tickets:**
- MBA-601: Add weather integration control parameter for --online CLI mode

## [0.13.24] - 2026-01-26

### Fixed
- **Online Mode Ground Impact** - Fixed `--ignore-ground-impact` flag not being passed to Flask API when using `--online` mode, causing trajectories to be truncated when bullet dropped below -100m threshold

### Added
- `ground_threshold` field in API client TrajectoryRequest

**Related Tickets:**
- Online solver returning partial results investigation

## [0.13.23] - 2026-01-25

### Added
- **BC Truing** - `--bc-adjustment` parameter for multiplying BC by a correction factor (e.g., 0.85)
- **Velocity Truing** - `--velocity-adjustment` parameter for adding offset to base velocity from chronograph data
- **CSV Profile Import** - `--profile` and `--profile-row` to load gun configurations from CSV files
- **CSV Location Import** - `--location` and `--site` to load shooting location environmental data from CSV files
- CSV parser handles Glenn's format (header with `#` prefix, comma-separated values)

### Changed
- `--velocity` and `--bc` are now optional when using `--profile` (loaded from CSV)
- Integration tests now use main `ballistics` binary instead of `ballistics-cli`

**Related Tickets:**
- MBA-591: CLI: Add --bc-adjustment parameter for BC truing
- MBA-592: CLI: Add --velocity-adjustment parameter
- MBA-593: CLI: Add CSV profile import (--profile)
- MBA-594: CLI: Add location presets (--location)

## [0.5.0] - 2025-01-31

### Added
- **Advanced Spin Drift Module** (`spin_drift_advanced.rs`) - Enhanced gyroscopic drift calculations with velocity-dependent effects (MBA-138)
- **Advanced Stability Module** (`stability_advanced.rs`) - Comprehensive stability analysis including transonic behavior (MBA-139)
- Both modules ported from private ballistics_rust repository as part of strategic reconciliation (MBA-144)

### Fixed
- **Critical Wind Vector Calculation** - Fixed swapped sin/cos in `wind.rs` `calc_vec()` function that caused 90° rotation of wind effects (MBA-146)
- **Wind Shear Bugs** - Multiple critical fixes in `wind_shear.rs` for proper altitude-dependent wind modeling (MBA-145)

### Verified
- **Transonic Drag Parity** - Verified parity with private repository implementation, all public functions available (MBA-140)

### Context
This release represents significant progress in the strategic reconciliation between the private ballistics_rust and public ballistics-engine repositories. The strategy established in MBA-144 v0.30.0 makes ballistics-engine the canonical physics source, with improvements flowing from private → public → wrapped in private.

**Module Reconciliation Progress:**
- 7 modules wrapped in ballistics_rust using ballistics-engine
- 2,164 lines of duplicate code eliminated (25.5% reduction)
- ballistics_rust now depends on ballistics-engine as authoritative source

**Related Tickets:**
- MBA-144: Strategic reconciliation with ballistics-engine as canonical source
- MBA-138: Port spin_drift_advanced.rs
- MBA-139: Port stability_advanced.rs
- MBA-140: Verify transonic_drag.rs parity
- MBA-142: Align Rust dependencies (nalgebra 0.33, rayon 1.10)
- MBA-145: Fix critical wind_shear.rs bugs
- MBA-146: Wind vector calculation bug
- MBA-147: Reynolds.rs missing FlowRegime enum

## [0.4.3] - 2025-01-20

### Fixed
- Fixed metric unit calculations
- Updated WASM to use RK45 integration method

## [0.4.0] - 2025-01-15

### Added
- RK45 integration method as default (improved accuracy)
- Muzzle height and target height parameters for realistic trajectory simulation

### Changed
- Major physics improvements across trajectory calculations

## [0.3.0] - 2024-12-01

### Added
- Initial public release
- Core trajectory physics
- RK4 and Euler integration methods
- Monte Carlo simulation support
- WASM support for web deployment

[0.5.0]: https://github.com/ajokela/ballistics-engine/compare/v0.4.3...v0.5.0
[0.4.3]: https://github.com/ajokela/ballistics-engine/compare/v0.4.0...v0.4.3
[0.4.0]: https://github.com/ajokela/ballistics-engine/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/ajokela/ballistics-engine/releases/tag/v0.3.0
