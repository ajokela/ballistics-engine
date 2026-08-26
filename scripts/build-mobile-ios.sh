#!/usr/bin/env bash
# Build BallisticsEngine.xcframework for the Ballistics Insight iOS app.
#
# Mobile feature set: bridge + ffi + pdf + profile-import, WITHOUT online (no
# network dependency, cleaner App Store privacy story) and WITHOUT cli.
# Keep panic = "unwind" (the bridge's catch_unwind depends on it) — do not add a
# panic="abort" profile for these builds.
#
# Output: target/mobile/BallisticsEngine.xcframework
set -euo pipefail

ENGINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ENGINE_DIR"

FEATURES="bridge,ffi,pdf,profile-import"
OUT="$ENGINE_DIR/target/mobile"
HEADER_DIR="$OUT/include"

# Simulator architectures.
#
# arm64 only by default. The x86_64 slice exists solely so an INTEL Mac can run
# the simulator, and Xcode 26 is arm64-only — its own xcodebuild, clang and
# swiftc carry no x86_64 slice — so no Intel Mac can host a simulator on a
# current toolchain at all. Building it cost a third of this script's runtime
# and 128 MB inside the fat library, for a slice nothing could load.
#
# Set BALLISTICS_IOS_SIM_X86_64=1 if you are on an older Xcode that still runs
# on Intel and genuinely need it.
SIM_X86_64="${BALLISTICS_IOS_SIM_X86_64:-0}"

TARGETS=(aarch64-apple-ios aarch64-apple-ios-sim)
[ "$SIM_X86_64" = "1" ] && TARGETS+=(x86_64-apple-ios)

rustup target add "${TARGETS[@]}" >/dev/null

for target in "${TARGETS[@]}"; do
  echo "==> $target"
  cargo build --release --no-default-features --features "$FEATURES" --target "$target"
done

mkdir -p "$OUT" "$HEADER_DIR"
cp "$ENGINE_DIR/include/ballistics_bridge.h" "$HEADER_DIR/"
# Legacy per-feature ABI header, if the app still uses it during migration.
[ -f "$ENGINE_DIR/BallisticsEngine.h" ] && cp "$ENGINE_DIR/BallisticsEngine.h" "$HEADER_DIR/" || true

cat > "$HEADER_DIR/module.modulemap" <<'EOF'
module BallisticsEngine {
    header "ballistics_bridge.h"
    export *
}
EOF

# Simulator library. `lipo -create` with one input still produces a valid
# single-architecture archive, so the xcframework layout below is unchanged
# either way — only the architectures inside the simulator slice differ.
SIM_FAT="$OUT/libballistics_engine_sim.a"
SIM_INPUTS=("target/aarch64-apple-ios-sim/release/libballistics_engine.a")
[ "$SIM_X86_64" = "1" ] && SIM_INPUTS+=("target/x86_64-apple-ios/release/libballistics_engine.a")
lipo -create "${SIM_INPUTS[@]}" -output "$SIM_FAT"

rm -rf "$OUT/BallisticsEngine.xcframework"
xcodebuild -create-xcframework \
  -library "target/aarch64-apple-ios/release/libballistics_engine.a" -headers "$HEADER_DIR" \
  -library "$SIM_FAT" -headers "$HEADER_DIR" \
  -output "$OUT/BallisticsEngine.xcframework"

echo "OK: $OUT/BallisticsEngine.xcframework (features: $FEATURES, engine $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2))"
echo "    simulator archs: $(lipo -archs "$SIM_FAT")"
