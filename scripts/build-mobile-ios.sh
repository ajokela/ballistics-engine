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

# arm64 only. The simulator used to carry an x86_64 slice as well, so an Intel
# Mac could run it — but Xcode 26 is arm64-only (its own xcodebuild, clang and
# swiftc have no x86_64 slice), so no Intel Mac can host a simulator on a
# current toolchain. Building it cost a third of this script's runtime and 64 MB
# of the shipped xcframework for something nothing could link against.
rustup target add aarch64-apple-ios aarch64-apple-ios-sim >/dev/null

for target in aarch64-apple-ios aarch64-apple-ios-sim; do
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

# Kept as a lipo step rather than a plain copy: -create with one input still
# produces a valid archive, and this is where a second architecture would go
# back if one were ever needed again.
SIM_FAT="$OUT/libballistics_engine_sim.a"
lipo -create "target/aarch64-apple-ios-sim/release/libballistics_engine.a" -output "$SIM_FAT"

rm -rf "$OUT/BallisticsEngine.xcframework"
xcodebuild -create-xcframework \
  -library "target/aarch64-apple-ios/release/libballistics_engine.a" -headers "$HEADER_DIR" \
  -library "$SIM_FAT" -headers "$HEADER_DIR" \
  -output "$OUT/BallisticsEngine.xcframework"

echo "OK: $OUT/BallisticsEngine.xcframework (features: $FEATURES, engine $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2))"
echo "    simulator archs: $(lipo -archs "$SIM_FAT")"
