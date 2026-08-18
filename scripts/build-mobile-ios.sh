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

rustup target add aarch64-apple-ios aarch64-apple-ios-sim x86_64-apple-ios >/dev/null

for target in aarch64-apple-ios aarch64-apple-ios-sim x86_64-apple-ios; do
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

# Fat simulator lib (arm64 + x86_64)
SIM_FAT="$OUT/libballistics_engine_sim.a"
lipo -create \
  "target/aarch64-apple-ios-sim/release/libballistics_engine.a" \
  "target/x86_64-apple-ios/release/libballistics_engine.a" \
  -output "$SIM_FAT"

rm -rf "$OUT/BallisticsEngine.xcframework"
xcodebuild -create-xcframework \
  -library "target/aarch64-apple-ios/release/libballistics_engine.a" -headers "$HEADER_DIR" \
  -library "$SIM_FAT" -headers "$HEADER_DIR" \
  -output "$OUT/BallisticsEngine.xcframework"

echo "OK: $OUT/BallisticsEngine.xcframework (features: $FEATURES, engine $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2))"
