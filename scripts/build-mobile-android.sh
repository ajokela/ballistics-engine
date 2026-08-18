#!/usr/bin/env bash
# Build libballistics_engine.so for the Ballistics Insight Android app.
#
# Uses cargo-ndk (install: `cargo install cargo-ndk`; needs ANDROID_NDK_HOME).
# ABIs: arm64-v8a + x86_64 (armeabi-v7a/x86 dropped — no measured user base).
# 16 KB page alignment is required for Google Play (Android 15+ devices).
# Mobile feature set matches iOS: bridge + ffi + pdf + profile-import, no online/cli.
#
# Output: target/mobile/jniLibs/<abi>/libballistics_engine.so
set -euo pipefail

ENGINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ENGINE_DIR"

FEATURES="bridge,ffi,pdf,profile-import"
OUT="$ENGINE_DIR/target/mobile/jniLibs"

command -v cargo-ndk >/dev/null || { echo "cargo-ndk not installed: cargo install cargo-ndk"; exit 1; }
: "${ANDROID_NDK_HOME:?set ANDROID_NDK_HOME to your NDK path}"

rustup target add aarch64-linux-android x86_64-linux-android >/dev/null

rm -rf "$OUT"
RUSTFLAGS="${RUSTFLAGS:-} -C link-arg=-Wl,-z,max-page-size=16384" \
  cargo ndk -t arm64-v8a -t x86_64 -o "$OUT" \
  build --release --no-default-features --features "$FEATURES"

echo "OK: $OUT (features: $FEATURES, engine $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2))"
find "$OUT" -name '*.so' -exec ls -la {} \;
