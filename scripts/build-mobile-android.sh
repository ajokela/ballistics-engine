#!/usr/bin/env bash
# Build libballistics_engine.so for the Ballistics Insight Android app.
#
# Uses cargo-ndk (install: `cargo install cargo-ndk`; needs ANDROID_NDK_HOME).
# ABIs: arm64-v8a + x86_64 (armeabi-v7a/x86 dropped — no measured user base).
# 16 KB page alignment is required for Google Play (Android 15+ devices).
# A DT_SONAME is required for a consumer to link the result portably (MBA-1541); this script
# sets one and then reads it back off the artifacts before declaring success.
# Mobile feature set matches iOS: bridge + ffi + pdf + profile-import, no online/cli.
#
# Output: target/mobile/jniLibs/<abi>/libballistics_engine.so
set -euo pipefail

ENGINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ENGINE_DIR"

FEATURES="bridge,ffi,pdf,profile-import"
OUT="$ENGINE_DIR/target/mobile/jniLibs"
SONAME="libballistics_engine.so"

command -v cargo-ndk >/dev/null || { echo "cargo-ndk not installed: cargo install cargo-ndk"; exit 1; }
: "${ANDROID_NDK_HOME:?set ANDROID_NDK_HOME to your NDK path}"

rustup target add aarch64-linux-android x86_64-linux-android >/dev/null

rm -rf "$OUT"
# MBA-1541: -soname is not optional for these artifacts. Without DT_SONAME the linker
# records whatever path the consumer linked against, so a CMake project using the ordinary
# `add_library(... SHARED IMPORTED)` + `IMPORTED_LOCATION /abs/path/libballistics_engine.so`
# pattern bakes that absolute HOST path into its own DT_NEEDED, and the app then fails to
# dlopen the library on device. An external integrator hit exactly that and had to fall back
# to plain -L/-l.
#
# It goes after `--`, on `rustc`, rather than into RUSTFLAGS next to max-page-size, because
# RUSTFLAGS reaches EVERY crate in the graph and a dependency here also declares a cdylib
# (printpdf), which cargo-ndk copies into the same jniLibs directory. A soname in RUSTFLAGS
# stamps `libballistics_engine.so` onto that file too, leaving two different libraries in one
# directory claiming one name. `cargo rustc`'s trailing arguments apply to this crate alone.
# max-page-size stays in RUSTFLAGS on purpose: Play's 16 KB requirement is about every .so
# shipped in the APK, not just ours.
RUSTFLAGS="${RUSTFLAGS:-} -C link-arg=-Wl,-z,max-page-size=16384" \
  cargo ndk -t arm64-v8a -t x86_64 -o "$OUT" \
  rustc --lib --release --no-default-features --features "$FEATURES" \
  -- -C link-arg=-Wl,-soname,"$SONAME"

# Prefer the NDK's own readelf: this script already requires ANDROID_NDK_HOME, so it is the
# one ELF reader we know is installed. macOS /usr/bin has none.
READELF=""
for candidate in "$ANDROID_NDK_HOME"/toolchains/llvm/prebuilt/*/bin/llvm-readelf; do
  [ -x "$candidate" ] && { READELF="$candidate"; break; }
done
[ -n "$READELF" ] || READELF="$(command -v llvm-readelf || command -v readelf || true)"
[ -n "$READELF" ] || { echo "no llvm-readelf/readelf found; cannot verify DT_SONAME" >&2; exit 1; }

# Two checks, because the fix has two ways to go wrong: the soname can go missing, and it can
# land on a file that is not ours (see the RUSTFLAGS note above). Both fail the build here
# rather than reaching a device.
checked=0
for so in "$OUT"/*/*.so; do
  [ -e "$so" ] || continue
  claimed="$("$READELF" -d "$so" | sed -n 's/.*SONAME.*\[\(.*\)\].*/\1/p')"
  if [ "$(basename "$so")" = "$SONAME" ]; then
    checked=$((checked + 1))
    if [ "$claimed" = "$SONAME" ]; then
      echo "    DT_SONAME ok: $so"
    else
      echo "FAIL: $so declares DT_SONAME '${claimed:-<none>}', expected '$SONAME'." >&2
      echo "      A consumer linking it by absolute path would bake that path into its own" >&2
      echo "      DT_NEEDED and fail to dlopen on device. Check the -soname arg above." >&2
      exit 1
    fi
  elif [ "$claimed" = "$SONAME" ]; then
    echo "FAIL: $so also claims DT_SONAME '$SONAME'." >&2
    echo "      Two libraries in one jniLibs directory answering to one name. The -soname" >&2
    echo "      arg has leaked out of this crate — check it is on 'cargo rustc --', not in" >&2
    echo "      RUSTFLAGS." >&2
    exit 1
  fi
done
[ "$checked" -gt 0 ] || { echo "FAIL: no $SONAME produced under $OUT" >&2; exit 1; }

echo "OK: $OUT (features: $FEATURES, engine $(grep -m1 '^version' Cargo.toml | cut -d'"' -f2))"
find "$OUT" -name '*.so' -exec ls -la {} \;
