# Calling the JSON bridge from Android (JNI)

The engine's mobile-facing surface is the JSON command bridge: three C symbols declared in
[`include/ballistics_bridge.h`](../include/ballistics_bridge.h), one request/response exchange
per call.

```c
char *ballistics_bridge_call(const char *request_json);            /* NUL-terminated */
char *ballistics_bridge_call_n(const uint8_t *request, size_t len); /* length-explicit */
void  ballistics_bridge_free(char *response);
```

On Android, use **`ballistics_bridge_call_n`** over a `ByteArray`, and return a `jbyteArray`.
Not `ballistics_bridge_call` over a `jstring`. The obvious `jstring` path compiles, passes an
ASCII smoke test, and ships a defect that surfaces later on a user's own text.

## Why a jstring is the wrong input

JNI does not hand out UTF-8. `GetStringUTFChars` and `GetStringUTFRegion` both produce
**modified UTF-8**, which is UTF-8 applied to UTF-16 code units rather than to Unicode scalar
values. Two things come out differently from standard UTF-8:

| Input | UTF-8 | Modified UTF-8 |
| --- | --- | --- |
| U+0000 | `00` | `C0 80` (so the byte stream never contains a NUL) |
| U+1F3AF 🎯 (any character above U+FFFF) | `F0 9F 8E AF` | `ED A0 BC ED BE AF` — its two surrogates, each in the 3-byte form |

Everything from U+0001 to U+FFFF encodes identically, which is what makes this a bad bug rather
than an obvious one. `"München"`, `"Ω"` and `"中文"` all survive the wrong path untouched. An
emoji in a profile name does not.

The engine decodes the request buffer with strict UTF-8 and refuses anything else, rather than
guessing at a repair. A modified-UTF-8 request therefore comes back as an ordinary in-band
failure envelope:

```json
{"ok":false,"api_version":1,"error":{"code":"invalid_json","message":"request is not valid UTF-8"}}
```

(This one is raised by the C ABI layer itself, before the request reaches the bridge, and its
envelope is the minimal shape: `ok`, `api_version`, `error`. Envelopes the bridge produces also
carry `engine_version`. Branch on `ok` and `error.code`, not on which fields are present.)

That refusal is deliberate. Accepting modified UTF-8 means accepting an overlong sequence
(`C0 80`) and an encoded surrogate, which is exactly what UTF-8 validation is for. The fix
belongs on the calling side.

**`GetStringUTFRegion` is not the fix.** The JNI specification has it write modified UTF-8 as
well; it differs from `GetStringUTFChars` in who owns the buffer, not in the encoding. Encoding
the string in Kotlin with `toByteArray(Charsets.UTF_8)` and passing the bytes is the short way
round. (`GetStringChars` plus your own UTF-16 → UTF-8 conversion in C is the long way round, and
gets you to the same place.)

What the engine does with each encoding, and where the two diverge, is driven by
`a_request_in_jni_modified_utf8_is_refused` and
`modified_utf8_and_utf8_agree_below_the_supplementary_planes` in
[`src/bridge/ffi.rs`](../src/bridge/ffi.rs) — so those two claims above fail a test if they stop
being true. What JNI does is the JNI specification's, not ours.

## The response has the same trap in reverse

`NewStringUTF` **expects** modified UTF-8. The bridge returns standard UTF-8 JSON, which can
contain a 4-byte sequence whenever the response echoes text a user supplied. Handing that to
`NewStringUTF` is the mirror image of the input bug.

Return a `jbyteArray` and decode it on the Kotlin side with `Charsets.UTF_8`.

## Worked example

Kotlin:

```kotlin
package com.example.ballistics

object BallisticsBridge {
    init { System.loadLibrary("ballistics_engine") }

    /** Bytes in, bytes out. Neither side ever sees a jstring. */
    private external fun nativeCall(request: ByteArray): ByteArray

    fun call(requestJson: String): String =
        String(nativeCall(requestJson.toByteArray(Charsets.UTF_8)), Charsets.UTF_8)
}
```

C:

```c
#include <jni.h>
#include <stdint.h>
#include <string.h>

#include "ballistics_bridge.h"

JNIEXPORT jbyteArray JNICALL
Java_com_example_ballistics_BallisticsBridge_nativeCall(JNIEnv *env,
                                                        jobject thiz,
                                                        jbyteArray request)
{
    (void)thiz;

    jsize request_len = (*env)->GetArrayLength(env, request);
    jbyte *body = (*env)->GetByteArrayElements(env, request, NULL);
    if (body == NULL) {
        return NULL; /* OutOfMemoryError is already pending */
    }

    char *response = ballistics_bridge_call_n((const uint8_t *)body, (size_t)request_len);

    /* JNI_ABORT: the engine does not write to the request buffer, so there is nothing to
     * copy back. */
    (*env)->ReleaseByteArrayElements(env, request, body, JNI_ABORT);

    /* The call never returns NULL, and never throws or aborts. Every failure — bad UTF-8,
     * bad JSON, unknown command, a command that failed, an internal panic — arrives as an
     * {"ok":false} envelope, so there is no error branch to write here. */
    jsize response_len = (jsize)strlen(response);
    jbyteArray out = (*env)->NewByteArray(env, response_len);
    if (out != NULL) {
        (*env)->SetByteArrayRegion(env, out, 0, response_len, (const jbyte *)response);
    }

    /* Exactly once, on every path, including the one where NewByteArray failed. */
    ballistics_bridge_free(response);
    return out;
}
```

In C++ the `JNIEnv` calls lose their first argument: `env->GetArrayLength(request)`. If you
would rather not spell out the mangled name, register the method with `RegisterNatives` from
`JNI_OnLoad` instead; nothing above changes.

### Things the example is relying on

- **`strlen` is safe on the response.** It is a NUL-terminated C string, and the engine builds
  it so that it has no interior NUL.
- **Free exactly once.** `ballistics_bridge_free` releases a pointer from either call; freeing
  NULL is a no-op. Nothing else may free it.
- **Calls are independent and thread-safe.** There is no shared mutable state, so a request may
  go on any thread and several may be in flight at once. Keep them off the main thread: a solve
  is real work.
- **Requests are capped at 1 MiB** and a larger one is refused with error code
  `resource_limit`.
- **Feature-detect with `meta.capabilities`** rather than assuming a command list. Mobile builds
  are `--no-default-features` with a chosen feature set, and what they answer to differs.

## Loading the library

`scripts/build-mobile-android.sh` writes `target/mobile/jniLibs/<abi>/libballistics_engine.so`
for `arm64-v8a` and `x86_64`. Copying those into the module's `jniLibs` source set is enough for
`System.loadLibrary("ballistics_engine")`.

Copy the engine library by name rather than the whole directory: cargo-ndk also drops any
`cdylib` a dependency happens to declare next to it (currently a `libprintpdf-<hash>.so`),
which nothing here loads and which is only APK weight.

If you link against the library from CMake instead, the ordinary imported-target pattern works:

```cmake
add_library(ballistics_engine SHARED IMPORTED)
set_target_properties(ballistics_engine PROPERTIES
    IMPORTED_LOCATION "${JNI_LIBS}/${ANDROID_ABI}/libballistics_engine.so"
    INTERFACE_INCLUDE_DIRECTORIES "${ENGINE_INCLUDE_DIR}")
target_link_libraries(your_jni_lib PRIVATE ballistics_engine)
```

This depends on the library declaring a `DT_SONAME`; without one, the linker records the path
you linked against — here an absolute host path — in your own `DT_NEEDED`, and the app fails to
`dlopen` it on device. Libraries built by that script carry
`DT_SONAME libballistics_engine.so`, and the script fails the build if they do not (MBA-1541).
If you are linking a library from elsewhere, `llvm-readelf -d libballistics_engine.so | grep
SONAME` says whether it has one.

## Sending a BC ladder with `solve`

A `solve` request can carry a velocity-banded BC ladder as `projectile.bc_segments`. Speeds are
in m/s, and the bands go fastest first without overlapping. The rules, and how the scalar
`ballistic_coefficient` still matters, are in
[`projectile.bc_segments`](SOLVE_JSON_V1.md#projectilebc_segments-mba-1598).

```kotlin
import org.json.JSONArray
import org.json.JSONObject

/** One band: [bc] flies while the bullet's airspeed is in [minMps, maxMps). */
data class BcBand(val minMps: Double, val maxMps: Double, val bc: Double)

/** Fastest band first. A ladder kept in fps converts with `fps * 0.3048`. */
fun bcSegments(bands: List<BcBand>): JSONArray =
    JSONArray(bands.sortedByDescending { it.minMps }.map { band ->
        JSONObject()
            .put("velocity_min_mps", band.minMps)
            .put("velocity_max_mps", band.maxMps)
            .put("ballistic_coefficient", band.bc)
    })

// projectile: the request's "projectile" JSONObject
projectile.put("bc_segments", bcSegments(listOf(
    BcBand(640.0, 1200.0, 0.462),
    BcBand(490.0, 640.0, 0.447),
    BcBand(0.0, 490.0, 0.424),
)))
```

A ladder the contract refuses comes back as an ordinary `{"ok":false}` envelope, and
`error.details.error.path` names the entry at fault, such as
`$.projectile.bc_segments[1].velocity_max_mps`. The `card.*` commands take a ladder too, but with
other member names, the request's `velocity_unit` and no order rule (they still refuse overlapping
bands, more than 64 bands and band BCs outside 0.0001 to 100), so a stored ladder has to be
written out for each command. The table at the end of that section lists the differences.

## Sight-in (zero-day) conditions with `solve`

A `solve` request has one atmosphere, so a zero made on another day takes two calls: one that
solves the zero in the sight-in air, and one that flies that angle in today's air. The rules,
including what the zero request should leave out, are in
[Sight-in (zero-day) conditions](SOLVE_JSON_V1.md#sight-in-zero-day-conditions).

```kotlin
import org.json.JSONObject

/**
 * [request] is a complete solve request for today's shot. [sightIn] is the atmosphere the rifle
 * was zeroed in: put "pressure_pa" in it only if a pressure was measured on that day. It replaces
 * the whole "atmosphere", so with effects.coriolis on it must carry the sight-in "latitude_rad"
 * too: Coriolis without a latitude is refused.
 */
fun solveWithSightIn(request: JSONObject, sightIn: JSONObject, zeroDistanceM: Double): JSONObject {
    // 1. The zero request: sight-in air, still air, level shot.
    val zero = JSONObject(request.toString())
    zero.put("atmosphere", sightIn)
    zero.put("wind", JSONObject())
    zero.getJSONObject("shot").apply {
        remove("muzzle_angle_rad")
        remove("shooting_angle_rad")
        remove("cant_angle_rad")
        put("zero_distance_m", zeroDistanceM)
    }
    val angle = call(zero).getJSONObject("resolved_request")
        .getJSONObject("shot").getDouble("muzzle_angle_rad")

    // 2. The shot request: today's air, wind and geometry, flying the pinned angle.
    val shot = JSONObject(request.toString())
    shot.getJSONObject("shot").apply {
        remove("zero_distance_m")
        put("muzzle_angle_rad", angle)
    }
    return call(shot)
}

/** One solve call; returns "result", or throws with the engine's error envelope. */
fun call(request: JSONObject): JSONObject {
    val envelope = JSONObject(BallisticsBridge.call(
        JSONObject().put("api_version", 1).put("command", "solve")
            .put("request", request).toString()))
    if (!envelope.getBoolean("ok")) throw IllegalStateException(envelope.toString())
    return envelope.getJSONObject("result")
}
```

Keep `zero_distance_m` in the second request only when it uses `rifle.sight_offset_lateral_m` or
`shot.zero_poi_right_m`, which need it.

## See also

- [`include/ballistics_bridge.h`](../include/ballistics_bridge.h) — the ABI contract:
  envelope shapes, ownership, and the guarantees the example leans on.
- [`docs/SOLVE_JSON_V1.md`](SOLVE_JSON_V1.md) — the request and result contract carried inside
  the envelope's `request` and `result`.
- [`docs/ATMOSPHERE_DENSITY_ALTITUDE.md`](ATMOSPHERE_DENSITY_ALTITUDE.md) — the
  `atmosphere.density_altitude` command: density altitude, pressure altitude and air density for
  an atmosphere, computed by the engine rather than by a copy of its formula, with a Kotlin
  example.
