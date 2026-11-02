/* ballistics-engine JSON command bridge — C ABI.
 *
 * One request/response exchange per call. The request is a UTF-8 JSON envelope:
 *
 *   { "api_version": 1, "command": "meta.capabilities" }
 *   { "api_version": 1, "command": "solve", "request": { ...solve-json v1... } }
 *
 * The response is always a UTF-8 JSON envelope. On success:
 *
 *   { "ok": true, "api_version": 1, "engine_version": "...",
 *     "command": "...", "result": { ... } }
 *
 * On ANY failure (invalid UTF-8, malformed JSON, unknown command, command
 * failure, internal error) the failure is reported IN-BAND:
 *
 *   { "ok": false, "api_version": 1, "engine_version": "...",
 *     "error": { "code": "...", "message": "...", "details": { ... } } }
 *
 * Branch on "ok" and "error"."code" rather than on which fields are present:
 * a request refused before it reaches the bridge — one whose bytes are not
 * UTF-8, or a NULL pointer — is answered by this ABI layer with the minimal
 * shape, "ok"/"api_version"/"error", and carries no "engine_version".
 *
 * "details" appears only when "code" is "command_failed", and not always then.
 * Its shape depends on the command:
 *
 *   - "solve" and "atmosphere.density_altitude": the command's whole solve-json
 *     error envelope, so the code to branch on is error.details.error.code:
 *
 *       "details": { "schema_version": 1, "status": "error",
 *                    "error": { "code": "invalid_value", "message": "...",
 *                               "path": "$.atmosphere.latitude_rad",
 *                               "line": null, "column": null } }
 *
 *     "path" is a string, or JSON null when no one field is at fault (as with
 *     "solve_failed"). "line" and "column" are null: malformed JSON fails the
 *     envelope itself first, as "invalid_json". The outer "message" is
 *     "solve request rejected" when the solve-json decoder refused the
 *     document, "solve failed" when the solve service refused or could not
 *     finish it, and "atmosphere.density_altitude request rejected" for every
 *     refusal of that command.
 *   - "reticle.describe", "reticle.hold", "reticle.holds", "reticle.import",
 *     "true.dsf", "true.plan", "true.dial_plan" and "profile.export_a7p": a
 *     flat object with a machine-readable "reason", plus the offending values
 *     when the error has any.
 *   - Every other command: no "details"; the reason is in "message" alone.
 *
 * Messages are text for a person; branch on codes. Per command, with a Kotlin
 * unwrap: docs/ANDROID_JNI_BRIDGE.md, "Reading a failure".
 *
 * Relative humidity, per command. One rule: a field named "relative_humidity"
 * (the SI wire of solve-json v1, beside temperature_k and pressure_pa) is a
 * FRACTION, 0 to 1 (0.5 = 50 %); every other humidity is a PERCENT, 0 to 100
 * (50 = 50 %). No field takes the other unit, and none guesses the unit from
 * the size of the number: 50 sent as a fraction is refused, and 0.5 sent as a
 * percent is 0.5 %. A reading up to 5 points outside the range (0.05 as a
 * fraction) is the saturated or dry air it reads and is used exactly as 100 %
 * or 0 %, so a sensor that reads 102 % in fog never fails a call. A value
 * farther out, or not finite, is refused, never clamped, and so is air that
 * cannot exist (a humidity needing more water-vapour pressure than the air's
 * total pressure), by one rule on every command. "if omitted" is what an
 * absent field means; "required" means its absence is refused.
 *
 *   command                      humidity field                unit      if omitted
 *   solve                        atmosphere.relative_humidity  fraction  0.5
 *   atmosphere.density_altitude  atmosphere.relative_humidity  fraction  0.5
 *   card.come_ups                humidity                      percent   50
 *   card.range_table             humidity                      percent   50
 *   card.wind                    humidity                      percent   50
 *   card.pdf                     humidity                      percent   50
 *   profile.validate             humidity                      percent   50
 *   profile.normalize            humidity                      percent   50
 *   profile.import_a7p           a7p:air_humidity_pct          percent   0
 *   profile.export_a7p           humidity                      percent   50
 *   true.fit                     model.humidity_pct            percent   required
 *   true.dsf                     model.humidity_pct            percent   required
 *   true.plan                    model.humidity_pct            percent   required
 *   true.wind                    humidity_pct                  percent   required
 *   true.tall_target             -                             none      -
 *   true.dial_plan               -                             none      -
 *   meta.capabilities            -                             none      -
 *   meta.version                 -                             none      -
 *   reticle.catalog              -                             none      -
 *   reticle.describe             -                             none      -
 *   reticle.hold                 -                             none      -
 *   reticle.holds                -                             none      -
 *   reticle.import               -                             none      -
 *   bc5d.info                    -                             none      -
 *
 * The profile.* humidity is the stored profile's own "humidity" key.
 * profile.validate reports one outside the band as a warning (valid: false)
 * rather than an error; profile.normalize echoes it unchanged, since it
 * re-serializes the document and reads nothing. profile.import_a7p reads the
 * file's whole-percent air humidity field: an absent one is 0 %, the format's
 * own meaning, because the file writes 0 by omitting the field.
 *
 * Contract:
 *   - The calls NEVER return NULL and NEVER throw/abort; check "ok" in the JSON.
 *   - The caller owns the input buffer; the engine does not retain it.
 *   - Every returned pointer must be released exactly once with
 *     ballistics_bridge_free(). Freeing NULL is a no-op.
 *   - Calls are thread-safe and independent; there is no shared mutable state.
 *   - Requests larger than 1 MiB are rejected with code "resource_limit".
 *   - Feature-detect with the "meta.capabilities" command instead of assuming a
 *     command list; builds differ (e.g. "pdf", "profile-import").
 *
 * ANDROID / JNI: call ballistics_bridge_call_n with the bytes of a Kotlin
 * String.toByteArray(Charsets.UTF_8), and return the response as a jbyteArray.
 * JNI's own string conversions — GetStringUTFChars AND GetStringUTFRegion —
 * produce MODIFIED UTF-8, which this ABI refuses as invalid UTF-8 the first
 * time a user types a character above U+FFFF, and NewStringUTF expects modified
 * UTF-8 on the way back. Worked Kotlin + C example, and why:
 * docs/ANDROID_JNI_BRIDGE.md.
 */

#ifndef BALLISTICS_BRIDGE_H
#define BALLISTICS_BRIDGE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Process one bridge request (NUL-terminated UTF-8 JSON envelope). */
char *ballistics_bridge_call(const char *request_json);

/* Length-explicit variant for buffers that are not NUL-terminated. */
char *ballistics_bridge_call_n(const uint8_t *request, size_t len);

/* Release a response returned by either call. NULL is a no-op. */
void ballistics_bridge_free(char *response);

#ifdef __cplusplus
}
#endif

#endif /* BALLISTICS_BRIDGE_H */
