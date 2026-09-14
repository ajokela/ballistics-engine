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
