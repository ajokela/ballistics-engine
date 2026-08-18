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
 * Contract:
 *   - The calls NEVER return NULL and NEVER throw/abort; check "ok" in the JSON.
 *   - The caller owns the input buffer; the engine does not retain it.
 *   - Every returned pointer must be released exactly once with
 *     ballistics_bridge_free(). Freeing NULL is a no-op.
 *   - Calls are thread-safe and independent; there is no shared mutable state.
 *   - Requests larger than 1 MiB are rejected with code "resource_limit".
 *   - Feature-detect with the "meta.capabilities" command instead of assuming a
 *     command list; builds differ (e.g. "pdf", "profile-import").
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
