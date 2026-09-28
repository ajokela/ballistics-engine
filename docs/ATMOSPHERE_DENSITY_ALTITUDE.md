# Density altitude over the bridge: `atmosphere.density_altitude`

One bridge command that answers "how thin is this air?" from the engine itself, so an app never
has to carry its own copy of a density-altitude formula. It takes the same `atmosphere` object a
[solve-json v1](SOLVE_JSON_V1.md#atmosphere) request does, decodes and resolves it with the
solve's own code, and reports:

- **pressure altitude**,
- **density altitude, under two names** — `faa_rule` and `density_matched` — because the field
  uses the phrase for two different quantities (see [Which one to show](#which-one-to-show)),
- **air density** in kg/m³, the CIPM-2007 humid-air density the solver itself uses at the muzzle.

It runs no trajectory. It is present in every build that has the bridge — mobile, desktop and
wasm32 — and is listed by `meta.capabilities` directly after `solve`.

## Request

```json
{
  "api_version": 1,
  "command": "atmosphere.density_altitude",
  "request": {
    "atmosphere": {
      "altitude_m": 0.0,
      "temperature_k": 288.15,
      "pressure_pa": 101325.0,
      "relative_humidity": 0.5
    }
  }
}
```

`request.atmosphere` is exactly solve-json v1's `atmosphere`, field for field. The units are SI
and are not negotiable: there is no `temperature_c`, and humidity is a fraction.

| Field | Default when omitted | Meaning |
| --- | --- | --- |
| `altitude_m` | `0` | Station altitude, meters. |
| `temperature_k` | ICAO standard at `altitude_m` | Station temperature, **kelvin** (°C + 273.15). |
| `pressure_pa` | ICAO standard at `altitude_m` | Pressure, **pascals** (hPa × 100). Station pressure unless `pressure_reference` says otherwise. |
| `pressure_reference` | `"absolute"` | `"absolute"`: `pressure_pa` is station pressure. `"qnh"`: it is a sea-level-corrected reading, reduced to station pressure at `altitude_m` first. |
| `relative_humidity` | `0.5` | **Fraction** 0–1. 50% is `0.5`; `50.0` is refused. |
| `latitude_rad` | absent | Validated as the solve validates it (radians, ±π/2 — a value in degrees is refused), echoed back, otherwise unused. |

Two practical rules:

- **Know which pressure you have.** It decides whether every figure below is right or off by
  about the station's elevation. In order of preference:
  - A phone's pressure sensor reads **station pressure**: send it as-is. This is the best
    source there is.
  - A weather meter can show either kind. Its *station pressure* reading is absolute. Its
    *barometric* ("baro") reading is corrected to sea level using the meter's reference
    altitude, so send it as `qnh` with the real `altitude_m` — unless that reference altitude
    is set to 0, in which case it is station pressure.
  - A weather service's *surface* or *station* pressure field is absolute: send it as-is, as
    long as the service's elevation for the point is close to the shooter's.
  - An altimeter setting — a METAR, ATIS or airport figure, "A3002" — is a QNH: send it as `qnh`
    with the real `altitude_m`.
  - A weather API's or widget's *sea-level pressure* (MSLP) is **not** a QNH, though it looks
    like one. It was reduced to sea level using the actual temperature, and `qnh` undoes it
    using the standard atmosphere, so at elevation the two disagree by roughly 1% of pressure
    per 1000 m per 20 °C away from standard. Sent as `qnh` it is approximate; sent as absolute
    it is wrong by the station's elevation. Prefer any of the sources above.
- **Send what you know.** Every omitted field is filled with an ICAO standard value and announced
  in `assumptions`, exactly as a solve would announce it. With an explicit temperature and an
  absolute pressure, `altitude_m` does not move any figure; it matters for a QNH reduction and
  for the defaults — but send it anyway, so the echoed atmosphere describes the real station.

## Response

The envelope's `result` for the request above — engine output, not hand-written:

```json
{
  "air_density_kg_m3": 1.2216312037539498,
  "assumptions": [],
  "atmosphere": {
    "altitude_m": 0.0,
    "pressure_pa": 101325.0,
    "relative_humidity": 0.5,
    "temperature_k": 288.15
  },
  "density_altitude": {
    "density_matched": { "ft": 108.61100435743677, "m": 33.10463412814673 },
    "faa_rule": { "ft": 0.0, "m": 0.0 }
  },
  "pressure_altitude": { "ft": 0.0, "m": 0.0 }
}
```

| Field | Meaning |
| --- | --- |
| `atmosphere` | The atmosphere these figures describe, AFTER resolution: defaults filled in, and `pressure_pa` already reduced to station pressure if a QNH was sent. |
| `pressure_altitude` | NWS pressure altitude of the station pressure. |
| `density_altitude.faa_rule` | Pressure altitude plus the FAA's 120 ft/°C rule of thumb. No humidity term. |
| `density_altitude.density_matched` | The ISA altitude whose standard density equals this air's actual density, humidity included. |
| `air_density_kg_m3` | CIPM-2007 humid-air density — what the solver flies the bullet through. |
| `assumptions` | Every default applied, plus a `qnh_reduced_to_station_pressure` notice whenever a QNH was reduced, in the solve's notice format (`code`, `message`, `path`). Informational: show them in a details view, never refuse a result over them, and key on `code` rather than on the list being empty. |

Every altitude is an object carrying both `m` and `ft`, so neither side has to convert.

**`result.atmosphere` is a report, not a request.** When a QNH was sent, it carries the reduced
station pressure *and* still echoes `"pressure_reference": "qnh"` — the same convention as the
solve's `resolved_request`. Sent back as-is, to this command or to `solve`, that pressure would
be reduced a second time. To reuse it, drop `pressure_reference`, or re-send the original.

## Which one to show

There is deliberately no field called just `density_altitude`. The two figures answer different
questions, and an app has to choose one by name.

**`faa_rule`** is the pilot's rule-of-thumb density altitude: NWS pressure altitude corrected by
120 ft for every °C the air is away from standard. It ignores humidity. It is the formula the
engine's DOPE card header uses, and it is what the engine's density-altitude *input* (the
CLI's `--density-altitude` and `DA` location column, and the C ABI's
`ballistics_density_altitude_*` conversions) inverts exactly. It is **not**
the National Weather Service's density altitude: the NWS calculator includes humidity.

**`density_matched`** is the textbook definition: the altitude in the standard atmosphere where
the air is as thin as this air actually is. Humidity is in it, because moist air is lighter than
dry air. It is what humidity-aware tools and weather meters mean by density altitude, the NWS's
own calculator included — those tools differ among themselves by a few meters depending on their
method (the NWS one reads about 5 m at ICAO standard dry air because its constants are rounded),
and this one is an exact match against the solver's own density.

**Entering a density altitude is a different matter.** The engine's density-altitude input
inverts `faa_rule` exactly, then applies the run's humidity on top. Any other figure —
`density_matched`, or a humidity-aware meter's reading — is off by the whole gap in the table
below: the humidity already inside it, counted again by the run's humidity, plus the FAA rule's
own approximation, which is there even in dry air. Warm, humid air comes out too thin (0.6% at
30 °C and 50% RH, 1.5% at 35 °C and 90%); cold or dry air comes out too dense (0.9% at −20 °C).
And enter the real temperature with it. A density altitude entered alone rebuilds the air at the
standard temperature for that altitude, which is different air, not merely a different pressure:
at 30 °C and 1013.25 hPa it comes back as 11.4 °C and 949 hPa, 0.3% denser, with a speed of sound
3% lower — and the trajectory moves with it.

The two figures coincide at ICAO standard sea-level dry air; anywhere else, a match is
coincidence. **Their difference is not a humidity correction**:

Each row sent with `altitude_m` 0 and an absolute pressure (so the altitude moves no figure); altitudes in meters; engine output.

| Temperature | Pressure | RH | Pressure altitude | `faa_rule` | `density_matched` | Difference | Air density (kg/m³) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 15 °C | 1013.25 hPa | 0% | 0.0 | 0.0 | 0.0 | 0.0 | 1.2255 |
| 15 °C | 1013.25 hPa | 50% | 0.0 | 0.0 | 33.1 | +33.1 | 1.2216 |
| 15 °C | 1013.25 hPa | 100% | 0.0 | 0.0 | 66.2 | +66.2 | 1.2178 |
| 30 °C | 1013.25 hPa | 0% | 0.0 | 548.6 | 526.9 | −21.7 | 1.1647 |
| 30 °C | 1013.25 hPa | 50% | 0.0 | 548.6 | 608.5 | +59.9 | 1.1555 |
| −10 °C | 1013.25 hPa | 50% | 0.0 | −914.4 | −953.2 | −38.8 | 1.3417 |
| 15 °C | 700 hPa | 0% | 3010.9 | 3727.5 | 3690.9 | −36.6 | 0.8465 |
| 5 °C | 850 hPa | 30% | 1456.7 | 1437.6 | 1449.6 | +12.0 | 1.0638 |
| 25 °C | 900 hPa | 80% | 988.1 | 1589.0 | 1670.2 | +81.2 | 1.0407 |

Look at the dry rows. With no water in the air at all, the two still differ — by 21.7 m at 30 °C
at sea level, and by 36.6 m on a 15 °C day at 700 hPa. The FAA rule is a straight line where the
real dependence curves, and 120 ft/°C is a little steeper than the exact slope at sea level
(119.0 ft/°C for the CIPM density `density_matched` uses; 118.6 for an ideal gas). Its error
grows with how far the temperature is from standard *for the pressure altitude*, not from
15 °C, so even a 15 °C day at altitude is enough. An app that shows `density_matched −
faa_rule` as "humidity effect" is showing the wrong thing.

Whichever figure is displayed, the trajectory is the same: the solver uses the full humid-air
density (`air_density_kg_m3`) either way. In this command density altitude is a readout,
not an input.

## Errors

Errors follow the bridge's usual envelope, `ok: false` with an `error` object, and are shaped
exactly like `solve`'s:

- **`invalid_request`** — there was no `request` at all.
- **`command_failed`** — anything wrong with the request. `error.message` is a fixed summary;
  the reason is in `error.details`, which is the solve-json error envelope, locating the
  offending field by path. Whenever `solve` would reject the same `atmosphere` object, this
  command rejects it with identical `details`, so an app that already highlights fields from
  `solve` errors needs no new code. This command also refuses some air `solve` still accepts
  (see [Scope](#scope)) — in the same envelope, at the path of the field to blame, or at
  `$.atmosphere` when only the combination of fields is out of range. Field-highlighting code
  should fall back to the whole atmosphere for that path.

Humidity sent as a percentage (`"relative_humidity": 50.0`) comes back as:

```json
{
  "ok": false,
  "api_version": 1,
  "engine_version": "…",
  "error": {
    "code": "command_failed",
    "message": "atmosphere.density_altitude request rejected",
    "details": {
      "schema_version": 1,
      "status": "error",
      "error": {
        "code": "invalid_value",
        "message": "value must be finite and in the inclusive range [0, 1]",
        "path": "$.atmosphere.relative_humidity",
        "line": null,
        "column": null
      }
    }
  }
}
```

The other mistakes worth knowing, each by its `details.error` (messages abridged with …):

| Sent | `code` | `path` | Message |
| --- | --- | --- | --- |
| `"temperature_c": 30.0` | `unknown_field` | `$.atmosphere.temperature_c` | unknown field `temperature_c` |
| `"temperature_k": null` | `invalid_value` | `$.atmosphere.temperature_k` | expected a number |
| `"temperature_k": 15.0` (°C in the kelvin field) | `invalid_value` | `$.atmosphere.temperature_k` | 15 K is outside 173.15–373.15 K … temperature_k is kelvin (°C + 273.15) |
| `"pressure_pa": 1013.25` (hPa in the pascal field) | `invalid_value` | `$.atmosphere.pressure_pa` | a station pressure of 1013.25 Pa is a pressure altitude of 25861 m … pressure_pa is pascals (hPa × 100) |
| `"temperature_k": 320, "pressure_pa": 25000` (each plausible, together beyond the troposphere) | `invalid_value` | `$.atmosphere` | the faa_rule density altitude of these conditions is 13989 m … |

To leave a field to its default, omit it; `null` is refused.

## Calling it from Kotlin

Using the `BallisticsBridge.call` wrapper from [ANDROID_JNI_BRIDGE.md](ANDROID_JNI_BRIDGE.md):

```kotlin
import org.json.JSONObject

/**
 * Station conditions in the units a phone usually holds them in. [stationPressureHpa] is
 * station pressure; for a sea-level-corrected reading add "pressure_reference": "qnh" (and
 * then [altitudeM] is required, not optional).
 */
fun densityAltitude(
    temperatureC: Double,
    stationPressureHpa: Double,
    humidityPercent: Double,
    altitudeM: Double? = null,
): JSONObject {
    val atmosphere = JSONObject()
        .put("temperature_k", temperatureC + 273.15)
        .put("pressure_pa", stationPressureHpa * 100.0)
        .put("relative_humidity", humidityPercent / 100.0)
    altitudeM?.let { atmosphere.put("altitude_m", it) }
    val request = JSONObject()
        .put("api_version", 1)
        .put("command", "atmosphere.density_altitude")
        .put("request", JSONObject().put("atmosphere", atmosphere))

    val response = JSONObject(BallisticsBridge.call(request.toString()))
    if (!response.getBoolean("ok")) {
        // The reason, and the field to blame, are in details; message is only a summary.
        val error = response.getJSONObject("error")
        val cause = error.optJSONObject("details")?.optJSONObject("error")
        throw IllegalArgumentException(
            if (cause != null) "${cause.optString("path")}: ${cause.optString("message")}"
            else error.getString("message")
        )
    }
    return response.getJSONObject("result")
}

val result = densityAltitude(temperatureC = 15.0, stationPressureHpa = 1013.25, humidityPercent = 50.0, altitudeM = 0.0)
val matchedMeters = result.getJSONObject("density_altitude").getJSONObject("density_matched").getDouble("m")
```

Build the JSON with `JSONObject` (or any JSON library), never with `String.format`. Kotlin's
`String.format` without an explicit `Locale` uses the device's, and on a phone set to a
comma-decimal locale `"%.2f".format(288.15)` produces `288,15` — which is not a JSON number, so
the request is refused on exactly the phones whose owners are least likely to report it.
`JSONObject` writes numbers the same way on every device. (It refuses a NaN or infinite
`Double` outright, so a sensor that has not produced a reading yet fails in your code, not
here.)

Feature-detect before calling: an engine built before this command answers `unknown_command`.
`meta.capabilities` lists `atmosphere.density_altitude` when it is there.

## Scope

- Standard-atmosphere arithmetic is the closed-form troposphere with geometric altitude; there is
  no geopotential correction. That covers anywhere a rifle is fired, and the command refuses
  anything outside it rather than extrapolating: a pressure altitude or either density altitude
  below −5,000 m or above 11,000 m is `invalid_value`.
- It also refuses air that cannot exist, which is also what unit mistakes produce: a temperature
  outside 173.15–373.15 K (−100 °C to 100 °C, where the engine's humid-air formulas are
  defined), and a relative humidity that would need more water-vapor pressure than the air's
  total pressure. `solve` does not yet make these checks; here they are what stands between a
  unit mistake and a confident wrong number.
- Each refusal is located at the field to blame — `temperature_k`, `pressure_pa` (or
  `altitude_m`, when the pressure was defaulted or a plausible QNH was reduced at an impossible
  altitude), `relative_humidity` — except a density altitude out of range with every field
  individually plausible, which is located at `$.atmosphere`.
- `density_matched` is anchored to the solver's own dry-air density at 15 °C and 1013.25 hPa, so
  ICAO standard dry air reads exactly 0 m rather than a few meters off from a rounding difference
  between two density formulas.

Implementation: `crate::atmosphere_service` (decoding and the service), `crate::atmosphere`
(`nws_pressure_altitude_ft`, `faa_rule_density_altitude_ft`, `density_matched_altitude_m`), and
`src/bridge/mod.rs` (the command).
