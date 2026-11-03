//! MBA-1519: the card request is denominated PER DIMENSION, with `units` as the preset
//! that fills in whatever a request does not state for itself.
//!
//! Two things have to be true at once and this file pins both. A request written before the
//! per-dimension fields existed must produce the card it always produced, down to the byte —
//! `CardRequestV1` backs the paid DOPE card, the saved-card reprint and the PDF export, and a
//! stored request is replayed years after it was written. And a request that DOES name a
//! dimension must be denominating the same shot differently, never solving a different one.

use ballistics_engine::card_service::{come_ups_v1, range_table_v1, wind_card_v1, CardRequestV1};
use serde_json::{json, Value};

/// The exact bytes for the old-shape requests below. They were first captured from pre-MBA-1519
/// `main` (bfc8e87) and pasted here unaltered, so at first nothing in the per-dimension work
/// computed them. They have been re-captured once since, from the 0.50.0 build (below), so they
/// are now that build's output, no longer an independent pre-1519 oracle. A conversion that
/// moves them by one ulp still fails here.
///
/// Re-captured once since, in 0.50.0 (the unified aloft vapour model): the
/// cards fly the default 50 % RH, and the air off the station's height now holds the station's
/// water-vapour mole fraction with a virtual-temperature lapse. Only row values moved, every one
/// by at most 6e-7 relative; the largest, old -> new, are the 600 yd energy (imperial come-ups
/// 1097.800875568174 -> 1097.8015258301768 ft-lb, +5.9e-7; metric 1375.5021828188696 ->
/// 1375.5028324586233 J, +4.7e-7) and the 600 yd 10 mph wind column (imperial wind card
/// -1.4188703891384 -> -1.4188699544098804 MOA, -3.1e-7). The per-dimension fields resolve
/// through the same arithmetic as before, which the tests below prove on these new bytes.
const OLD_SHAPE_IMPERIAL_COME_UPS: &str = r#"{"schema_version":1,"kind":"come_ups","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_adj":-0.0002383115567345896,"come_up":-0.0002383115567345896,"velocity":2471.487048765613,"energy":2373.1265064596223,"time":0.11723436605171646},{"range":200.0,"drop_adj":0.5820871779613539,"come_up":0.5823254895180885,"velocity":2299.9539267097994,"energy":2055.14507487883,"time":0.2430640594744383},{"range":300.0,"drop_adj":1.386201591216661,"come_up":0.8041144132553071,"velocity":2135.7154928928917,"energy":1772.1113552947409,"time":0.378434154516867},{"range":400.0,"drop_adj":2.318304715233842,"come_up":0.9321031240171809,"velocity":1978.346157757052,"energy":1520.5782701358603,"time":0.5243858682936872},{"range":500.0,"drop_adj":3.371952280872122,"come_up":1.05364756563828,"velocity":1826.945972476667,"energy":1296.7481015646051,"time":0.6821843656291716},{"range":600.0,"drop_adj":4.557991163092783,"come_up":1.1860388822206613,"velocity":1680.9694047434455,"energy":1097.8015258301768,"time":0.8533797629410692}]}"#;

const OLD_SHAPE_IMPERIAL_RANGE_TABLE: &str = r#"{"schema_version":1,"kind":"range_table","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MOA"},"rows":[{"range":100.0,"drop_linear":-0.0008589363928091085,"drop_adj":-0.0002385934424469746,"wind_linear":-0.7087072147063299,"wind_adj":-0.676815390044545,"velocity":2471.4883185643516,"energy":2373.128944986517,"time":0.1172343314927059},{"range":200.0,"drop_linear":4.191018410500851,"drop_adj":0.5820858903473404,"wind_linear":-2.9301920467095806,"wind_adj":-1.3991667023038248,"velocity":2299.956009350782,"energy":2055.1487968066895,"time":0.24306390827015528},{"range":300.0,"drop_linear":14.97094343299353,"drop_adj":1.3861984660179194,"wind_linear":-6.830787563781233,"wind_adj":-2.174467374470359,"velocity":2135.717838332994,"energy":1772.1152475578579,"time":0.37843378363303637},{"range":400.0,"drop_linear":33.38349773580004,"drop_adj":2.318298453875003,"wind_linear":-12.593747978105789,"wind_adj":-3.006757329772757,"velocity":1978.347995929339,"energy":1520.5810958154054,"time":0.5243851536870159},{"range":500.0,"drop_linear":60.69494565889362,"drop_adj":3.3719414254940903,"wind_linear":-20.441742298755262,"wind_adj":-3.9043727790622547,"velocity":1826.9464915586927,"energy":1296.7488384432397,"time":0.6821831686069002},{"range":600.0,"drop_linear":98.45223836537404,"drop_adj":4.557973998396947,"wind_linear":-30.647591015253415,"wind_adj":-4.878074903261169,"velocity":1680.9677658066714,"energy":1097.799385129211,"time":0.8533779224134287}]}"#;

const OLD_SHAPE_IMPERIAL_WIND_CARD: &str = r#"{"schema_version":1,"kind":"wind_card","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"wind_speeds":[5.0,10.0],"wind_angles_deg":[90.0],"rows":[{"range":100.0,"wind_columns":[-0.09843092401894904,-0.19686311519620275]},{"range":200.0,"wind_columns":[-0.2034841727657323,-0.40697111759855287]},{"range":300.0,"wind_columns":[-0.3162378984861873,-0.6324803299797438]},{"range":400.0,"wind_columns":[-0.43727964063634106,-0.8745658318129019]},{"range":500.0,"wind_columns":[-0.567821785895692,-1.1356523499308477]},{"range":600.0,"wind_columns":[-0.7094293533839818,-1.4188699544098804]}]}"#;

const OLD_SHAPE_METRIC_COME_UPS: &str = r#"{"schema_version":1,"kind":"come_ups","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_adj":-0.00043294391023204873,"come_up":-0.00043294391023204873,"velocity":749.2214282012793,"energy":3182.7570435040066,"time":0.12855680088144575},{"range":200.0,"drop_adj":0.6846897349190779,"come_up":0.68512267882931,"velocity":693.1635696521198,"energy":2724.297698812993,"time":0.26732691722816576},{"range":300.0,"drop_adj":1.595074848156121,"come_up":0.9103851132370431,"velocity":639.6678285781318,"energy":2320.021957816611,"time":0.41751355410555474},{"range":400.0,"drop_adj":2.6509268767283105,"come_up":1.0558520285721895,"velocity":588.5998039514741,"energy":1964.3702575305435,"time":0.5804883582429626},{"range":500.0,"drop_adj":3.8525518172889583,"come_up":1.2016249405606478,"velocity":539.6244473704896,"energy":1651.0732360324052,"time":0.7579317727858298},{"range":600.0,"drop_adj":5.216641740421487,"come_up":1.3640899231325285,"velocity":492.53737650535805,"energy":1375.5028324586233,"time":0.9519021598335182}]}"#;

const OLD_SHAPE_METRIC_RANGE_TABLE: &str = r#"{"schema_version":1,"kind":"range_table","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_linear":-0.04332759417122845,"drop_adj":-0.00043327594171228454,"wind_linear":-21.381465662817384,"wind_adj":-0.21381465662817384,"velocity":749.2218385791115,"energy":3182.760530141095,"time":0.1285567597210484},{"range":200.0,"drop_linear":136.93764594308783,"drop_adj":0.6846882297154391,"wind_linear":-88.72285091936837,"wind_adj":-0.44361425459684184,"velocity":693.1642252183849,"energy":2724.3028518758824,"time":0.26732673673372326},{"range":300.0,"drop_linear":478.52131533125873,"drop_adj":1.5950710511041957,"wind_linear":-207.43857856392995,"wind_adj":-0.6914619285464332,"velocity":639.6685329345956,"energy":2320.027067104013,"time":0.41751310715369977},{"range":400.0,"drop_linear":1060.3677459787593,"drop_adj":2.650919364946898,"wind_linear":-383.70105887857613,"wind_adj":-0.9592526471964404,"velocity":588.6002789113701,"energy":1964.3734277613219,"time":0.5804874934385215},{"range":500.0,"drop_linear":1926.2693036639475,"drop_adj":3.852538607327895,"wind_linear":-625.0722860177899,"wind_adj":-1.2501445720355797,"velocity":539.6243969482913,"energy":1651.0729274864902,"time":0.7579303127080137},{"range":600.0,"drop_linear":3129.972369772147,"drop_adj":5.216620616286912,"wind_linear":-940.8148894286986,"wind_adj":-1.5680248157144978,"velocity":492.53648900592646,"energy":1375.4978754536571,"time":0.9518998980382989}]}"#;

const OLD_SHAPE_METRIC_WIND_CARD: &str = r#"{"schema_version":1,"kind":"wind_card","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"wind_speeds":[5.0,10.0],"wind_angles_deg":[90.0],"rows":[{"range":100.0,"wind_columns":[-0.23757232343576448,-0.47515989045655027]},{"range":200.0,"wind_columns":[-0.4929057857660132,-0.9858449939317048]},{"range":300.0,"wind_columns":[-0.768292779030205,-1.5366407332402738]},{"range":400.0,"wind_columns":[-1.0658388103164342,-2.1317576900487962]},{"range":500.0,"wind_columns":[-1.3890529493133676,-2.778214046650609]},{"range":600.0,"wind_columns":[-1.7422542163917805,-3.48464804299762]}]}"#;

/// The imperial request whose bytes the goldens above hold. Written the way a card request
/// was written before MBA-1519: one `units` scalar and no per-dimension field anywhere.
fn pre_1519_imperial() -> Value {
    json!({
        "units": "imperial",
        "muzzle_velocity": 2650.0,
        "ballistic_coefficient": 0.243,
        "mass": 175.0,
        "diameter": 0.308,
        "drag_model": "g7",
        "zero_distance": 100.0,
        "altitude": 1000.0,
        "temperature": 72.0,
        "pressure": 29.5,
        "wind_speed": 10.0,
        "wind_direction_deg": 90.0,
        "start": 100.0,
        "end": 600.0,
        "step": 100.0,
        "adjustment_unit": "mil",
        "windage_unit": "moa"
    })
}

/// The metric half of the same pin.
fn pre_1519_metric() -> Value {
    json!({
        "units": "metric",
        "muzzle_velocity": 807.72,
        "ballistic_coefficient": 0.243,
        "mass": 11.34,
        "diameter": 7.82,
        "drag_model": "g7",
        "zero_distance": 100.0,
        "altitude": 500.0,
        "wind_speed": 4.5,
        "wind_direction_deg": 90.0,
        "start": 100.0,
        "end": 600.0,
        "step": 100.0,
        "adjustment_unit": "mil"
    })
}

fn parse(v: &Value) -> CardRequestV1 {
    serde_json::from_value(v.clone()).expect("card request must parse")
}

fn come_ups_json(v: &Value) -> String {
    serde_json::to_string(&come_ups_v1(&parse(v)).expect("come-ups card")).unwrap()
}
fn range_table_json(v: &Value) -> String {
    serde_json::to_string(&range_table_v1(&parse(v)).expect("range-table card")).unwrap()
}
fn wind_card_json(v: &Value) -> String {
    serde_json::to_string(&wind_card_v1(&parse(v)).expect("wind card")).unwrap()
}

/// `wind_speeds` is what makes a request a wind card; the goldens were captured with these.
fn with_wind_sweep(mut v: Value) -> Value {
    v["wind_speeds"] = json!([5.0, 10.0]);
    v
}

fn rows(card: &str) -> Vec<Value> {
    let parsed: Value = serde_json::from_str(card).unwrap();
    parsed["rows"].as_array().unwrap().clone()
}

fn units_block(card: &str) -> Value {
    let parsed: Value = serde_json::from_str(card).unwrap();
    parsed["units"].clone()
}

/// Relative agreement, for the cross-unit relations where one side is a conversion of the
/// other. The tolerance is far tighter than any units defect could hide under: mistaking
/// millimetres for centimetres is a factor of ten, the closest pair of units on this request
/// (true MOA and SMOA) differ by 4.7%, and the closest pair of NUMBERS is 1013.25 hPa vs
/// 29.92 inHg at 3.4%.
#[track_caller]
fn assert_close(actual: f64, expected: f64, what: &str) {
    let tol = expected.abs().max(1.0) * 1e-12;
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: expected {expected}, got {actual}"
    );
}

// -----------------------------------------------------------------------------------------
// 1. A request that predates the per-dimension fields
// -----------------------------------------------------------------------------------------

/// The whole additive claim, tested rather than asserted in prose: the six documents a
/// pre-MBA-1519 request produced, byte for byte.
///
/// `CardRequestV1` is `deny_unknown_fields`, so the new fields could only ever be additive by
/// being optional — but optional is not the same as inert. Every one of them had to resolve
/// to the value the `units` scalar used to hard-code, through the same multiplications in the
/// same order, or a row moves in the last digits and a saved card reprints differently from
/// the screen it was saved off.
#[test]
fn an_old_shape_request_still_produces_its_pinned_bytes() {
    let imperial = pre_1519_imperial();
    assert_eq!(come_ups_json(&imperial), OLD_SHAPE_IMPERIAL_COME_UPS);
    assert_eq!(range_table_json(&imperial), OLD_SHAPE_IMPERIAL_RANGE_TABLE);
    assert_eq!(
        wind_card_json(&with_wind_sweep(imperial)),
        OLD_SHAPE_IMPERIAL_WIND_CARD
    );

    let metric = pre_1519_metric();
    assert_eq!(come_ups_json(&metric), OLD_SHAPE_METRIC_COME_UPS);
    assert_eq!(range_table_json(&metric), OLD_SHAPE_METRIC_RANGE_TABLE);
    assert_eq!(
        wind_card_json(&with_wind_sweep(metric)),
        OLD_SHAPE_METRIC_WIND_CARD
    );
}

/// The preset IS the per-dimension record, not a parallel mode: spelling out every dimension
/// the imperial preset implies produces the identical document, and so does the metric one.
///
/// This is the half the goldens cannot show. They prove the old shape did not move; this
/// proves the two shapes meet — that "apply a preset and it writes every per-dimension value"
/// is what the engine actually does, rather than the preset keeping a private path of its own.
#[test]
fn naming_every_dimension_explicitly_matches_the_preset_that_implies_it() {
    let mut imperial = pre_1519_imperial();
    let spelled_out = imperial.as_object_mut().unwrap();
    spelled_out.insert("distance_unit".into(), json!("yards"));
    spelled_out.insert("velocity_unit".into(), json!("fps"));
    spelled_out.insert("mass_unit".into(), json!("grains"));
    spelled_out.insert("diameter_unit".into(), json!("inches"));
    spelled_out.insert("sight_height_unit".into(), json!("inches"));
    spelled_out.insert("drop_unit".into(), json!("inches"));
    spelled_out.insert("wind_speed_unit".into(), json!("mph"));
    spelled_out.insert("temperature_unit".into(), json!("fahrenheit"));
    spelled_out.insert("pressure_unit".into(), json!("inhg"));
    spelled_out.insert("energy_unit".into(), json!("ftlb"));
    spelled_out.insert("altitude_unit".into(), json!("meters"));
    assert_eq!(come_ups_json(&imperial), OLD_SHAPE_IMPERIAL_COME_UPS);
    assert_eq!(range_table_json(&imperial), OLD_SHAPE_IMPERIAL_RANGE_TABLE);

    let mut metric = pre_1519_metric();
    let spelled_out = metric.as_object_mut().unwrap();
    spelled_out.insert("distance_unit".into(), json!("meters"));
    spelled_out.insert("velocity_unit".into(), json!("mps"));
    spelled_out.insert("mass_unit".into(), json!("grams"));
    spelled_out.insert("diameter_unit".into(), json!("mm"));
    spelled_out.insert("sight_height_unit".into(), json!("mm"));
    spelled_out.insert("drop_unit".into(), json!("mm"));
    spelled_out.insert("wind_speed_unit".into(), json!("mps"));
    spelled_out.insert("temperature_unit".into(), json!("celsius"));
    spelled_out.insert("pressure_unit".into(), json!("hpa"));
    spelled_out.insert("energy_unit".into(), json!("joules"));
    spelled_out.insert("altitude_unit".into(), json!("meters"));
    assert_eq!(come_ups_json(&metric), OLD_SHAPE_METRIC_COME_UPS);
    assert_eq!(range_table_json(&metric), OLD_SHAPE_METRIC_RANGE_TABLE);
}

// -----------------------------------------------------------------------------------------
// 2. A dimension that is named denominates, it does not re-solve
// -----------------------------------------------------------------------------------------

/// Output dimensions are presentation and nothing else: naming ft-lb and inches on an
/// otherwise metric card leaves the range, the dial and the time of flight bit-identical and
/// converts only the two columns that were asked for.
#[test]
fn naming_an_output_dimension_converts_the_column_and_nothing_else() {
    let metric = pre_1519_metric();
    let mut mixed = metric.clone();
    mixed["energy_unit"] = json!("ftlb");
    mixed["drop_unit"] = json!("inches");

    let base = range_table_json(&metric);
    let mixed_card = range_table_json(&mixed);
    assert_eq!(units_block(&base)["energy"], json!("J"));
    assert_eq!(units_block(&mixed_card)["energy"], json!("ft-lb"));
    assert_eq!(units_block(&mixed_card)["drop"], json!("in"));
    // Untouched dimensions keep the metric preset's labels.
    assert_eq!(units_block(&mixed_card)["distance"], json!("m"));
    assert_eq!(units_block(&mixed_card)["velocity"], json!("m/s"));

    let (base_rows, mixed_rows) = (rows(&base), rows(&mixed_card));
    assert_eq!(base_rows.len(), mixed_rows.len());
    assert!(!base_rows.is_empty(), "the fixture must produce rows");
    for (b, m) in base_rows.iter().zip(&mixed_rows) {
        // Bit-identical: the physics, and every column neither field renamed.
        assert_eq!(b["range"], m["range"], "range moved");
        assert_eq!(b["drop_adj"], m["drop_adj"], "elevation dial moved");
        assert_eq!(b["wind_adj"], m["wind_adj"], "windage dial moved");
        assert_eq!(b["velocity"], m["velocity"], "velocity column moved");
        assert_eq!(b["time"], m["time"], "time of flight moved");
        // Converted: joules -> ft-lb, and millimetres -> inches.
        let j = b["energy"].as_f64().unwrap();
        assert_close(m["energy"].as_f64().unwrap(), j * 0.737562, "energy");
        for column in ["drop_linear", "wind_linear"] {
            let mm = b[column].as_f64().unwrap();
            assert_close(m[column].as_f64().unwrap(), mm / 25.4, column);
        }
    }
}

/// The Finnish case from the design: metres downrange, a muzzle velocity off an American box
/// in fps. Stating the velocity in fps and the distance in metres describes ONE shot, and it
/// is the same shot the all-metric card describes — every row bit-identical except the
/// velocity column, which is the one thing that changed denomination.
#[test]
fn naming_an_input_dimension_states_the_same_shot_in_another_unit() {
    // 2650 fps expressed in m/s using the engine's own factor, so the two requests reach the
    // solver with the identical f64 and any difference downstream is a real defect rather
    // than a rounding artifact of the fixture.
    let mps = 2650.0 * 0.3048;
    let mut all_metric = pre_1519_metric();
    all_metric["muzzle_velocity"] = json!(mps);
    let mut fps_velocity = all_metric.clone();
    fps_velocity["muzzle_velocity"] = json!(2650.0);
    fps_velocity["velocity_unit"] = json!("fps");

    let base = range_table_json(&all_metric);
    let mixed = range_table_json(&fps_velocity);
    assert_eq!(units_block(&base)["velocity"], json!("m/s"));
    assert_eq!(units_block(&mixed)["velocity"], json!("fps"));
    assert_eq!(units_block(&mixed)["distance"], json!("m"));

    for (b, m) in rows(&base).iter().zip(&rows(&mixed)) {
        assert_eq!(b["range"], m["range"], "range moved");
        assert_eq!(b["drop_adj"], m["drop_adj"], "elevation dial moved");
        assert_eq!(b["drop_linear"], m["drop_linear"], "drop column moved");
        assert_eq!(b["energy"], m["energy"], "energy moved");
        assert_eq!(b["time"], m["time"], "time of flight moved");
        let mps = b["velocity"].as_f64().unwrap();
        assert_close(m["velocity"].as_f64().unwrap(), mps / 0.3048, "velocity");
    }
}

/// A wind speed in km/h or knots is a wind speed, not a second wind: 36 km/h, 10 m/s and
/// 19.438... knots put the same air on the bullet and produce the same drift.
#[test]
fn the_wind_axis_accepts_the_two_units_no_preset_can_reach() {
    let mut mps = pre_1519_metric();
    mps["wind_speed"] = json!(10.0);
    let mut kph = mps.clone();
    kph["wind_speed"] = json!(36.0);
    kph["wind_speed_unit"] = json!("kph");
    let mut knots = mps.clone();
    knots["wind_speed"] = json!(10.0 / (1852.0 / 3600.0));
    knots["wind_speed_unit"] = json!("knots");

    let base = range_table_json(&mps);
    assert_eq!(units_block(&base)["wind_speed"], json!("m/s"));
    for (name, card) in [
        ("km/h", range_table_json(&kph)),
        ("kn", range_table_json(&knots)),
    ] {
        assert_eq!(units_block(&card)["wind_speed"], json!(name));
        for (b, m) in rows(&base).iter().zip(&rows(&card)) {
            assert_close(
                m["wind_adj"].as_f64().unwrap(),
                b["wind_adj"].as_f64().unwrap(),
                "windage dial",
            );
        }
    }
}

// -----------------------------------------------------------------------------------------
// 3. The dimensions that do NOT take a unit from the preset, or at all
// -----------------------------------------------------------------------------------------

/// `altitude` has ignored the `units` scalar since this module was written — it goes to the
/// atmosphere unconverted, and that field is metres. So `altitude_unit` must default to
/// metres on an IMPERIAL card too; if the preset were allowed to fill it with feet, every
/// stored imperial card carrying `altitude: 1000` would quietly become a card shot at 304.8 m.
///
/// The three-way comparison is the point: the unstated card must equal the metres card and
/// must NOT equal the feet card. Equality alone would pass on a build where every altitude
/// was ignored.
#[test]
fn altitude_defaults_to_metres_even_on_an_imperial_card() {
    let unstated = pre_1519_imperial();
    let mut metres = unstated.clone();
    metres["altitude_unit"] = json!("meters");
    let mut feet = unstated.clone();
    feet["altitude_unit"] = json!("feet");

    assert_eq!(
        range_table_json(&unstated),
        range_table_json(&metres),
        "an imperial card with no altitude_unit must already be in metres"
    );
    assert_ne!(
        range_table_json(&unstated),
        range_table_json(&feet),
        "1000 ft and 1000 m are different air; the two cards cannot be the same document"
    );

    // And feet mean feet: 1000 ft is the 304.8 m card.
    let mut as_metres = unstated.clone();
    as_metres["altitude"] = json!(1000.0 * 0.3048);
    for (f, m) in rows(&range_table_json(&feet))
        .iter()
        .zip(&rows(&range_table_json(&as_metres)))
    {
        assert_close(
            f["drop_adj"].as_f64().unwrap(),
            m["drop_adj"].as_f64().unwrap(),
            "elevation dial at 1000 ft vs 304.8 m",
        );
    }
}

/// The drop COLUMN's unit is not the sight height's unit. Asking for a card that prints drop
/// in centimetres must not reinterpret a 1.5-inch scope as a 1.5-centimetre one — a 62%
/// error in the sight height, which shifts every dial on the card.
#[test]
fn the_drop_column_unit_does_not_reinterpret_the_sight_height() {
    let imperial = pre_1519_imperial();
    let mut cm_drop = imperial.clone();
    cm_drop["drop_unit"] = json!("cm");

    let base = range_table_json(&imperial);
    let cm = range_table_json(&cm_drop);
    assert_eq!(units_block(&cm)["drop"], json!("cm"));
    for (b, c) in rows(&base).iter().zip(&rows(&cm)) {
        // The solve is untouched, so the dial is bit-identical.
        assert_eq!(b["drop_adj"], c["drop_adj"], "elevation dial moved");
        // Only the linear column changed denomination: inches -> cm.
        assert_close(
            c["drop_linear"].as_f64().unwrap(),
            b["drop_linear"].as_f64().unwrap() * 2.54,
            "drop_linear",
        );
    }

    // Saying so explicitly is the same card, which is what "the sight height has its own
    // dimension" means.
    let mut explicit = cm_drop.clone();
    explicit["sight_height_unit"] = json!("inches");
    assert_eq!(cm, range_table_json(&explicit));
}

/// The angular half of the drop dimension was already per-dimension and stays where it is:
/// SMOA is a value of `adjustment_unit`, the linear-at-distance spellings are header text
/// rather than units, and clicks still need a graduation rather than a picker row.
#[test]
fn the_angular_axis_is_unchanged_and_the_pseudo_units_are_still_refused() {
    // SMOA is expressible on the angular axis, `iphy` is its alias, and the two produce one
    // card under two headings.
    let mut smoa = pre_1519_imperial();
    smoa["adjustment_unit"] = json!("smoa");
    smoa["windage_unit"] = json!("smoa");
    let mut iphy = smoa.clone();
    iphy["adjustment_unit"] = json!("iphy");
    iphy["windage_unit"] = json!("iphy");
    let smoa_card = range_table_json(&smoa);
    assert_eq!(
        units_block(&smoa_card)["elevation_adjustment"],
        json!("SMOA")
    );
    assert_eq!(
        units_block(&range_table_json(&iphy))["elevation_adjustment"],
        json!("IPHY")
    );
    for (s, i) in rows(&smoa_card).iter().zip(&rows(&range_table_json(&iphy))) {
        assert_eq!(
            s["drop_adj"], i["drop_adj"],
            "SMOA and IPHY must print one number"
        );
    }

    // The linear-at-distance spellings are NOT values of the linear drop dimension. Accepting
    // them would give a caller two ways to ask for a column that prints identical numbers.
    //
    // Each refusal is checked to be a refusal of the VALUE (`unknown variant`), not of the
    // FIELD (`unknown field`): `deny_unknown_fields` would reject every one of these on a
    // build where `drop_unit` does not exist at all, so an `is_err()` assertion here proves
    // nothing about the prohibition and would survive the field being renamed.
    for spelling in ["smoa", "moa", "mil", "inches@100yd", "cm/100m"] {
        let mut bogus = pre_1519_imperial();
        bogus["drop_unit"] = json!(spelling);
        let err = serde_json::from_value::<CardRequestV1>(bogus)
            .expect_err(&format!(
                "drop_unit must not accept the angular/at-distance spelling {spelling:?}"
            ))
            .to_string();
        assert!(
            err.contains("unknown variant"),
            "drop_unit must refuse {spelling:?} as a value, not as an unknown field: {err}"
        );
    }

    // The positive counterpart, which is what makes the loop above mean something: the field
    // is real, and the linear value the design DOES give it parses and prints.
    let mut cm = pre_1519_imperial();
    cm["drop_unit"] = json!("cm");
    assert_eq!(
        units_block(&range_table_json(&cm))["drop"],
        json!("cm"),
        "drop_unit must accept the linear value it exists to carry"
    );

    // Clicks remain a graduation, not a unit: the refusal is unchanged.
    let mut clicks = pre_1519_imperial();
    clicks["adjustment_unit"] = json!("clicks");
    clicks["windage_unit"] = json!("clicks");
    let err = range_table_v1(&parse(&clicks)).unwrap_err().to_string();
    assert!(
        err.contains("elevation_click_value"),
        "clicks without a graduation must still be refused: {err}"
    );

    // And there is no clicks unit on any of the new dimensions — again as a refusal of the
    // value rather than of the field.
    for field in [
        "distance_unit",
        "velocity_unit",
        "drop_unit",
        "wind_speed_unit",
    ] {
        let mut bogus = pre_1519_imperial();
        bogus[field] = json!("clicks");
        let err = serde_json::from_value::<CardRequestV1>(bogus)
            .expect_err(&format!("{field} must not accept 'clicks'"))
            .to_string();
        assert!(
            err.contains("unknown variant"),
            "{field} must refuse 'clicks' as a value, not as an unknown field: {err}"
        );
    }
}

/// A caller that reads a card's `units` block and writes those strings back into the next
/// request is doing the obvious thing. Every label the response prints is therefore an
/// accepted spelling of the unit it names — ALL SEVEN fields of the block, including the two
/// dial labels, which print uppercase (`MIL`/`MOA`/`SMOA`/`IPHY`/`CLICKS`) while
/// `AdjustmentUnit` used to deserialize lowercase only. A test that echoed five of seven
/// would pass while the property it is named for did not hold, so the coverage is the point
/// of the test rather than an incidental detail of it.
#[test]
fn every_label_the_response_prints_is_a_spelling_the_request_accepts() {
    /// The block's seven fields, paired with the request field that denominates each.
    const BLOCK: [(&str, &str); 7] = [
        ("distance_unit", "distance"),
        ("velocity_unit", "velocity"),
        ("drop_unit", "drop"),
        ("wind_speed_unit", "wind_speed"),
        ("energy_unit", "energy"),
        ("adjustment_unit", "elevation_adjustment"),
        ("windage_unit", "windage_adjustment"),
    ];

    // Every dial spelling the response can print, not only the default one: each card is
    // echoed back into a request and must come out as the same document.
    for dial in ["mil", "moa", "smoa", "iphy"] {
        let mut base = pre_1519_metric();
        base["adjustment_unit"] = json!(dial);
        base["windage_unit"] = json!(dial);
        let card = range_table_json(&base);
        let labels = units_block(&card);

        let mut echoed = base.clone();
        for (field, label) in BLOCK {
            let printed = labels[label].clone();
            assert!(printed.is_string(), "{label} must be a printed label");
            echoed[field] = printed;
        }
        // The dial labels really are the uppercase spellings, so this is not a test that
        // would still pass if the response started printing lowercase.
        assert_eq!(echoed["adjustment_unit"], json!(dial.to_uppercase()));
        assert_eq!(echoed["windage_unit"], json!(dial.to_uppercase()));

        assert_eq!(
            range_table_json(&echoed),
            card,
            "echoing a card's own labels back must reproduce that card ({dial})"
        );
    }

    // Clicks is the seventh label and needs its graduation, so it is echoed with one.
    let mut clicks = pre_1519_metric();
    clicks["adjustment_unit"] = json!("clicks");
    clicks["windage_unit"] = json!("clicks");
    clicks["elevation_click_value"] = json!("0.1mil");
    clicks["windage_click_value"] = json!("0.1mil");
    let card = range_table_json(&clicks);
    assert_eq!(units_block(&card)["elevation_adjustment"], json!("CLICKS"));
    let mut echoed = clicks.clone();
    echoed["adjustment_unit"] = units_block(&card)["elevation_adjustment"].clone();
    echoed["windage_unit"] = units_block(&card)["windage_adjustment"].clone();
    assert_eq!(
        range_table_json(&echoed),
        card,
        "echoing CLICKS back must reproduce the clicks card"
    );
}

/// The vocabularies that are deliberately NARROWER than the units they could have taken, both
/// of them because a wrong value here is invisible in the numbers a shooter reads.
///
/// `mass_unit` refuses the one-character tokens: `"g"` for `"grains"` is a 15.43x mass error
/// that leaves dial, drift, velocity and time bit-identical and moves only the energy column,
/// and no plausibility guard can catch it because the grain and gram ranges overlap. And
/// `diameter_unit` refuses `cm`, which the other linear dimensions do accept.
#[test]
fn the_ambiguous_short_spellings_are_refused_and_the_spelled_out_ones_are_not() {
    for (field, refused) in [
        ("mass_unit", "g"),
        ("mass_unit", "gr"),
        ("diameter_unit", "cm"),
        ("diameter_unit", "centimetres"),
    ] {
        let mut bogus = pre_1519_imperial();
        bogus[field] = json!(refused);
        let err = serde_json::from_value::<CardRequestV1>(bogus)
            .expect_err(&format!("{field} must not accept {refused:?}"))
            .to_string();
        assert!(
            err.contains("unknown variant"),
            "{field}: {refused:?} must be refused as a VALUE, not as an unknown field: {err}"
        );
    }

    // The counterpart, in the same test: the spellings that replace them do parse, and a
    // grains card and the same load in grams are two different cards.
    let mut grains = pre_1519_imperial();
    grains["mass_unit"] = json!("grains");
    grains["diameter_unit"] = json!("inches");
    let mut grams = grains.clone();
    grams["mass_unit"] = json!("grams");
    assert_eq!(
        range_table_json(&grains),
        OLD_SHAPE_IMPERIAL_RANGE_TABLE,
        "spelling out the imperial mass and diameter must be the imperial card"
    );
    assert_ne!(
        range_table_json(&grams),
        range_table_json(&grains),
        "175 grams is not 175 grains and must not print the same card"
    );
}

// -----------------------------------------------------------------------------------------
// 4. The reprint path, which is where a units mismatch would reach paper
// -----------------------------------------------------------------------------------------

/// `card.pdf`'s stored-card check compares the stored column labels against the request's own
/// axes and refuses a document whose rows mean something else. That check has to read the
/// REQUEST'S RESOLVED distance unit, not its preset: a card stored in metres reprints from an
/// otherwise-imperial request that says `distance_unit: "meters"`, and still does not reprint
/// from one that does not.
#[cfg(feature = "pdf")]
#[test]
fn a_reprint_matches_the_stored_labels_against_the_resolved_axes() {
    use ballistics_engine::card_service::{pdf_card_v1, StoredCardV1};

    // Metres downrange, mph of wind, grains of bullet — a card no preset can express.
    let mut mixed = pre_1519_imperial();
    mixed["distance_unit"] = json!("meters");
    mixed["zero_distance"] = json!(100.0);
    mixed["start"] = json!(100.0);
    mixed["end"] = json!(500.0);
    mixed["windage_unit"] = json!("mil");

    let card: Value = serde_json::from_str(&range_table_json(&mixed)).unwrap();
    assert_eq!(card["units"]["distance"], json!("m"));
    let stored: StoredCardV1 =
        serde_json::from_value(json!({ "card": card })).expect("stored card must parse");

    let printed = pdf_card_v1(&parse(&mixed), Some(&stored)).expect("the reprint must be allowed");
    assert!(printed.pdf_bytes.starts_with(b"%PDF-"));
    assert_eq!(printed.row_count, 5);

    // Drop the per-dimension field and the very same stored rows become a foreign document:
    // the request is back on yards and the rows are metres.
    let mut yards_again = mixed.clone();
    yards_again.as_object_mut().unwrap().remove("distance_unit");
    let err = pdf_card_v1(&parse(&yards_again), Some(&stored))
        .expect_err("metre rows must not print under a yard heading")
        .to_string();
    assert!(
        err.contains("units.distance") && err.contains("'m'") && err.contains("'yd'"),
        "the refusal must name both labels: {err}"
    );
}
