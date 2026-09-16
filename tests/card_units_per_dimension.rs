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

/// The exact bytes pre-MBA-1519 `main` (bfc8e87) emitted for the old-shape requests below,
/// captured from that build and pasted here unaltered. These are golden in the strict sense:
/// nothing in this branch computed them, so a conversion that moved by one ulp fails here.
const PRE_1519_IMPERIAL_COME_UPS: &str = r#"{"schema_version":1,"kind":"come_ups","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_adj":-0.0002383121383369078,"come_up":-0.0002383121383369078,"velocity":2471.487053309274,"energy":2373.126515185286,"time":0.11723436596556229},{"range":200.0,"drop_adj":0.5820871745599303,"come_up":0.5823254866982672,"velocity":2299.9539308062826,"energy":2055.1450821997314,"time":0.24306405909286793},{"range":300.0,"drop_adj":1.386201586122091,"come_up":0.8041144115621608,"velocity":2135.7154714024578,"energy":1772.111319631338,"time":0.3784341545052153},{"range":400.0,"drop_adj":2.318304723103024,"come_up":0.932103136980933,"velocity":1978.346062249845,"energy":1520.578123320118,"time":0.5243858721453875},{"range":500.0,"drop_adj":3.3719523511523235,"come_up":1.0536476280492995,"velocity":1826.9457275283778,"energy":1296.7477538409037,"time":0.6821843832267613},{"range":600.0,"drop_adj":4.557991414586017,"come_up":1.186039063433694,"velocity":1680.9689068981302,"energy":1097.800875568174,"time":0.8533798162941146}]}"#;

const PRE_1519_IMPERIAL_RANGE_TABLE: &str = r#"{"schema_version":1,"kind":"range_table","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MOA"},"rows":[{"range":100.0,"drop_linear":-0.0008589384865487697,"drop_adj":-0.00023859402404132488,"wind_linear":-0.7087071995432003,"wind_adj":-0.6768153755637563,"velocity":2471.4883231079493,"energy":2373.1289537120633,"time":0.11723433140655293},{"range":200.0,"drop_linear":4.191018386010942,"drop_adj":0.5820858869459642,"wind_linear":-2.9301919795532134,"wind_adj":-1.3991666702366594,"velocity":2299.9560134472567,"energy":2055.1488041275807,"time":0.2430639078885901},{"range":300.0,"drop_linear":14.970943377972961,"drop_adj":1.3861984609234224,"wind_linear":-6.830787561730556,"wind_adj":-2.1744673738175604,"velocity":2135.717816842412,"energy":1772.1152118941698,"time":0.3784337836213863},{"range":400.0,"drop_linear":33.38349784911674,"drop_adj":2.3182984617442184,"wind_linear":-12.593748656005069,"wind_adj":-3.00675749162121,"velocity":1978.3479004209103,"energy":1520.580948997648,"time":0.5243851575386943},{"range":500.0,"drop_linear":60.69494692393412,"drop_adj":3.371941495774118,"wind_linear":-20.44174539593106,"wind_adj":-3.9043733706228325,"velocity":1826.9462466069162,"energy":1296.7484907144892,"time":0.682183186204365},{"range":600.0,"drop_linear":98.4522437976048,"drop_adj":4.557974249889111,"wind_linear":-30.64760040538944,"wind_adj":-4.878076397857819,"velocity":1680.967267951786,"energy":1097.798734855342,"time":0.8533779757660104}]}"#;

const PRE_1519_IMPERIAL_WIND_CARD: &str = r#"{"schema_version":1,"kind":"wind_card","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"yd","velocity":"fps","energy":"ft-lb","drop":"in","wind_speed":"mph","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"wind_speeds":[5.0,10.0],"wind_angles_deg":[90.0],"rows":[{"range":100.0,"wind_columns":[-0.0984309219129806,-0.1968631109842223]},{"range":200.0,"wind_columns":[-0.20348416810214265,-0.40697110827127964]},{"range":300.0,"wind_columns":[-0.31623789839125876,-0.6324803297898663]},{"range":400.0,"wind_columns":[-0.43727966417441166,-0.8745658788892408]},{"range":500.0,"wind_columns":[-0.567821871927897,-1.13565252199617]},{"range":600.0,"wind_columns":[-0.7094295707468258,-1.4188703891384]}]}"#;

const PRE_1519_METRIC_COME_UPS: &str = r#"{"schema_version":1,"kind":"come_ups","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_adj":-0.0004329443375408339,"come_up":-0.0004329443375408339,"velocity":749.2214291165774,"energy":3182.7570512805232,"time":0.12855680081771045},{"range":200.0,"drop_adj":0.6846897325451241,"come_up":0.6851226768826649,"velocity":693.1635701488178,"energy":2724.2977027173024,"time":0.26732691697149746},{"range":300.0,"drop_adj":1.5950748453135082,"come_up":0.9103851127683841,"velocity":639.6678228528231,"energy":2320.021916286218,"time":0.41751355430375675},{"range":400.0,"drop_adj":2.650926886730171,"come_up":1.055852041416663,"velocity":588.5997808435905,"energy":1964.3701032921458,"time":0.5804883620121942},{"range":500.0,"drop_adj":3.8525518849018194,"come_up":1.2016249981716483,"velocity":539.6243895183324,"energy":1651.072882015658,"time":0.7579317889455178},{"range":600.0,"drop_adj":5.216641973928331,"come_up":1.3640900890265115,"velocity":492.5372601943863,"energy":1375.5021828188696,"time":0.9519022082225888}]}"#;

const PRE_1519_METRIC_RANGE_TABLE: &str = r#"{"schema_version":1,"kind":"range_table","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"rows":[{"range":100.0,"drop_linear":-0.043327636901674266,"drop_adj":-0.00043327636901674264,"wind_linear":-21.381465376008535,"wind_adj":-0.21381465376008535,"velocity":749.2218394944025,"energy":3182.760537917555,"time":0.12855675965731383},{"range":200.0,"drop_linear":136.93764546830292,"drop_adj":0.6846882273415146,"wind_linear":-88.7228497643612,"wind_adj":-0.44361424882180606,"velocity":693.1642257150872,"energy":2724.3028557802286,"time":0.2673267364770581},{"range":300.0,"drop_linear":478.52131447848353,"drop_adj":1.5950710482616117,"wind_linear":-207.43857945583872,"wind_adj":-0.6914619315194623,"velocity":639.6685272092561,"energy":2320.027025573351,"time":0.417513107351898},{"range":400.0,"drop_linear":1060.3677499794865,"drop_adj":2.6509193749487165,"wind_linear":-383.7010758401171,"wind_adj":-0.9592526896002928,"velocity":588.6002558031348,"energy":1964.3732735204521,"time":0.5804874972077267},{"range":500.0,"drop_linear":1926.2693374702399,"drop_adj":3.8525386749404795,"wind_linear":-625.0723587363835,"wind_adj":-1.2501447174727671,"velocity":539.6243390951057,"energy":1651.0725734634839,"time":0.757930328867579},{"range":600.0,"drop_linear":3129.9725098755666,"drop_adj":5.216620849792611,"wind_linear":-940.8151071795139,"wind_adj":-1.568025178632523,"velocity":492.5363726924396,"energy":1375.4972258010262,"time":0.9518999464269216}]}"#;

const PRE_1519_METRIC_WIND_CARD: &str = r#"{"schema_version":1,"kind":"wind_card","zero_distance":100.0,"bc_for_solve":0.243,"units":{"distance":"m","velocity":"m/s","energy":"J","drop":"mm","wind_speed":"m/s","elevation_adjustment":"MIL","windage_adjustment":"MIL"},"wind_speeds":[5.0,10.0],"wind_angles_deg":[90.0],"rows":[{"range":100.0,"wind_columns":[-0.2375723202489912,-0.475159884082741]},{"range":200.0,"wind_columns":[-0.492905779349289,-0.985844981097689]},{"range":300.0,"wind_columns":[-0.7682927823335854,-1.5366407398475035]},{"range":400.0,"wind_columns":[-1.0658388574319029,-2.1317577842821906]},{"range":500.0,"wind_columns":[-1.3890531109105293,-2.7782143698540036]},{"range":600.0,"wind_columns":[-1.7422546196349102,-3.4846488495115473]}]}"#;


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
fn an_old_shape_request_still_produces_the_exact_bytes_it_always_did() {
    let imperial = pre_1519_imperial();
    assert_eq!(come_ups_json(&imperial), PRE_1519_IMPERIAL_COME_UPS);
    assert_eq!(range_table_json(&imperial), PRE_1519_IMPERIAL_RANGE_TABLE);
    assert_eq!(
        wind_card_json(&with_wind_sweep(imperial)),
        PRE_1519_IMPERIAL_WIND_CARD
    );

    let metric = pre_1519_metric();
    assert_eq!(come_ups_json(&metric), PRE_1519_METRIC_COME_UPS);
    assert_eq!(range_table_json(&metric), PRE_1519_METRIC_RANGE_TABLE);
    assert_eq!(
        wind_card_json(&with_wind_sweep(metric)),
        PRE_1519_METRIC_WIND_CARD
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
    assert_eq!(come_ups_json(&imperial), PRE_1519_IMPERIAL_COME_UPS);
    assert_eq!(range_table_json(&imperial), PRE_1519_IMPERIAL_RANGE_TABLE);

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
    assert_eq!(come_ups_json(&metric), PRE_1519_METRIC_COME_UPS);
    assert_eq!(range_table_json(&metric), PRE_1519_METRIC_RANGE_TABLE);
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
    for (name, card) in [("km/h", range_table_json(&kph)), ("kn", range_table_json(&knots))] {
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
    for (f, m) in rows(&range_table_json(&feet)).iter().zip(&rows(&range_table_json(&as_metres))) {
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
    assert_eq!(units_block(&smoa_card)["elevation_adjustment"], json!("SMOA"));
    assert_eq!(units_block(&range_table_json(&iphy))["elevation_adjustment"], json!("IPHY"));
    for (s, i) in rows(&smoa_card).iter().zip(&rows(&range_table_json(&iphy))) {
        assert_eq!(s["drop_adj"], i["drop_adj"], "SMOA and IPHY must print one number");
    }

    // The linear-at-distance spellings are NOT values of the linear drop dimension. Accepting
    // them would give a caller two ways to ask for a column that prints identical numbers.
    for spelling in ["smoa", "moa", "mil", "inches@100yd", "cm/100m"] {
        let mut bogus = pre_1519_imperial();
        bogus["drop_unit"] = json!(spelling);
        assert!(
            serde_json::from_value::<CardRequestV1>(bogus).is_err(),
            "drop_unit must not accept the angular/at-distance spelling {spelling:?}"
        );
    }

    // Clicks remain a graduation, not a unit: the refusal is unchanged.
    let mut clicks = pre_1519_imperial();
    clicks["adjustment_unit"] = json!("clicks");
    clicks["windage_unit"] = json!("clicks");
    let err = range_table_v1(&parse(&clicks)).unwrap_err().to_string();
    assert!(
        err.contains("elevation_click_value"),
        "clicks without a graduation must still be refused: {err}"
    );

    // And there is no clicks unit on any of the new dimensions.
    for field in ["distance_unit", "velocity_unit", "drop_unit", "wind_speed_unit"] {
        let mut bogus = pre_1519_imperial();
        bogus[field] = json!("clicks");
        assert!(
            serde_json::from_value::<CardRequestV1>(bogus).is_err(),
            "{field} must not accept 'clicks'"
        );
    }
}

/// A caller that reads a card's `units` block and writes those strings back into the next
/// request is doing the obvious thing. Every label the response prints is therefore an
/// accepted spelling of the unit it names.
#[test]
fn every_label_the_response_prints_is_a_spelling_the_request_accepts() {
    let card = range_table_json(&pre_1519_metric());
    let labels = units_block(&card);
    let mut echoed = pre_1519_metric();
    for (field, label) in [
        ("distance_unit", "distance"),
        ("velocity_unit", "velocity"),
        ("drop_unit", "drop"),
        ("wind_speed_unit", "wind_speed"),
        ("energy_unit", "energy"),
    ] {
        echoed[field] = labels[label].clone();
    }
    assert_eq!(
        range_table_json(&echoed),
        card,
        "echoing a card's own labels back must reproduce that card"
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
