//! Сервис приложения: расчёт примера по умолчанию и экспорт.

use femgeo::service::Service;
use serde_json::{json, Value};

#[test]
fn default_example_runs_and_exports() {
    let mut s = Service::new();
    let input = s.dispatch("defaults", Value::Null, &mut |_| {}).unwrap();
    let t = std::time::Instant::now();
    let out = s.dispatch("run", input, &mut |_| {}).unwrap();
    println!("расчёт {:?}", t.elapsed());
    let vs = out["variants"].as_array().unwrap();
    assert!(!vs.is_empty());
    for v in vs {
        println!("свай {}: {:.1} кг/м³, Кисп {:.3}/{:.3}, s {:.3}, ds {:.5}, невозм {}", v["count"], v["kg_m3"].as_f64().unwrap(), v["util_ground"].as_f64().unwrap(), v["util_material"].as_f64().unwrap(), v["s_max"].as_f64().unwrap(), v["ds_rel"].as_f64().unwrap(), v["design_impossible"]);
    }
    assert!(vs.iter().all(|v| v["feasible"].as_bool().unwrap()));
    let dir = std::env::temp_dir();
    let (csv, dxf) = (dir.join("femgeo_piles.csv"), dir.join("femgeo_piles.dxf"));
    s.dispatch("export_csv", json!({ "path": csv, "variant": 0 }), &mut |_| {}).unwrap();
    s.dispatch("export_dxf", json!({ "path": dxf, "variant": 0 }), &mut |_| {}).unwrap();
    let c = std::fs::read_to_string(&csv).unwrap();
    assert_eq!(c.lines().count(), 1 + vs[0]["count"].as_u64().unwrap() as usize);
    assert!(std::fs::read_to_string(&dxf).unwrap().contains("CIRCLE"));
    assert!(s.dispatch("nonsense", Value::Null, &mut |_| {}).is_err());
}

fn small(s: &mut Service) -> Value {
    let mut i = s.dispatch("defaults", Value::Null, &mut |_| {}).unwrap();
    i["raft"]["outline"] = json!([[0, 0], [12, 0], [12, 12], [0, 12]]);
    i
}

fn run(s: &mut Service, i: Value) -> Value {
    s.dispatch("run", i, &mut |_| {}).unwrap()
}

/// Продавливание: тонкая плита не проходит; Кисп вычислен по F/(Rbt·u·h0), у сваи
/// у края контур урезан (Кисп больше, чем был бы у внутренней сваи с той же N).
#[test]
fn punching_governs_thin_raft() {
    let mut s = Service::new();
    let mut i = small(&mut s);
    i["raft"]["h"] = json!(0.3);
    i["raft"]["cover"] = json!(0.05);
    i["acrc_lim"] = json!(0.0);
    let out = run(&mut s, i.clone());
    let v = &out["variants"][0];
    println!("h = 0,3: Кисп продавливания {}", v["util_punch"]);
    assert!(!v["feasible"].as_bool().unwrap());
    assert!(v["util_punch"].as_f64().unwrap() > 1.0);
    // Проверка формулы на наиболее нагруженной свае: u = π(d + h0) × доля внутри плиты.
    let (rbt, h0, d) = (1150.0, 0.25, 0.6);
    let n_max = v["piles"].as_array().unwrap().iter().map(|p| p["n"].as_f64().unwrap()).fold(0.0, f64::max);
    let full = n_max / (rbt * std::f64::consts::PI * (d + h0) * h0);
    assert!(v["util_punch"].as_f64().unwrap() >= full * (1.0 - 1e-9));
    i["punching"] = json!(false);
    let off = run(&mut s, i);
    assert_eq!(off["variants"][0]["util_punch"].as_f64().unwrap(), 0.0);
}

/// Трещины: строже предел — больше расход.
#[test]
fn crack_limit_raises_rebar() {
    let mut s = Service::new();
    let mut i = small(&mut s);
    i["raft"]["h"] = json!(0.6);
    i["acrc_lim"] = json!(0.0);
    let a = run(&mut s, i.clone())["variants"][0]["kg_m3"].as_f64().unwrap();
    i["acrc_lim"] = json!(0.1);
    let b = run(&mut s, i)["variants"][0]["kg_m3"].as_f64().unwrap();
    println!("кг/м³: без ПС2 {a:.1}, a_crc ≤ 0,1 мм {b:.1}");
    assert!(b > a);
}

/// Существующие сваи: первая строка — существующая расстановка; свободный режим от неё
/// стартует и держится ближе к ней, чем без ориентира.
#[test]
fn existing_piles_are_a_guide() {
    let mut s = Service::new();
    let mut i = small(&mut s);
    i["layout"] = json!("free");
    i["acrc_lim"] = json!(0.0);
    // Существующие: сетка 2,5 м со сдвигом (не в узлах сетки плиты).
    let ex: Vec<[f64; 2]> = (0..5).flat_map(|a| (0..5).map(move |b| [1.2 + 2.4 * a as f64, 1.1 + 2.45 * b as f64])).collect();
    i["existing"] = json!(ex);
    let out = run(&mut s, i.clone());
    let vs = out["variants"].as_array().unwrap();
    assert!(vs[0]["reference"].as_bool().unwrap());
    assert_eq!(vs[0]["count"].as_u64().unwrap(), 25);
    assert!(vs[0]["shift"].as_f64().unwrap() < 0.5);
    let with: Vec<f64> = vs[1..].iter().map(|v| v["shift"].as_f64().unwrap()).collect();
    i["existing"] = json!([]);
    let free = run(&mut s, i.clone());
    assert!(free["variants"][0]["reference"].as_bool() == Some(false));
    println!("смещение от существующих (с ориентиром): {with:?}");
    assert!(with[0] < 1.0, "старт свободного режима — существующие сваи");
}
