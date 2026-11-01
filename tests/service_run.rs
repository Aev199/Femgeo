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
