//! Замер расчёта через сервис: `cargo run --release --example bench_service -- 48 36`
use femgeo::service::Service;
use serde_json::{json, Value};
fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (w, h) = (a.first().copied().unwrap_or(24.0), a.get(1).copied().unwrap_or(18.0));
    let mut s = Service::new();
    let mut i = s.dispatch("defaults", Value::Null, &mut |_| {}).unwrap();
    i["raft"]["outline"] = json!([[0, 0], [w, 0], [w, h], [0, h]]);
    let t = std::time::Instant::now();
    let out = s.dispatch("run", i, &mut |st| eprintln!("  {st}")).unwrap();
    println!("{w}×{h} м: вариантов {}, {:?}", out["variants"].as_array().unwrap().len(), t.elapsed());
}
