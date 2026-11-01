//! Замер сборки и разложения: плита n × n Q4 на пружинах. `cargo run --release --example bench_fem -- 200`
use femgeo::fem::{LoadCase, Model, Shell, ShellMaterial, Spring};
use std::time::Instant;

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(100);
    let mut m = Model::default();
    let id = |i: usize, j: usize| i * (n + 1) + j;
    for i in 0..=n {
        for j in 0..=n {
            m.nodes.push([j as f64 * 0.5, i as f64 * 0.5, 0.0]);
        }
    }
    let mat = ShellMaterial { e: 3e7, nu: 0.2, h: 1.0 };
    for i in 0..n {
        for j in 0..n {
            m.shells.push(Shell { nodes: vec![id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)], mat });
        }
    }
    for k in 0..m.nodes.len() {
        m.springs.push(Spring { node: k, k: [1e2, 1e2, 1e4, 0.0, 0.0, 0.0] });
    }
    let t = Instant::now();
    let lin = m.factorize().unwrap();
    let tf = t.elapsed();
    let lc = LoadCase { pressure: (0..m.shells.len()).map(|s| (s, -10.0)).collect(), ..Default::default() };
    let t = Instant::now();
    let u = lin.solve(&m, &lc).unwrap();
    println!("{} уравнений: сборка+разложение {tf:?}, решение {:?}, w0 = {:.5}", lin.equations(), t.elapsed(), u[2]);
}
