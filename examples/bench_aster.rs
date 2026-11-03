//! Замер 3D-основания (Code_Aster): поле nx × ny свай с шагом 3d.
//! `cargo run --release --example bench_aster -- 9 10 [h] [hz]`
use femgeo::aster::{self, Options};
use femgeo::piles::{Layer, Pile, SoilColumn};

fn soil(_: f64, _: f64) -> SoilColumn {
    SoilColumn { top: -10.85, layers: vec![Layer { bottom: -20.0, e: 20_000.0, nu: 0.3, gamma: 0.0, phi: 0.0, c: 0.0 }, Layer { bottom: -200.0, e: 45_000.0, nu: 0.3, gamma: 0.0, phi: 0.0, c: 0.0 }] }
}

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (nx, ny) = (a.first().copied().unwrap_or(4.0) as usize, a.get(1).copied().unwrap_or(4.0) as usize);
    let (d, l) = (0.88, 22.1);
    let ps: Vec<Pile> = (0..ny).flat_map(|j| (0..nx).map(move |i| Pile { x: 3.0 * d * i as f64, y: 3.0 * d * j as f64, head: -10.85, length: l, d, ea: 3.25e7 * std::f64::consts::PI * d * d / 4.0 })).collect();
    let mut o = Options::new(Options::find_run_aster().expect("run_aster"), std::env::temp_dir().join("femgeo_bench_aster"), &ps);
    if let Some(&h) = a.get(2) { o.h = h; }
    if let Some(&hz) = a.get(3) { o.hz = hz; }
    let g = aster::grid(&ps, &soil, &o).unwrap();
    println!("свай {}, узлов {} (линейных), элементов {}", ps.len(), g.n_nodes(), g.n_elems());
    let f = aster::flexibility(&ps, &soil, &o).unwrap();
    let n = ps.len();
    let diag: Vec<f64> = (0..n).map(|i| f.flexibility[i * n + i]).collect();
    let row: f64 = (0..n).map(|j| f.flexibility[j]).sum();
    println!("{:.1} с; F11 {:.3e} м/кН, Fii {:.3e}..{:.3e}, Σ строки 1 {:.3e}", f.seconds, diag[0], diag.iter().cloned().fold(f64::MAX, f64::min), diag.iter().cloned().fold(0.0, f64::max), row);
}
