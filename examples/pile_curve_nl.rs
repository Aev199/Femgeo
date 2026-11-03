//! Кривая «нагрузка — осадка» одиночной сваи на упругопластическом основании (Code_Aster,
//! Кулон–Мор): `cargo run --release --example pile_curve_nl -- [N, кН]`.
use femgeo::aster::{self, Options};
use femgeo::piles::{Layer, Pile, SoilColumn};
fn soil(_: f64, _: f64) -> SoilColumn { SoilColumn { top: 0.0, layers: vec![Layer { bottom: -200.0, e: 20_000.0, nu: 0.3, gamma: 19.0, phi: 28.0, c: 10.0 }] } }
fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let d = 0.6;
    let ps = [Pile { x: 0.0, y: 0.0, head: -0.8, length: 15.0, d, ea: 3.0e7 * std::f64::consts::PI * d * d / 4.0 }];
    let o = Options::new(Options::find_run_aster().unwrap(), std::env::temp_dir().join("femgeo_pile_curve"), &ps);
    let f = aster::flexibility(&ps, &soil, &o).unwrap();
    println!("F11 {:.4e}, {:.1}s", f.flexibility[0], f.seconds);
    let n = a.first().copied().unwrap_or(2000.0);
    match aster::settlements_nl(&ps, &soil, &[n], 10, &o) {
        Ok(r) => { println!("полный {}", r.complete); for (t, s) in &r.curve { println!("{:.2} N {:6.0}: s {:.4} м, лин {:.4}", t, t * n, s[0], t * n * f.flexibility[0]); } println!("{:.1}s", r.seconds); }
        Err(e) => println!("ERR {e}"),
    }
}
