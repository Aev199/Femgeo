//! Замер горячих мест подбора (не тест). `cargo run --release --example profile_sp63`
use femgeo::sp63::capra_maury::{self, PlateForces, PlateSection, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass, StripLimits};
use femgeo::sp63::strip_design;
use std::time::Instant;

fn main() {
    let (c, r) = (ConcreteClass::B25.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let sec = PlateSection { h: 0.2, a_top: 0.03, a_bot: 0.03, ds_x: 0.012, ds_y: 0.012, concrete_diagram: ConcreteDiagramType::L3 };
    let table = strip_design::uls_table(0.2, 0.07, -0.07, &c.diagram(ConcreteDiagramType::L3), &r.diagram(), &StripLimits::new(&c, &r));
    let t = Instant::now();
    let mut acc = 0.0;
    for i in 0..100_000 {
        acc += table.eta(&[1e-4, 2e-4], -50.0 + i as f64 * 1e-3, 20.0);
    }
    println!("eta: {:?} на вызов ({acc:.1})", t.elapsed() / 100_000);
    let mut calls = [0usize; 4];
    let t = Instant::now();
    for i in 0..18_000 {
        let i = i as f64 / 18.0;
        let res = strip_design::required_uls(-100.0 + 0.2 * i, 10.0 + 0.03 * i, &table);
        calls[res.status as usize] += 1;
    }
    println!("required_uls: {:?} на фасетку, статусы {calls:?}", t.elapsed() / 18_000);
    let ctx = UlsContext::new(&sec, &c, &r);
    let t = Instant::now();
    for i in 0..1000 {
        let i = i as f64;
        capra_maury::design_uls(&ctx, &PlateForces { nx: -100.0 + 0.2 * i, ny: 50.0 - 0.1 * i, nxy: 20.0, mx: 10.0 + 0.03 * i, my: -5.0 + 0.02 * i, mxy: 8.0 }, 10.0);
    }
    println!("design_uls: {:?} на строку", t.elapsed() / 1000);
    let rows: Vec<PlateForces> = (0..4000)
        .map(|i| {
            let i = i as f64 / 4.0;
            PlateForces { nx: -100.0 + 0.2 * i, ny: 50.0 - 0.1 * i, nxy: 20.0, mx: 10.0 + 0.03 * i, my: -5.0 + 0.02 * i, mxy: 8.0 }
        })
        .collect();
    let t = Instant::now();
    let seq: Vec<_> = rows.iter().map(|f| capra_maury::design_uls(&ctx, f, 10.0)).collect();
    let ts = t.elapsed();
    let t = Instant::now();
    let par = capra_maury::design_uls_batch(&ctx, &rows, 10.0);
    let tp = t.elapsed();
    let same = seq.iter().zip(&par).all(|(a, b)| a.total().to_bits() == b.total().to_bits());
    println!("4000 строк: последовательно {ts:?}, параллельно {tp:?}, совпадение {same}");
}
