//! Замер 3D-основания по файлу свай «x y z_головы» (сваи d 0,88, L 22,1, грунт ЖК Остров):
//! `cargo run --release --example bench_aster_file -- piles.txt [grid|flex|curve]`.
use femgeo::aster::{self, Options};
use femgeo::piles::{Layer, Pile, SoilColumn};

fn soil(_: f64, _: f64) -> SoilColumn {
    let l = |bottom, e, nu, gamma, phi, c| Layer { bottom, e, nu, gamma, phi, c };
    SoilColumn { top: 0.65, layers: vec![l(-3.95, 17000.0, 0.4, 16.6, 17.0, 40.0), l(-8.85, 19800.0, 0.35, 17.4, 18.0, 45.0), l(-22.95, 20300.0, 0.24, 17.3, 17.0, 50.0), l(-25.75, 50000.0, 0.25, 18.0, 40.0, 1.0), l(-41.05, 39700.0, 0.2, 20.4, 20.0, 60.0), l(-56.05, 50000.0, 0.25, 18.0, 40.0, 1.0)] }
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let (d, len) = (0.88, 22.1);
    let ps: Vec<Pile> = std::fs::read_to_string(&a[0]).unwrap().lines().map(|l| {
        let v: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
        Pile { x: v[0], y: v[1], head: -10.85, length: len, d, ea: 3.25e7 * std::f64::consts::PI * d * d / 4.0 }
    }).collect();
    let mut o = Options::new(Options::find_run_aster().expect("run_aster"), std::env::temp_dir().join("femgeo_bench_file"), &ps);
    let env = |k: &str| std::env::var(k).ok().map(|v| v.parse::<f64>().unwrap());
    if let Some(v) = env("H") { o.h_zone = v; }
    if let Some(v) = env("HZ") { o.hz = v; }
    if let Some(v) = env("LAT") { o.lateral = v; }
    if let Some(v) = env("BELOW") { o.below = v; }
    if let Some(v) = env("SNAP") { o.snap = v; }
    if let Some(v) = env("HP") { o.h_pile = v; }
    if let Some(v) = env("GROW") { o.grow = v; }
    println!("h {} hz {} lateral {:.1} below {:.1} snap {}", o.h, o.hz, o.lateral, o.below, o.snap);
    let g = aster::model(&ps, &soil, &o).unwrap();
    println!("свай {}, узлов {}, элементов {}", ps.len(), g.nodes, g.elements);
    let mode = a.get(1).map(String::as_str).unwrap_or("grid");
    if mode == "flex" {
        let f = aster::flexibility(&ps, &soil, &o).unwrap();
        println!("F: {:.0} с, F11 {:.3e}, узлов {}", f.seconds, f.flexibility[0], f.nodes);
    }
    if mode == "curve" {
        let mut oc = o.clone();
        oc.lateral = len;
        oc.below = len;
        let c = aster::pile_curve(&ps[0], &soil, 11_000.0, 10, &oc).unwrap();
        println!("кривая: {:.0} с, до {:.0} кН, {:?}", c.seconds, c.n_max, c.points.iter().map(|(n, d)| (n.round(), (d * 1000.0).round() / 1000.0)).collect::<Vec<_>>());
    }
}
