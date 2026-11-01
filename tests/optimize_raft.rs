//! Подбор расстановки свай: плита 12×12 м, h = 0,8 м, сетка 1 м, равномерная нагрузка ПС1.

use femgeo::fem::{LoadCase, Model, Shell, ShellMaterial};
use femgeo::optimize::{PileType, Problem, SettlementLimits};
use femgeo::piles::layout::{PileKind, Region, Rules};
use femgeo::piles::{Layer, SoilColumn};
use femgeo::sp63::capra_maury::{PlateSection, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass};
use femgeo::sp63::{pile, section::Section};

const A: f64 = 12.0;
const N: usize = 12;
const Q: f64 = 300.0;

fn raft() -> Model {
    let mut m = Model::default();
    let id = |i: usize, j: usize| i * (N + 1) + j;
    for i in 0..=N {
        for j in 0..=N {
            m.nodes.push([A * j as f64 / N as f64, A * i as f64 / N as f64, 0.0]);
        }
    }
    let mat = ShellMaterial { e: 3.0e7, nu: 0.2, h: 0.8 };
    for i in 0..N {
        for j in 0..N {
            m.shells.push(Shell { nodes: vec![id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)], mat });
        }
    }
    for k in 0..m.nodes.len() {
        m.fixed.push((k, [false, false, false, false, false, true]));
    }
    m.fixed.push((0, [true, true, false, false, false, false]));
    m.fixed.push((N, [false, true, false, false, false, false]));
    m
}

fn soil(_: f64, _: f64) -> SoilColumn {
    SoilColumn { top: 0.0, layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3 }, Layer { bottom: -60.0, e: 40_000.0, nu: 0.3 }] }
}

fn material_util_per_kn(d: f64, l0: f64) -> f64 {
    let (c, r) = (ConcreteClass::B25.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let s = Section::circular(d / 2.0, 96, 20, 8, 0.016, d / 2.0 - 0.06, 0.0, c.diagram(ConcreteDiagramType::L2), r.diagram(), c.eps_b2);
    let (n, mx, my) = pile::with_accidental_eccentricity(-1000.0, 0.0, 0.0, l0, d, Some((1.0, 0.0))).unwrap();
    pile::all_factor(&s, n, mx, my, &Default::default()).utilization / 1000.0
}

fn run() -> Vec<femgeo::optimize::Evaluation> {
    let d = 0.6;
    let (c, r) = (ConcreteClass::B25.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let sec = PlateSection { h: 0.8, a_top: 0.05, a_bot: 0.05, ds_x: 0.016, ds_y: 0.016, concrete_diagram: ConcreteDiagramType::L3 };
    let design = UlsContext::new(&sec, &c, &r);
    let m = raft();
    let lc = LoadCase { pressure: (0..m.shells.len()).map(|s| (s, -Q)).collect(), ..Default::default() };
    let p = Problem {
        raft: m,
        pile_head: -0.8,
        pile: PileType { length: 15.0, d, ea: 3.0e7 * std::f64::consts::PI * d * d / 4.0, fd_ground: 2_500.0, material_util_per_kn: material_util_per_kn(d, 3.0) },
        region: Region { outer: vec![[0.0, 0.0], [A, 0.0], [A, A], [0.0, A]], holes: vec![] },
        rules: Rules::sp24(PileKind::Bored, d),
        soil: &soil,
        combinations: vec![lc],
        limits: SettlementLimits { s_max: 0.15, ds_rel_max: 0.002 },
        design: &design,
        rebar_axis: [1.0, 0.0, 0.0],
        raft_h: 0.8,
    };
    let t = std::time::Instant::now();
    let f = p.frontier(49, 8, &mut |_| {}).unwrap();
    println!("фронт за {:?}", t.elapsed());
    for e in &f {
        println!("свай {:2}: кг/м³ {:6.1}, Кисп грунт {:.3}, материал {:.3}, s {:.4} м, Δs/L {:.5}, допустим {}", e.piles.len(), e.kg_m3, e.util_ground, e.util_material, e.s_max, e.ds_rel, e.feasible);
        assert!(p.layout_violations(&e.piles).is_empty(), "нарушены правила расстановки");
        // Равновесие: сумма реакций свай равна нагрузке.
        let sum: f64 = e.pile_n.iter().sum();
        assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A, "Σ реакций {sum}");
    }
    f
}

#[test]
fn frontier_on_uniform_raft() {
    let f = run();
    assert!(f.len() >= 2, "фронт из одного варианта");
    assert!(f.iter().all(|e| e.feasible));
    for w in f.windows(2) {
        assert!(w[1].piles.len() < w[0].piles.len());
    }
    // Нижняя граница числа свай по грунту: Q·A²/Fd = 17,3.
    assert!(f.last().unwrap().piles.len() >= 18);
    // Детерминизм.
    let g = run();
    assert_eq!(f.iter().map(|e| e.piles.clone()).collect::<Vec<_>>(), g.iter().map(|e| e.piles.clone()).collect::<Vec<_>>());
}
