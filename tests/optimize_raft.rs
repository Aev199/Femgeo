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
    run_with(None, 49).0
}

fn run_with(building: Option<&femgeo::fem::condense::Building>, n: usize) -> (Vec<femgeo::optimize::Evaluation>, femgeo::optimize::Evaluation) {
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
        cases: vec![lc],
        uls: vec![vec![1.0]],
        sls: vec![],
        sls_long: None,
        // Для проверки механизма рядов предел Δs/L ослаблен: ограничивает несущая способность.
        limits: SettlementLimits { s_max: 0.15, ds_rel_max: if n == usize::MAX { 0.004 } else { 0.002 } },
        design: &design,
        rebar_axis: [1.0, 0.0, 0.0],
        raft_h: 0.8,
        building,
        crack_short: None,
        crack_long: None,
        punching: None,
        existing: vec![],
        cache: Default::default(),
    };
    if n == usize::MAX {
        let t = std::time::Instant::now();
        let f = p.regular(&mut |_| {}).unwrap();
        println!("ряды за {:?}", t.elapsed());
        for e in &f {
            println!("свай {:2}: кг/м³ {:6.1}, Кисп {:.3}, Δs/L {:.5}", e.piles.len(), e.kg_m3, e.util_ground, e.ds_rel);
            assert!(p.layout_violations(&e.piles).is_empty());
            let sum: f64 = e.pile_n.iter().sum();
            assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A);
            // Регулярность: координаты свай образуют сетку (мало различных X и Y).
            let mut xs: Vec<i64> = e.piles.iter().map(|&k| (p.raft.nodes[k][0] * 1000.0).round() as i64).collect();
            let mut ys: Vec<i64> = e.piles.iter().map(|&k| (p.raft.nodes[k][1] * 1000.0).round() as i64).collect();
            xs.sort(); xs.dedup(); ys.sort(); ys.dedup();
            assert_eq!(xs.len() * ys.len(), e.piles.len(), "не сетка");
        }
        let e = f[0].clone();
        return (f, e);
    }
    if n == 1 {
        // Прямой расчёт для сравнения: свайное поле в модели, разложение заново.
        let grid: Vec<usize> = p.candidates().into_iter().filter(|&c| { let q = p.raft.nodes[c]; (q[0] as i64 - 1) % 3 == 0 && (q[1] as i64 - 1) % 3 == 0 }).collect();
        let e = p.evaluate(&grid).unwrap();
        let ps: Vec<femgeo::piles::Pile> = grid.iter().map(|&k| femgeo::piles::Pile { x: p.raft.nodes[k][0], y: p.raft.nodes[k][1], head: p.pile_head, length: p.pile.length, d: p.pile.d, ea: p.pile.ea }).collect();
        let kk = femgeo::piles::Field::new(&ps, &soil).stiffness().unwrap();
        let mut m = p.raft.clone();
        m.spring_groups.push(femgeo::fem::SpringGroup { nodes: grid.clone(), k: kk.clone() });
        let u = m.factorize().unwrap().solve(&m, &p.cases[0]).unwrap();
        let np = grid.len();
        for i in 0..np {
            let r: f64 = -(0..np).map(|j| kk[i * np + j] * u[6 * grid[j] + 2]).sum::<f64>();
            assert!((r - e.pile_n[i]).abs() < 1e-4 * r.abs().max(1.0), "свая {i}: прямой {r}, Вудбери {}", e.pile_n[i]);
        }
        return (vec![], e);
    }
    if n == 0 {
        let grid: Vec<usize> = p.candidates().into_iter().filter(|&c| { let q = p.raft.nodes[c]; (q[0] as i64 - 1) % 2 == 0 && (q[1] as i64 - 1) % 2 == 0 }).collect();
        return (vec![], p.evaluate(&grid).unwrap());
    }
    let t = std::time::Instant::now();
    let f = p.frontier(n, 8, &mut |_| {}).unwrap();
    println!("фронт за {:?}", t.elapsed());
    for e in &f {
        println!("свай {:2}: кг/м³ {:6.1}, Кисп грунт {:.3}, материал {:.3}, s {:.4} м, Δs/L {:.5}, допустим {}", e.piles.len(), e.kg_m3, e.util_ground, e.util_material, e.s_max, e.ds_rel, e.feasible);
        assert!(p.layout_violations(&e.piles).is_empty(), "нарушены правила расстановки");
        // Равновесие: сумма реакций свай равна нагрузке.
        let sum: f64 = e.pile_n.iter().sum();
        assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A, "Σ реакций {sum}");
    }
    let e = f[0].clone();
    (f, e)
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

/// Каркас 4×4 колонны (шаг 4 м, высота 3 м) с перекрытием 0,25 м над плитой; нагрузка на
/// перекрытие 15 кПа. Равновесие с конденсированным зданием и меньшая неравномерность осадок.
#[test]
fn building_stiffness_in_evaluation() {
    use femgeo::fem::{condense, Beam};
    let mut s = Model::default();
    let mut iface = vec![];
    let top0 = 16;
    for i in 0..4 {
        for j in 0..4 {
            s.nodes.push([1.0 + 3.0 * j as f64, 1.0 + 3.0 * i as f64, 0.0]);
        }
    }
    for k in 0..16 {
        let p = s.nodes[k];
        iface.push((k, (p[1] as usize) * (N + 1) + p[0] as usize));
    }
    for i in 0..4 {
        for j in 0..4 {
            s.nodes.push([1.0 + 3.0 * j as f64, 1.0 + 3.0 * i as f64, 3.0]);
        }
    }
    let top = |i: usize, j: usize| top0 + i * 4 + j;
    let slab = ShellMaterial { e: 3.0e7, nu: 0.2, h: 0.25 };
    for i in 0..3 {
        for j in 0..3 {
            s.shells.push(Shell { nodes: vec![top(i, j), top(i, j + 1), top(i + 1, j + 1), top(i + 1, j)], mat: slab });
        }
    }
    let (a, iz) = (0.25, 0.5f64.powi(4) / 12.0);
    for k in 0..16 {
        s.beams.push(Beam { n1: k, n2: top0 + k, e: 3.0e7, g: 1.25e7, a, iy: iz, iz, j: 2.0 * iz, ref_dir: None });
    }
    let slab_load = LoadCase { pressure: (0..9).map(|e| (e, -15.0)).collect(), ..Default::default() };
    let b = condense::building(&s, &iface, &[slab_load]).unwrap();
    let total_building = 15.0 * 9.0 * 9.0;
    let fz: f64 = b.loads[0].iter().map(|(_, f)| -f[2]).sum();
    assert!((fz - total_building).abs() < 1e-6 * total_building, "приведённая вертикальная нагрузка {fz}");

    let (_, with) = run_with(Some(&b), 0);
    let sum: f64 = with.pile_n.iter().sum();
    assert!((sum - (Q * A * A + total_building)).abs() < 1e-3 * sum, "Σ реакций {sum}");
    // Без здания — те же сваи и только нагрузка на плиту (меньше): даже так Δs/L больше.
    let (_, bare) = run_with(None, 0);
    println!("Δs/L: со зданием {:.5}, без {:.5} (без веса здания)", with.ds_rel, bare.ds_rel);
    assert!(with.ds_rel < bare.ds_rel, "здание должно уменьшать неравномерность осадок");
}

#[test]
fn regular_rows_frontier() {
    let (f, _) = run_with(None, usize::MAX);
    assert!(f.len() >= 2 && f.iter().all(|e| e.feasible));
    for w in f.windows(2) {
        assert!(w[1].piles.len() < w[0].piles.len() && w[1].kg_m3 > w[0].kg_m3, "фронт не монотонен");
    }
}

#[test]
fn woodbury_matches_direct_solution() {
    run_with(None, 1);
}
