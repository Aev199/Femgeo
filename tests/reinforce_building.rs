//! Армирование здания на синтетической модели: плита 8×8 м на пружинах, стены по периметру
//! высотой 3 м, перекрытие; симметрия относительно диагонали x = y.

use femgeo::fem::{LoadCase, Model, Shell, ShellMaterial, Spring};
use femgeo::reinforce::{self, Group};
use femgeo::sp63::capra_maury::{PlateSection, SlsContext, SlsOptions, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass};

const N: usize = 8;

/// Узел (i, j) уровня z (0 — плита, 1 — перекрытие); элементы: 0..64 плита, 64..128
/// перекрытие, далее стены. Возвращает модель, число элементов плиты и перекрытия.
fn building(roof_h: f64) -> Model {
    let mut m = Model::default();
    let id = |lvl: usize, i: usize, j: usize| lvl * (N + 1) * (N + 1) + i * (N + 1) + j;
    for lvl in 0..2 {
        for i in 0..=N {
            for j in 0..=N {
                m.nodes.push([j as f64, i as f64, 3.0 * lvl as f64]);
            }
        }
    }
    let raft = ShellMaterial { e: 3.25e7, nu: 0.2, h: 0.6 };
    let roof = ShellMaterial { e: 3.25e7, nu: 0.2, h: roof_h };
    let wall = ShellMaterial { e: 3.25e7, nu: 0.2, h: 0.25 };
    for (lvl, mat) in [(0, raft), (1, roof)] {
        for i in 0..N {
            for j in 0..N {
                m.shells.push(Shell { nodes: vec![id(lvl, i, j), id(lvl, i, j + 1), id(lvl, i + 1, j + 1), id(lvl, i + 1, j)], mat });
            }
        }
    }
    // Стены: четыре стороны, один элемент по высоте.
    for k in 0..N {
        for (a, b) in [((0, k), (0, k + 1)), ((k, N), (k + 1, N)), ((N, k + 1), (N, k)), ((k + 1, 0), (k, 0))] {
            m.shells.push(Shell { nodes: vec![id(0, a.0, a.1), id(0, b.0, b.1), id(1, b.0, b.1), id(1, a.0, a.1)], mat: wall });
        }
    }
    for k in 0..(N + 1) * (N + 1) {
        m.springs.push(Spring { node: k, k: [0.0, 0.0, 3e4, 0.0, 0.0, 0.0] });
    }
    for k in 0..m.nodes.len() {
        m.fixed.push((k, [false, false, false, false, false, true]));
    }
    m.fixed.push((0, [true, true, false, false, false, false]));
    m.fixed.push((N, [false, true, false, false, false, false]));
    m
}

fn run(roof_h: f64, background: [f64; 4]) -> (Model, Vec<reinforce::Element>) {
    run_q(roof_h, background, 8.0)
}

fn run_q(roof_h: f64, background: [f64; 4], q_roof: f64) -> (Model, Vec<reinforce::Element>) {
    let m = building(roof_h);
    let (c, r) = (ConcreteClass::B30.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let mk = |h: f64| PlateSection { h, a_top: 0.04, a_bot: 0.04, ds_x: 0.016, ds_y: 0.016, concrete_diagram: ConcreteDiagramType::L3 };
    let (sr, sw, sroof) = (mk(0.6), mk(0.25), mk(roof_h));
    let (ur, uw, uroof) = (UlsContext::new(&sr, &c, &r), UlsContext::new(&sw, &c, &r), UlsContext::new(&sroof, &c, &r));
    let (cn, rn) = (ConcreteClass::B30.chars(CalcType::Sls), RebarClass::A500.chars(CalcType::Sls));
    let crack = SlsContext::new(&sr, &cn, &rn, Some((&c, &r)), SlsOptions { acrc_lim_mm: 0.3, phi1: 1.4, ..Default::default() });
    let g = |section: PlateSection, uls| Group { section, concrete: c.clone(), rebar: r.clone(), rsw: RebarClass::A500.rsw(), uls, crack_short: None, crack_long: None, background };
    let mut groups = vec![g(sr, &ur), g(sroof, &uroof), g(sw, &uw)];
    groups[0].crack_long = Some(&crack);
    let group_of: Vec<Option<usize>> = (0..m.shells.len()).map(|e| Some(if e < 64 { 0 } else if e < 128 { 1 } else { 2 })).collect();
    // Загружения: перекрытие q_roof, плита 5 кПа (постоянные).
    let roof = LoadCase { pressure: (64..128).map(|e| (e, -q_roof)).collect(), ..Default::default() };
    let raft = LoadCase { pressure: (0..64).map(|e| (e, -5.0)).collect(), ..Default::default() };
    let out = reinforce::design(&m, &groups, &group_of, &[roof, raft], &[vec![1.1, 1.1]], &[vec![1.0, 1.0]], Some(&vec![1.0, 1.0])).expect("армирование");
    (m, out)
}

#[test]
fn fields_symmetric_and_background() {
    let (m, out) = run(0.2, [0.0; 4]);
    // Симметрия x ↔ y: элемент плиты (i, j) ↔ (j, i), верх x ↔ верх y, низ x ↔ низ y.
    let tol = |a: f64| 1e-3 * a.max(1e-5);
    for i in 0..N {
        for j in 0..N {
            let (a, b) = (&out[i * N + j], &out[j * N + i]);
            for (p, q) in [(0, 1), (2, 3)] {
                assert!((a.required[p] - b.required[q]).abs() <= tol(a.required[p]), "({i},{j}) {:?} vs {:?}", a.required, b.required);
            }
        }
    }
    let total = |o: &[reinforce::Element]| o.iter().map(|e| e.required.iter().sum::<f64>()).sum::<f64>();
    assert!(total(&out) > 0.0);
    assert!(out.iter().all(|e| e.addition == e.required));
    // Фон больше требуемого — добавки нет.
    let (_, out2) = run(0.2, [0.01; 4]);
    assert!(out2.iter().all(|e| e.addition == [0.0; 4]));
    // Оси армирования стен: горизонталь в плоскости.
    let ax = reinforce::rebar_axis([0.0, -1.0, 0.0]);
    assert!(ax[2] == 0.0 && (ax[0].abs() - 1.0).abs() < 1e-12);
    // DXF: по грани на каждый элемент с ненулевым значением.
    let v: Vec<f64> = out.iter().map(|e| e.required[2]).collect();
    let d = reinforce::dxf_fields(&m, &[("AS_BOT_X", &v)], "см²/м");
    assert_eq!(d.matches("3DFACE").count(), v.iter().filter(|&&x| x > 1e-12).count());
}

#[test]
fn thin_roof_reported_impossible() {
    let (_, out) = run_q(0.1, [0.0; 4], 80.0);
    assert!(out[64..128].iter().any(|e| e.impossible));
    assert!(!out[..64].iter().any(|e| e.impossible));
}

#[test]
fn rebar_outside_section_is_an_error() {
    let r = std::panic::catch_unwind(|| run(0.06, [0.0; 4]));
    assert!(r.is_err(), "сечение с арматурой вне толщины должно отклоняться");
}
