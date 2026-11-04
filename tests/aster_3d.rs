//! Этап 3: 3D-основание в Code_Aster. Тесты с решателем пропускаются, если `run_aster` не
//! найден (`FEMGEO_ASTER` или /opt/aster/bin/run_aster).

use femgeo::aster::{self, Options};
use femgeo::fem::{LoadCase, Model, Shell, ShellMaterial};
use femgeo::optimize::{PileType, Problem, SettlementLimits};
use femgeo::piles::layout::{PileKind, Region, Rules};
use femgeo::piles::{Field, Layer, Pile, SinglePile, SoilColumn};
use femgeo::sp63::capra_maury::{PlateSection, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass};

fn uniform(_: f64, _: f64) -> SoilColumn {
    SoilColumn { top: 0.0, layers: vec![Layer { bottom: -100.0, e: 20_000.0, nu: 0.3, gamma: 0.0, phi: 0.0, c: 0.0 }] }
}

fn layered(_: f64, _: f64) -> SoilColumn {
    SoilColumn { top: 0.0, layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3, gamma: 0.0, phi: 0.0, c: 0.0 }, Layer { bottom: -60.0, e: 40_000.0, nu: 0.3, gamma: 0.0, phi: 0.0, c: 0.0 }] }
}

fn pile(x: f64, y: f64, d: f64, l: f64) -> Pile {
    Pile { x, y, head: -0.8, length: l, d, ea: 3.0e7 * std::f64::consts::PI * d * d / 4.0 }
}

fn workdir(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("femgeo_aster_{name}_{}", std::process::id()))
}

#[test]
fn grid_geometry() {
    let ps = [pile(0.0, 0.0, 0.6, 10.0), pile(1.8, 0.0, 0.6, 10.0)];
    let o = Options::new("run_aster".into(), workdir("grid"), &ps);
    let g = aster::grid(&ps, &layered, &o).unwrap();
    let b = aster::square_side(0.6);
    // Объём свай — n·b²·L.
    let mut v = 0.0;
    let (nx, ny) = (g.x.len() - 1, g.y.len() - 1);
    for (e, &m) in g.mat.iter().enumerate() {
        let (i, j, k) = (e % nx, (e / nx) % ny, e / (nx * ny));
        if m != aster::VOID && g.materials[m].pile {
            v += (g.x[i + 1] - g.x[i]) * (g.y[j + 1] - g.y[j]) * (g.z[k + 1] - g.z[k]);
        }
    }
    assert!((v - 2.0 * b * b * 10.0).abs() < 1e-6, "объём свай {v}");
    assert!(g.heads.iter().all(|h| h.len() == 4));
    assert_eq!(g.z.last().copied(), Some(-0.8));
    assert!(g.z.iter().any(|z| (z + 8.0).abs() < 1e-9), "граница слоя не в сетке");
    // Боковые границы и низ — на max(L, B) от свай.
    assert!((g.x[0] - (-b / 2.0 - 10.0)).abs() < 1e-9 && (g.z[0] - (-10.8 - 10.0)).abs() < 1e-9);
    let m = aster::mail(&g);
    assert_eq!(m.lines().filter(|l| l.starts_with('M') && l.contains(" N")).count(), g.n_elems());
}

#[test]
fn single_pile_vs_sp24() {
    let Some(exe) = Options::find_run_aster() else {
        println!("Code_Aster не найден — пропуск");
        return;
    };
    let ps = [pile(0.0, 0.0, 0.6, 15.0)];
    let o = Options::new(exe, workdir("single"), &ps);
    let f = aster::flexibility(&ps, &uniform, &o).unwrap();
    let sp = SinglePile::new(&ps[0], &uniform(0.0, 0.0));
    let r = f.flexibility[0] / sp.compliance;
    println!("одиночная свая: 3D {:.4e} м/кН, СП 24 {:.4e}, отношение {r:.3}; узлов {}, {:.1} с", f.flexibility[0], sp.compliance, f.nodes, f.seconds);
    // Сравнение разных моделей (сплошная среда и формула (7.33)): порог задан заранее.
    assert!(r > 0.7 && r < 1.4, "отношение {r}");
}

fn raft(a: f64, n: usize) -> Model {
    let mut m = Model::default();
    let id = |i: usize, j: usize| i * (n + 1) + j;
    for i in 0..=n {
        for j in 0..=n {
            m.nodes.push([a * j as f64 / n as f64, a * i as f64 / n as f64, 0.0]);
        }
    }
    let mat = ShellMaterial { e: 3.0e7, nu: 0.2, h: 0.8 };
    for i in 0..n {
        for j in 0..n {
            m.shells.push(Shell { nodes: vec![id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)], mat });
        }
    }
    for k in 0..m.nodes.len() {
        m.fixed.push((k, [false, false, false, false, false, true]));
    }
    m.fixed.push((0, [true, true, false, false, false, false]));
    m.fixed.push((n, [false, true, false, false, false, false]));
    m
}

const A: f64 = 12.0;
const Q: f64 = 300.0;

/// Плита 12 × 12 м, q = 300 кПа, сваи d 0,6, L 15, и сетка свай 4 × 4 с шагом 3 м.
fn problem(soil: &'static (dyn Fn(f64, f64) -> SoilColumn + Sync)) -> (Problem<'static>, Vec<usize>) {
    let d = 0.6;
    let (c, r) = (ConcreteClass::B25.chars(CalcType::Uls), RebarClass::A500.chars(CalcType::Uls));
    let sec = PlateSection { h: 0.8, a_top: 0.05, a_bot: 0.05, ds_x: 0.016, ds_y: 0.016, concrete_diagram: ConcreteDiagramType::L3 };
    let design: &'static UlsContext = Box::leak(Box::new(UlsContext::new(&sec, &c, &r)));
    let m = raft(A, 12);
    let lc = LoadCase { pressure: (0..m.shells.len()).map(|s| (s, -Q)).collect(), ..Default::default() };
    let p = Problem {
        raft: m,
        pile_head: -0.8,
        pile: PileType { length: 15.0, d, ea: 3.0e7 * std::f64::consts::PI * d * d / 4.0, fd_ground: 2_500.0, material_util_per_kn: 0.0 },
        region: Region { outer: vec![[0.0, 0.0], [A, 0.0], [A, A], [0.0, A]], holes: vec![] },
        rules: Rules::sp24(PileKind::Bored, d),
        soil,
        cases: vec![lc],
        uls: vec![vec![1.0]],
        sls: vec![],
        sls_long: None,
        limits: SettlementLimits { s_max: 0.15, ds_rel_max: 0.002 },
        design,
        rebar_axis: [1.0, 0.0, 0.0],
        raft_h: 0.8,
        building: None,
        crack_short: None,
        crack_long: None,
        punching: None,
        existing: vec![],
        pile_fu: None,
        conditional: None,
        structure_punch: vec![],
        cache: Default::default(),
    };
    let grid: Vec<usize> = p.candidates().into_iter().filter(|&k| { let q = p.raft.nodes[k]; (q[0] as i64 - 1) % 3 == 0 && (q[1] as i64 - 1) % 3 == 0 }).collect();
    assert_eq!(grid.len(), 16);
    (p, grid)
}

fn corner_ratio(p: &Problem, grid: &[usize], e: &femgeo::optimize::Evaluation) -> f64 {
    let i = grid.iter().position(|&k| p.raft.nodes[k][0] == 1.0 && p.raft.nodes[k][1] == 1.0).unwrap();
    let j = grid.iter().position(|&k| p.raft.nodes[k][0] == 4.0 && p.raft.nodes[k][1] == 4.0).unwrap();
    e.pile_n[i] / e.pile_n[j]
}

#[test]
fn raft_on_3d_soil() {
    let Some(exe) = Options::find_run_aster() else {
        println!("Code_Aster не найден — пропуск");
        return;
    };
    let (p, grid) = problem(&layered);
    let ps = p.piles_of(&grid);
    let o = Options::new(exe, workdir("raft"), &ps);
    let f = aster::flexibility(&ps, &layered, &o).unwrap();
    println!("3D: узлов {}, элементов {}, {:.1} с, несимметрия {:.2e}", f.nodes, f.elements, f.seconds, f.asymmetry);
    // Одна факторизация, симметричная матрица: несимметрия — ошибка округления.
    assert!(f.asymmetry < 1e-6);
    let e3 = p.evaluate_flex(&grid, &f.flexibility).unwrap();
    let e1 = p.evaluate(&grid).unwrap();
    let sum: f64 = e3.pile_n.iter().sum();
    assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A, "Σ реакций {sum}");
    let corner = |e: &femgeo::optimize::Evaluation| corner_ratio(&p, &grid, e);
    let f24 = Field::new(&ps, &layered).flexibility;
    println!("СП 24: s {:.4} м, угол/центр {:.3}, F11 {:.3e}; 3D: s {:.4} м, угол/центр {:.3}, F11 {:.3e}", e1.s_max, corner(&e1), f24[0], e3.s_max, corner(&e3), f.flexibility[0]);
    // Упругое полупространство: угловые сваи нагружены больше центральных.
    assert!(corner(&e3) > 1.0);
    assert!(e3.s_max > 0.0 && e3.s_max.is_finite());
}

#[test]
fn service_refine_3d() {
    if Options::find_run_aster().is_none() {
        println!("Code_Aster не найден — пропуск");
        return;
    }
    let mut s = femgeo::service::Service::new();
    let mut input = s.dispatch("defaults", serde_json::Value::Null, &mut |_| {}).unwrap();
    input["refine_3d"] = serde_json::json!(1);
    // Малая плита (меньше свай — короче нелинейный расчёт); грунт с γ, φ, c (Кулон–Мор).
    input["raft"]["outline"] = serde_json::json!([[0.0, 0.0], [9.0, 0.0], [9.0, 9.0], [0.0, 9.0]]);
    input["soil"]["layers"] = serde_json::json!([
        { "bottom": -8.0, "e": 15000.0, "nu": 0.3, "gamma": 18.0, "phi": 20.0, "c": 20.0 },
        { "bottom": -60.0, "e": 40000.0, "nu": 0.3, "gamma": 19.5, "phi": 30.0, "c": 5.0 }
    ]);
    let out = s.dispatch("run", input, &mut |m| if m.starts_with("3D") { println!("{m}") }).unwrap();
    let vs = out["variants"].as_array().unwrap();
    let v3: Vec<&serde_json::Value> = vs.iter().filter(|v| v["soil_3d"] == true).collect();
    assert_eq!(v3.len(), 1);
    // Та же расстановка, что у финалиста; равновесие сохраняется.
    let n3 = v3[0]["piles"].as_array().unwrap();
    let base = vs.iter().find(|v| v["soil_3d"] == false && v["piles"].as_array().unwrap().len() == n3.len() && v["feasible"] == true).unwrap();
    let sum = |v: &serde_json::Value| v["piles"].as_array().unwrap().iter().map(|p| p["n"].as_f64().unwrap()).sum::<f64>();
    println!("Σ N: {:.0} / {:.0}", sum(base), sum(v3[0]));
    assert!((sum(base) - sum(v3[0])).abs() < 1e-3 * sum(base));
    let rp = std::env::temp_dir().join("femgeo_aster_report.html");
    let i = vs.len() - 1;
    assert_eq!(vs[i]["soil_3d"], true);
    s.dispatch("export_report", serde_json::json!({ "path": rp, "variant": i }), &mut |_| {}).unwrap();
    let html = std::fs::read_to_string(&rp).unwrap();
    assert!(html.contains("3D-основании") && html.contains("Кулон–Мор"));
    assert!(!vs[i]["nl_factor"].as_array().unwrap().is_empty());
}

/// Связка «плита — нелинейное основание» на модельном законе s_i = (F·N)_i·(1 + N_i/Nu):
/// равновесие, сходимость, перераспределение с наиболее нагруженных (угловых) свай.
#[test]
fn nonlinear_coupling_model_law() {
    let (p, grid) = problem(&layered);
    let f = Field::new(&p.piles_of(&grid), &layered).flexibility;
    let n = grid.len();
    let nu: f64 = std::env::var("NU").map(|v| v.parse().unwrap()).unwrap_or(5000.0);
    let mut calls = 0;
    let mut settle = |load: &[f64]| -> Result<Vec<f64>, String> {
        calls += 1;
        Ok((0..n).map(|i| 1.0 + load[i] / nu).collect())
    };
    let (e, it) = p.refine_nonlinear(&grid, &f, &mut settle, 0.01, 10).unwrap();
    let e1 = p.evaluate_flex(&grid, &f).unwrap();
    if std::env::var("FEMGEO_DBG").is_ok() {
        eprintln!("N_sls {:?}\nN {:?}", e1.pile_n_sls.iter().map(|x| x.round()).collect::<Vec<_>>(), e1.pile_n.iter().map(|x| x.round()).collect::<Vec<_>>());
    }
    let sum: f64 = e.pile_n.iter().sum();
    println!("итераций {it}, угол/центр: линейно {:.3}, нелинейно {:.3}; s {:.4} → {:.4}", corner_ratio(&p, &grid, &e1), corner_ratio(&p, &grid, &e), e1.s_max, e.s_max);
    assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A);
    assert!(corner_ratio(&p, &grid, &e) < corner_ratio(&p, &grid, &e1));
    assert!(e.s_max > e1.s_max);
    // Неподвижная точка: осадки свай в модели равны закону при своих нагрузках.
    let r = settle(&e.pile_n_sls).unwrap();
    let s_law: Vec<f64> = (0..n).map(|i| (0..n).map(|j| f[i * n + j] * e.pile_n_sls[j]).sum::<f64>() * r[i]).collect();
    for (k, &node) in grid.iter().enumerate() {
        assert!((e.w_nodes[node] - s_law[k]).abs() < 0.02 * s_law[k], "свая {k}: {} / {}", e.w_nodes[node], s_law[k]);
    }
}

fn soil_nl(_: f64, _: f64) -> SoilColumn {
    SoilColumn { top: 0.0, layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3, gamma: 19.0, phi: 22.0, c: 15.0 }, Layer { bottom: -60.0, e: 40_000.0, nu: 0.3, gamma: 19.5, phi: 30.0, c: 5.0 }] }
}

/// Гибридный метод (кривая одиночной сваи + упругое взаимовлияние) против полного группового
/// нелинейного расчёта той же плиты (`raft_on_nonlinear_soil`: s 46,2 мм, угол/центр 1,336;
/// 11 мин). Допуски заданы до расчёта: s ±15 %, угол/центр ±5 %.
#[test]
fn raft_hybrid_vs_group() {
    let Some(exe) = Options::find_run_aster() else {
        println!("Code_Aster не найден — пропуск");
        return;
    };
    let (p, grid) = problem(&soil_nl);
    let ps = p.piles_of(&grid);
    let o = Options::new(exe, workdir("raft_hy"), &ps);
    let t = std::time::Instant::now();
    let f = aster::flexibility(&ps, &soil_nl, &o).unwrap();
    let e0 = p.evaluate_flex(&grid, &f.flexibility).unwrap();
    let n_top = 1.5 * e0.pile_n_sls.iter().cloned().fold(0.0, f64::max);
    let mut oc = o.clone();
    oc.workdir = o.workdir.join("pile");
    oc.lateral = 15.0;
    oc.below = 15.0;
    let cv = aster::pile_curve(&ps[0], &soil_nl, n_top, 10, &oc).unwrap();
    println!("кривая до {:.0} кН (нужно {n_top:.0}), полная {}", cv.n_max, cv.complete);
    let mut settle = |n: &[f64]| -> Result<Vec<f64>, String> { n.iter().map(|&x| cv.ratio(x).ok_or(format!("срыв: {x:.0} > {:.0}", cv.n_max))).collect() };
    let (e, it) = p.refine_nonlinear(&grid, &f.flexibility, &mut settle, 0.02, 12).unwrap();
    let r = corner_ratio(&p, &grid, &e);
    println!("гибрид: {:.0} с (кривая {:.0} с), итераций {it}; кривая {:?}; s {:.4} м, угол/центр {r:.3}", t.elapsed().as_secs_f64(), cv.seconds, cv.points.iter().map(|(n, d)| (n.round(), (d * 1000.0).round() / 1000.0)).collect::<Vec<_>>(), e.s_max);
    assert!((e.s_max / 0.0462 - 1.0).abs() < 0.15, "s {}", e.s_max);
    assert!((r / 1.336 - 1.0).abs() < 0.05, "угол/центр {r}");
}

/// Нелинейное основание в Code_Aster: та же плита, грунт с γ, φ, c (полный групповой расчёт —
/// эталон гибридного метода; ≈ 11 мин, в обычном прогоне пропускается).
#[test]
#[ignore]
fn raft_on_nonlinear_soil() {
    let Some(exe) = Options::find_run_aster() else {
        println!("Code_Aster не найден — пропуск");
        return;
    };
    let (p, grid) = problem(&soil_nl);
    let ps = p.piles_of(&grid);
    let o = Options::new(exe, workdir("raft_nl"), &ps);
    let f = aster::flexibility(&ps, &soil_nl, &o).unwrap();
    let e1 = p.evaluate_flex(&grid, &f.flexibility).unwrap();
    let mut t = 0.0;
    let mut settle = |load: &[f64]| -> Result<Vec<f64>, String> {
        let r = aster::settlements_nl(&ps, &soil_nl, load, 10, &o)?;
        t += r.seconds;
        assert!(r.complete, "срыв");
        Ok(r.s.iter().zip(&r.s_el).map(|(a, b)| a / b).collect())
    };
    let (e, it) = p.refine_nonlinear(&grid, &f.flexibility, &mut settle, 0.02, 12).unwrap();
    let (lo, hi) = e.nl_factor.iter().fold((f64::MAX, 0.0f64), |a, &d| (a.0.min(d), a.1.max(d)));
    println!("итераций {it}, {t:.0} с; D {lo:.3}–{hi:.3}; угол/центр {:.3} → {:.3}; s {:.4} → {:.4}", corner_ratio(&p, &grid, &e1), corner_ratio(&p, &grid, &e), e1.s_max, e.s_max);
    let sum: f64 = e.pile_n.iter().sum();
    assert!((sum - Q * A * A).abs() < 1e-3 * Q * A * A);
    // Пластичность только увеличивает податливость (с точностью шага).
    assert!(lo > 0.98, "D {lo}");
}
