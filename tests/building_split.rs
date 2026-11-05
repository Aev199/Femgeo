//! Здание целиком → плита + сконденсированная надземная часть.

use femgeo::building::{split, Mesh};
use femgeo::fem::{Beam, Shell, ShellMaterial};

fn grid(m: &mut Mesh, z: f64, h: f64, hole: bool) -> Vec<usize> {
    let base = m.nodes.len();
    for i in 0..=10 {
        for j in 0..=10 {
            m.nodes.push([j as f64, i as f64, z]);
        }
    }
    let mut ids = vec![];
    for i in 0..10 {
        for j in 0..10 {
            if hole && (4..6).contains(&i) && (4..6).contains(&j) {
                continue;
            }
            let n = |a: usize, b: usize| base + a * 11 + b;
            ids.push(m.shells.len());
            m.shells.push(Shell { nodes: vec![n(i, j), n(i, j + 1), n(i + 1, j + 1), n(i + 1, j)], mat: ShellMaterial { e: 3e7, nu: 0.2, h } });
        }
    }
    ids
}

#[test]
fn raft_with_columns_and_slab() {
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 1.0, true);
    grid(&mut m, 3.0, 0.2, false);
    for (x, y) in [(2, 2), (8, 2), (2, 8), (8, 8)] {
        let (a, b) = (y * 11 + x, 121 + y * 11 + x);
        m.beams.push(Beam { n1: a, n2: b, e: 3e7, g: 1.25e7, a: 0.16, iy: 2.1e-3, iz: 2.1e-3, j: 3.6e-3, ref_dir: None });
    }
    // Загружение: 1 кН вниз в каждом узле перекрытия и 2 кН в узле плиты.
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -1.0, 0.0, 0.0, 0.0])).chain([(0, [0.0, 0.0, -2.0, 0.0, 0.0, 0.0])]).collect());
    m.piles.push(([2.0, 2.0, 0.0], [0.0, 0.0, 5e4, 0.0, 0.0, 0.0]));
    m.piles.push(([20.0, 2.0, 0.0], [0.0, 0.0, 5e4, 0.0, 0.0, 0.0]));
    let s = split(&m).unwrap();
    assert_eq!(s.raft.shells.len(), 96);
    assert_eq!(s.region.holes.len(), 1);
    let area = |p: &[[f64; 2]]| (0..p.len()).map(|i| p[i][0] * p[(i + 1) % p.len()][1] - p[(i + 1) % p.len()][0] * p[i][1]).sum::<f64>().abs() / 2.0;
    assert!((area(&s.region.outer) - 100.0).abs() < 1e-9 && (area(&s.region.holes[0]) - 4.0).abs() < 1e-9);
    assert_eq!((s.z_mid, s.h), (0.0, 1.0));
    assert_eq!(s.piles.len(), 1);
    assert_eq!(s.piles[0].0, [2.0, 2.0]);
    assert_eq!(s.notes, vec!["сваи КЭ 56 вне плиты: 1 из 2".to_string()]);
    let b = s.building.unwrap();
    assert_eq!(b.block.dofs.len(), 4 * 3);
    // Перекрытие на колоннах: вся нагрузка (121 кН) приходит в 4 узла сопряжения.
    let fz: f64 = b.loads[0].iter().map(|l| l.1[2]).sum();
    assert!((fz + 121.0).abs() < 1e-6, "{fz}");
    assert_eq!(s.raft_cases[0].nodal.iter().map(|l| l.1[2]).sum::<f64>(), -2.0);
}

/// Плита 10 × 10 м (отверстие 2 × 2), 4 колонны, перекрытие: подбор по сетке здания.
#[test]
fn run_on_building_mesh() {
    use femgeo::service::{default_input, Service};
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 0.8, true);
    grid(&mut m, 3.0, 0.2, false);
    for (x, y) in [(2, 2), (8, 2), (2, 8), (8, 8)] {
        m.beams.push(Beam { n1: y * 11 + x, n2: 121 + y * 11 + x, e: 3e7, g: 1.25e7, a: 0.16, iy: 2.1e-3, iz: 2.1e-3, j: 3.6e-3, ref_dir: None });
    }
    // Постоянное: 30 кН в узлах перекрытия и 20 кН в узлах плиты; кратковременное — 10 кН в узлах перекрытия.
    let raft_nodes: Vec<usize> = { let mut v: Vec<usize> = m.raft.iter().flat_map(|&s| m.shells[s].nodes.clone()).collect(); v.sort(); v.dedup(); v };
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -30.0, 0.0, 0.0, 0.0])).chain(raft_nodes.iter().map(|&n| (n, [0.0, 0.0, -20.0, 0.0, 0.0, 0.0]))).collect());
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -10.0, 0.0, 0.0, 0.0])).collect());
    for (x, y) in [(2.0, 2.0), (8.0, 2.0), (2.0, 8.0), (8.0, 8.0)] {
        m.piles.push(([x, y, 0.0], [0.0, 0.0, 1e5, 0.0, 0.0, 0.0]));
    }
    let mut input = default_input();
    assert_eq!(input.cases.len(), 2);
    input.pile.fd = 1500.0;
    let out = Service::new().run_mesh(input, &m, &mut |_| {}).unwrap();
    let vs = out["variants"].as_array().unwrap();
    assert!(vs[0]["reference"].as_bool().unwrap());
    assert_eq!(vs[0]["count"], 4);
    assert!(vs.len() > 1);
    let best = vs.iter().skip(1).find(|v| v["feasible"].as_bool().unwrap()).expect("нет допустимого варианта");
    // Нормативная нагрузка 121·40 + 117·20 ≈ 7180 кН; по ПС1 при Fd = 1500 кН — не меньше 5 свай.
    assert!(best["count"].as_u64().unwrap() >= 5, "{}", best["count"]);
    assert_eq!(out["holes"].as_array().unwrap().len(), 1);
    // Продавливание колонной 0,4 × 0,4: F ≤ ПС1 четверти перекрытия (121·(30·1,1 + 10·1,2)/4 ≈ 1361 кН),
    // несущая способность Rbt·u·h0 = 1150·(1,6 + 4·0,74)·0,74 ≈ 3880 кН (B30, h0 = 0,8 − 0,06).
    // С моментами (СП 63, п. 8.1.46) каждая моментная часть ≤ 0,5·F/Fb, по двум осям — не больше 2·F/Fb.
    let ups = best["util_punch_struct"].as_f64().unwrap();
    assert!(ups > 0.0 && ups <= 2.0 * 1361.5 / (1150.0 * 4.56 * 0.74) + 1e-3, "{ups}");
}

/// Предел Fu: сваи не больше Fu, сумма сохраняется; условный фундамент — одинаковая жёсткость.
#[test]
fn pile_cap_and_conditional() {
    use femgeo::service::{default_input, Service};
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 0.8, false);
    let raft_nodes: Vec<usize> = (0..121).collect();
    m.cases.push(raft_nodes.iter().map(|&n| (n, [0.0, 0.0, -60.0, 0.0, 0.0, 0.0])).collect());
    m.cases.push(raft_nodes.iter().map(|&n| (n, [0.0, 0.0, -10.0, 0.0, 0.0, 0.0])).collect());
    for i in 0..5 {
        for j in 0..5 {
            m.piles.push(([1.0 + 2.0 * i as f64, 1.0 + 2.0 * j as f64, 0.0], [0.0; 6]));
        }
    }
    let run = |fu: Option<f64>, model: &str| {
        let mut input = default_input();
        input.pile.fd = 1e5;
        input.pile_fu = fu;
        input.pile_model = model.into();
        if let femgeo::service::SoilInput::Layers { layers, .. } = &mut input.soil {
            layers.iter_mut().for_each(|l| { l.gamma = 19.0; l.phi = 22.0; });
        }
        let out = Service::new().run_mesh(input, &m, &mut |_| {}).unwrap();
        let v = out["variants"][0].clone();
        assert!(v["reference"].as_bool().unwrap());
        let n: Vec<f64> = v["piles"].as_array().unwrap().iter().map(|p| p["n"].as_f64().unwrap()).collect();
        (n.iter().cloned().fold(0.0, f64::max), n.iter().sum::<f64>(), v)
    };
    let (lin_max, lin_sum, _) = run(None, "sp24");
    let fu = 0.5 * (lin_max + lin_sum / 25.0);
    let (cap_max, cap_sum, v) = run(Some(fu), "sp24");
    assert!(v["capped"].as_bool().unwrap());
    assert!(cap_max <= fu * (1.0 + 1e-6) && cap_max > 0.99 * fu, "{cap_max} {fu}");
    assert!((cap_sum / lin_sum - 1.0).abs() < 0.02, "{cap_sum} {lin_sum}");
    // Упругий спрос — до перераспределения: наибольший равен линейному, сваи сверх Fu отмечены.
    let ps = v["piles"].as_array().unwrap();
    let demand = ps.iter().map(|p| p["demand"].as_f64().unwrap()).fold(0.0, f64::max);
    assert!((demand / lin_max - 1.0).abs() < 1e-6, "{demand} {lin_max}");
    let capped = ps.iter().filter(|p| p["capped"].as_bool().unwrap()).count();
    assert!(capped > 0 && capped == v["capped_count"].as_u64().unwrap() as usize);
    assert!(ps.iter().all(|p| p["capped"].as_bool().unwrap() == (p["demand"].as_f64().unwrap() > fu)));
    let (_, cond_sum, v) = run(None, "conditional");
    assert!(v["conditional"]["s"].as_f64().unwrap() > 0.0);
    assert!((cond_sum / lin_sum - 1.0).abs() < 0.02);
}

/// Конструкции: колонна 0,4 × 0,6 и Г-образная стена (2 торца и угол).
#[test]
fn structures_detected() {
    use femgeo::building::StructureKind;
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 1.0, false);
    // Колонна в узле (5, 8): A = 0,24, Iy = 0,4·0,6³/12.
    let top = m.nodes.len();
    m.nodes.push([5.0, 8.0, 3.0]);
    m.beams.push(Beam { n1: 8 * 11 + 5, n2: top, e: 3e7, g: 1.25e7, a: 0.24, iy: 0.4 * 0.216 / 12.0, iz: 0.6 * 0.064 / 12.0, j: 1e-3, ref_dir: None });
    // Стена: (2,2) → (6,2) → (6,6), пластины 1 × 3 м, t = 0,2.
    let wall = |a: [usize; 2], b: [usize; 2], m: &mut Mesh| {
        let (na, nb) = (a[1] * 11 + a[0], b[1] * 11 + b[0]);
        let (ta, tb) = (m.nodes.len(), m.nodes.len() + 1);
        m.nodes.push([a[0] as f64, a[1] as f64, 3.0]);
        m.nodes.push([b[0] as f64, b[1] as f64, 3.0]);
        m.shells.push(Shell { nodes: vec![na, nb, tb, ta], mat: ShellMaterial { e: 3e7, nu: 0.2, h: 0.2 } });
    };
    for x in 2..6 {
        wall([x, 2], [x + 1, 2], &mut m);
    }
    for y in 2..6 {
        wall([6, y], [6, y + 1], &mut m);
    }
    m.cases.push(vec![(top, [0.0, 0.0, -100.0, 0.0, 0.0, 0.0])]);
    let s = split(&m).unwrap();
    let at = |k: usize| [s.raft.nodes[s.structures[k].node][0], s.raft.nodes[s.structures[k].node][1]];
    let col: Vec<usize> = (0..s.structures.len()).filter(|&k| matches!(s.structures[k].kind, StructureKind::Column { .. })).collect();
    assert_eq!(col.len(), 1);
    if let StructureKind::Column { b, h } = s.structures[col[0]].kind {
        assert!((b - 0.4).abs() < 1e-9 && (h - 0.6).abs() < 1e-9, "{b} {h}");
    }
    let mut ends: Vec<[f64; 2]> = (0..s.structures.len()).filter(|&k| matches!(s.structures[k].kind, StructureKind::WallEnd { .. })).map(at).collect();
    ends.sort_by(|a, b| a[0].total_cmp(&b[0]));
    assert_eq!(ends, vec![[2.0, 2.0], [6.0, 6.0]]);
    let corners: Vec<[f64; 2]> = (0..s.structures.len()).filter(|&k| matches!(s.structures[k].kind, StructureKind::WallCorner { .. })).map(at).collect();
    assert_eq!(corners, vec![[6.0, 2.0]]);
}

/// Армирование подземной части: стены подвала по контуру и перекрытие ниже поверхности грунта.
#[test]
fn underground_reinforcement() {
    use femgeo::service::{default_input, Service, SoilInput};
    use serde_json::json;
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 0.8, false);
    let top = grid(&mut m, 3.0, 0.25, false);
    // Стены по контуру: пластины 1 × 3 м между узлами плиты и перекрытия.
    let mut walls = 0;
    for k in 0..10 {
        for (a, b) in [((0, k), (0, k + 1)), ((10, k), (10, k + 1)), ((k, 0), (k + 1, 0)), ((k, 10), (k + 1, 10))] {
            let (na, nb) = (a.1 * 11 + a.0, b.1 * 11 + b.0);
            m.shells.push(Shell { nodes: vec![na, nb, 121 + nb, 121 + na], mat: ShellMaterial { e: 3e7, nu: 0.2, h: 0.3 } });
            walls += 1;
        }
    }
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -40.0, 0.0, 0.0, 0.0])).collect());
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -10.0, 0.0, 0.0, 0.0])).collect());
    for (x, y) in [(2.0, 2.0), (8.0, 2.0), (2.0, 8.0), (8.0, 8.0), (5.0, 5.0)] {
        m.piles.push(([x, y, 0.0], [0.0; 6]));
    }
    let mut input = default_input();
    input.pile.fd = 1e5;
    if let SoilInput::Layers { top, .. } = &mut input.soil {
        *top = 3.5;
    }
    let mut s = Service::new();
    s.run_mesh(input, &m, &mut |_| {}).unwrap();
    let r = s.dispatch("design_underground", json!({ "variant": 0 }), &mut |_| {}).unwrap();
    assert_eq!(r["elements"].as_u64().unwrap() as usize, walls + top.len(), "{r}");
    assert!(r["kg"].as_f64().unwrap() > 0.0 && r["impossible"].as_u64().unwrap() == 0, "{r}");
    let path = std::env::temp_dir().join("femgeo_underground.dxf");
    s.dispatch("export_underground_dxf", json!({ "path": path }), &mut |_| {}).unwrap();
    assert!(std::fs::read_to_string(&path).unwrap().contains("AS_TOP_X"));
}

/// Проверка заданной расстановки: тот же результат, что у существующих свай в расчёте;
/// повторная проверка (другая расстановка) — на сохранённой подготовке, без расхождений.
#[test]
fn check_piles_matches_reference() {
    use femgeo::service::{default_input, Service};
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 0.8, false);
    m.cases.push((0..121).map(|n| (n, [0.0, 0.0, -60.0, 0.0, 0.0, 0.0])).collect());
    m.cases.push((0..121).map(|n| (n, [0.0, 0.0, -10.0, 0.0, 0.0, 0.0])).collect());
    let pts: Vec<[f64; 2]> = (0..25).map(|k| [1.0 + 2.0 * (k / 5) as f64, 1.0 + 2.0 * (k % 5) as f64]).collect();
    m.piles = pts.iter().map(|p| ([p[0], p[1], 0.0], [0.0; 6])).collect();
    let n_of = |v: &serde_json::Value| -> Vec<f64> { v["piles"].as_array().unwrap().iter().map(|p| p["n"].as_f64().unwrap()).collect() };
    let mut s = Service::new();
    let reference = s.run_mesh(default_input(), &m, &mut |_| {}).unwrap()["variants"][0].clone();
    assert!(reference["reference"].as_bool().unwrap());
    let mut input = default_input();
    input.check = Some(pts.clone());
    let first = s.run_mesh(default_input(), &m, &mut |_| {}).unwrap()["variants"].as_array().unwrap().len();
    let out = s.run_mesh(input.clone(), &m, &mut |_| {}).unwrap();
    // Добавлен к вариантам расчёта, отмечен как проверка.
    let vs = out["variants"].as_array().unwrap();
    assert_eq!(vs.len(), first + 1);
    assert!(vs[..first].iter().all(|v| !v["checked"].as_bool().unwrap()));
    let v = &vs[first];
    assert!(v["checked"].as_bool().unwrap());
    let (a, b) = (n_of(v), n_of(&reference));
    assert!(a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-6 * y.abs().max(1.0)), "{a:?} {b:?}");
    // Без одной сваи, точка вне плиты и слишком близкая пара — замечания.
    let mut fewer = pts[1..].to_vec();
    fewer.push([50.0, 50.0]);
    input.check = Some(fewer.clone());
    let out = s.run_mesh(input.clone(), &m, &mut |_| {}).unwrap();
    assert_eq!(out["variants"].as_array().unwrap().len(), first + 2);
    let last = |o: &serde_json::Value| o["variants"].as_array().unwrap().last().unwrap().clone();
    assert_eq!(last(&out)["count"], 24);
    assert_eq!(out["check_notes"].as_array().unwrap().len(), 1);
    // Та же расстановка с новым сервисом (без подготовки) — тот же результат.
    let fresh = Service::new().run_mesh(input.clone(), &m, &mut |_| {}).unwrap();
    assert_eq!(fresh["variants"].as_array().unwrap().len(), 1);
    assert_eq!(n_of(&fresh["variants"][0]), n_of(&last(&out)));
    fewer.push([2.0, 1.0]);
    input.check = Some(fewer);
    let out = s.run_mesh(input, &m, &mut |_| {}).unwrap();
    assert!(!last(&out)["feasible"].as_bool().unwrap());
}
