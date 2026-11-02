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
    let (_, cond_sum, v) = run(None, "conditional");
    assert!(v["conditional"]["s"].as_f64().unwrap() > 0.0);
    assert!((cond_sum / lin_sum - 1.0).abs() < 0.02);
}
