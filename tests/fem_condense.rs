//! Конденсация здания на плиту: совпадение с полной моделью «здание + плита».

use femgeo::fem::condense::building;
use femgeo::fem::{Beam, LoadCase, Model, Shell, ShellMaterial, Spring};

const E: f64 = 3.0e7;

/// Плита 8×8 м, сетка 1 м, на винклеровских пружинах; узлы — сетка 9×9.
fn raft() -> Model {
    let mut m = Model::default();
    let id = |i: usize, j: usize| i * 9 + j;
    for i in 0..9 {
        for j in 0..9 {
            m.nodes.push([j as f64, i as f64, 0.0]);
        }
    }
    let mat = ShellMaterial { e: E, nu: 0.2, h: 0.6 };
    for i in 0..8 {
        for j in 0..8 {
            m.shells.push(Shell { nodes: vec![id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)], mat });
        }
    }
    for k in 0..m.nodes.len() {
        m.springs.push(Spring { node: k, k: [0.0, 0.0, 2e4, 0.0, 0.0, 0.0] });
        m.fixed.push((k, [false, false, false, false, false, true]));
    }
    m.fixed.push((0, [true, true, false, false, false, false]));
    m.fixed.push((8, [false, true, false, false, false, false]));
    m
}

/// Здание: 4 колонны 0,4×0,4 высотой 3 м и плита перекрытия 8×8 (сетка 2 м) над ними,
/// стена толщиной 0,2 по краю y = 0 (элементы 2×3 м). Узлы основания — первые в списке.
/// Возвращает модель и пары (узел здания, узел плиты по координатам).
fn superstructure() -> (Model, Vec<(usize, usize)>) {
    let mut m = Model::default();
    let mut iface = vec![];
    let raft_id = |x: f64, y: f64| (y as usize) * 9 + x as usize;
    let mut base = |m: &mut Model, x: f64, y: f64| {
        m.nodes.push([x, y, 0.0]);
        iface.push((m.nodes.len() - 1, raft_id(x, y)));
        m.nodes.len() - 1
    };
    let cols: Vec<usize> = [(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0)].iter().map(|&(x, y)| base(&mut m, x, y)).collect();
    let wall: Vec<usize> = (0..5).map(|k| base(&mut m, 2.0 * k as f64, 0.0)).collect();
    let top0 = m.nodes.len();
    for i in 0..5 {
        for j in 0..5 {
            m.nodes.push([2.0 * j as f64, 2.0 * i as f64, 3.0]);
        }
    }
    let top = |i: usize, j: usize| top0 + i * 5 + j;
    let slab = ShellMaterial { e: E, nu: 0.2, h: 0.2 };
    for i in 0..4 {
        for j in 0..4 {
            m.shells.push(Shell { nodes: vec![top(i, j), top(i, j + 1), top(i + 1, j + 1), top(i + 1, j)], mat: slab });
        }
    }
    let (a, i) = (0.16, 0.4f64.powi(4) / 12.0);
    for (k, &(ci, cj)) in [(1usize, 1usize), (1, 3), (3, 3), (3, 1)].iter().enumerate() {
        m.beams.push(Beam { n1: cols[k], n2: top(ci, cj), e: E, g: E / 2.4, a, iy: i, iz: i, j: 2.0 * i, ref_dir: None });
    }
    let w = ShellMaterial { e: E, nu: 0.2, h: 0.2 };
    for k in 0..4 {
        m.shells.push(Shell { nodes: vec![wall[k], wall[k + 1], top(0, k + 1), top(0, k)], mat: w });
    }
    (m, iface)
}

#[test]
fn condensed_building_matches_full_model() {
    let r = raft();
    let (s, iface) = superstructure();
    let slab_load = LoadCase { pressure: (0..16).map(|e| (e, -10.0)).collect(), ..Default::default() };
    let raft_load = LoadCase { pressure: (0..64).map(|e| (e, -20.0)).collect(), ..Default::default() };

    // Полная модель: узлы здания, кроме основания, добавлены к плите; основание — узлы плиты.
    let mut full = r.clone();
    let mut map: Vec<usize> = vec![usize::MAX; s.nodes.len()];
    for &(sn, rn) in &iface {
        map[sn] = rn;
    }
    for (k, p) in s.nodes.iter().enumerate() {
        if map[k] == usize::MAX {
            map[k] = full.nodes.len();
            full.nodes.push(*p);
        }
    }
    let shell0 = full.shells.len();
    for sh in &s.shells {
        full.shells.push(Shell { nodes: sh.nodes.iter().map(|&n| map[n]).collect(), mat: sh.mat });
    }
    for b in &s.beams {
        full.beams.push(Beam { n1: map[b.n1], n2: map[b.n2], ..b.clone() });
    }
    for &(_, rn) in &iface {
        full.fixed.push((rn, [true, true, false, false, false, true]));
    }
    let mut lc = raft_load.clone();
    lc.pressure.extend(slab_load.pressure.iter().map(|&(e, q)| (shell0 + e, q)));
    let u_full = full.factorize().unwrap().solve(&full, &lc).unwrap();

    // Конденсированная: плита + суперэлемент + приведённые нагрузки.
    let b = building(&s, &iface, &[slab_load]).unwrap();
    let mut cm = r.clone();
    cm.blocks.push(b.block.clone());
    for &(_, rn) in &iface {
        cm.fixed.push((rn, [true, true, false, false, false, true]));
    }
    let mut lc2 = raft_load.clone();
    lc2.nodal.extend(b.loads[0].iter().copied());
    let u_c = cm.factorize().unwrap().solve(&cm, &lc2).unwrap();

    let scale = (0..81).map(|n| u_full[6 * n + 2].abs()).fold(0.0, f64::max);
    for n in 0..81 {
        for k in [2, 3, 4] {
            let (a, c) = (u_full[6 * n + k], u_c[6 * n + k]);
            assert!((a - c).abs() < 1e-8 * scale.max(1e-12) * if k == 2 { 1.0 } else { 10.0 }, "узел {n}, ст. св. {k}: {a} vs {c}");
        }
    }

    // Здание ужесточает плиту: неравномерность осадок меньше, чем без него.
    let mut bare = r.clone();
    let mut lc3 = raft_load;
    lc3.nodal.extend(b.loads[0].iter().map(|&(n, f)| (n, [0.0, 0.0, f[2], 0.0, 0.0, 0.0])));
    bare.fixed.extend(iface.iter().map(|&(_, rn)| (rn, [true, true, false, false, false, true])));
    let u_b = bare.factorize().unwrap().solve(&bare, &lc3).unwrap();
    let spread = |u: &[f64]| {
        let w: Vec<f64> = (0..81).map(|n| u[6 * n + 2]).collect();
        w.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - w.iter().cloned().fold(f64::INFINITY, f64::min)
    };
    assert!(spread(&u_c) < spread(&u_b), "со зданием {} без {}", spread(&u_c), spread(&u_b));
}
