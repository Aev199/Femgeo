//! Приёмка линейного МКЭ: аналитические эталоны и инвариантность.

use femgeo::fem::{self, Beam, LoadCase, Model, Shell, ShellMaterial, Spring, V3};

const E: f64 = 3.0e7;
const NU: f64 = 0.3;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// Поворот вокруг оси (1, 2, 3)/|…| на 37° и перенос — проверка инвариантности.
fn rot(p: V3) -> V3 {
    let k = [1.0, 2.0, 3.0].map(|x: f64| x / 14f64.sqrt());
    let (s, c) = 37f64.to_radians().sin_cos();
    let kp = fem::shell::cross(k, p);
    let kd = k[0] * p[0] + k[1] * p[1] + k[2] * p[2];
    [0, 1, 2].map(|i| p[i] * c + kp[i] * s + k[i] * kd * (1.0 - c) + [5.0, -3.0, 2.0][i])
}

fn rot_vec(v: V3) -> V3 {
    let o = rot([0.0; 3]);
    let p = rot(v);
    [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
}

// ── стержни ──

fn cantilever(transform: bool) -> f64 {
    let n = 10;
    let l = 3.0;
    let tr = |p: V3| if transform { rot(p) } else { p };
    let mut m = Model::default();
    for i in 0..=n {
        m.nodes.push(tr([l * i as f64 / n as f64, 0.0, 0.0]));
    }
    for i in 0..n {
        m.beams.push(Beam { n1: i, n2: i + 1, e: E, g: E / 2.6, a: 0.06, iy: 4.5e-4, iz: 2e-4, j: 3e-4, ref_dir: Some(if transform { rot_vec([0.0, 0.0, 1.0]) } else { [0.0, 0.0, 1.0] }) });
    }
    m.fixed.push((0, [true; 6]));
    let p = if transform { rot_vec([0.0, 0.0, -10.0]) } else { [0.0, 0.0, -10.0] };
    let lc = LoadCase { nodal: vec![(n, [p[0], p[1], p[2], 0.0, 0.0, 0.0])], ..Default::default() };
    let u = m.factorize().unwrap().solve(&m, &lc).unwrap();
    let t = &u[6 * n..6 * n + 3];
    (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt()
}

#[test]
fn cantilever_tip_deflection_and_invariance() {
    let exact = 10.0 * 27.0 / (3.0 * E * 4.5e-4);
    let w = cantilever(false);
    assert!(rel(w, exact) < 1e-9, "{w} vs {exact}");
    assert!(rel(cantilever(true), exact) < 1e-8);
}

#[test]
fn beam_end_forces() {
    let mut m = Model::default();
    m.nodes = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]];
    let b = Beam { n1: 0, n2: 1, e: E, g: E / 2.6, a: 0.5, iy: 0.02, iz: 0.02, j: 0.04, ref_dir: None };
    m.beams.push(b.clone());
    m.fixed.push((0, [true; 6]));
    let lc = LoadCase { nodal: vec![(1, [5.0, 0.0, -100.0, 0.0, 0.0, 0.0])], ..Default::default() };
    let u = m.factorize().unwrap().solve(&m, &lc).unwrap();
    let (s, e) = fem::beam_forces(&m.nodes, &b, &u).unwrap();
    assert!((s[0] + 100.0).abs() < 1e-6 && (e[0] + 100.0).abs() < 1e-6, "сжатие {s:?}");
    let ms = (s[4] * s[4] + s[5] * s[5]).sqrt();
    assert!((ms - 20.0).abs() < 1e-6, "момент у заделки {ms}");
    assert!((e[4] * e[4] + e[5] * e[5]).sqrt() < 1e-6);
}

// ── пластины ──

/// Квадратная пластина a × a, сетка n × n, шарнирное опирание по контуру; `skew` — сдвиг
/// внутренних узлов (доля шага), `tri` — каждый квадрат делится на два треугольника.
fn plate(n: usize, h: f64, skew: f64, tri: bool, transform: bool) -> (Model, LoadCase, usize) {
    let a = 1.0;
    let mut m = Model::default();
    let id = |i: usize, j: usize| i * (n + 1) + j;
    for i in 0..=n {
        for j in 0..=n {
            let (mut x, mut y) = (a * j as f64 / n as f64, a * i as f64 / n as f64);
            if i > 0 && i < n && j > 0 && j < n && (i, j) != (n / 2, n / 2) {
                let s = ((i * 7 + j * 13) % 5) as f64 / 4.0 - 0.5;
                let t = ((i * 11 + j * 3) % 5) as f64 / 4.0 - 0.5;
                x += skew * s * a / n as f64;
                y += skew * t * a / n as f64;
            }
            m.nodes.push(if transform { rot([x, y, 0.0]) } else { [x, y, 0.0] });
        }
    }
    let mat = ShellMaterial { e: E, nu: NU, h };
    for i in 0..n {
        for j in 0..n {
            let q = [id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)];
            if tri {
                m.shells.push(Shell { nodes: vec![q[0], q[1], q[2]], mat });
                m.shells.push(Shell { nodes: vec![q[0], q[2], q[3]], mat });
            } else {
                m.shells.push(Shell { nodes: q.to_vec(), mat });
            }
        }
    }
    for i in 0..=n {
        for j in 0..=n {
            if i == 0 || j == 0 || i == n || j == n {
                // Прогиб закреплён; мембранные перемещения — в одном узле (против сдвига).
                let fix_all = (i, j) == (0, 0);
                m.fixed.push((id(i, j), if transform { [true, true, true, false, false, false] } else { [fix_all, fix_all || (i, j) == (0, n), true, false, false, fix_all] }));
            }
        }
    }
    if !transform {
        // Поворот вокруг нормали без жёсткости в узлах — закрепить.
        for k in 0..m.nodes.len() {
            m.fixed.push((k, [false, false, false, false, false, true]));
        }
    }
    let lc = LoadCase { pressure: (0..m.shells.len()).map(|s| (s, -1.0)).collect(), ..Default::default() };
    (m, lc, id(n / 2, n / 2))
}

fn d_plate(h: f64) -> f64 {
    E * h.powi(3) / (12.0 * (1.0 - NU * NU))
}

fn center_w(m: &Model, lc: &LoadCase, c: usize) -> (f64, Vec<f64>) {
    let u = m.factorize().unwrap().solve(m, lc).unwrap();
    let t = &u[6 * c..6 * c + 3];
    ((t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt(), u)
}

#[test]
fn navier_plate_q4_deflection_and_moment() {
    let h = 0.01;
    let (m, lc, c) = plate(16, h, 0.0, false, false);
    let (w, u) = center_w(&m, &lc, c);
    let exact = 0.00406235 / d_plate(h);
    assert!(rel(w, exact) < 0.01, "w = {w}, Навье {exact}");
    // Момент Mx в центре: 0,0479·q·a² (ν = 0,3). Среднее по четырём элементам у центра.
    let n = 16;
    let mut mx = 0.0;
    for (i, j) in [(7, 7), (7, 8), (8, 7), (8, 8)] {
        let s = &m.shells[i * n + j];
        let cs: Vec<V3> = s.nodes.iter().map(|&k| m.nodes[k]).collect();
        let ue: Vec<f64> = s.nodes.iter().flat_map(|&k| u[6 * k..6 * k + 6].to_vec()).collect();
        let (_, f) = fem::shell::forces(&cs, &s.mat, &ue).unwrap();
        mx += f[3] / 4.0;
    }
    // Давление вниз (−z): изгиб растягивает низ — M < 0.
    assert!(rel(-mx, 0.0479) < 0.03, "Mx = {mx}");
}

#[test]
fn navier_plate_skewed_mesh() {
    let h = 0.01;
    let exact = 0.00406235 / d_plate(h);
    let (m, lc, c) = plate(16, h, 0.4, false, false);
    let (w, _) = center_w(&m, &lc, c);
    assert!(rel(w, exact) < 0.02, "искажённая сетка: w = {w}, Навье {exact}");
}

#[test]
fn navier_plate_triangles() {
    let h = 0.01;
    let exact = 0.00406235 / d_plate(h);
    let (m, lc, c) = plate(16, h, 0.0, true, false);
    let (w, _) = center_w(&m, &lc, c);
    assert!(rel(w, exact) < 0.05, "T3: w = {w}, Навье {exact}");
}

#[test]
fn plate_rotated_in_space_gives_same_deflection() {
    let h = 0.01;
    let (m0, lc0, c0) = plate(8, h, 0.3, false, false);
    let (m1, lc1, c1) = plate(8, h, 0.3, false, true);
    let (w0, _) = center_w(&m0, &lc0, c0);
    let (w1, _) = center_w(&m1, &lc1, c1);
    assert!(rel(w1, w0) < 1e-4, "{w0} vs {w1}");
}

/// Мембранный патч-тест: искажённая сетка, равномерное растяжение по краю x = a →
/// Nx постоянна во всех элементах.
#[test]
fn membrane_patch_test_distorted() {
    let (mut m, _, _) = plate(4, 0.1, 0.5, false, false);
    let n = 4;
    m.fixed.clear();
    let id = |i: usize, j: usize| i * (n + 1) + j;
    for i in 0..=n {
        m.fixed.push((id(i, 0), [true, i == 0, true, true, true, true]));
        for j in 1..=n {
            m.fixed.push((id(i, j), [false, false, true, true, true, true]));
        }
    }
    // Узловые силы на краю x = 1 эквивалентны Nx = 10 кН/м (шаг 0,25).
    let nodal = (0..=n).map(|i| (id(i, n), [10.0 * if i == 0 || i == n { 0.125 } else { 0.25 }, 0.0, 0.0, 0.0, 0.0, 0.0])).collect();
    let lc = LoadCase { nodal, ..Default::default() };
    let u = m.factorize().unwrap().solve(&m, &lc).unwrap();
    for s in &m.shells {
        let cs: Vec<V3> = s.nodes.iter().map(|&k| m.nodes[k]).collect();
        let ue: Vec<f64> = s.nodes.iter().flat_map(|&k| u[6 * k..6 * k + 6].to_vec()).collect();
        let (r, f) = fem::shell::forces(&cs, &s.mat, &ue).unwrap();
        let f = fem::shell::rotate_forces(&r, &f, [1.0, 0.0, 0.0]);
        assert!((f[0] - 10.0).abs() < 1e-6 && f[1].abs() < 1e-6 && f[2].abs() < 1e-6, "{f:?}");
    }
}

/// Плита на винклеровских пружинах под равномерной нагрузкой: w = q·A/(Σk).
#[test]
fn plate_on_springs_uniform_settlement() {
    let (mut m, lc, c) = plate(6, 0.3, 0.0, false, false);
    m.fixed.retain(|&(_, f)| f == [false, false, false, false, false, true]);
    m.fixed.push((0, [true, true, false, false, false, true]));
    m.fixed.push((6, [false, true, false, false, false, false]));
    let kn = 1e4;
    for k in 0..m.nodes.len() {
        m.springs.push(Spring { node: k, k: [0.0, 0.0, kn, 0.0, 0.0, 0.0] });
    }
    let (w, _) = center_w(&m, &lc, c);
    let exact = 1.0 / (kn * m.nodes.len() as f64);
    // Узловые нагрузки на краях меньше, чем в середине — плита жёсткая, осадка почти равномерная.
    assert!(rel(w, exact) < 0.02, "{w} vs {exact}");
}

/// Изменяемая схема распознаётся, а не решается молча.
#[test]
fn mechanism_is_reported() {
    let mut m = Model::default();
    m.nodes = vec![[0.0; 3], [1.0, 0.0, 0.0]];
    m.beams.push(Beam { n1: 0, n2: 1, e: E, g: E / 2.6, a: 0.1, iy: 1e-3, iz: 1e-3, j: 1e-3, ref_dir: None });
    assert!(m.factorize().is_err());
}
