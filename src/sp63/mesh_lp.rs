//! Линейная задача подбора ортогональной сетки по требованиям фасеток:
//! min Asx + Asy при Asx·cos²θ + Asy·sin²θ ≥ r(θ) для всех θ, Asx, Asy ≥ 0.
//!
//! Две переменные — решается точно перебором вершин допустимой области.
//!
//! Источник: OpenCS `CScore/ShellSimpl/OrthogonalMeshLp.cs`, коммит `387e157`.

/// Оптимум ортогональной сетки, м²/м, и углы активных фасеток, градусы.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshSolution {
    pub asx: f64,
    pub asy: f64,
    pub active_thetas: Vec<f64>,
}

const REL_TOL: f64 = 1e-9;

/// `demands` — пары (θ, градусы; r(θ), м²/м).
pub fn solve(demands: &[(f64, f64)]) -> MeshSolution {
    struct Con {
        theta: f64,
        c2: f64,
        s2: f64,
        r: f64,
    }
    let mut cons = Vec::with_capacity(demands.len());
    let mut scale: f64 = 0.0;
    for &(theta, r) in demands {
        if !(r > 0.0) {
            continue;
        }
        let a = theta.to_radians();
        let (c, s) = (a.cos(), a.sin());
        cons.push(Con { theta, c2: c * c, s2: s * s, r });
        scale = scale.max(r);
    }
    if cons.is_empty() {
        return MeshSolution { asx: 0.0, asy: 0.0, active_thetas: vec![] };
    }
    let tol = REL_TOL * scale;

    let feasible = |x: f64, y: f64| -> bool {
        if x < -tol || y < -tol {
            return false;
        }
        cons.iter().all(|k| k.c2 * x + k.s2 * y >= k.r - tol)
    };

    let (mut best_x, mut best_y, mut best_sum) = (f64::NAN, f64::NAN, f64::INFINITY);
    let mut consider = |x: f64, y: f64| {
        if !x.is_finite() || !y.is_finite() {
            return;
        }
        let (x, y) = (x.max(0.0), y.max(0.0));
        let sum = x + y;
        // При равной сумме (ребро оптимума параллельно целевой) — более равномерная сетка.
        if sum > best_sum + tol {
            return;
        }
        if sum > best_sum - tol && (x - y).abs() >= (best_x - best_y).abs() {
            return;
        }
        if !feasible(x, y) {
            return;
        }
        best_x = x;
        best_y = y;
        best_sum = sum;
    };

    consider(scale, scale);
    for k in &cons {
        if k.c2 > 1e-12 {
            consider(k.r / k.c2, 0.0);
        }
        if k.s2 > 1e-12 {
            consider(0.0, k.r / k.s2);
        }
    }
    for i in 0..cons.len() {
        for j in i + 1..cons.len() {
            let (a, b) = (&cons[i], &cons[j]);
            let det = a.c2 * b.s2 - a.s2 * b.c2;
            if det.abs() < 1e-12 {
                continue;
            }
            consider((a.r * b.s2 - a.s2 * b.r) / det, (a.c2 * b.r - a.r * b.c2) / det);
        }
    }

    assert!(!best_x.is_nan(), "ортогональная сетка не найдена: допустимая область пуста");
    let active = cons
        .iter()
        .filter(|k| k.c2 * best_x + k.s2 * best_y <= k.r + 1e-6 * scale)
        .map(|k| k.theta)
        .collect();
    MeshSolution { asx: best_x, asy: best_y, active_thetas: active }
}
