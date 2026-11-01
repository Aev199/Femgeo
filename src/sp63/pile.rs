//! Прочность нормального сечения сваи: предельный множитель нагрузки по НДМ и случайный
//! эксцентриситет.
//!
//! Источник: OpenCS `CScore/LimitForceSolver.cs` (`AllFactor`, `MomentFactor`, `Bisect`,
//! `IsFeasible`), коммит `387e157`; случайный эксцентриситет — как в `cscore-pile-poc`
//! (`pile_rebar.py`, режим `sp63_stat_determinate`).

use super::section::{self, Plane, Section};

/// Параметры поиска (значения по умолчанию — как в OpenCS).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LimitOptions {
    pub solver_tol: f64,
    pub solver_max_iter: usize,
    pub bisect_tol: f64,
    pub bisect_max_iter: usize,
    /// Учитывать растянутый бетон (в ULS `cscore-pile-poc` — нет).
    pub ten: bool,
}

impl Default for LimitOptions {
    fn default() -> Self {
        Self { solver_tol: 0.5, solver_max_iter: 60, bisect_tol: 1e-4, bisect_max_iter: 60, ten: false }
    }
}

/// Что исчерпывает несущую способность.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Governing {
    Concrete,
    Rebar,
    Both,
    None,
}

/// Итог проверки.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// Множитель найден с точностью `bisect_tol`.
    Ok,
    /// Множитель не ограничен сверху до 1e12 (нагрузки практически нет).
    Unbounded,
    /// Уже при k = 0 нет допустимого состояния или поиск не сошёлся — результат не определён.
    Indeterminate,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LimitResult {
    pub factor: f64,
    /// Кисп = 1/k.
    pub utilization: f64,
    pub status: Status,
    pub governing: Governing,
    pub plane: Option<Plane>,
    pub n_limit: f64,
    pub mx_limit: f64,
    pub my_limit: f64,
}

fn feasible(s: &Section, n: f64, mx: f64, my: f64, o: &LimitOptions) -> Option<Plane> {
    let e = section::solve(s, n, mx, my, o.ten, o.solver_tol, o.solver_max_iter);
    if !e.converged || e.residual > o.solver_tol {
        return None;
    }
    let p = e.plane;
    if s.eps_contour_min(&p) < s.eps_cu || s.eps_bars_max(&p) > s.eps_su {
        return None;
    }
    Some(p)
}

fn bisect(s: &Section, o: &LimitOptions, load: impl Fn(f64) -> (f64, f64, f64)) -> LimitResult {
    let at = |k: f64| {
        let (n, mx, my) = load(k);
        feasible(s, n, mx, my, o)
    };
    let mut grid = vec![0.0];
    let mut k = 1e-6;
    while k < 1e12 {
        grid.push(k);
        k *= 2.0;
    }
    grid.push(1e12);

    let (mut a, mut b) = (None, None);
    for &k in &grid {
        if at(k).is_some() {
            a = Some(k);
        } else if a.is_some() {
            b = Some(k);
            break;
        }
    }
    let result = |k: f64, status: Status, plane: Option<Plane>| {
        let (n, mx, my) = load(k);
        let governing = match plane {
            None => Governing::None,
            Some(p) => {
                let c = (s.eps_contour_min(&p) - s.eps_cu).abs() < s.eps_cu.abs() * 0.02;
                let r = (s.eps_bars_max(&p) - s.eps_su).abs() < s.eps_su.abs() * 0.02;
                match (c, r) {
                    (true, true) => Governing::Both,
                    (false, true) => Governing::Rebar,
                    _ => Governing::Concrete,
                }
            }
        };
        LimitResult {
            factor: k,
            utilization: if k > 1e-15 { 1.0 / k } else { f64::INFINITY },
            status,
            governing,
            plane,
            n_limit: n,
            mx_limit: mx,
            my_limit: my,
        }
    };
    let Some(a) = a else { return result(0.0, Status::Indeterminate, None) };
    let Some(b) = b else { return result(a, Status::Unbounded, None) };

    let (mut left, mut right) = (a, b);
    let mut best = at(left);
    for _ in 0..o.bisect_max_iter {
        let mid = 0.5 * (left + right);
        match at(mid) {
            Some(p) => {
                left = mid;
                best = Some(p);
            }
            None => right = mid,
        }
        if right - left <= o.bisect_tol {
            break;
        }
    }
    let status = if right - left <= o.bisect_tol * 10.0 { Status::Ok } else { Status::Indeterminate };
    result(left, status, best)
}

/// Множитель ко всем усилиям: k·(N, Mx, My) на границе (как `AllFactor`).
pub fn all_factor(s: &Section, n: f64, mx: f64, my: f64, o: &LimitOptions) -> LimitResult {
    bisect(s, o, |k| (k * n, k * mx, k * my))
}

/// Множитель к моментам при неизменной N (как `MomentFactor` без η).
pub fn moment_factor(s: &Section, n: f64, mx: f64, my: f64, o: &LimitOptions) -> LimitResult {
    bisect(s, o, |k| (n, k * mx, k * my))
}

/// Случайный эксцентриситет e_a = max(l0/600, D/30, 10 мм) (п. 8.1.7) для сжатого элемента
/// статически определимой схемы: к моменту добавляется |N|·e_a по направлению вектора
/// (Mx, My). При N ≥ 0 усилия не меняются. При нулевом моменте направление задаётся
/// `zero_moment_dir` (единичный вектор); None — ошибка.
pub fn with_accidental_eccentricity(
    n: f64,
    mx: f64,
    my: f64,
    l0: f64,
    d: f64,
    zero_moment_dir: Option<(f64, f64)>,
) -> Result<(f64, f64, f64), String> {
    if n >= 0.0 {
        return Ok((n, mx, my));
    }
    let ea = (l0 / 600.0).max(d / 30.0).max(0.01);
    let m = (mx * mx + my * my).sqrt();
    let (ux, uy) = if m > 1e-12 {
        (mx / m, my / m)
    } else {
        zero_moment_dir.ok_or("центральное сжатие: задайте направление случайного эксцентриситета")?
    };
    let dm = n.abs() * ea;
    Ok((n, mx + dm * ux, my + dm * uy))
}
