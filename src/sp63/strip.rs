//! Одноосная плоская НДМ полосы плиты 1 м × h (пп. 8.1.20–8.1.30, 8.2.14–8.2.18 СП 63).
//!
//! Источник: OpenCS `CScore/ShellSimpl/StripNdmSection.cs`, `StripNdm.cs`,
//! `StripNdmUlsTable.cs`, коммит `387e157`. Добавлен вариант ψs = 1 (п. 8.2.15) —
//! консервативная проверка, выбранная в `cscore-pile-poc`.
//!
//! Знаки: «+» N — растяжение, «+» M растягивает грань +Z (верх); ε(z) = ε₀ + κ·z.

use std::sync::Arc;

use super::diagram::Diagram;
use super::materials::StripLimits;

/// Ряд арматуры полосы: координата центра от срединной плоскости (+ вверх), м, и площадь, м²/м.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StripRow {
    pub z: f64,
    pub as_: f64,
}

/// Отклик полосы на плоскость деформаций: усилия и касательная матрица ∂(N, M)/∂(ε₀, κ).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StripResponse {
    pub n: f64,
    pub m: f64,
    pub k00: f64,
    pub k01: f64,
    pub k11: f64,
}

impl StripResponse {
    pub fn det(&self) -> f64 {
        self.k00 * self.k11 - self.k01 * self.k01
    }
}

const GAUSS_W0: f64 = 8.0 / 9.0;
const GAUSS_W1: f64 = 5.0 / 9.0;

/// Полоса 1 м × h: бетон по толщине и ряды арматуры.
///
/// Бетон интегрируется кусочной квадратурой: толщина делится в точках, где ε(z) попадает в
/// узлы диаграммы и в ноль; на каждом куске — Гаусс–Лежандр 3 точки (для кусочно-линейных
/// диаграмм интеграл точный). Площадь бетона под арматурой не вычитается.
#[derive(Clone, Debug)]
pub struct StripSection {
    pub h: f64,
    pub rows: Vec<StripRow>,
    /// Диаграммы разделяются между копиями полосы (подбор меняет только площади рядов).
    pub concrete: Arc<Diagram>,
    pub rebar: Arc<Diagram>,
    breaks: Arc<[f64]>,
}

/// Наибольшее число узлов разбиения толщины (узлы диаграммы бетона и ноль).
const MAX_BREAKS: usize = 30;

impl StripSection {
    pub fn new(h: f64, rows: Vec<StripRow>, concrete: Diagram, rebar: Diagram) -> Self {
        Self::shared(h, rows, Arc::new(concrete), Arc::new(rebar))
    }

    /// Полоса с уже разделяемыми диаграммами.
    pub fn shared(h: f64, rows: Vec<StripRow>, concrete: Arc<Diagram>, rebar: Arc<Diagram>) -> Self {
        assert!(h > 0.0, "толщина полосы должна быть положительной");
        let mut breaks: Vec<f64> = vec![0.0];
        breaks.extend_from_slice(concrete.compression.xs());
        breaks.extend_from_slice(concrete.tension.xs());
        breaks.sort_by(f64::total_cmp);
        breaks.dedup();
        assert!(breaks.len() <= MAX_BREAKS, "слишком много узлов диаграммы бетона");
        Self { h, rows, concrete, rebar, breaks: breaks.into() }
    }

    /// Та же полоса с другими площадями рядов.
    pub fn with_areas(&self, areas: &[f64]) -> Self {
        let mut s = self.clone();
        s.set_areas(areas);
        s
    }

    /// Задать площади рядов на месте.
    pub fn set_areas(&mut self, areas: &[f64]) {
        assert_eq!(areas.len(), self.rows.len(), "число площадей не совпадает с числом рядов");
        for (r, &a) in self.rows.iter_mut().zip(areas) {
            r.as_ = a;
        }
    }

    /// Усилия и касательная матрица при плоскости деформаций (ε₀, κ).
    ///
    /// `tension_concrete` — учитывать растянутый бетон: ПС1 и сечение с трещиной — нет,
    /// поиск M_crc (п. 8.2.14) — да.
    pub fn integrate(&self, eps0: f64, kappa: f64, tension_concrete: bool) -> StripResponse {
        let h = self.h;
        let half = h / 2.0;
        let gx = (3.0f64 / 5.0).sqrt();
        let mut r = StripResponse::default();

        let mut buf = [0.0f64; MAX_BREAKS + 2];
        let mut count = 0;
        buf[count] = -half;
        count += 1;
        if kappa.abs() * h > 1e-14 {
            for &e in self.breaks.iter() {
                let z = (e - eps0) / kappa;
                if z > -half && z < half {
                    buf[count] = z;
                    count += 1;
                }
            }
        }
        buf[count] = half;
        count += 1;
        let zs = &mut buf[..count];
        zs.sort_unstable_by(f64::total_cmp);

        let add_concrete = |z: f64, w: f64, r: &mut StripResponse| {
            let (sig, e2) = self.concrete.sig(eps0 + kappa * z, tension_concrete, true);
            r.n += sig * w;
            r.m += sig * w * z;
            r.k00 += e2 * w;
            r.k01 += e2 * w * z;
            r.k11 += e2 * w * z * z;
        };

        for w in zs.windows(2) {
            let (a, b) = (w[0], w[1]);
            let hl = 0.5 * (b - a);
            if hl < 1e-15 {
                continue;
            }
            let mid = 0.5 * (a + b);
            add_concrete(mid, hl * GAUSS_W0, &mut r);
            add_concrete(mid - hl * gx, hl * GAUSS_W1, &mut r);
            add_concrete(mid + hl * gx, hl * GAUSS_W1, &mut r);
        }

        for row in &self.rows {
            if row.as_ <= 0.0 {
                continue;
            }
            let (sig, e2) = self.rebar.sig(eps0 + kappa * row.z, true, true);
            let (z, w) = (row.z, row.as_);
            r.n += sig * w;
            r.m += sig * w * z;
            r.k00 += e2 * w;
            r.k01 += e2 * w * z;
            r.k11 += e2 * w * z * z;
        }
        r
    }

    /// Начальная (упругая) касательная матрица полного сечения.
    pub fn elastic_stiffness(&self) -> StripResponse {
        let probe = 1e-9;
        let (_, eb) = self.concrete.sig(-probe, false, true);
        let (_, es) = self.rebar.sig(probe, true, true);
        let h = self.h;
        let mut r = StripResponse { k00: eb * h, k11: eb * h * h * h / 12.0, ..Default::default() };
        for row in &self.rows {
            if row.as_ <= 0.0 {
                continue;
            }
            r.k00 += es * row.as_;
            r.k01 += es * row.as_ * row.z;
            r.k11 += es * row.as_ * row.z * row.z;
        }
        r
    }

    /// Масштаб усилий для допусков: |σb,min|·h, кН/м.
    pub fn force_scale(&self) -> f64 {
        let mut s_min: f64 = 0.0;
        for &x in self.concrete.compression.xs() {
            s_min = s_min.min(self.concrete.stress(x, false, true));
        }
        (s_min.abs() * self.h).max(1.0)
    }
}

/// Найденное НДС полосы.
#[derive(Clone, Debug, PartialEq)]
pub struct StripState {
    pub eps0: f64,
    pub kappa: f64,
    pub converged: bool,
    pub iterations: usize,
    pub reason: String,
    pub n: f64,
    pub m: f64,
}

impl StripState {
    pub fn eps_at(&self, z: f64) -> f64 {
        self.eps0 + self.kappa * z
    }

    /// Координата нейтральной оси от срединной плоскости; NaN при κ ≈ 0.
    pub fn neutral_axis(&self) -> f64 {
        if self.kappa.abs() > 1e-14 {
            -self.eps0 / self.kappa
        } else {
            f64::NAN
        }
    }

    /// Высота сжатой зоны от сжатой грани: h — всё сечение сжато, 0 — сжатия нет.
    pub fn compression_height(&self, h: f64) -> f64 {
        let top = self.eps_at(h / 2.0);
        let bot = self.eps_at(-h / 2.0);
        if top >= 0.0 && bot >= 0.0 {
            return 0.0;
        }
        if top <= 0.0 && bot <= 0.0 {
            return h;
        }
        let comp = top.min(bot);
        let tens = top.max(bot);
        h * -comp / (tens - comp)
    }

    fn failed(eps0: f64, kappa: f64, iterations: usize, reason: &str) -> Self {
        Self { eps0, kappa, converged: false, iterations, reason: reason.into(), n: 0.0, m: 0.0 }
    }
}

/// Шарнир предельной плоскости деформаций (п. 8.1.30).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pivot {
    None,
    /// Растянутый ряд на εs,ult.
    A,
    /// Сжатая грань на εb2 (двузначная эпюра).
    B,
    /// Однозначная эпюра сжатия, εb,ult интерполяцией между εb2 и εb0.
    C,
}

/// Результат ПС1 полосы по лучу усилий.
#[derive(Clone, Debug)]
pub struct StripUls {
    /// Предельный множитель нагрузки (∞ при нулевой нагрузке).
    pub lambda: f64,
    /// Кисп = 1/λ; `f64::MAX`, если луч не пересекает границу.
    pub eta: f64,
    pub pivot: Pivot,
    /// Шарнир на верхней грани или ряду.
    pub top_face: bool,
    pub ultimate: Option<StripState>,
    pub reason: String,
}

/// Способ определения σs,crc для ψs (п. 8.2.15, 8.2.18).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SigmaSCrcMethod {
    /// σs,crc — σs сечения с трещиной при M = M_crc и той же N (общее определение п. 8.2.18).
    #[default]
    ReleasedConcrete8137,
    /// Ф. (8.138): ψs = 1 − 0,8·M_crc/M (изгибаемые элементы).
    CrackingMoment8138,
    /// ψs = 1 (допускается п. 8.2.15) — консервативная проверка.
    PsiOne,
}

/// ПС2 одного ряда полосы.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StripSlsRow {
    pub z: f64,
    pub as_: f64,
    pub sigma_s: f64,
    pub sigma_s_crc: f64,
    pub psi_s: f64,
    pub h_bt: f64,
    pub ls_m: f64,
    pub phi3: f64,
    pub acrc_mm: f64,
}

/// Результат ПС2 полосы.
#[derive(Clone, Debug)]
pub struct StripSls {
    /// Момент трещинообразования грани, растянутой моментом (модуль), кН·м/м.
    pub mcrc: f64,
    pub crack_face_top: bool,
    pub cracked: bool,
    pub h_bt: f64,
    pub rows: Vec<StripSlsRow>,
    pub acrc_mm: f64,
    pub state: Option<StripState>,
    pub converged: bool,
    pub reason: String,
}

fn norm(rn: f64, rm: f64, tol_n: f64, tol_m: f64) -> f64 {
    rn.abs() / tol_n + rm.abs() / tol_m
}

/// НДС при заданных (N, M): Ньютон 2×2 с демпфированием, при неудаче — вложенный одномерный
/// поиск. Нет решения — `converged = false` с причиной.
pub fn solve(s: &StripSection, n: f64, m: f64, tension_concrete: bool) -> StripState {
    let fs = s.force_scale();
    let tol_n = 1e-9 * fs;
    let tol_m = 1e-9 * fs * s.h;
    let newton = solve_newton(s, n, m, tension_concrete, tol_n, tol_m);
    if newton.converged {
        return newton;
    }
    let mut nested = solve_nested(s, n, m, tension_concrete, 1e3 * tol_n, 1e3 * tol_m);
    nested.iterations += newton.iterations;
    nested
}

fn solve_newton(s: &StripSection, n: f64, m: f64, ten: bool, tol_n: f64, tol_m: f64) -> StripState {
    let el = s.elastic_stiffness();
    let det0 = el.det();
    let mut e0 = (el.k11 * n - el.k01 * m) / det0;
    let mut k = (el.k00 * m - el.k01 * n) / det0;

    let mut f = s.integrate(e0, k, ten);
    let mut rn = n - f.n;
    let mut rm = m - f.m;
    let mut nrm = norm(rn, rm, tol_n, tol_m);
    let mut it = 0;
    while it < 60 {
        if rn.abs() <= tol_n && rm.abs() <= tol_m {
            return StripState { eps0: e0, kappa: k, converged: true, iterations: it, reason: String::new(), n: f.n, m: f.m };
        }
        let (mut k00, mut k01, mut k11, mut det) = (f.k00, f.k01, f.k11, f.det());
        if !(k00 > 0.0) || !(det > 1e-10 * k00 * k11) || !(k11 > 0.0) {
            // Вырожденная касательная — шаг по упругой матрице.
            k00 = el.k00;
            k01 = el.k01;
            k11 = el.k11;
            det = det0;
        }
        let de0 = (k11 * rn - k01 * rm) / det;
        let dk = (k00 * rm - k01 * rn) / det;

        let mut alpha = 1.0;
        let mut f_new = StripResponse::default();
        let (mut rn_new, mut rm_new, mut norm_new) = (0.0, 0.0, f64::INFINITY);
        for _ in 0..12 {
            f_new = s.integrate(e0 + alpha * de0, k + alpha * dk, ten);
            rn_new = n - f_new.n;
            rm_new = m - f_new.m;
            norm_new = norm(rn_new, rm_new, tol_n, tol_m);
            if norm_new < nrm {
                break;
            }
            alpha *= 0.5;
        }
        if !(norm_new < nrm) {
            break;
        }
        e0 += alpha * de0;
        k += alpha * dk;
        f = f_new;
        rn = rn_new;
        rm = rm_new;
        nrm = norm_new;
        it += 1;
    }
    StripState { eps0: e0, kappa: k, converged: false, iterations: it, reason: "Ньютон не сошёлся".into(), n: f.n, m: f.m }
}

fn solve_nested(s: &StripSection, n: f64, m: f64, ten: bool, tol_n: f64, tol_m: f64) -> StripState {
    let h = s.h;
    let mut evals = 0usize;

    // ε₀ при фиксированной κ: N(ε₀) неубывающая для монотонных диаграмм.
    let eps0_for = |kappa: f64, evals: &mut usize| -> Option<f64> {
        let mut lo = -0.1 - kappa.abs() * h;
        let mut hi = 0.1 + kappa.abs() * h;
        let n_lo = s.integrate(lo, kappa, ten).n - n;
        let n_hi = s.integrate(hi, kappa, ten).n - n;
        *evals += 2;
        if n_lo > tol_n || n_hi < -tol_n {
            return None;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            let r = s.integrate(mid, kappa, ten).n - n;
            *evals += 1;
            if r.abs() <= tol_n || hi - lo < 1e-15 {
                return Some(mid);
            }
            if r < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    };

    let g = |kappa: f64, evals: &mut usize| -> Option<(f64, f64)> {
        let e0 = eps0_for(kappa, evals)?;
        *evals += 1;
        Some((s.integrate(e0, kappa, ten).m - m, e0))
    };

    let finalize = |e0: f64, kappa: f64, evals: usize| -> StripState {
        let f = s.integrate(e0, kappa, ten);
        let ok = (f.n - n).abs() <= tol_n && (f.m - m).abs() <= tol_m;
        StripState {
            eps0: e0,
            kappa,
            converged: ok,
            iterations: evals,
            reason: if ok {
                String::new()
            } else {
                "Равновесие не найдено: отклик полосы разрывный (срыв растянутого бетона) или за пределом".into()
            },
            n: f.n,
            m: f.m,
        }
    };

    let Some((g0, e00)) = g(0.0, &mut evals) else {
        return StripState::failed(0.0, 0.0, evals, "Продольная сила вне несущей способности полосы");
    };
    if g0.abs() <= tol_m {
        return finalize(e00, 0.0, evals);
    }

    let dir = if g0 < 0.0 { 1.0 } else { -1.0 }; // M(κ) возрастает с κ
    let (mut k_lo, mut g_lo, mut e_lo) = (0.0, g0, e00);
    let mut step = 1e-4 / h;
    let (mut k_hi, mut g_hi, mut e_hi) = (0.0, f64::NAN, 0.0);
    let mut bracket = false;
    while step < 2.0 / h {
        k_hi = dir * step;
        let Some((gv, ev)) = g(k_hi, &mut evals) else { break };
        g_hi = gv;
        e_hi = ev;
        if g_hi.signum() != g_lo.signum() || g_hi.abs() <= tol_m {
            bracket = true;
            break;
        }
        k_lo = k_hi;
        g_lo = g_hi;
        e_lo = e_hi;
        step *= 2.0;
    }
    if !bracket {
        return StripState::failed(e_lo, k_lo, evals, "Момент вне несущей способности полосы при данной продольной силе");
    }
    if g_hi.abs() <= tol_m {
        return finalize(e_hi, k_hi, evals);
    }

    let (mut ka, mut kb, mut ga) = (k_lo, k_hi, g_lo);
    for _ in 0..200 {
        let mid = 0.5 * (ka + kb);
        let Some((gv, e_mid)) = g(mid, &mut evals) else {
            kb = mid;
            continue;
        };
        if gv.abs() <= tol_m || (kb - ka).abs() * h < 1e-15 {
            return finalize(e_mid, mid, evals);
        }
        if gv.signum() == ga.signum() {
            ka = mid;
            ga = gv;
        } else {
            kb = mid;
        }
    }
    match g(0.5 * (ka + kb), &mut evals) {
        Some((_, e_last)) => finalize(e_last, 0.5 * (ka + kb), evals),
        None => StripState::failed(e_lo, k_lo, evals, "Поиск κ не сошёлся"),
    }
}

/// Использование деформаций по п. 8.1.30: max(|εb|/εb,ult, εs/εs,ult).
pub fn strain_utilization(s: &StripSection, st: &StripState, lim: &StripLimits) -> f64 {
    let top = st.eps_at(s.h / 2.0);
    let bot = st.eps_at(-s.h / 2.0);
    let e2 = top.min(bot);
    let e1 = top.max(bot);
    let mut uc = 0.0;
    if e2 < 0.0 {
        let ult = if e1 >= 0.0 {
            lim.eps_b2
        } else {
            lim.eps_b2 - (lim.eps_b2 - lim.eps_b0) * (e1 / e2).clamp(0.0, 1.0)
        };
        uc = -e2 / ult;
    }
    let mut us: f64 = 0.0;
    for row in &s.rows {
        if row.as_ > 0.0 {
            us = us.max(st.eps_at(row.z) / lim.eps_s_ult);
        }
    }
    f64::max(uc, us)
}

/// Отрезок контура предельных плоскостей.
#[derive(Clone, Copy, Debug)]
enum SegKind {
    ABot,
    BTop,
    CTop,
    CBot,
    BBot,
    ATop,
}

const SEGMENTS: [(SegKind, Pivot, bool); 6] = [
    (SegKind::ABot, Pivot::A, false),
    (SegKind::BTop, Pivot::B, true),
    (SegKind::CTop, Pivot::C, true),
    (SegKind::CBot, Pivot::C, false),
    (SegKind::BBot, Pivot::B, false),
    (SegKind::ATop, Pivot::A, true),
];

/// Замкнутый контур предельных плоскостей (п. 8.1.30) в деформациях граней: шарнир A у нижнего
/// ряда, B у верхней грани, C сверху и снизу, B у нижней грани, A у верхнего ряда. Ряд со
/// стороны без арматуры заменяется условным рядом на грани — это лишь сужает контур (в запас).
struct Contour {
    h: f64,
    half: f64,
    zt: f64,
    zb: f64,
    eu: f64,
    b2: f64,
    b0: f64,
    u_top: f64,
    u_bot: f64,
    x_start_top: f64,
    x_start_bot: f64,
}

/// Плоскость через две точки (z, ε) → деформации граней (верх, низ).
fn through(z1: f64, e1: f64, z2: f64, e2: f64, half: f64) -> (f64, f64) {
    let k = (e2 - e1) / (z2 - z1);
    (e1 + k * (half - z1), e1 + k * (-half - z1))
}

impl Contour {
    fn new(s: &StripSection, lim: &StripLimits) -> Self {
        let h = s.h;
        let half = h / 2.0;
        let (mut zt, mut zb) = (half, -half);
        let (mut has_top, mut has_bot) = (false, false);
        for r in &s.rows {
            if r.z > 0.0 {
                if !has_top || r.z > zt {
                    zt = r.z;
                }
                has_top = true;
            } else {
                if !has_bot || r.z < zb {
                    zb = r.z;
                }
                has_bot = true;
            }
        }
        if !has_top {
            zt = half;
        }
        if !has_bot {
            zb = -half;
        }
        let (eu, b2, b0) = (lim.eps_s_ult, lim.eps_b2, lim.eps_b0);
        let (u_top, u_bot) = through(zb, eu, zt, eu, half);
        let mut c = Self { h, half, zt, zb, eu, b2, b0, u_top, u_bot, x_start_top: 0.0, x_start_bot: 0.0 };
        let (_, b_start_bot) = c.a_bot(1.0);
        c.x_start_top = h * b2 / (b2 + b_start_bot);
        let (a_top_start_top, _) = c.a_top(0.0);
        c.x_start_bot = h * b2 / (b2 + a_top_start_top);
        c
    }

    fn a_bot(&self, t: f64) -> (f64, f64) {
        let et = if t < 0.5 { self.u_top * (1.0 - 2.0 * t) } else { -self.b2 * (2.0 * t - 1.0) };
        through(self.zb, self.eu, self.half, et, self.half)
    }

    fn a_top(&self, t: f64) -> (f64, f64) {
        let eb = if t < 0.5 { -self.b2 * (1.0 - 2.0 * t) } else { self.u_bot * (2.0 * t - 1.0) };
        through(self.zt, self.eu, -self.half, eb, self.half)
    }

    fn c_top(&self, t: f64) -> (f64, f64) {
        let (b2, b0) = (self.b2, self.b0);
        let a = b2 - (b2 - b0) * t;
        (-a, -a * (b2 - a) / (b2 - b0))
    }

    fn plane(&self, kind: SegKind, t: f64) -> (f64, f64) {
        let (h, b2) = (self.h, self.b2);
        match kind {
            SegKind::ABot => self.a_bot(t),
            SegKind::BTop => {
                let x = self.x_start_top + (h - self.x_start_top) * t;
                (-b2, -b2 + b2 * h / x)
            }
            SegKind::CTop => self.c_top(t),
            SegKind::CBot => {
                let (x, y) = self.c_top(1.0 - t);
                (y, x)
            }
            SegKind::BBot => {
                let x = h - (h - self.x_start_bot) * t;
                (-b2 + b2 * h / x, -b2)
            }
            SegKind::ATop => self.a_top(t),
        }
    }
}

fn plane_state(top: f64, bot: f64, h: f64, f: &StripResponse) -> StripState {
    StripState {
        eps0: (top + bot) / 2.0,
        kappa: (top - bot) / h,
        converged: true,
        iterations: 0,
        reason: String::new(),
        n: f.n,
        m: f.m,
    }
}

/// ПС1: предельный множитель λ на луче (N, M) и Кисп = 1/λ. Пересечение контура с лучом
/// ищется по углу в нормированных координатах (N·h/6, M); при нескольких — наименьший λ.
pub fn ultimate_factor(s: &StripSection, n: f64, m: f64, lim: &StripLimits) -> StripUls {
    ultimate_factor_samples(s, n, m, lim, 24)
}

pub fn ultimate_factor_samples(s: &StripSection, n: f64, m: f64, lim: &StripLimits, samples: usize) -> StripUls {
    use std::f64::consts::PI;
    let h = s.h;
    let sc = h / 6.0;
    let (tu, tv) = (n * sc, m);
    let tr = (tu * tu + tv * tv).sqrt();
    if tr < 1e-12 {
        return StripUls {
            lambda: f64::INFINITY,
            eta: 0.0,
            pivot: Pivot::None,
            top_face: false,
            ultimate: None,
            reason: "нулевая нагрузка".into(),
        };
    }
    let phi = tv.atan2(tu);
    let contour = Contour::new(s, lim);

    let angle = |kind: SegKind, t: f64| -> (f64, f64, StripState) {
        let (top, bot) = contour.plane(kind, t);
        let f = s.integrate((top + bot) / 2.0, (top - bot) / h, false);
        let st = plane_state(top, bot, h, &f);
        let (u, v) = (f.n * sc, f.m);
        let radius = (u * u + v * v).sqrt();
        let mut d = v.atan2(u) - phi;
        while d <= -PI {
            d += 2.0 * PI;
        }
        while d > PI {
            d -= 2.0 * PI;
        }
        (d, radius, st)
    };

    let mut best_lambda = f64::INFINITY;
    let mut best_pivot = Pivot::None;
    let mut best_top = false;
    let mut best_state: Option<StripState> = None;

    for &(kind, pivot, top_face) in &SEGMENTS {
        let mut t_prev = 0.0;
        let (mut d_prev, mut r_prev, _) = angle(kind, 0.0);
        for i in 1..=samples {
            let t = i as f64 / samples as f64;
            let (d, r, _) = angle(kind, t);
            let cross = ((d_prev <= 0.0 && d >= 0.0) || (d_prev >= 0.0 && d <= 0.0))
                && (d - d_prev).abs() < PI
                && (r_prev > 1e-12 || r > 1e-12);
            if cross {
                let (mut a, mut b, mut da) = (t_prev, t, d_prev);
                let mut tm;
                let mut found: Option<(f64, StripState)> = None;
                for _ in 0..50 {
                    tm = 0.5 * (a + b);
                    let (dm, rm, stm) = angle(kind, tm);
                    found = Some((rm, stm));
                    if dm.abs() < 1e-12 {
                        break;
                    }
                    if dm.signum() == da.signum() {
                        a = tm;
                        da = dm;
                    } else {
                        b = tm;
                    }
                }
                let (rm, stm) = found.expect("бисекция выполняет хотя бы один шаг");
                let lambda = rm / tr;
                if lambda < best_lambda {
                    best_lambda = lambda;
                    best_pivot = pivot;
                    best_top = top_face;
                    best_state = Some(stm);
                }
            }
            t_prev = t;
            d_prev = d;
            r_prev = r;
        }
    }

    match best_state {
        None => StripUls {
            lambda: 0.0,
            eta: f64::MAX,
            pivot: Pivot::None,
            top_face: false,
            ultimate: None,
            reason: "луч усилий не пересекает предельную границу полосы".into(),
        },
        Some(st) => StripUls {
            lambda: best_lambda,
            eta: if best_lambda > 1e-12 { 1.0 / best_lambda } else { f64::MAX },
            pivot: best_pivot,
            top_face: best_top,
            ultimate: Some(st),
            reason: String::new(),
        },
    }
}

/// Предельный контур полосы, затабулированный по плоскостям: бетонные усилия и напряжения
/// рядов в каждой точке. Для заданных площадей контур — замкнутая ломаная в координатах
/// (N·h/6, M); Кисп = 1/λ по наименьшему пересечению луча с ломаной.
#[derive(Clone, Debug)]
pub struct UlsTable {
    pub h: f64,
    sc: f64,
    z: Vec<f64>,
    nc: Vec<f64>,
    mc: Vec<f64>,
    sig: Vec<f64>,
    pivots: Vec<Pivot>,
}

impl UlsTable {
    pub fn row_count(&self) -> usize {
        self.z.len()
    }

    /// Кисп и шарнир при площадях рядов `areas`; 0 — нулевая нагрузка, `f64::MAX` — луч не
    /// пересекает контур.
    pub fn eta_pivot(&self, areas: &[f64], n: f64, m: f64) -> (f64, Pivot) {
        let mut pivot = Pivot::None;
        let (tu, tv) = (n * self.sc, m);
        let tr = (tu * tu + tv * tv).sqrt();
        if tr < 1e-12 {
            return (0.0, pivot);
        }
        let (du, dv) = (tu / tr, tv / tr);
        if self.z.len() == 2 && areas.len() == 2 {
            return self.eta_two_rows(areas[0], areas[1], du, dv, tr);
        }
        let n_rows = self.z.len();
        let count = self.nc.len();
        let mut best = f64::INFINITY;
        let (mut pu0, mut pv0) = (0.0, 0.0);
        for k in 0..=count {
            let i = k % count;
            let (mut nn, mut mm) = (self.nc[i], self.mc[i]);
            for r in 0..n_rows {
                let f = self.sig[i * n_rows + r] * areas[r];
                nn += f;
                mm += f * self.z[r];
            }
            let (u, v) = (nn * self.sc, mm);
            if k == 0 {
                pu0 = u;
                pv0 = v;
                continue;
            }
            let (eu, ev) = (u - pu0, v - pv0);
            let cde = du * ev - dv * eu;
            if cde.abs() > 1e-300 {
                let s = -(du * pv0 - dv * pu0) / cde;
                if (0.0..=1.0).contains(&s) {
                    let t = (pu0 + s * eu) * du + (pv0 + s * ev) * dv;
                    if t > 0.0 && t < best {
                        best = t;
                        pivot = self.pivots[(k - 1) % count];
                    }
                }
            }
            pu0 = u;
            pv0 = v;
        }
        if !best.is_finite() {
            return (f64::MAX, pivot);
        }
        let lambda = best / tr;
        (if lambda > 1e-12 { 1.0 / lambda } else { f64::MAX }, pivot)
    }

    /// То же для двух рядов: те же операции в том же порядке, без `%` и проверок границ.
    fn eta_two_rows(&self, a0: f64, a1: f64, du: f64, dv: f64, tr: f64) -> (f64, Pivot) {
        let (z0, z1) = (self.z[0], self.z[1]);
        let sc = self.sc;
        let count = self.nc.len();
        let point = |i: usize| -> (f64, f64) {
            let f0 = self.sig[2 * i] * a0;
            let f1 = self.sig[2 * i + 1] * a1;
            let mut nn = self.nc[i];
            let mut mm = self.mc[i];
            nn += f0;
            mm += f0 * z0;
            nn += f1;
            mm += f1 * z1;
            (nn * sc, mm)
        };
        let (mut pu0, mut pv0) = point(0);
        let mut best = f64::INFINITY;
        let mut best_k = usize::MAX;
        for k in 1..=count {
            let (u, v) = if k == count { point(0) } else { point(k) };
            let (eu, ev) = (u - pu0, v - pv0);
            let cde = du * ev - dv * eu;
            if cde.abs() > 1e-300 {
                let s = -(du * pv0 - dv * pu0) / cde;
                if (0.0..=1.0).contains(&s) {
                    let t = (pu0 + s * eu) * du + (pv0 + s * ev) * dv;
                    if t > 0.0 && t < best {
                        best = t;
                        best_k = k - 1;
                    }
                }
            }
            pu0 = u;
            pv0 = v;
        }
        let pivot = if best_k == usize::MAX { Pivot::None } else { self.pivots[best_k] };
        if !best.is_finite() {
            return (f64::MAX, pivot);
        }
        let lambda = best / tr;
        (if lambda > 1e-12 { 1.0 / lambda } else { f64::MAX }, pivot)
    }

    pub fn eta(&self, areas: &[f64], n: f64, m: f64) -> f64 {
        self.eta_pivot(areas, n, m).0
    }
}

/// Таблица предельного контура для многократной оценки Кисп при разных площадях рядов:
/// плоскости контура от площадей не зависят, арматура входит линейно.
pub fn tabulate(geometry: &StripSection, lim: &StripLimits) -> UlsTable {
    tabulate_samples(geometry, lim, 128)
}

pub fn tabulate_samples(geometry: &StripSection, lim: &StripLimits, samples: usize) -> UlsTable {
    let contour = Contour::new(geometry, lim);
    let zero = geometry.with_areas(&vec![0.0; geometry.rows.len()]);
    let h = geometry.h;
    let n_rows = geometry.rows.len();
    let count = SEGMENTS.len() * samples;
    let mut nc = Vec::with_capacity(count);
    let mut mc = Vec::with_capacity(count);
    let mut sig = Vec::with_capacity(count * n_rows);
    let mut pivots = Vec::with_capacity(count);
    for &(kind, pivot, _) in &SEGMENTS {
        for i in 0..samples {
            let (top, bot) = contour.plane(kind, i as f64 / samples as f64);
            let (e0, kappa) = ((top + bot) / 2.0, (top - bot) / h);
            let f = zero.integrate(e0, kappa, false);
            nc.push(f.n);
            mc.push(f.m);
            pivots.push(pivot);
            for r in &geometry.rows {
                sig.push(geometry.rebar.stress(e0 + kappa * r.z, true, true));
            }
        }
    }
    UlsTable {
        h,
        sc: h / 6.0,
        z: geometry.rows.iter().map(|r| r.z).collect(),
        nc,
        mc,
        sig,
        pivots,
    }
}

/// Высота растянутой зоны от грани (верх/низ) по деформациям граней.
fn tension_zone(top: f64, bot: f64, from_top: bool, h: f64) -> f64 {
    let (face, other) = if from_top { (top, bot) } else { (bot, top) };
    if face <= 0.0 {
        return 0.0;
    }
    if other >= 0.0 {
        return h;
    }
    h * face / (face - other)
}

/// Параметры ширины раскрытия трещин (п. 8.2.15).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrackParams {
    /// Модуль упругости арматуры, кПа.
    pub es: f64,
    /// φ1: 1,0 — непродолжительное, 1,4 — продолжительное действие.
    pub phi1: f64,
    /// φ2: 0,5 — периодический профиль.
    pub phi2: f64,
    pub method: SigmaSCrcMethod,
}

/// ПС2 по деформационной модели: M_crc (п. 8.2.14) — плоскость с растянутой гранью на
/// εbt,ult при той же N; σs — из НДС сечения с трещиной (п. 8.2.16); σs,crc — по
/// [`SigmaSCrcMethod`]; ψs, l_s, φ3, a_crc — пп. 8.2.15–8.2.17. Сквозное растяжение даёт a_crc
/// для каждого растянутого ряда.
///
/// `s` — полоса с диаграммами для ПС2; `row_diameters` — диаметры стержней рядов, м.
pub fn crack_width(s: &StripSection, n: f64, m: f64, row_diameters: &[f64], p: &CrackParams) -> StripSls {
    let h = s.h;
    let half = h / 2.0;
    let ebt = s.concrete.tension.last_x();

    let sign = if m > 0.0 {
        1.0
    } else if m < 0.0 {
        -1.0
    } else {
        let z0 = solve(s, n, 0.0, true);
        if !z0.converged || z0.eps_at(half) >= z0.eps_at(-half) {
            1.0
        } else {
            -1.0
        }
    };
    let face_top = sign > 0.0;

    // Семейство плоскостей с гранью на εbt,ult: ε₀ = εbt − q·h/2, κ = sign·q, q ≥ 0.
    let at = |q: f64| s.integrate(ebt - q * half, sign * q, true);

    let mcrc;
    let (plane_top, plane_bot);
    let mut cracked_at_zero = false;
    let mut reason = String::new();
    let f0 = at(0.0);
    if n >= f0.n {
        cracked_at_zero = true;
        mcrc = 0.0;
        plane_top = ebt;
        plane_bot = ebt;
    } else {
        let mut lo = 0.0;
        let mut hi = ebt / h;
        let mut grow = 0;
        while at(hi).n > n && grow < 60 {
            lo = hi;
            hi *= 2.0;
            grow += 1;
        }
        if at(hi).n > n {
            reason = "трещина не образуется: сечение исчерпывает сжатие раньше".into();
            plane_top = 0.0;
            plane_bot = 0.0;
            mcrc = f64::MAX;
        } else {
            let mut i = 0;
            while i < 100 && hi - lo > 1e-14 * f64::max(1.0, hi) {
                let mid = 0.5 * (lo + hi);
                if at(mid).n > n {
                    lo = mid;
                } else {
                    hi = mid;
                }
                i += 1;
            }
            let q = 0.5 * (lo + hi);
            let fq = at(q);
            let e0 = ebt - q * half;
            let k = sign * q;
            plane_top = e0 + k * half;
            plane_bot = e0 - k * half;
            if sign * fq.m <= 0.0 {
                cracked_at_zero = true;
                mcrc = 0.0;
            } else {
                mcrc = sign * fq.m;
            }
        }
    }
    let h_bt_state = tension_zone(plane_top, plane_bot, face_top, h);

    let cracked = if cracked_at_zero { n > 0.0 || sign * m > 0.0 } else { sign * m > mcrc };
    let empty_rows = || -> Vec<StripSlsRow> {
        s.rows.iter().map(|r| StripSlsRow { z: r.z, as_: r.as_, ..Default::default() }).collect()
    };
    if !cracked {
        return StripSls {
            mcrc: if mcrc == f64::MAX { 0.0 } else { mcrc },
            crack_face_top: face_top,
            cracked: false,
            h_bt: h_bt_state,
            rows: empty_rows(),
            acrc_mm: 0.0,
            state: None,
            converged: true,
            reason,
        };
    }

    let st = solve(s, n, m, false);
    if !st.converged {
        let r = format!("НДС сечения с трещиной не найдено: {}", st.reason);
        return StripSls {
            mcrc,
            crack_face_top: face_top,
            cracked: true,
            h_bt: h_bt_state,
            rows: empty_rows(),
            acrc_mm: 0.0,
            state: Some(st),
            converged: false,
            reason: r,
        };
    }

    let mut st_crc: Option<StripState> = None;
    if p.method == SigmaSCrcMethod::ReleasedConcrete8137 {
        let c = solve(s, n, sign * mcrc, false);
        if c.converged {
            st_crc = Some(c);
        } else {
            reason = "σs,crc — по (8.138): НДС при M = M_crc не найдено".into();
        }
    }
    let ratio = if m.abs() > 1e-9 { (mcrc / m.abs()).clamp(0.0, 1.0) } else { 0.0 };
    let phi3 = if n > 1e-3 { 1.2 } else { 1.0 };

    let mut worst: f64 = 0.0;
    let mut rows_out = Vec::with_capacity(s.rows.len());
    for (i, r) in s.rows.iter().enumerate() {
        let eps = st.eps_at(r.z);
        let blank = StripSlsRow { z: r.z, as_: r.as_, ..Default::default() };
        if r.as_ <= 0.0 || eps <= 0.0 {
            rows_out.push(blank);
            continue;
        }
        let sigma = s.rebar.stress(eps, true, true);
        if sigma < 1e-3 {
            rows_out.push(blank);
            continue;
        }
        let sigma_crc = match &st_crc {
            Some(c) => s.rebar.stress(c.eps_at(r.z), true, true).max(0.0).min(sigma),
            None => sigma * ratio,
        };
        let psi = if p.method == SigmaSCrcMethod::PsiOne {
            1.0
        } else {
            (1.0 - 0.8 * sigma_crc / sigma).clamp(0.1, 1.0)
        };

        // П. 8.2.17: высота растянутой зоны — по состоянию трещинообразования, от грани ряда;
        // 2a ≤ h_bt ≤ 0,5·h0.
        let row_top = r.z > 0.0;
        let a = half - r.z.abs();
        let h0 = half + r.z.abs();
        let xt = tension_zone(plane_top, plane_bot, row_top, h);
        let h_bt = xt.max(2.0 * a).min(h0 / 2.0);
        let ds = match row_diameters.get(i) {
            Some(&d) if d > 0.0 => d,
            _ => 0.012,
        };
        // l_s по (8.144): 10·ds ≥ 100 мм ≤ l_s ≤ 40·ds ≤ 400 мм. При ds < 2,5 мм границы
        // меняются местами — берётся нижняя (в запас); в OpenCS это исключение.
        let (ls_lo, ls_hi) = ((10.0 * ds).max(0.10), (40.0 * ds).min(0.40));
        let ls_raw = 0.5 * h_bt / r.as_ * ds;
        let ls = if ls_lo <= ls_hi { ls_raw.clamp(ls_lo, ls_hi) } else { ls_lo };
        let acrc = p.phi1 * p.phi2 * phi3 * psi * sigma / p.es * ls * 1000.0;
        worst = worst.max(acrc);
        rows_out.push(StripSlsRow { z: r.z, as_: r.as_, sigma_s: sigma, sigma_s_crc: sigma_crc, psi_s: psi, h_bt, ls_m: ls, phi3, acrc_mm: acrc });
    }
    StripSls {
        mcrc,
        crack_face_top: face_top,
        cracked: true,
        h_bt: h_bt_state,
        rows: rows_out,
        acrc_mm: worst,
        state: Some(st),
        converged: true,
        reason,
    }
}
