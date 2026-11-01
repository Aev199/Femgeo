//! Подбор площадей полосы по НДМ: ПС1 — минимальная сумма As_t + As_c при Кисп ≤ 1,
//! ПС2 — As_t при a_crc ≤ a_crc,ult. Точность — 0,5 % площади.
//!
//! Источник: OpenCS `CScore/ShellSimpl/StripNdmDesign.cs`, коммит `387e157`.

use super::diagram::Diagram;
use super::materials::StripLimits;
use super::strip::{self, CrackParams, StripRow, StripSection, StripSls, UlsTable};

/// Исход подбора площадей полосы.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DesignStatus {
    /// Арматура не нужна: бетон воспринимает усилия (ПС1) или трещины нет (ПС2).
    NotRequired,
    /// Достаточно ряда у растянутой моментом грани.
    TensionOnly,
    /// Нужны оба ряда (сжатая арматура, сквозное растяжение).
    WithCompression,
    /// Нет решения и при наибольшей площади — недостаточна толщина.
    Impossible,
}

/// Требуемые площади полосы, м²/м.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignResult {
    /// Ряд у грани, растянутой моментом (при M ≥ 0 — верхний).
    pub as_tension: f64,
    /// Ряд у противоположной грани.
    pub as_compression: f64,
    pub tension_top: bool,
    pub status: DesignStatus,
    /// ПС1 — Кисп при найденных площадях, ПС2 — a_crc, мм.
    pub value: f64,
    pub reason: String,
}

impl DesignResult {
    pub fn top(&self) -> f64 {
        if self.tension_top { self.as_tension } else { self.as_compression }
    }
    pub fn bot(&self) -> f64 {
        if self.tension_top { self.as_compression } else { self.as_tension }
    }
    fn new(at: f64, ac: f64, tension_top: bool, status: DesignStatus, value: f64) -> Self {
        Self { as_tension: at, as_compression: ac, tension_top, status, value, reason: String::new() }
    }
}

const REL_TOL: f64 = 0.005;

/// Наибольшая рассматриваемая площадь ряда — 0,1·h·1 м.
pub fn max_area(h: f64) -> f64 {
    0.1 * h
}

/// Таблица предельного контура полосы подбора: ряды [верх, низ].
pub fn uls_table(h: f64, z_top: f64, z_bot: f64, concrete: &Diagram, rebar: &Diagram, lim: &StripLimits) -> UlsTable {
    let s = StripSection::new(
        h,
        vec![StripRow { z: z_top, as_: 0.0 }, StripRow { z: z_bot, as_: 0.0 }],
        concrete.clone(),
        rebar.clone(),
    );
    strip::tabulate(&s, lim)
}

/// Наименьшая x ∈ (0, as_max], при которой ok(x); None — нет и при as_max.
fn bisect(mut ok: impl FnMut(f64) -> bool, as_max: f64) -> Option<f64> {
    if !ok(as_max) {
        return None;
    }
    let (mut lo, mut hi) = (0.0, as_max);
    while hi - lo > REL_TOL * hi {
        let mid = 0.5 * (lo + hi);
        if ok(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// ПС1 по таблице контура (ряды [верх, низ]): As_c = 0 и бисекция по As_t; не хватает —
/// наименьшая As_c, при которой решение есть, и золотое сечение по As_c для минимума суммы.
pub fn required_uls(n: f64, m: f64, table: &UlsTable) -> DesignResult {
    let tension_top = m >= 0.0;
    let as_max = max_area(table.h);
    let eta = |at: f64, ac: f64| -> f64 {
        let areas = if tension_top { [at, ac] } else { [ac, at] };
        table.eta(&areas, n, m)
    };

    let eta0 = eta(0.0, 0.0);
    if eta0 <= 1.0 {
        return DesignResult::new(0.0, 0.0, tension_top, DesignStatus::NotRequired, eta0);
    }

    if let Some(a1) = bisect(|at| eta(at, 0.0) <= 1.0, as_max) {
        return DesignResult::new(a1, 0.0, tension_top, DesignStatus::TensionOnly, eta(a1, 0.0));
    }

    let Some(ac_min) = bisect(|ac| eta(as_max, ac) <= 1.0, as_max) else {
        let mut r = DesignResult::new(as_max, as_max, tension_top, DesignStatus::Impossible, eta(as_max, as_max));
        r.reason = "сечение не проходит: прочность не обеспечена и при наибольшей площади обоих рядов — недостаточна толщина".into();
        return r;
    };

    let f = |ac: f64| -> (f64, f64) {
        match bisect(|x| eta(x, ac) <= 1.0, as_max) {
            Some(v) => (ac + v, v),
            None => (f64::INFINITY, as_max),
        }
    };

    let mut best = ac_min;
    let mut best_f = f(best);
    let consider = |ac: f64, fv: (f64, f64), best: &mut f64, best_f: &mut (f64, f64)| {
        if fv.0 < best_f.0 {
            *best = ac;
            *best_f = fv;
        }
    };

    const GR: f64 = 0.6180339887498949;
    let (mut a, mut b) = (ac_min, as_max);
    let mut c = b - GR * (b - a);
    let mut d = a + GR * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    consider(c, fc, &mut best, &mut best_f);
    consider(d, fd, &mut best, &mut best_f);
    while b - a > REL_TOL * best_f.0.max(1e-6) {
        if fc.0 <= fd.0 {
            b = d;
            d = c;
            fd = fc;
            c = b - GR * (b - a);
            fc = f(c);
            consider(c, fc, &mut best, &mut best_f);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + GR * (b - a);
            fd = f(d);
            consider(d, fd, &mut best, &mut best_f);
        }
    }
    DesignResult::new(best_f.1, best, tension_top, DesignStatus::WithCompression, eta(best_f.1, best))
}

/// ПС1 при заданной площади ряда у грани, сжатой моментом: наименьшая As_t при Кисп ≤ 1.
/// Если не хватает и наибольшей As_t — подбор обоих рядов ([`required_uls`]).
pub fn required_tension_uls(n: f64, m: f64, table: &UlsTable, as_other: f64) -> DesignResult {
    let tension_top = m >= 0.0;
    let as_max = max_area(table.h);
    let eta = |at: f64| -> f64 {
        let areas = if tension_top { [at, as_other] } else { [as_other, at] };
        table.eta(&areas, n, m)
    };
    let eta0 = eta(0.0);
    if eta0 <= 1.0 {
        return DesignResult::new(0.0, as_other, tension_top, DesignStatus::NotRequired, eta0);
    }
    if eta(as_max) > 1.0 {
        return required_uls(n, m, table);
    }
    let (mut lo, mut hi) = (0.0, as_max);
    while hi - lo > REL_TOL * hi {
        let mid = 0.5 * (lo + hi);
        if eta(mid) <= 1.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    DesignResult::new(hi, as_other, tension_top, DesignStatus::TensionOnly, eta(hi))
}

/// Полоса и параметры подбора по ПС2.
#[derive(Clone, Debug)]
pub struct SlsStrip<'a> {
    pub h: f64,
    pub z_top: f64,
    pub z_bot: f64,
    /// Диаграмма бетона вида расчёта для ПС2 (п. 6.1.26).
    pub concrete: &'a Diagram,
    pub rebar: &'a Diagram,
    pub ds_top: f64,
    pub ds_bot: f64,
    pub crack: CrackParams,
}

/// ПС2: бисекция по As_t до a_crc ≤ a_crc,ult. При сквозном растяжении второй ряд
/// наращивается пропорционально своей доле усилия, но не меньше `as_compression`.
pub fn required_sls(n: f64, m: f64, p: &SlsStrip, as_compression: f64, acrc_lim_mm: f64) -> DesignResult {
    let tension_top = m >= 0.0;
    let (zt, zc) = if tension_top { (p.z_top, p.z_bot) } else { (p.z_bot, p.z_top) };
    let ds = if tension_top { [p.ds_top, p.ds_bot] } else { [p.ds_bot, p.ds_top] };
    let sec = StripSection::new(
        p.h,
        vec![StripRow { z: zt, as_: 0.0 }, StripRow { z: zc, as_: 0.0 }],
        p.concrete.clone(),
        p.rebar.clone(),
    );
    let as_max = max_area(p.h);

    let ft = (m - n * zc) / (zt - zc);
    let fcomp = n - ft;
    let share = if ft > 0.0 && fcomp > 0.0 { fcomp / ft } else { 0.0 };
    let ac = |at: f64| as_compression.max(at * share);
    let sls = |at: f64| -> StripSls { strip::crack_width(&sec.with_areas(&[at, ac(at)]), n, m, &ds, &p.crack) };

    const TINY: f64 = 1e-7;
    let first = sls(TINY);
    if !first.cracked {
        return DesignResult::new(0.0, as_compression, tension_top, DesignStatus::NotRequired, 0.0);
    }
    let ok = |at: f64| {
        let r = sls(at);
        r.converged && r.acrc_mm <= acrc_lim_mm
    };
    if !ok(as_max) {
        let mut r = DesignResult::new(as_max, ac(as_max), tension_top, DesignStatus::Impossible, sls(as_max).acrc_mm);
        r.reason = "ширина раскрытия трещин не обеспечивается и при наибольшей площади".into();
        return r;
    }
    let (mut lo, mut hi) = (TINY, as_max);
    while hi - lo > REL_TOL * hi {
        let mid = 0.5 * (lo + hi);
        if ok(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let a_c = ac(hi);
    let status = if a_c > as_compression + 1e-12 { DesignStatus::WithCompression } else { DesignStatus::TensionOnly };
    DesignResult::new(hi, a_c, tension_top, status, sls(hi).acrc_mm)
}
