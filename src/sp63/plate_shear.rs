//! Поперечная сила в пластине по направлению главной поперечной силы Q₀ = √(Qx² + Qy²)
//! (пп. 8.1.32–8.1.34 СП 63): бетонная полоса (8.55), бетон без хомутов Qb1 = 0,5·φn·Rbt·b·h0,
//! требуемая интенсивность поперечной арматуры по (8.60).
//!
//! Источник: OpenCS `CScore/ShellSimpl/PlateShearCheck.cs`, `CScore/Sp63Shear/PhiNCalculator.cs`,
//! `CScore/Sp63Shear/ShearFormulas.cs` (StripCapacity, MinConcreteShear), коммит `387e157`.

use super::capra_maury::{PlateForces, PlateSection};
use super::materials::{ConcreteChars, RebarChars};

/// φb1 в (8.55).
pub const PHI_B1: f64 = 0.3;
const INFINITE: f64 = 1e9;

/// Результат проверки на поперечную силу.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateShear {
    pub q0: f64,
    pub theta_deg: f64,
    /// N, M на фасетке θ_Q.
    pub n: f64,
    pub m: f64,
    pub h0: f64,
    pub phi_n: f64,
    pub phi_n_to_strip: bool,
    /// Несущая способность бетонной полосы (8.55), кН/м.
    pub q_strip: f64,
    /// Qb1 — бетон без поперечной арматуры, кН/м.
    pub qb1: f64,
    pub eta_strip: f64,
    pub eta_concrete: f64,
    /// Требуемая qsw, кН/м²; 0 — поперечная арматура не нужна.
    pub qsw_required: f64,
    pub qsw_min: f64,
    /// Требуемая Asw/sw, м²/м² (площадь хомутов на 1 м² плиты).
    pub asw_per_sw: f64,
}

impl PlateShear {
    pub fn eta(&self) -> f64 {
        self.eta_strip.max(self.eta_concrete)
    }
    pub fn stirrups_required(&self) -> bool {
        self.eta_concrete > 1.0
    }
}

/// φn по п. 8.1.34 для пластины («остальные элементы»): (значение, применяется к полосе).
pub fn phi_n(n: f64, h: f64, as_total: f64, c: &ConcreteChars, es: f64) -> (f64, bool) {
    if n.abs() < 1e-12 {
        return (1.0, false);
    }
    let compression = n < 0.0;
    let nu = if compression { c.rb / (c.eps_b0 * c.eb) } else { c.rbt / (c.eps_bt0 * c.eb) };
    let a_red = h + es / c.eb / nu * as_total;
    if !(a_red > 0.0) || !a_red.is_finite() {
        return (1.0, false);
    }
    let sigma = n.abs() / a_red;
    if compression {
        let ratio = sigma / c.rb;
        let v = if ratio <= 0.25 {
            1.0 + ratio
        } else if ratio <= 0.5 {
            1.25
        } else {
            (2.5 * (1.0 - ratio)).max(0.0)
        };
        (v, true)
    } else {
        ((1.0 - sigma / (2.0 * c.rbt)).max(0.0), false)
    }
}

fn ratio(applied: f64, capacity: f64) -> f64 {
    if capacity > 1e-12 { applied / capacity } else { INFINITE }
}

/// Проверка при поперечных силах `qx`, `qy` (кН/м). `as_n_total` — продольная арматура обеих
/// граней, приведённая к направлению Q₀, м²/м; `rsw` — расчётное сопротивление поперечной
/// арматуры (табл. 6.15), кПа. None — поперечной силы нет.
pub fn check(
    qx: f64,
    qy: f64,
    f: &PlateForces,
    section: &PlateSection,
    concrete: &ConcreteChars,
    rebar: &RebarChars,
    rsw: f64,
    as_n_total: f64,
) -> Option<PlateShear> {
    let q0 = (qx * qx + qy * qy).sqrt();
    if q0 <= 1e-9 {
        return None;
    }
    let theta = qy.atan2(qx);
    let (c, s) = (theta.cos(), theta.sin());
    let (c2, s2, cs) = (c * c, s * s, c * s);
    let n = f.nx * c2 + f.ny * s2 + 2.0 * f.nxy * cs;
    let m = f.mx * c2 + f.my * s2 + 2.0 * f.mxy * cs;
    let h = section.h;
    let (h0_top, h0_bot) = (h - section.a_top, h - section.a_bot);
    let eps = 1e-9;
    let h0 = if m > eps {
        h0_top
    } else if m < -eps {
        h0_bot
    } else {
        h0_top.min(h0_bot)
    };
    let (rb, rbt) = (concrete.rb, concrete.rbt);
    let (phi, to_strip) = phi_n(n, h, as_n_total, concrete, rebar.es);
    let q_strip = (if to_strip { phi } else { 1.0 }) * PHI_B1 * rb * h0;
    let qb1 = phi * (0.5 * rbt * h0).min(2.5 * rbt * h0);
    let (eta_strip, eta_concrete) = (ratio(q0, q_strip), ratio(q0, qb1));
    let qsw_min = 0.25 * rbt;
    let (mut qsw_required, mut asw_per_sw) = (0.0, 0.0);
    if eta_concrete > 1.0 {
        // (8.60): Q₀ ≤ Qb1 + qsw·h0; меньше 0,25·Rbt·b хомуты не учитываются.
        qsw_required = ((q0 - qb1) / h0).max(qsw_min);
        asw_per_sw = if rsw > 0.0 { qsw_required / rsw } else { INFINITE };
    }
    Some(PlateShear {
        q0,
        theta_deg: theta.to_degrees(),
        n,
        m,
        h0,
        phi_n: phi,
        phi_n_to_strip: to_strip,
        q_strip,
        qb1,
        eta_strip,
        eta_concrete,
        qsw_required,
        qsw_min,
        asw_per_sw,
    })
}
