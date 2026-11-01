//! Полоса плиты h = 200 мм, B25, A500, Ø12 шаг 200 у каждой грани, центр арматуры в 30 мм от
//! грани — фикстура OpenCS `CScore.Tests/ShellSimpl/StripNdmFixture.cs` (коммит `387e157`).
#![allow(dead_code)]

use femgeo::sp63::capra_maury::PlateSection;
use femgeo::sp63::diagram::Diagram;
use femgeo::sp63::materials::{CalcType, ConcreteChars, ConcreteDiagramType, RebarChars, StripLimits};
use femgeo::sp63::strip::{CrackParams, SigmaSCrcMethod, StripRow, StripSection};

pub const H: f64 = 0.2;
pub const Z: f64 = 0.07;
pub const DS: f64 = 0.012;
pub const ES: f64 = 200_000_000.0;
pub const EB: f64 = 30_000_000.0;
pub const RS: f64 = 435_000.0;

pub fn as_per_meter() -> f64 {
    std::f64::consts::PI * DS * DS / 4.0 / 0.2
}

pub fn concrete(ct: CalcType) -> ConcreteChars {
    match ct {
        CalcType::Uls => ConcreteChars::short_term(EB, 14_500.0, 1_050.0),
        CalcType::Sls => ConcreteChars::short_term(EB, 18_500.0, 1_550.0),
    }
}

/// В фикстуре OpenCS Rsc не задан — сжатая ветвь по Rs.
pub fn rebar(ct: CalcType) -> RebarChars {
    match ct {
        CalcType::Uls => RebarChars::physical_yield(ES, RS, RS),
        CalcType::Sls => RebarChars::physical_yield(ES, 500_000.0, 500_000.0),
    }
}

pub fn concrete_diagram(ct: CalcType, t: ConcreteDiagramType) -> Diagram {
    concrete(ct).diagram(t)
}

pub fn rebar_diagram(ct: CalcType) -> Diagram {
    rebar(ct).diagram()
}

pub fn strip_t(ct: CalcType, t: ConcreteDiagramType, as_top: f64, as_bot: f64) -> StripSection {
    StripSection::new(
        H,
        vec![StripRow { z: Z, as_: as_top }, StripRow { z: -Z, as_: as_bot }],
        concrete_diagram(ct, t),
        rebar_diagram(ct),
    )
}

pub fn strip(ct: CalcType) -> StripSection {
    strip_t(ct, ConcreteDiagramType::L3, as_per_meter(), as_per_meter())
}

pub fn strip_as(ct: CalcType, as_top: f64, as_bot: f64) -> StripSection {
    strip_t(ct, ConcreteDiagramType::L3, as_top, as_bot)
}

pub fn limits() -> StripLimits {
    StripLimits::new(&concrete(CalcType::Uls), &rebar(CalcType::Uls))
}

pub fn crack(method: SigmaSCrcMethod) -> CrackParams {
    CrackParams { es: ES, phi1: 1.0, phi2: 0.5, method }
}

pub fn plate() -> PlateSection {
    PlateSection {
        h: H,
        a_top: H / 2.0 - Z,
        a_bot: H / 2.0 - Z,
        ds_x: DS,
        ds_y: DS,
        concrete_diagram: ConcreteDiagramType::L3,
    }
}

pub fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Детерминированный генератор для «случайных» тестов (System.Random в C# не воспроизводим здесь).
pub struct Lcg(u64);

impl Lcg {
    pub fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407))
    }
    pub fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}
