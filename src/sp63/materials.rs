//! Характеристики бетона и арматуры, диаграммы СП 63 и предельные деформации.
//!
//! Построение диаграмм — по OpenCS `CScore/MaterialChars.cs` (`D2L`, `D3L`,
//! `TensionTrilinear`, срезка сжатой ветви арматуры по Rsc) и
//! `CScore/ShellSimpl/StripNdmLimits.cs`, коммит `387e157`.
//!
//! Все величины — положительные модули, кПа и безразмерные деформации; знаки ставятся при
//! построении диаграмм.

use super::diagram::{Branch, Diagram, MaterialKind};

/// Вид расчёта (п. 6.1 СП 63). Первая версия — только непродолжительное действие нагрузки.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CalcType {
    /// Первая группа предельных состояний, непродолжительное действие (расчётные значения).
    Uls,
    /// Вторая группа, непродолжительное действие (нормативные значения, п. 6.1.26).
    Sls,
}

/// Тип диаграммы бетона.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConcreteDiagramType {
    /// Двухлинейная (п. 6.1.20, εb1,red).
    L2,
    /// Трёхлинейная (п. 6.1.20, σb1 = 0,6·Rb).
    L3,
}

/// Характеристики бетона для одного вида расчёта.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcreteChars {
    /// Начальный модуль упругости Eb, кПа.
    pub eb: f64,
    /// Сопротивление сжатию (Rb или Rb,n), кПа.
    pub rb: f64,
    /// Сопротивление растяжению (Rbt или Rbt,n), кПа.
    pub rbt: f64,
    /// εb0.
    pub eps_b0: f64,
    /// εb2.
    pub eps_b2: f64,
    /// εb1,red (двухлинейная диаграмма).
    pub eps_b1_red: f64,
    /// εbt0.
    pub eps_bt0: f64,
    /// εbt2.
    pub eps_bt2: f64,
    /// εbt1,red (двухлинейная диаграмма).
    pub eps_bt1_red: f64,
}

impl ConcreteChars {
    /// Тяжёлый бетон при непродолжительном действии нагрузки: деформации по п. 6.1.20 СП 63
    /// (εb0 = 0,002, εb2 = 0,0035, εb1,red = 0,0015, εbt0 = 0,0001, εbt2 = 0,00015,
    /// εbt1,red = 0,00008).
    pub fn short_term(eb: f64, rb: f64, rbt: f64) -> Self {
        Self {
            eb,
            rb,
            rbt,
            eps_b0: 0.002,
            eps_b2: 0.0035,
            eps_b1_red: 0.0015,
            eps_bt0: 0.0001,
            eps_bt2: 0.00015,
            eps_bt1_red: 0.00008,
        }
    }

    /// Диаграмма бетона указанного типа.
    pub fn diagram(&self, kind: ConcreteDiagramType) -> Diagram {
        let (rb, rbt, e) = (self.rb, self.rbt, self.eb);
        let (compression, tension) = match kind {
            ConcreteDiagramType::L2 => (
                Branch::new(&[(-self.eps_b2, -rb), (-self.eps_b1_red, -rb), (0.0, 0.0)]),
                Branch::new(&[(0.0, 0.0), (self.eps_bt1_red, rbt), (self.eps_bt2, rbt)]),
            ),
            ConcreteDiagramType::L3 => (
                Branch::new(&[
                    (-self.eps_b2, -rb),
                    (-self.eps_b0, -rb),
                    (-0.6 * rb / e, -0.6 * rb),
                    (0.0, 0.0),
                ]),
                Branch::new(&[
                    (0.0, 0.0),
                    (0.6 * rbt / e, 0.6 * rbt),
                    (self.eps_bt0, rbt),
                    (self.eps_bt2, rbt),
                ]),
            ),
        };
        Diagram {
            kind: MaterialKind::Concrete,
            compression: compression.expect("узлы диаграммы бетона различны"),
            tension: tension.expect("узлы диаграммы бетона различны"),
        }
    }
}

/// Характеристики арматуры с физическим пределом текучести для одного вида расчёта.
#[derive(Clone, Debug, PartialEq)]
pub struct RebarChars {
    /// Модуль упругости Es, кПа.
    pub es: f64,
    /// Сопротивление растяжению (Rs или Rs,n), кПа.
    pub rs: f64,
    /// Сопротивление сжатию Rsc, кПа (п. 6.2.13: сжатая ветвь срезается по Rsc).
    pub rsc: f64,
    /// Предельная деформация εs,ult (п. 8.1.30): 0,025 при физическом пределе текучести.
    pub eps_s_ult: f64,
}

impl RebarChars {
    /// Арматура с физическим пределом текучести: εs,ult = 0,025.
    pub fn physical_yield(es: f64, rs: f64, rsc: f64) -> Self {
        Self { es, rs, rsc, eps_s_ult: 0.025 }
    }

    /// Двухлинейная диаграмма (п. 6.2.14); сжатая ветвь — по Rs, срезанная по Rsc.
    pub fn diagram(&self) -> Diagram {
        let rs = self.rs;
        let rc = self.rsc.min(rs);
        let e2 = self.eps_s_ult;
        Diagram {
            kind: MaterialKind::Rebar,
            compression: Branch::new(&[(-e2, -rc), (-rc / self.es, -rc), (0.0, 0.0)])
                .expect("узлы диаграммы арматуры различны"),
            tension: Branch::new(&[(0.0, 0.0), (rs / self.es, rs), (e2, rs)])
                .expect("узлы диаграммы арматуры различны"),
        }
    }
}

/// Предельные деформации полосы (п. 8.1.30): εb0, εb2 и εs,ult.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StripLimits {
    pub eps_b0: f64,
    pub eps_b2: f64,
    pub eps_s_ult: f64,
}

impl StripLimits {
    pub fn new(concrete: &ConcreteChars, rebar: &RebarChars) -> Self {
        Self {
            eps_b0: concrete.eps_b0,
            eps_b2: concrete.eps_b2,
            eps_s_ult: rebar.eps_s_ult,
        }
    }
}

/// Класс тяжёлого бетона по прочности на сжатие.
///
/// Значения — СП 63.13330.2018, табл. 6.7 (нормативные Rb,n, Rbt,n), табл. 6.8 (расчётные
/// Rb, Rbt) и табл. 6.11 (Eb), МПа. **Требуют сверки с текстом нормы** перед выдачей
/// результатов (`docs/DEV_STATE.md`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConcreteClass {
    B15,
    B20,
    B25,
    B30,
    B35,
    B40,
    B45,
    B50,
    B55,
    B60,
}

impl ConcreteClass {
    /// (Rb, Rbt, Rb,n, Rbt,n, Eb), МПа.
    fn table(self) -> (f64, f64, f64, f64, f64) {
        use ConcreteClass::*;
        match self {
            B15 => (8.5, 0.75, 11.0, 1.10, 24_000.0),
            B20 => (11.5, 0.90, 15.0, 1.35, 27_500.0),
            B25 => (14.5, 1.05, 18.5, 1.55, 30_000.0),
            B30 => (17.0, 1.15, 22.0, 1.75, 32_500.0),
            B35 => (19.5, 1.30, 25.5, 1.95, 34_500.0),
            B40 => (22.0, 1.40, 29.0, 2.10, 36_000.0),
            B45 => (25.0, 1.50, 32.0, 2.25, 37_000.0),
            B50 => (27.5, 1.60, 36.0, 2.45, 38_000.0),
            B55 => (30.0, 1.70, 39.5, 2.60, 39_000.0),
            B60 => (33.0, 1.80, 43.0, 2.75, 39_500.0),
        }
    }

    /// Характеристики при непродолжительном действии нагрузки, γb1 = 1.
    pub fn chars(self, calc: CalcType) -> ConcreteChars {
        let (rb, rbt, rbn, rbtn, eb) = self.table();
        let k = 1000.0;
        match calc {
            CalcType::Uls => ConcreteChars::short_term(eb * k, rb * k, rbt * k),
            CalcType::Sls => ConcreteChars::short_term(eb * k, rbn * k, rbtn * k),
        }
    }
}

/// Класс арматуры с физическим пределом текучести.
///
/// Значения — СП 63.13330.2018, табл. 6.13 (Rs,n) и 6.14 (Rs, Rsc), Es = 200 000 МПа.
/// **Требуют сверки с текстом нормы** (в частности Rsc для А500).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RebarClass {
    A240,
    A400,
    A500,
}

impl RebarClass {
    /// (Rs, Rsc, Rs,n), МПа.
    fn table(self) -> (f64, f64, f64) {
        match self {
            RebarClass::A240 => (210.0, 210.0, 240.0),
            RebarClass::A400 => (350.0, 350.0, 400.0),
            RebarClass::A500 => (435.0, 400.0, 500.0),
        }
    }

    /// Rsw поперечной арматуры, кПа (табл. 6.15).
    pub fn rsw(self) -> f64 {
        1000.0
            * match self {
                RebarClass::A240 => 170.0,
                RebarClass::A400 => 280.0,
                RebarClass::A500 => 300.0,
            }
    }

    pub fn chars(self, calc: CalcType) -> RebarChars {
        let (rs, rsc, rsn) = self.table();
        let k = 1000.0;
        let es = 200_000.0 * k;
        match calc {
            CalcType::Uls => RebarChars::physical_yield(es, rs * k, rsc * k),
            CalcType::Sls => RebarChars::physical_yield(es, rsn * k, rsn * k),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concrete_l3_nodes() {
        let c = ConcreteChars::short_term(30e6, 14_500.0, 1_050.0);
        let d = c.diagram(ConcreteDiagramType::L3);
        assert_eq!(d.stress(-0.003, false, true), -14_500.0);
        let e1 = 0.6 * 14_500.0 / 30e6;
        assert!((d.stress(-e1 / 2.0, false, true) + 0.3 * 14_500.0).abs() < 1e-6);
        assert_eq!(d.stress(0.0001, false, true), 0.0);
        assert!((d.stress(0.0001, true, true) - 1_050.0).abs() < 1e-9);
        // За εbt2 растянутый бетон не работает.
        assert_eq!(d.stress(0.0002, true, true), 0.0);
    }

    #[test]
    fn rebar_compression_clipped_by_rsc() {
        let r = RebarChars::physical_yield(2e8, 435_000.0, 400_000.0);
        let d = r.diagram();
        assert_eq!(d.stress(0.01, true, true), 435_000.0);
        assert_eq!(d.stress(-0.01, true, true), -400_000.0);
        assert!((d.stress(-0.001, true, true) + 200_000.0).abs() < 1e-6);
        assert_eq!(d.stress(-0.01, true, false), 0.0);
    }
}
