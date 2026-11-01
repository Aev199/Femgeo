//! Свайное поле: осадка одиночной сваи и взаимовлияние свай по СП 24.13330.2021, п. 7.4
//! (формулы (7.33)–(7.41), таблица 7.18), матрица податливости и жёсткости поля.
//!
//! Сверено с текстом СП 24.13330.2021 (2026-11-01): выбор формулы по κ = G1·l/(G2·d),
//! β с множителем 0,3 (7.35), λ1 (7.36), kν (7.37), короткие сваи (7.38), взаимовлияние
//! (7.40)–(7.41); G1, ν1 — по глубине погружения, G2, ν2 — в пределах 0,5·l под нижним концом.
//! Принято: ν в ζ0 и mν (таблица 7.18) — (ν1 + ν2)/2, как для kν. Единицы: кН, м, кПа.

pub mod layout;

/// Слой грунта: подошва слоя (отметка, м, ось z вверх), модуль деформации E и ν.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Layer {
    pub bottom: f64,
    pub e: f64,
    pub nu: f64,
}

/// Грунтовая колонка в точке: отметка поверхности и слои сверху вниз (последний продолжается
/// неограниченно вниз).
#[derive(Clone, Debug, PartialEq)]
pub struct SoilColumn {
    pub top: f64,
    pub layers: Vec<Layer>,
}

impl SoilColumn {
    /// Средние по толщине G = E/(2(1+ν)) и ν в интервале отметок [z_lo, z_hi].
    pub fn average(&self, z_lo: f64, z_hi: f64) -> (f64, f64) {
        let (mut g, mut nu, mut t) = (0.0, 0.0, 0.0);
        let mut upper = self.top;
        for (i, l) in self.layers.iter().enumerate() {
            let lower = if i + 1 == self.layers.len() { f64::NEG_INFINITY } else { l.bottom };
            let (a, b) = (lower.max(z_lo), upper.min(z_hi));
            if b > a {
                let h = b - a;
                g += l.e / (2.0 * (1.0 + l.nu)) * h;
                nu += l.nu * h;
                t += h;
            }
            upper = l.bottom;
        }
        if t <= 0.0 {
            return (f64::NAN, f64::NAN);
        }
        (g / t, nu / t)
    }
}

/// Свая: отметка головы, длина, диаметр (для квадратной — 1,13·b), жёсткость ствола E·A.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pile {
    pub x: f64,
    pub y: f64,
    pub head: f64,
    pub length: f64,
    pub d: f64,
    pub ea: f64,
}

/// Параметры одиночной сваи по п. 7.4.3.
#[derive(Clone, Debug, PartialEq)]
pub struct SinglePile {
    /// G1, ν1 — по длине сваи; G2, ν2 — под нижним концом (в пределах 0,5·l).
    pub g1: f64,
    pub nu1: f64,
    pub g2: f64,
    pub nu2: f64,
    /// κ = G1·l/(G2·d): ≥ 7,5 — длинная висячая свая (7.34), иначе короткая (7.38).
    pub kappa: f64,
    /// β по (7.35); NaN для короткой сваи.
    pub beta: f64,
    /// Податливость s/N, м/кН.
    pub compliance: f64,
    /// Нарушенные условия применимости п. 7.4.2 (l/d > 5, κ > 1).
    pub warnings: Vec<String>,
}

fn k_nu(nu: f64) -> f64 {
    2.82 - 3.78 * nu + 2.18 * nu * nu
}

/// Таблица 7.18: ν, ζ0, mν (линейная интерполяция по ν).
const T718_NU: [f64; 11] = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5];
const T718_ZETA0: [f64; 11] = [0.455, 0.437, 0.419, 0.4, 0.38, 0.361, 0.34, 0.319, 0.297, 0.274, 0.25];
const T718_M: [f64; 11] = [1.345, 1.373, 1.405, 1.446, 1.491, 1.54, 1.607, 1.685, 1.786, 1.916, 2.01];

fn table_718(row: &[f64; 11], nu: f64) -> f64 {
    let nu = nu.clamp(0.0, 0.5);
    let i = ((nu / 0.05).floor() as usize).min(9);
    let t = (nu - T718_NU[i]) / 0.05;
    row[i] + t * (row[i + 1] - row[i])
}

impl SinglePile {
    /// Осадка одиночной сваи без уширения (п. 7.4.2), κ = G1·l/(G2·d):
    /// κ ≥ 7,5 — s = β·N/(G1·l), β = β′/λ1 + 0,3·(1 − β′/α′)/χ, β′ = 0,17·ln(kν·G1·l/(G2·d)),
    /// α′ = 0,17·ln(kν1·l/d), χ = E·A/(G1·l²), λ1 = 2,12·χ^¾/(1 + 2,12·χ^¾);
    /// κ < 7,5 — s = ζ′·N/(G2·d), ζ′ = ζ0/(1 + κ/mν) (таблица 7.18).
    pub fn new(p: &Pile, soil: &SoilColumn) -> Self {
        let tip = p.head - p.length;
        let (g1, nu1) = soil.average(tip, p.head.min(soil.top));
        let (g2, nu2) = soil.average(tip - 0.5 * p.length, tip);
        let l = p.length;
        let nu = 0.5 * (nu1 + nu2);
        let kappa = g1 * l / (g2 * p.d);
        let (beta, compliance) = if kappa >= 7.5 {
            let beta_r = 0.17 * (k_nu(nu) * g1 * l / (g2 * p.d)).ln();
            let alpha_r = 0.17 * (k_nu(nu1) * l / p.d).ln();
            let chi = p.ea / (g1 * l * l);
            let c = 2.12 * chi.powf(0.75);
            let lambda1 = c / (1.0 + c);
            let beta = beta_r / lambda1 + 0.3 * (1.0 - beta_r / alpha_r) / chi;
            (beta, beta / (g1 * l))
        } else {
            let zeta = table_718(&T718_ZETA0, nu) / (1.0 + kappa / table_718(&T718_M, nu));
            (f64::NAN, zeta / (g2 * p.d))
        };
        let mut warnings = vec![];
        if l / p.d <= 5.0 {
            warnings.push(format!("l/d = {:.2} ≤ 5: п. 7.4.2 неприменим", l / p.d));
        }
        if kappa <= 1.0 {
            warnings.push("G1·l/(G2·d) ≤ 1: п. 7.4.2 неприменим".into());
        }
        Self { g1, nu1, g2, nu2, kappa, beta, compliance, warnings }
    }
}

/// Коэффициент дополнительной осадки от соседней сваи на расстоянии `a` (п. 7.4.4):
/// δ = 0,17·ln(kν·G1·l/(2·G2·a)), 0 — при аргументе ≤ 1.
pub fn interaction(sp: &SinglePile, l: f64, a: f64) -> f64 {
    let arg = k_nu(0.5 * (sp.nu1 + sp.nu2)) * sp.g1 * l / (2.0 * sp.g2 * a);
    if arg > 1.0 { 0.17 * arg.ln() } else { 0.0 }
}

/// Поле свай: параметры каждой сваи и матрица податливости F (м/кН, n×n, построчно):
/// F_ii = β_i/(G1_i·l_i), F_ij = δ_ij/(G1_j·l_j) — осадка сваи i от единичной силы на свае j;
/// F симметризуется (средним), как и подразумевает п. 7.4.4 для свай одного типа.
pub struct Field {
    pub piles: Vec<SinglePile>,
    pub flexibility: Vec<f64>,
}

impl Field {
    pub fn new(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn) -> Self {
        let sp: Vec<SinglePile> = piles.iter().map(|p| SinglePile::new(p, &soil(p.x, p.y))).collect();
        let n = piles.len();
        let mut f = vec![0.0; n * n];
        for i in 0..n {
            f[i * n + i] = sp[i].compliance;
            for j in 0..n {
                if i != j {
                    let a = ((piles[i].x - piles[j].x).powi(2) + (piles[i].y - piles[j].y).powi(2)).sqrt();
                    f[i * n + j] = interaction(&sp[j], piles[j].length, a) / (sp[j].g1 * piles[j].length);
                }
            }
        }
        for i in 0..n {
            for j in i + 1..n {
                let m = 0.5 * (f[i * n + j] + f[j * n + i]);
                f[i * n + j] = m;
                f[j * n + i] = m;
            }
        }
        Self { piles: sp, flexibility: f }
    }

    /// Матрица жёсткости поля K = F⁻¹ (кН/м). Ошибка — F не положительно определена.
    pub fn stiffness(&self) -> Result<Vec<f64>, String> {
        let n = self.piles.len();
        // Холецкий F = L·Lᵀ, затем K = L⁻ᵀ·L⁻¹.
        let mut l = self.flexibility.clone();
        for j in 0..n {
            let mut d = l[j * n + j];
            for k in 0..j {
                d -= l[j * n + k] * l[j * n + k];
            }
            if !(d > 0.0) {
                return Err("матрица податливости свайного поля не положительно определена".into());
            }
            let d = d.sqrt();
            l[j * n + j] = d;
            for i in j + 1..n {
                let mut s = l[i * n + j];
                for k in 0..j {
                    s -= l[i * n + k] * l[j * n + k];
                }
                l[i * n + j] = s / d;
            }
        }
        let mut k = vec![0.0; n * n];
        let mut col = vec![0.0; n];
        for c in 0..n {
            // L·y = e_c, Lᵀ·x = y.
            for i in 0..n {
                let mut s = if i == c { 1.0 } else { 0.0 };
                for p in 0..i {
                    s -= l[i * n + p] * col[p];
                }
                col[i] = s / l[i * n + i];
            }
            for i in (0..n).rev() {
                let mut s = col[i];
                for p in i + 1..n {
                    s -= l[p * n + i] * col[p];
                }
                col[i] = s / l[i * n + i];
            }
            for i in 0..n {
                k[i * n + c] = col[i];
            }
        }
        Ok(k)
    }
}
