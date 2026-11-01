//! Свайное поле: осадка одиночной сваи и взаимовлияние свай по СП 24.13330.2021, п. 7.4
//! (формулы пп. 7.4.2–7.4.4), матрица податливости и жёсткости поля для быстрой модели.
//!
//! **Формулы и правила осреднения параметров грунта требуют сверки с текстом СП 24**
//! (`docs/DEV_STATE.md`). Единицы: кН, м, кПа.

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
    pub beta: f64,
    /// Податливость s/N = β/(G1·l), м/кН.
    pub compliance: f64,
    /// Нарушенные условия применимости п. 7.4.2 (l/d > 5, G1·l/(G2·d) > 1).
    pub warnings: Vec<String>,
}

fn k_nu(nu: f64) -> f64 {
    2.82 - 3.78 * nu + 2.18 * nu * nu
}

impl SinglePile {
    /// Осадка одиночной сваи без уширения: s = β·N/(G1·l),
    /// β = β'/λ1 + (1 − β'/α')/χ, β' = 0,17·ln(kν·G1·l/(G2·d)), α' = 0,17·ln(kν1·l/d),
    /// χ = E·A/(G1·l²), λ1 = 2,12·χ^¾/(1 + 2,12·χ^¾), kν = 2,82 − 3,78ν + 2,18ν².
    pub fn new(p: &Pile, soil: &SoilColumn) -> Self {
        let tip = p.head - p.length;
        let (g1, nu1) = soil.average(tip, p.head.min(soil.top));
        let (g2, nu2) = soil.average(tip - 0.5 * p.length, tip);
        let l = p.length;
        let kv = k_nu(0.5 * (nu1 + nu2));
        let kv1 = k_nu(nu1);
        let beta_r = 0.17 * (kv * g1 * l / (g2 * p.d)).ln();
        let alpha_r = 0.17 * (kv1 * l / p.d).ln();
        let chi = p.ea / (g1 * l * l);
        let c = 2.12 * chi.powf(0.75);
        let lambda1 = c / (1.0 + c);
        let beta = beta_r / lambda1 + (1.0 - beta_r / alpha_r) / chi;
        let mut warnings = vec![];
        if l / p.d <= 5.0 {
            warnings.push(format!("l/d = {:.2} ≤ 5: п. 7.4.2 неприменим", l / p.d));
        }
        if g1 * l / (g2 * p.d) <= 1.0 {
            warnings.push("G1·l/(G2·d) ≤ 1: п. 7.4.2 неприменим".into());
        }
        Self { g1, nu1, g2, nu2, beta, compliance: beta / (g1 * l), warnings }
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
