//! Осадка свайного фундамента как условного (СП 24.13330.2021, п. 7.4.6 — контур по
//! наружным граням свай, расширенный на L·tg(φ_II,mt/4); осадка — по СП 22.13330.2016,
//! п. 5.6.31, метод послойного суммирования).
//!
//! s = β Σ (σzp,i − σzγ,i) h_i / E_i + β Σ σzγ,i h_i / Ee,i, β = 0,8; σzp = α p, σzγ = α σzg,0
//! (σzg,0 — природное давление на отметке подошвы условного фундамента); α — по Буссинеску
//! для фактического контура в его центре тяжести (для прямоугольника совпадает с табл. 5.8
//! СП 22). Ee = 5E (п. 5.6.31, при отсутствии данных — допущение). Нижняя граница сжимаемой
//! толщи — σzp = 0,5 σzg (п. 5.6.41). Подземные воды не учитываются (γ без взвешивания —
//! допущение). Номера пунктов СП 24 и СП 22 — по памяти, сверить с текстом норм.

use super::SoilColumn;

/// Результат расчёта условного фундамента.
#[derive(Clone, Debug)]
pub struct Conditional {
    /// Осадка, м.
    pub s: f64,
    /// Среднее давление по подошве, кПа; площадь подошвы, м²; φ_mt, град.
    pub p: f64,
    pub area: f64,
    pub phi_mt: f64,
    /// Мощность сжимаемой толщи, м.
    pub hc: f64,
}

fn area_centroid(p: &[[f64; 2]]) -> (f64, [f64; 2]) {
    let n = p.len();
    let (mut a, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (q, r) = (p[i], p[(i + 1) % n]);
        let c = q[0] * r[1] - r[0] * q[1];
        a += c;
        cx += (q[0] + r[0]) * c;
        cy += (q[1] + r[1]) * c;
    }
    (a / 2.0, [cx / (3.0 * a), cy / (3.0 * a)])
}

/// Контур, смещённый наружу на `d` (стыки рёбер — по пересечению смещённых прямых).
pub fn offset(p: &[[f64; 2]], d: f64) -> Vec<[f64; 2]> {
    let n = p.len();
    let sign = if area_centroid(p).0 > 0.0 { 1.0 } else { -1.0 };
    let line = |i: usize| {
        let (a, b) = (p[i], p[(i + 1) % n]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l = (dx * dx + dy * dy).sqrt();
        // Внешняя нормаль при обходе против часовой стрелки — (dy, −dx).
        let nrm = [sign * dy / l, -sign * dx / l];
        ([a[0] + d * nrm[0], a[1] + d * nrm[1]], [dx, dy])
    };
    (0..n)
        .map(|i| {
            let (p1, d1) = line((i + n - 1) % n);
            let (p2, d2) = line(i);
            let den = d1[0] * d2[1] - d1[1] * d2[0];
            if den.abs() < 1e-12 {
                return p2;
            }
            let t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / den;
            [p1[0] + t * d1[0], p1[1] + t * d1[1]]
        })
        .collect()
}

/// α = σz/p под точкой `c` на глубине `z` от равномерной нагрузки на многоугольнике
/// (Буссинеск, интегрирование по лучам).
pub fn alpha(p: &[[f64; 2]], c: [f64; 2], z: f64) -> f64 {
    let rays = 720;
    let f = |r: f64| z.powi(3) / (r * r + z * z).powf(1.5);
    let mut sum = 0.0;
    for k in 0..rays {
        let th = (k as f64 + 0.5) * std::f64::consts::TAU / rays as f64;
        let dir = [th.cos(), th.sin()];
        let mut hits: Vec<f64> = vec![];
        for i in 0..p.len() {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            let e = [b[0] - a[0], b[1] - a[1]];
            let den = dir[0] * e[1] - dir[1] * e[0];
            if den.abs() < 1e-14 {
                continue;
            }
            let w = [a[0] - c[0], a[1] - c[1]];
            let t = (w[0] * e[1] - w[1] * e[0]) / den;
            let u = (w[0] * dir[1] - w[1] * dir[0]) / den;
            if t > 0.0 && (0.0..1.0).contains(&u) {
                hits.push(t);
            }
        }
        hits.sort_by(|a, b| a.total_cmp(b));
        // Точка внутри — нечётное число пересечений: интервалы [0, h0], [h1, h2], …
        let mut r = if hits.len() % 2 == 1 { vec![0.0] } else { vec![] };
        r.extend(hits);
        for pair in r.chunks(2) {
            if let [ra, rb] = pair {
                sum += f(*ra) - f(*rb);
            }
        }
    }
    sum / rays as f64
}

/// `outline` — контур по наружным граням крайних свай; `z_head`, `length` — отметка голов
/// и длина свай; `n` — нагрузка на уровне голов, кН.
pub fn settlement(col: &SoilColumn, outline: &[[f64; 2]], z_head: f64, length: f64, n: f64) -> Result<Conditional, String> {
    let z_tip = z_head - length;
    let layer_at = |z: f64| col.layers.iter().find(|l| l.bottom < z).or(col.layers.last()).unwrap();
    // Средневзвешенный φ по длине сваи.
    let dz = 0.05;
    let steps = (length / dz).ceil() as usize;
    let mut phi_sum = 0.0;
    for k in 0..steps {
        let z = z_head - (k as f64 + 0.5) * length / steps as f64;
        let l = layer_at(z);
        if l.phi <= 0.0 {
            return Err(format!("условный фундамент: нет φ слоя на отметке {z:.2}"));
        }
        phi_sum += l.phi;
    }
    let phi_mt = phi_sum / steps as f64;
    let poly = offset(outline, length * (phi_mt / 4.0).to_radians().tan());
    let (a, c) = area_centroid(&poly);
    let area = a.abs();
    // Природное давление σzg на отметке z (от поверхности).
    let sigma_g = |z: f64| -> Result<f64, String> {
        let k = ((col.top - z) / dz).ceil().max(1.0) as usize;
        let mut s = 0.0;
        for i in 0..k {
            let zz = col.top - (i as f64 + 0.5) * (col.top - z) / k as f64;
            let l = layer_at(zz);
            if l.gamma <= 0.0 {
                return Err(format!("условный фундамент: нет γ слоя на отметке {zz:.2}"));
            }
            s += l.gamma * (col.top - z) / k as f64;
        }
        Ok(s)
    };
    let s_tip = sigma_g(z_tip)?;
    let s_head = sigma_g(z_head)?;
    // Вес блока (грунт и сваи) от голов до острия.
    let g_block = (s_tip - s_head) * area;
    let p = (n + g_block) / area;
    let (beta, h) = (0.8, 0.2);
    let mut s = 0.0;
    let mut z = 0.0;
    loop {
        let zc = z + h / 2.0;
        let al = alpha(&poly, c, zc);
        let (szp, szg_nat) = (al * p, sigma_g(z_tip - zc)?);
        if szp <= 0.5 * szg_nat {
            break;
        }
        let szy = al * s_tip;
        let e = layer_at(z_tip - zc).e;
        s += beta * ((szp - szy).max(0.0) * h / e + szy * h / (5.0 * e));
        z += h;
        if z > 200.0 {
            return Err("условный фундамент: сжимаемая толща больше 200 м".into());
        }
    }
    Ok(Conditional { s, p, area, phi_mt, hc: z })
}
