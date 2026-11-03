//! Этап 3 (точный): матрица податливости свайного поля на 3D-основании, решатель Code_Aster
//! (GPL, отдельный процесс `run_aster`).
//!
//! Модель: упругий слоистый массив (гексаэдры HEXA8, E и ν слоёв — из грунтовой колонки в
//! каждой точке плана), сваи — объёмные стойки квадратного сечения равной площади
//! (b = d·√π/2, E = EA/b²), полное сцепление с грунтом (без проскальзывания — допущение).
//! Верх модели — отметка голов свай (грунт выше и выемка котлована не моделируются,
//! контакт плиты с грунтом не учитывается — как в этапах 1–2); низ закреплён, боковые грани —
//! по нормали. На голову каждой сваи — единичная сила (поровну по узлам головы), все
//! загружения — одной факторизацией (MUMPS). F_ij — средняя осадка головы i от
//! силы 1 кН на свае j, м/кН; F симметризуется. Используется вместо матрицы СП 24 п. 7.4 или
//! условного фундамента в `optimize::Problem::evaluate_flex`.

use crate::piles::{Pile, SoilColumn};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

/// Параметры модели и запуска.
#[derive(Clone, Debug)]
pub struct Options {
    /// Путь к `run_aster`.
    pub run_aster: PathBuf,
    /// Рабочий каталог (создаётся).
    pub workdir: PathBuf,
    /// Наибольший размер элемента в зоне свай (план и глубина), м.
    pub h: f64,
    /// Наибольший размер элемента по глубине в пределах длины свай, м.
    pub hz: f64,
    /// Рост размера элемента вне зоны свай: h + (g − 1)·расстояние.
    pub grade: f64,
    /// Расстояние от крайних свай до боковых границ, м.
    pub lateral: f64,
    /// Глубина модели ниже острия свай, м.
    pub below: f64,
    /// Квадратичные элементы (HEXA20): сходимость по сетке заметно лучше HEXA8.
    pub quadratic: bool,
    /// Квадратичные элементы в нелинейном расчёте (по умолчанию да: HEXA8 «запирается» при
    /// пластическом течении без изменения объёма — D одиночной сваи 1,10 против 1,49).
    pub quadratic_nl: bool,
    pub ncpus: usize,
    pub memory_mb: usize,
    pub time_limit_s: usize,
}

impl Options {
    /// Путь к `run_aster`: переменная `FEMGEO_ASTER`, иначе `/opt/aster/bin/run_aster`.
    pub fn find_run_aster() -> Option<PathBuf> {
        let p = std::env::var_os("FEMGEO_ASTER").map(PathBuf::from).unwrap_or_else(|| PathBuf::from("/opt/aster/bin/run_aster"));
        p.is_file().then_some(p)
    }

    /// Значения по умолчанию: боковые границы и глубина ниже острия — max(L, B), B — наибольший
    /// размер свайного поля в плане (допущение, проверять чувствительность).
    pub fn new(run_aster: PathBuf, workdir: PathBuf, piles: &[Pile]) -> Self {
        let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
        for p in piles {
            lo = [lo[0].min(p.x), lo[1].min(p.y)];
            hi = [hi[0].max(p.x), hi[1].max(p.y)];
        }
        let b = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(0.0);
        let l = piles.iter().map(|p| p.length).fold(0.0, f64::max);
        let r = l.max(b);
        Self { run_aster, workdir, h: 1.0, hz: 2.0, grade: 1.3, lateral: r, below: r, quadratic: true, quadratic_nl: true, ncpus: 4, memory_mb: 12_000, time_limit_s: 7200 }
    }
}

/// Структурированная сетка: линии по осям, материалы элементов, узлы голов свай.
#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// По возрастанию (снизу вверх); последняя — отметка голов.
    pub z: Vec<f64>,
    /// Материал элемента (i, j, k) — индекс в `materials`, порядок i быстрее j быстрее k.
    pub mat: Vec<usize>,
    /// Материалы; `VOID` в `mat` — нет элемента (слой оголовков вне свай).
    pub materials: Vec<Material>,
    /// Над отметкой голов — слой оголовков (только ячейки свай, нагрузка — к их верху).
    pub stub: bool,
    /// Номера узлов головы каждой сваи (с нуля, порядок записи сетки).
    pub heads: Vec<Vec<usize>>,
}

/// Материал элемента: E, ν (кПа), γ (кН/м³), φ (град), c (кПа); `pile` — бетон сваи.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Material {
    pub e: f64,
    pub nu: f64,
    pub gamma: f64,
    pub phi: f64,
    pub c: f64,
    pub pile: bool,
    /// Упругие E, ν шага собственного веса (природные напряжения): грунт — ν0 = K0/(1 + K0);
    /// свая — грунт своего слоя («свая без возмущения»), оголовок — сама свая.
    pub e0: f64,
    pub nu0: f64,
}

/// ν0 = K0/(1 + K0), K0 = 1 − sin φ (Яки); φ = 0 — ν слоя.
fn nu_geostatic(l: &crate::piles::Layer) -> f64 {
    if l.phi > 0.0 {
        let k0 = 1.0 - l.phi.to_radians().sin();
        k0 / (1.0 + k0)
    } else {
        l.nu
    }
}

pub const VOID: usize = usize::MAX;

impl Grid {
    pub fn node(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.y.len() + j) * self.x.len() + i
    }
    pub fn n_nodes(&self) -> usize {
        self.x.len() * self.y.len() * self.z.len()
    }
    pub fn n_elems(&self) -> usize {
        self.mat.iter().filter(|&&m| m != VOID).count()
    }
    /// Узел принадлежит хотя бы одному элементу.
    fn used(&self, i: usize, j: usize, k: usize) -> bool {
        let (nx, ny, nz) = (self.x.len() - 1, self.y.len() - 1, self.z.len() - 1);
        for kk in k.saturating_sub(1)..=k.min(nz - 1) {
            for jj in j.saturating_sub(1)..=j.min(ny - 1) {
                for ii in i.saturating_sub(1)..=i.min(nx - 1) {
                    if self.mat[(kk * ny + jj) * nx + ii] != VOID {
                        return true;
                    }
                }
            }
        }
        false
    }
}

/// Сторона квадрата равной площади для круглой сваи диаметром d.
pub fn square_side(d: f64) -> f64 {
    d * std::f64::consts::PI.sqrt() / 2.0
}

/// Линии сетки по оси: жёсткие точки `hard` (включаются), зона [a, b] с шагом ≤ h, вне зоны —
/// шаг h + (g − 1)·расстояние до зоны, до границ [lo, hi].
fn lines(mut hard: Vec<f64>, a: f64, b: f64, lo: f64, hi: f64, h: f64, g: f64) -> Vec<f64> {
    hard.push(lo);
    hard.push(hi);
    hard.retain(|v| *v >= lo - 1e-9 && *v <= hi + 1e-9);
    hard.sort_by(f64::total_cmp);
    hard.dedup_by(|p, q| (*p - *q).abs() < 0.02);
    let mut out = vec![hard[0]];
    for w in hard.windows(2) {
        let (p, q) = (w[0], w[1]);
        let mid = 0.5 * (p + q);
        if mid >= a && mid <= b {
            let n = ((q - p) / h - 1e-9).ceil().max(1.0) as usize;
            out.extend((1..=n).map(|s| if s == n { q } else { p + (q - p) * s as f64 / n as f64 }));
            continue;
        }
        // Вне зоны: шаги растут от ближнего к зоне конца, последний подгоняется.
        let (near, d0) = if mid > b { (p, p - b) } else { (q, a - q) };
        let len = q - p;
        let mut t = vec![0.0];
        while *t.last().unwrap() < len - 1e-9 {
            let tl = *t.last().unwrap();
            t.push(tl + h + (g - 1.0) * (d0 + tl));
        }
        if t.len() > 2 && t[t.len() - 1] - len > 0.5 * (t[t.len() - 1] - t[t.len() - 2]) {
            t.pop();
        }
        let sc = len / t.last().unwrap();
        let u: Vec<f64> = t[1..].iter().map(|v| v * sc).collect();
        let n = u.len();
        let pts: Vec<f64> = if near == p {
            (0..n).map(|k| if k + 1 == n { q } else { p + u[k] }).collect()
        } else {
            (0..n - 1).rev().map(|k| q - u[k]).chain(std::iter::once(q)).collect()
        };
        out.extend(pts);
    }
    out
}

fn layer_at(col: &SoilColumn, z: f64) -> Option<crate::piles::Layer> {
    col.layers.iter().enumerate().find(|(i, l)| z > l.bottom || i + 1 == col.layers.len()).map(|(_, l)| *l)
}



/// Сетка для свай `piles` (одна отметка голов) и грунта `soil`.
pub fn grid(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn, o: &Options) -> Result<Grid, String> {
    grid_with(piles, soil, o, false)
}

/// Сетка; `stub` — слой оголовков высотой b над отметкой голов (нагрузка не в узлах, общих с грунтом).
pub fn grid_with(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn, o: &Options, stub: bool) -> Result<Grid, String> {
    if piles.is_empty() {
        return Err("нет свай".into());
    }
    let head = piles[0].head;
    if piles.iter().any(|p| (p.head - head).abs() > 1e-6) {
        return Err("3D-модель: головы свай на разных отметках не поддерживаются".into());
    }
    let sq: Vec<[f64; 5]> = piles
        .iter()
        .map(|p| {
            let b = square_side(p.d);
            [p.x - b / 2.0, p.x + b / 2.0, p.y - b / 2.0, p.y + b / 2.0, p.head - p.length]
        })
        .collect();
    let (x0, x1) = (sq.iter().map(|s| s[0]).fold(f64::MAX, f64::min), sq.iter().map(|s| s[1]).fold(f64::MIN, f64::max));
    let (y0, y1) = (sq.iter().map(|s| s[2]).fold(f64::MAX, f64::min), sq.iter().map(|s| s[3]).fold(f64::MIN, f64::max));
    let tip = sq.iter().map(|s| s[4]).fold(f64::MAX, f64::min);
    let x = lines(sq.iter().flat_map(|s| [s[0], s[1]]).collect(), x0, x1, x0 - o.lateral, x1 + o.lateral, o.h, o.grade);
    let y = lines(sq.iter().flat_map(|s| [s[2], s[3]]).collect(), y0, y1, y0 - o.lateral, y1 + o.lateral, o.h, o.grade);
    let bottom = tip - o.below;
    let cx = 0.5 * (x0 + x1);
    let cy = 0.5 * (y0 + y1);
    let mut zh: Vec<f64> = sq.iter().map(|s| s[4]).collect();
    zh.extend(soil(cx, cy).layers.iter().map(|l| l.bottom));
    let mut z = lines(zh, tip, head, bottom, head, o.hz, o.grade);
    if stub {
        z.push(head + piles.iter().map(|p| square_side(p.d)).fold(0.0, f64::max));
    }
    let (nx, ny, nz) = (x.len() - 1, y.len() - 1, z.len() - 1);
    // Колонки грунта — по центрам ячеек плана (кэш).
    let cols: Vec<SoilColumn> = (0..ny).flat_map(|j| (0..nx).map(move |i| (i, j))).map(|(i, j)| soil(0.5 * (x[i] + x[i + 1]), 0.5 * (y[j] + y[j + 1]))).collect();
    let mut materials: Vec<Material> = vec![];
    let mut key = |m: Material| -> usize {
        if let Some(p) = materials.iter().position(|q| *q == m) {
            return p;
        }
        materials.push(m);
        materials.len() - 1
    };
    let pile_e: Vec<f64> = piles.iter().map(|p| p.ea / square_side(p.d).powi(2)).collect();
    let mut mat = Vec::with_capacity(nx * ny * nz);
    for k in 0..nz {
        let zc = 0.5 * (z[k] + z[k + 1]);
        for j in 0..ny {
            let yc = 0.5 * (y[j] + y[j + 1]);
            for i in 0..nx {
                let xc = 0.5 * (x[i] + x[i + 1]);
                let p = sq.iter().position(|s| xc > s[0] && xc < s[1] && yc > s[2] && yc < s[3] && zc > s[4]);
                if p.is_none() && zc > head {
                    mat.push(VOID);
                    continue;
                }
                if p.is_some() && zc > head {
                    // Оголовок: без собственного веса.
                    let e = pile_e[p.unwrap()];
                    mat.push(key(Material { e, nu: 0.2, gamma: 0.0, phi: 0.0, c: 0.0, pile: true, e0: e, nu0: 0.2 }));
                    continue;
                }
                let l = layer_at(&cols[j * nx + i], zc).ok_or("грунт: нет слоёв")?;
                if !(l.e > 0.0) || !(l.nu >= 0.0 && l.nu < 0.5) {
                    return Err(format!("грунт в точке ({xc:.2}, {yc:.2}, {zc:.2}): E {}, ν {}", l.e, l.nu));
                }
                let (e0, nu0) = (l.e, nu_geostatic(&l));
                // Свая: γ слоя (вес сваи сверх грунта не учитывается — допущение, мало).
                mat.push(match p {
                    Some(p) => key(Material { e: pile_e[p], nu: 0.2, gamma: l.gamma, phi: 0.0, c: 0.0, pile: true, e0, nu0 }),
                    None => key(Material { e: l.e, nu: l.nu, gamma: l.gamma, phi: l.phi, c: l.c, pile: false, e0, nu0 }),
                });
            }
        }
    }
    let kt = z.len() - 1;
    let mut g = Grid { x, y, z, mat, materials, stub, heads: vec![] };
    g.heads = sq
        .iter()
        .map(|s| {
            let mut v = vec![];
            for j in 0..g.y.len() {
                for i in 0..g.x.len() {
                    if g.x[i] > s[0] - 0.03 && g.x[i] < s[1] + 0.03 && g.y[j] > s[2] - 0.03 && g.y[j] < s[3] + 0.03 {
                        v.push(g.node(i, j, kt));
                    }
                }
            }
            v
        })
        .collect();
    if g.heads.iter().any(|h| h.len() < 4) {
        return Err("3D-модель: голова сваи не попала в узлы сетки".into());
    }
    Ok(g)
}

/// Сетка в формате ASTER (текст `.mail`).
pub fn mail(g: &Grid) -> String {
    let mut s = String::with_capacity(64 * g.n_nodes());
    s.push_str("COOR_3D\n");
    for k in 0..g.z.len() {
        for j in 0..g.y.len() {
            for i in 0..g.x.len() {
                if g.stub && k + 1 == g.z.len() && !g.used(i, j, k) {
                    continue;
                }
                let _ = writeln!(s, "N{} {:.6} {:.6} {:.6}", g.node(i, j, k) + 1, g.x[i], g.y[j], g.z[k]);
            }
        }
    }
    s.push_str("FINSF\nHEXA8\n");
    let (nx, ny) = (g.x.len() - 1, g.y.len() - 1);
    let mut e = 0;
    let mut groups: Vec<Vec<usize>> = vec![vec![]; g.materials.len()];
    for k in 0..g.z.len() - 1 {
        for j in 0..ny {
            for i in 0..nx {
                let m = g.mat[(k * ny + j) * nx + i];
                if m == VOID {
                    continue;
                }
                let n = [g.node(i, j, k), g.node(i + 1, j, k), g.node(i + 1, j + 1, k), g.node(i, j + 1, k), g.node(i, j, k + 1), g.node(i + 1, j, k + 1), g.node(i + 1, j + 1, k + 1), g.node(i, j + 1, k + 1)];
                let _ = writeln!(s, "M{} N{} N{} N{} N{} N{} N{} N{} N{}", e + 1, n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1, n[4] + 1, n[5] + 1, n[6] + 1, n[7] + 1);
                groups[m].push(e);
                e += 1;
            }
        }
    }
    s.push_str("FINSF\n");
    let list = |s: &mut String, name: &str, kind: &str, ids: &mut dyn Iterator<Item = usize>| {
        let _ = writeln!(s, "{kind}\n{name}");
        for (c, id) in ids.enumerate() {
            let _ = write!(s, "{}{}{}", if kind == "GROUP_MA" { "M" } else { "N" }, id + 1, if c % 10 == 9 { "\n" } else { " " });
        }
        s.push_str("\nFINSF\n");
    };
    for (m, ids) in groups.iter().enumerate() {
        if !ids.is_empty() {
            list(&mut s, &format!("MAT{m}"), "GROUP_MA", &mut ids.iter().copied());
        }
    }
    for (p, h) in g.heads.iter().enumerate() {
        list(&mut s, &format!("H{p}"), "GROUP_NO", &mut h.iter().copied());
    }
    s.push_str("FIN\n");
    s
}

/// Общая часть командного файла: сетка, модель, материалы, закрепления, узлы голов. `nl` —
/// нелинейный расчёт: плотность ρ = γ/g и Кулон–Мор для грунта.
fn comm_head(g: &Grid, quadratic: bool, nl: bool) -> Result<String, String> {
    let mut s = String::from("import numpy as np\nDEBUT(LANG='EN')\nmesh = LIRE_MAILLAGE(FORMAT='ASTER', UNITE=20)\n");
    if quadratic {
        s.push_str("mesh = CREA_MAILLAGE(MAILLAGE=mesh, LINE_QUAD=_F(TOUT='OUI'))\n");
    }
    // Граничные узлы (в том числе средние узлы рёбер) — по плоскостям.
    let (x0, x1, y0, y1, z0) = (g.x[0], *g.x.last().unwrap(), g.y[0], *g.y.last().unwrap(), g.z[0]);
    let plane = |name: &str, p: [f64; 3], n: [f64; 3]| format!("_F(NOM='{name}', OPTION='PLAN', POINT=({}, {}, {}), VECT_NORMALE=({}, {}, {}), PRECISION=1e-4)", p[0], p[1], p[2], n[0], n[1], n[2]);
    let _ = writeln!(
        s,
        "mesh = DEFI_GROUP(reuse=mesh, MAILLAGE=mesh, CREA_GROUP_NO=({}, {}, {}, {}, {}))",
        plane("BOTTOM", [x0, y0, z0], [0.0, 0.0, 1.0]),
        plane("XLO", [x0, y0, z0], [1.0, 0.0, 0.0]),
        plane("XHI", [x1, y0, z0], [1.0, 0.0, 0.0]),
        plane("YLO", [x0, y0, z0], [0.0, 1.0, 0.0]),
        plane("YHI", [x0, y1, z0], [0.0, 1.0, 0.0])
    );
    s.push_str("mod = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION='3D'))\n");
    let used: Vec<usize> = (0..g.materials.len()).filter(|m| g.mat.contains(m)).collect();
    s.push_str("aff = []\naff0 = []\nsoil_groups = []\npile_groups = []\n");
    for &m in &used {
        let Material { e, nu, gamma, phi, c, pile, e0, nu0 } = g.materials[m];
        if !nl {
            let _ = writeln!(s, "m{m} = DEFI_MATERIAU(ELAS=_F(E={e:e}, NU={nu}))");
        } else if pile {
            let _ = writeln!(s, "m{m} = DEFI_MATERIAU(ELAS=_F(E={e:e}, NU={nu}, RHO={}))\npile_groups.append('MAT{m}')", gamma / G);
        } else {
            if !(gamma > 0.0) || !(phi > 0.0 || c > 0.0) {
                return Err(format!("нелинейный грунт: у слоя E {e} не заданы γ и φ/c"));
            }
            // Кулон–Мор, неассоциированный (ψ = 0 — допущение); c не меньше `MC_C_MIN`
            // (вершина конуса на поверхности при c = 0 — расходимость).
            let _ = writeln!(
                s,
                "m{m} = DEFI_MATERIAU(ELAS=_F(E={e:e}, NU={nu}, RHO={}), MOHR_COULOMB=_F(PHI={phi}, ANGDIL=0.0, COHESION={}))\nsoil_groups.append('MAT{m}')",
                gamma / G,
                c.max(MC_C_MIN)
            );
        }
        let _ = writeln!(s, "aff.append(_F(GROUP_MA='MAT{m}', MATER=m{m}))");
        if nl {
            let _ = writeln!(s, "g{m} = DEFI_MATERIAU(ELAS=_F(E={e0:e}, NU={nu0}, RHO={}))\naff0.append(_F(GROUP_MA='MAT{m}', MATER=g{m}))", gamma / G);
        }
    }
    s.push_str("cm = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=aff)\n");
    if nl {
        s.push_str("cm0 = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=aff0)\n");
    }
    // Закрепления — исключением (AFFE_CHAR_CINE): быстрее множителей Лагранжа, результат тот же.
    s.push_str("bc = AFFE_CHAR_CINE(MODELE=mod, MECA_IMPO=(_F(GROUP_NO='BOTTOM', DX=0, DY=0, DZ=0), _F(GROUP_NO=('XLO', 'XHI'), DX=0), _F(GROUP_NO=('YLO', 'YHI'), DY=0)))\n");
    let _ = writeln!(s, "heads = [list(mesh.getNodes('H%d' % p)) for p in range({})]", g.heads.len());
    Ok(s)
}

const G: f64 = 9.81;
/// Допуск невязки — наибольшая узловая сила, доля наибольшей силы на узел головы. Не
/// относительный: RESI_GLOB_RELA нормируется весом массива (≫ нагрузок свай), D «шумит».
const NL_RESI: f64 = 1e-4;
/// Наименьшее сцепление, кПа (допущение для сходимости у поверхности).
const MC_C_MIN: f64 = 1.0;

/// Командный файл линейного расчёта: единичные силы на головы (одна факторизация), F в `out`.
pub fn comm(g: &Grid, out: &Path, quadratic: bool) -> String {
    let mut s = comm_head(g, quadratic, false).unwrap();
    s.push_str(
        "ke = CALC_MATR_ELEM(OPTION='RIGI_MECA', MODELE=mod, CHAM_MATER=cm)\n\
         nu = NUME_DDL(MATR_RIGI=ke)\n\
         K = ASSE_MATRICE(MATR_ELEM=ke, NUME_DDL=nu, CHAR_CINE=bc)\n\
         K = FACTORISER(reuse=K, MATR_ASSE=K, METHODE='MUMPS')\n\
         n = len(heads)\n\
         f = np.zeros((n, n))\n\
         for j, h in enumerate(heads):\n    \
             c = AFFE_CHAR_MECA(MODELE=mod, FORCE_NODALE=_F(GROUP_NO='H%d' % j, FZ=-1.0 / len(h)))\n    \
             fv = ASSE_VECTEUR(VECT_ELEM=CALC_VECT_ELEM(OPTION='CHAR_MECA', CHARGE=c), NUME_DDL=nu)\n    \
             u = RESOUDRE(MATR=K, CHAM_NO=fv, RESI_RELA=-1.0)\n    \
             sf = u.toSimpleFieldOnNodes()\n    \
             iz = list(sf.getComponents()).index('DZ')\n    \
             v = np.array(sf.toNumpy()[0], copy=True)\n    \
             for i, hh in enumerate(heads):\n        \
                 f[i, j] = -np.mean(v[hh, iz])\n",
    );
    let _ = writeln!(s, "np.savetxt('{}', f, fmt='%.10e')\nFIN()", out.display());
    s
}

/// Командный файл нелинейного расчёта: 1) собственный вес, упруго (природные напряжения,
/// K0 = 1 − sin φ); 2) от этого напряжённого состояния (перемещения — с нуля) —
/// силы `loads` (кН, сжатие > 0) на головы, `steps` шагов, грунт — Кулон–Мор. В `out` —
/// строки «время, осадки голов» по шагам (м).
pub fn comm_nl(g: &Grid, loads: &[f64], steps: usize, out: &Path, quadratic: bool) -> Result<String, String> {
    let mut s = comm_head(g, quadratic, true)?;
    let n: Vec<String> = loads.iter().map(|v| format!("{v:e}")).collect();
    let _ = writeln!(s, "N = [{}]", n.join(", "));
    let _ = writeln!(s, "G = {G}\nSTEPS = {steps}\nRESI = {NL_RESI:e}\nOUT = '{}'", out.display());
    s.push_str(
        r#"grav = AFFE_CHAR_MECA(MODELE=mod, PESANTEUR=_F(GRAVITE=G, DIRECTION=(0.0, 0.0, -1.0)))
load = AFFE_CHAR_MECA(MODELE=mod, FORCE_NODALE=[_F(GROUP_NO='H%d' % p, FZ=-N[p] / len(h)) for p, h in enumerate(heads)])
L0 = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))
res0 = STAT_NON_LINE(MODELE=mod, CHAM_MATER=cm0, EXCIT=(_F(CHARGE=bc), _F(CHARGE=grav)),
    COMPORTEMENT=_F(TOUT='OUI', RELATION='ELAS'), INCREMENT=_F(LIST_INST=L0), SOLVEUR=_F(METHODE='MUMPS'))
sig0 = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_SIEF_R', RESULTAT=res0, NOM_CHAM='SIEF_ELGA', INST=1.0)
ramp = DEFI_FONCTION(NOM_PARA='INST', VALE=(1.0, 0.0, 2.0, 1.0))
L = DEFI_LIST_REEL(DEBUT=1.0, INTERVALLE=_F(JUSQU_A=2.0, NOMBRE=STEPS))
LI = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=L), ECHEC=_F(ACTION='DECOUPE', SUBD_NIVEAU=5, SUBD_PAS=4))
comp = [_F(GROUP_MA=tuple(soil_groups), RELATION='MOHR_COULOMB')]
if pile_groups:
    comp.append(_F(GROUP_MA=tuple(pile_groups), RELATION='ELAS'))
# Упругое решение на той же сетке (для отношения s_нл/s_упр, не зависящего от сетки).
rl = MECA_STATIQUE(MODELE=mod, CHAM_MATER=cm, EXCIT=(_F(CHARGE=bc), _F(CHARGE=load)), SOLVEUR=_F(METHODE='MUMPS'))
sfl = rl.getField('DEPL', rl.getLastIndex()).toSimpleFieldOnNodes()
vl = np.array(sfl.toNumpy()[0], copy=True)
izl = list(sfl.getComponents()).index('DZ')
rows = [[-1.0] + [-np.mean(vl[h, izl]) for h in heads]]
from libaster import AsterError
kw = dict(MODELE=mod, CHAM_MATER=cm, EXCIT=(_F(CHARGE=bc), _F(CHARGE=grav), _F(CHARGE=load, FONC_MULT=ramp)),
    COMPORTEMENT=comp, NEWTON=_F(MATRICE='ELASTIQUE', REAC_INCR=0, REAC_ITER=0),
    CONVERGENCE=_F(RESI_GLOB_MAXI=RESI * max(N) / 4.0, ITER_GLOB_MAXI=500), SOLVEUR=_F(METHODE='MUMPS'))
# Первый шаг, затем остальные одним вызовом (reuse): при несходимости пройденные шаги
# остаются в res (результат дальше не определён). Матрица — начальная упругая (одно
# разложение на вызов, симметричная): у срыва сходится надёжнее касательной.
res = None
try:
    res = STAT_NON_LINE(ETAT_INIT=_F(SIGM=sig0), INCREMENT=_F(LIST_INST=LI, INST_INIT=1.0, INST_FIN=1.0 + 1.0 / STEPS), **kw)
    if STEPS > 1:
        res = STAT_NON_LINE(reuse=res, RESULTAT=res, ETAT_INIT=_F(EVOL_NOLI=res), INCREMENT=_F(LIST_INST=LI, INST_FIN=2.0), **kw)
except AsterError:
    pass
for idx in (res.getIndexes() if res is not None else []):
    sf = res.getField('DEPL', idx).toSimpleFieldOnNodes()
    iz = list(sf.getComponents()).index('DZ')
    v = np.array(sf.toNumpy()[0], copy=True)
    rows.append([res.getTime(idx)] + [-np.mean(v[h, iz]) for h in heads])
np.savetxt(OUT, np.array(rows) if rows else np.zeros((0, len(heads) + 1)), fmt='%.10e')
FIN()
"#,
    );
    Ok(s)
}

/// Запуск Code_Aster в `o.workdir`: сетка, командный файл; возвращает текст `out` и время, с.
fn execute(o: &Options, g: &Grid, comm_text: &dyn Fn(&Path) -> Result<String, String>, out_name: &str, strict: bool) -> Result<(String, f64), String> {
    std::fs::create_dir_all(&o.workdir).map_err(|e| format!("{}: {e}", o.workdir.display()))?;
    let wd = o.workdir.canonicalize().map_err(|e| e.to_string())?;
    let out = wd.join(out_name);
    let _ = std::fs::remove_file(&out);
    let w = |name: &str, text: &str| std::fs::write(wd.join(name), text).map_err(|e| format!("{name}: {e}"));
    w("model.mail", &mail(g))?;
    w("model.comm", &comm_text(&out)?)?;
    w(
        "model.export",
        &format!(
            "P time_limit {}\nP memory_limit {}\nP ncpus {}\nF comm {} D 1\nF mail {} D 20\nF mess {} R 6\n",
            o.time_limit_s,
            o.memory_mb,
            o.ncpus,
            wd.join("model.comm").display(),
            wd.join("model.mail").display(),
            wd.join("model.mess").display()
        ),
    )?;
    let t = std::time::Instant::now();
    let r = std::process::Command::new(&o.run_aster)
        .arg(wd.join("model.export"))
        .current_dir(&wd)
        .stdout(std::fs::File::create(wd.join("run.log")).map_err(|e| e.to_string())?)
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(|e| format!("{}: {e}", o.run_aster.display()))?;
    let seconds = t.elapsed().as_secs_f64();
    let log = std::fs::read_to_string(wd.join("run.log")).unwrap_or_default();
    let diag = log.lines().find(|l| l.contains("DIAGNOSTIC JOB")).unwrap_or("").to_string();
    if !(diag.contains(": OK") || diag.contains("<A>") || (!strict && out.is_file())) {
        return Err(format!("Code_Aster: {} (код {:?}), журнал {}", diag.trim(), r.code(), wd.join("run.log").display()));
    }
    let text = std::fs::read_to_string(&out).map_err(|e| format!("{}: {e}", out.display()))?;
    Ok((text, seconds))
}

fn numbers(text: &str) -> Result<Vec<f64>, String> {
    text.split_whitespace().map(|t| t.parse::<f64>().map_err(|e| e.to_string())).collect()
}

/// Нелинейный расчёт: осадки голов свай под силами `loads` (кН, сжатие > 0).
#[derive(Clone, Debug)]
pub struct Nonlinear {
    /// Осадки на последнем сошедшемся шаге, м; `complete` — дошёл до полной нагрузки (иначе
    /// результат при полной нагрузке не определён: срыв или несходимость).
    pub s: Vec<f64>,
    /// Упругие осадки под теми же силами на той же сетке, м.
    pub s_el: Vec<f64>,
    pub complete: bool,
    /// Кривая: доля нагрузки и осадки по шагам.
    pub curve: Vec<(f64, Vec<f64>)>,
    pub seconds: f64,
}

/// Осадки свай на упругопластическом основании (Code_Aster, Кулон–Мор). Несходимость —
/// ошибка (результат не определён).
pub fn settlements_nl(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn, loads: &[f64], steps: usize, o: &Options) -> Result<Nonlinear, String> {
    if loads.len() != piles.len() {
        return Err("число сил не равно числу свай".into());
    }
    let g = grid_with(piles, soil, o, true)?;
    let (text, seconds) = execute(o, &g, &|out| comm_nl(&g, loads, steps.max(1), out, o.quadratic_nl), "settle.txt", false)?;
    let v = numbers(&text)?;
    let n = piles.len();
    if v.len() % (n + 1) != 0 || v.len() < 3 * (n + 1) {
        return Err("нелинейный расчёт не сошёлся на первом шаге — результат не определён".into());
    }
    let s_el = v[1..n + 1].to_vec();
    let curve: Vec<(f64, Vec<f64>)> = v[n + 1..].chunks(n + 1).map(|r| (r[0] - 1.0, r[1..].to_vec())).collect();
    let last = curve.last().unwrap();
    let complete = (last.0 - 1.0).abs() < 1e-6 && last.1.iter().all(|x| x.is_finite());
    Ok(Nonlinear { s: last.1.clone(), s_el, complete, curve, seconds })
}

/// Результат: матрица податливости и сведения о модели.
#[derive(Clone, Debug)]
pub struct Flex3d {
    /// F, м/кН, n × n построчно, симметризованная.
    pub flexibility: Vec<f64>,
    /// Наибольшая несимметрия |F_ij − F_ji| / max F_ii до симметризации.
    pub asymmetry: f64,
    pub nodes: usize,
    pub elements: usize,
    pub seconds: f64,
}

/// Матрица податливости свай на 3D-основании (запуск Code_Aster).
pub fn flexibility(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn, o: &Options) -> Result<Flex3d, String> {
    let g = grid_with(piles, soil, o, true)?;
    let (text, seconds) = execute(o, &g, &|out| Ok(comm(&g, out, o.quadratic)), "flex.txt", true)?;
    let v = numbers(&text)?;
    let n = piles.len();
    if v.len() != n * n {
        return Err(format!("Code_Aster: {} значений вместо {}", v.len(), n * n));
    }
    let fmax = (0..n).map(|i| v[i * n + i]).fold(0.0, f64::max);
    if !(fmax > 0.0) || v.iter().any(|x| !x.is_finite()) {
        return Err("Code_Aster: недопустимая матрица податливости".into());
    }
    let mut f = v.clone();
    let mut asym: f64 = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            asym = asym.max((v[i * n + j] - v[j * n + i]).abs() / fmax);
            let m = 0.5 * (v[i * n + j] + v[j * n + i]);
            f[i * n + j] = m;
            f[j * n + i] = m;
        }
    }
    Ok(Flex3d { flexibility: f, asymmetry: asym, nodes: g.n_nodes(), elements: g.n_elems(), seconds })
}
