//! Армирование всего здания: поля требуемой площади по элементам оболочек (плиты и стены).
//!
//! Для каждого элемента: огибающая ПС1 (Капра–Мори) по сочетаниям ПС1, трещины по полным
//! (непродолжительное) и длительному (продолжительное) сочетаниям ПС2, поперечная арматура
//! по Q₀, добавка сверх фоновой сетки. Оси армирования: у плит x — проекция глобальной X; у
//! стен (нормаль почти горизонтальна) x — горизонталь в плоскости стены, y — вверх.

use crate::fem::{self, LoadCase, Model, V3};
use crate::sp63::capra_maury::{self, PlateForces, PlateSection, SlsContext, UlsContext};
use crate::sp63::materials::{ConcreteChars, RebarChars};
use crate::sp63::plate_shear;

/// Группа элементов с одним сечением, материалами и фоновым армированием.
pub struct Group<'a> {
    pub section: PlateSection,
    pub concrete: ConcreteChars,
    pub rebar: RebarChars,
    pub rsw: f64,
    pub uls: &'a UlsContext,
    pub crack_short: Option<&'a SlsContext>,
    pub crack_long: Option<&'a SlsContext>,
    /// Фоновая сетка, м²/м: [верх x, верх y, низ x, низ y].
    pub background: [f64; 4],
}

/// Результат по элементу.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Element {
    /// Требуемая площадь, м²/м: [верх x, верх y, низ x, низ y].
    pub required: [f64; 4],
    /// Добавка сверх фоновой сетки, м²/м.
    pub addition: [f64; 4],
    /// Требуемая поперечная арматура Asw/sw, м²/м²; 0 — не требуется.
    pub asw: f64,
    /// Подбор продольной арматуры невозможен (мала толщина) хотя бы в одном сочетании.
    pub impossible: bool,
    /// Не проходит бетонная полоса (8.55) — поперечная арматура не помогает.
    pub strip_fails: bool,
}

/// Ось x армирования элемента с нормалью `n`.
pub fn rebar_axis(n: V3) -> V3 {
    if n[2].abs() < 0.7 {
        // Стена: горизонталь в плоскости, y тогда направлена по вертикали (вверх или вниз).
        let h = fem::shell::cross([0.0, 0.0, 1.0], n);
        let l = fem::shell::norm(h);
        [h[0] / l, h[1] / l, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    }
}

/// Поля армирования. `groups[g]` — группы, `group_of[e]` — группа элемента (None — не
/// армируется); `cases` — загружения, `uls`, `sls`, `sls_long` — коэффициенты сочетаний.
pub fn design(
    m: &Model,
    groups: &[Group],
    group_of: &[Option<usize>],
    cases: &[LoadCase],
    uls: &[Vec<f64>],
    sls: &[Vec<f64>],
    sls_long: Option<&Vec<f64>>,
) -> Result<Vec<Element>, String> {
    if group_of.len() != m.shells.len() {
        return Err("группы заданы не для всех элементов".into());
    }
    for (g, grp) in groups.iter().enumerate() {
        if !(grp.section.z_top() > 0.0 && grp.section.z_bot() < 0.0) {
            return Err(format!("группа {}: арматура вне сечения (расстояние до центра стержней ≥ h/2)", g + 1));
        }
    }
    let lin = m.factorize()?;
    let u_cases: Vec<Vec<f64>> = cases.iter().map(|c| lin.solve(m, c)).collect::<Result<_, _>>()?;
    let combine = |c: &[f64]| -> Vec<f64> {
        let mut u = vec![0.0; m.ndof()];
        for (f, uc) in c.iter().zip(&u_cases) {
            if *f != 0.0 {
                for (a, b) in u.iter_mut().zip(uc) {
                    *a += f * b;
                }
            }
        }
        u
    };
    // Усилия в осях армирования: [Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy].
    let forces = |u: &[f64]| -> Vec<[f64; 8]> {
        m.shells
            .iter()
            .map(|s| {
                let c: Vec<V3> = s.nodes.iter().map(|&n| m.nodes[n]).collect();
                let ue: Vec<f64> = s.nodes.iter().flat_map(|&n| u[6 * n..6 * n + 6].to_vec()).collect();
                let (r, f) = fem::shell::forces(&c, &s.mat, &ue).expect("элемент проверен при сборке");
                fem::shell::rotate_forces(&r, &f, rebar_axis(r[2]))
            })
            .collect()
    };
    let pf = |f: &[f64; 8]| PlateForces { nx: f[0], ny: f[1], nxy: f[2], mx: f[3], my: f[4], mxy: f[5] };
    let mut out = vec![Element::default(); m.shells.len()];
    let members: Vec<Vec<usize>> = (0..groups.len()).map(|g| (0..m.shells.len()).filter(|&e| group_of[e] == Some(g)).collect()).collect();
    let take = |e: usize, d: &capra_maury::Design, out: &mut Vec<Element>| {
        let a = [d.top.asx, d.top.asy, d.bottom.asx, d.bottom.asy];
        for k in 0..4 {
            out[e].required[k] = out[e].required[k].max(a[k]);
        }
        out[e].impossible |= !d.feasible();
    };
    for c in uls {
        let f = forces(&combine(c));
        for (g, grp) in groups.iter().enumerate() {
            let pfs: Vec<PlateForces> = members[g].iter().map(|&e| pf(&f[e])).collect();
            for (k, d) in capra_maury::design_uls_batch(grp.uls, &pfs, 10.0).iter().enumerate() {
                take(members[g][k], d, &mut out);
            }
        }
        // Поперечная сила — по фоновой и уже подобранной продольной арматуре (в запас — меньшей).
        for (g, grp) in groups.iter().enumerate() {
            for &e in &members[g] {
                let fe = &f[e];
                let theta = fe[7].atan2(fe[6]);
                let (c2, s2) = (theta.cos().powi(2), theta.sin().powi(2));
                let b = grp.background;
                let as_n = (b[0] + b[2]) * c2 + (b[1] + b[3]) * s2;
                if let Some(r) = plate_shear::check(fe[6], fe[7], &pf(fe), &grp.section, &grp.concrete, &grp.rebar, grp.rsw, as_n) {
                    out[e].asw = out[e].asw.max(r.asw_per_sw);
                    out[e].strip_fails |= r.eta_strip > 1.0;
                }
            }
        }
    }
    let cracks = |c: &[f64], long: bool, out: &mut Vec<Element>| {
        let f = forces(&combine(c));
        for (g, grp) in groups.iter().enumerate() {
            if let Some(ctx) = if long { grp.crack_long } else { grp.crack_short } {
                let pfs: Vec<PlateForces> = members[g].iter().map(|&e| pf(&f[e])).collect();
                for (k, d) in capra_maury::design_sls_batch(ctx, &pfs, 10.0).iter().enumerate() {
                    take(members[g][k], d, out);
                }
            }
        }
    };
    for c in sls {
        cracks(c, false, &mut out);
    }
    if let Some(c) = sls_long {
        cracks(c, true, &mut out);
    }
    for (e, el) in out.iter_mut().enumerate() {
        if let Some(g) = group_of[e] {
            for k in 0..4 {
                el.addition[k] = (el.required[k] - groups[g].background[k]).max(0.0);
            }
        }
    }
    Ok(out)
}

/// Изополя в DXF R12: элементы — грани 3DFACE в своих координатах (стены — в пространстве),
/// по слою на поле; цвет ACI по 8 интервалам общей шкалы от 0 до max; шкала — слой `LEGEND`.
pub fn dxf_fields(m: &Model, fields: &[(&str, &[f64])], unit: &str) -> String {
    const ACI: [u8; 8] = [151, 141, 131, 121, 111, 51, 41, 11];
    let hi = fields.iter().flat_map(|(_, v)| v.iter().copied()).fold(0.0, f64::max).max(1e-12);
    let mut s = String::from("0\nSECTION\n2\nENTITIES\n");
    for (layer, values) in fields {
        for (e, sh) in m.shells.iter().enumerate() {
            let v = values.get(e).copied().unwrap_or(0.0);
            if v <= 1e-12 {
                continue;
            }
            let bin = (((v / hi) * 8.0).ceil() as usize).clamp(1, 8) - 1;
            let p: Vec<V3> = sh.nodes.iter().map(|&n| m.nodes[n]).collect();
            let q = [p[0], p[1], p[2], *p.get(3).unwrap_or(&p[2])];
            s += &format!("0\n3DFACE\n8\n{layer}\n62\n{}\n", ACI[bin]);
            for (k, pt) in q.iter().enumerate() {
                s += &format!("{}\n{}\n{}\n{}\n{}\n{}\n", 10 + k, pt[0], 20 + k, pt[1], 30 + k, pt[2]);
            }
        }
    }
    let (x0, y0) = m.nodes.iter().fold((f64::INFINITY, f64::INFINITY), |a, p| (a.0.min(p[0]), a.1.min(p[1])));
    for b in 0..8 {
        let t = format!("{:.2}–{:.2} {unit}", hi * b as f64 / 8.0, hi * (b + 1) as f64 / 8.0);
        s += &format!("0\nTEXT\n8\nLEGEND\n62\n{}\n10\n{}\n20\n{}\n30\n0\n40\n0.4\n1\n{t}\n", ACI[b], x0, y0 - 1.5 - 0.6 * b as f64);
    }
    s + "0\nENDSEC\n0\nEOF\n"
}
