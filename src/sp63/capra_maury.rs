//! Подбор арматуры пластины по Капра–Мори на НДМ полосы.
//!
//! Для каждой фасетки θ (0…180°, шаг `step_deg`) — N(θ), M(θ), требуемые площади рядов полосы,
//! затем для верхней и нижней сеток — линейная задача [`mesh_lp`](super::mesh_lp).
//!
//! Источник: OpenCS `CScore/ShellSimpl/CapraMauryDesignSolver.cs` и `StripNdmAdapter.cs`,
//! коммит `387e157`.

use std::sync::Arc;

use super::diagram::Diagram;
use super::materials::{ConcreteChars, ConcreteDiagramType, RebarChars, StripLimits};
use super::mesh_lp;
use super::strip::{CrackParams, SigmaSCrcMethod, UlsTable};
use super::strip_design::{self, DesignResult, DesignStatus, SlsStrip};

/// Усилия пластины в точке: Nx, Ny, Nxy — кН/м; Mx, My, Mxy — кН·м/м.
/// «+» N — растяжение, «+» M растягивает верхнюю грань.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PlateForces {
    pub nx: f64,
    pub ny: f64,
    pub nxy: f64,
    pub mx: f64,
    pub my: f64,
    pub mxy: f64,
}

/// Геометрия плитного сечения для подбора.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateSection {
    /// Толщина, м.
    pub h: f64,
    /// Расстояние от верхней грани до центра верхнего ряда, м.
    pub a_top: f64,
    /// Расстояние от нижней грани до центра нижнего ряда, м.
    pub a_bot: f64,
    /// Диаметры стержней направлений x и y, м (для ПС2 берётся наибольший — в запас).
    pub ds_x: f64,
    pub ds_y: f64,
    pub concrete_diagram: ConcreteDiagramType,
}

impl PlateSection {
    pub fn z_top(&self) -> f64 {
        self.h / 2.0 - self.a_top
    }
    pub fn z_bot(&self) -> f64 {
        -(self.h / 2.0 - self.a_bot)
    }
}

/// Подобранная ортогональная сетка одной грани, м²/м.
#[derive(Clone, Debug, PartialEq)]
pub struct Mesh {
    pub asx: f64,
    pub asy: f64,
    /// Управляющие фасетки — активные ограничения ЛП, градусы.
    pub active_thetas: Vec<f64>,
}

/// Фасетка подбора: усилия, требования к рядам и обеспечиваемые площади, м²/м.
#[derive(Clone, Debug, PartialEq)]
pub struct DesignFacet {
    pub theta_deg: f64,
    pub n: f64,
    pub m: f64,
    pub as_top: f64,
    pub as_bot: f64,
    pub as_n_top: f64,
    pub as_n_bot: f64,
    pub status: DesignStatus,
    /// Кисп полосы при требуемых площадях (ПС1; в подборе по ПС2 — 0).
    pub eta: f64,
    /// a_crc при требуемых площадях, мм (ПС2; в подборе по ПС1 — 0).
    pub acrc_mm: f64,
    /// ПС2: As ряда у грани, сжатой моментом, из ПС1 той же фасетки.
    pub as_c_uls: f64,
}

/// Предельное состояние подбора.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Limit {
    Uls,
    Sls,
}

/// Результат подбора сеток.
#[derive(Clone, Debug)]
pub struct Design {
    pub top: Mesh,
    pub bottom: Mesh,
    pub facets: Vec<DesignFacet>,
    pub diagnostics: Vec<String>,
    pub step_deg: f64,
    pub limit: Limit,
    pub acrc_lim_mm: f64,
}

impl Design {
    /// Подбор возможен на всех фасетках.
    pub fn feasible(&self) -> bool {
        self.facets.iter().all(|f| f.status != DesignStatus::Impossible)
    }
    /// Сумма площадей четырёх направлений, м²/м.
    pub fn total(&self) -> f64 {
        self.top.asx + self.top.asy + self.bottom.asx + self.bottom.asy
    }
}

/// Контекст подбора ПС1 — общий для всех строк усилий одного сечения.
#[derive(Clone, Debug)]
pub struct UlsContext {
    pub h: f64,
    pub z_top: f64,
    pub z_bot: f64,
    table: UlsTable,
}

impl UlsContext {
    /// `concrete`, `rebar` — расчётные характеристики (ПС1).
    pub fn new(section: &PlateSection, concrete: &ConcreteChars, rebar: &RebarChars) -> Self {
        let cd = concrete.diagram(section.concrete_diagram);
        let rd = rebar.diagram();
        let lim = StripLimits::new(concrete, rebar);
        let (z_top, z_bot) = (section.z_top(), section.z_bot());
        Self { h: section.h, z_top, z_bot, table: strip_design::uls_table(section.h, z_top, z_bot, &cd, &rd, &lim) }
    }
}

/// Контекст подбора ПС2.
#[derive(Clone, Debug)]
pub struct SlsContext {
    pub h: f64,
    pub z_top: f64,
    pub z_bot: f64,
    /// Диаметр стержней обоих рядов в расчёте a_crc, м.
    pub ds: f64,
    pub acrc_lim_mm: f64,
    concrete: Arc<Diagram>,
    rebar: Arc<Diagram>,
    crack: CrackParams,
    uls_table: Option<UlsTable>,
}

/// Параметры подбора по трещинам.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlsOptions {
    pub acrc_lim_mm: f64,
    pub phi1: f64,
    pub phi2: f64,
    pub method: SigmaSCrcMethod,
}

impl Default for SlsOptions {
    fn default() -> Self {
        Self { acrc_lim_mm: 0.3, phi1: 1.0, phi2: 0.5, method: SigmaSCrcMethod::ReleasedConcrete8137 }
    }
}

impl SlsContext {
    /// `concrete_sls`, `rebar_sls` — нормативные характеристики (п. 6.1.26). `uls` —
    /// расчётные характеристики для отправной As_c из ПС1 той же фасетки; None — As_c = 0.
    pub fn new(
        section: &PlateSection,
        concrete_sls: &ConcreteChars,
        rebar_sls: &RebarChars,
        uls: Option<(&ConcreteChars, &RebarChars)>,
        opts: SlsOptions,
    ) -> Self {
        let (z_top, z_bot) = (section.z_top(), section.z_bot());
        let uls_table = uls.map(|(c, r)| {
            let lim = StripLimits::new(c, r);
            strip_design::uls_table(section.h, z_top, z_bot, &c.diagram(section.concrete_diagram), &r.diagram(), &lim)
        });
        Self {
            h: section.h,
            z_top,
            z_bot,
            ds: section.ds_x.max(section.ds_y),
            acrc_lim_mm: opts.acrc_lim_mm,
            concrete: Arc::new(concrete_sls.diagram(section.concrete_diagram)),
            rebar: Arc::new(rebar_sls.diagram()),
            crack: CrackParams { es: rebar_sls.es, phi1: opts.phi1, phi2: opts.phi2, method: opts.method },
            uls_table,
        }
    }

    pub fn with_uls(&self) -> bool {
        self.uls_table.is_some()
    }

    fn strip(&self, ds: f64) -> SlsStrip {
        SlsStrip {
            h: self.h,
            z_top: self.z_top,
            z_bot: self.z_bot,
            concrete: self.concrete.clone(),
            rebar: self.rebar.clone(),
            ds_top: ds,
            ds_bot: ds,
            crack: self.crack,
        }
    }
}

/// Фасетки 0…180° с шагом `step_deg`: (θ, N(θ), M(θ)).
fn facets(f: &PlateForces, step_deg: f64) -> Vec<(f64, f64, f64)> {
    let count = ((180.0 / step_deg).round() as usize).max(1);
    (0..count)
        .map(|i| {
            let theta = i as f64 * step_deg;
            let a = theta.to_radians();
            let (c, s) = (a.cos(), a.sin());
            let (c2, s2, cs) = (c * c, s * s, c * s);
            (theta, f.nx * c2 + f.ny * s2 + 2.0 * f.nxy * cs, f.mx * c2 + f.my * s2 + 2.0 * f.mxy * cs)
        })
        .collect()
}

fn c2s2(theta: f64) -> (f64, f64) {
    let a = theta.to_radians();
    (a.cos() * a.cos(), a.sin() * a.sin())
}

struct Raw {
    theta: f64,
    n: f64,
    m: f64,
    r: DesignResult,
    as_c_uls: f64,
}

fn meshes(raw: Vec<Raw>, diagnostics: Vec<String>, step_deg: f64, limit: Limit, acrc_lim_mm: f64) -> Design {
    let top = mesh_lp::solve(&raw.iter().map(|f| (f.theta, f.r.top())).collect::<Vec<_>>());
    let bot = mesh_lp::solve(&raw.iter().map(|f| (f.theta, f.r.bot())).collect::<Vec<_>>());
    let sls = limit == Limit::Sls;
    let facets = raw
        .iter()
        .map(|f| {
            let (c2, s2) = c2s2(f.theta);
            DesignFacet {
                theta_deg: f.theta,
                n: f.n,
                m: f.m,
                as_top: f.r.top(),
                as_bot: f.r.bot(),
                as_n_top: top.asx * c2 + top.asy * s2,
                as_n_bot: bot.asx * c2 + bot.asy * s2,
                status: f.r.status,
                eta: if sls { 0.0 } else { f.r.value },
                acrc_mm: if sls { f.r.value } else { 0.0 },
                as_c_uls: f.as_c_uls,
            }
        })
        .collect();
    Design {
        top: Mesh { asx: top.asx, asy: top.asy, active_thetas: top.active_thetas },
        bottom: Mesh { asx: bot.asx, asy: bot.asy, active_thetas: bot.active_thetas },
        facets,
        diagnostics,
        step_deg,
        limit,
        acrc_lim_mm: if sls { acrc_lim_mm } else { 0.0 },
    }
}

fn diagnostic(theta: f64, r: &DesignResult, n: f64, m: f64) -> String {
    format!("θ = {theta}°: {} (N = {n:.1} кН/м, M = {m:.2} кН·м/м)", r.reason)
}

/// Подбор по ПС1.
pub fn design_uls(ctx: &UlsContext, f: &PlateForces, step_deg: f64) -> Design {
    let mut raw = Vec::new();
    let mut diagnostics = Vec::new();
    for (theta, n, m) in facets(f, step_deg) {
        let r = strip_design::required_uls(n, m, &ctx.table);
        if r.status == DesignStatus::Impossible {
            diagnostics.push(diagnostic(theta, &r, n, m));
        }
        raw.push(Raw { theta, n, m, r, as_c_uls: 0.0 });
    }
    meshes(raw, diagnostics, step_deg, Limit::Uls, 0.0)
}

/// Подбор по ПС2: As ряда у растянутой моментом грани до a_crc ≤ a_crc,ult; ряд у второй
/// грани — As_c из ПС1 той же фасетки (если задано), при сквозном растяжении наращивается.
pub fn design_sls(ctx: &SlsContext, f: &PlateForces, step_deg: f64) -> Design {
    let mut raw = Vec::new();
    let mut diagnostics = Vec::new();
    let strip = ctx.strip(ctx.ds);
    for (theta, n, m) in facets(f, step_deg) {
        let mut as_c = 0.0;
        if let Some(t) = &ctx.uls_table {
            let u = strip_design::required_uls(n, m, t);
            if u.status != DesignStatus::Impossible {
                as_c = u.as_compression;
            }
        }
        let r = strip_design::required_sls(n, m, &strip, as_c, ctx.acrc_lim_mm);
        if r.status == DesignStatus::Impossible {
            diagnostics.push(diagnostic(theta, &r, n, m));
        }
        raw.push(Raw { theta, n, m, r, as_c_uls: as_c });
    }
    meshes(raw, diagnostics, step_deg, Limit::Sls, ctx.acrc_lim_mm)
}

/// Заданное армирование пластины, м²/м, и диаметры направлений, м.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ProvidedMeshes {
    pub asx_top: f64,
    pub asy_top: f64,
    pub asx_bot: f64,
    pub asy_bot: f64,
    pub ds_x: f64,
    pub ds_y: f64,
}

impl ProvidedMeshes {
    fn at(&self, theta: f64) -> (f64, f64, f64, f64) {
        let (c2, s2) = c2s2(theta);
        (self.asx_top * c2 + self.asy_top * s2, self.asx_bot * c2 + self.asy_bot * s2, c2, s2)
    }
}

/// Фасетка проверки через ограничения: требуемые площади рядов и обеспечиваемые сетками.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstraintFacet {
    pub theta_deg: f64,
    pub n: f64,
    pub m: f64,
    pub r_top: f64,
    pub r_bot: f64,
    pub as_n_top: f64,
    pub as_n_bot: f64,
    pub status: DesignStatus,
}

impl ConstraintFacet {
    fn ratio(r: f64, as_n: f64) -> f64 {
        if r <= 1e-12 {
            0.0
        } else if as_n <= 1e-12 {
            1e9
        } else {
            r / as_n
        }
    }
    pub fn eta_top(&self) -> f64 {
        Self::ratio(self.r_top, self.as_n_top)
    }
    pub fn eta_bot(&self) -> f64 {
        Self::ratio(self.r_bot, self.as_n_bot)
    }
    pub fn eta(&self) -> f64 {
        self.eta_top().max(self.eta_bot())
    }
}

/// Проверка ПС1 через ограничения: на каждой фасетке — наименьшая площадь растянутого ряда
/// при заданной площади второго; второй ряд — требование, только если не хватает и
/// наибольшей площади первого.
pub fn check_uls(ctx: &UlsContext, f: &PlateForces, p: &ProvidedMeshes, step_deg: f64) -> Vec<ConstraintFacet> {
    facets(f, step_deg)
        .into_iter()
        .map(|(theta, n, m)| {
            let (p_top, p_bot, _, _) = p.at(theta);
            let tension_top = m >= 0.0;
            let r = strip_design::required_tension_uls(n, m, &ctx.table, if tension_top { p_bot } else { p_top });
            let other = match r.status {
                DesignStatus::WithCompression | DesignStatus::Impossible => r.as_compression,
                _ => 0.0,
            };
            ConstraintFacet {
                theta_deg: theta,
                n,
                m,
                r_top: if tension_top { r.as_tension } else { other },
                r_bot: if tension_top { other } else { r.as_tension },
                as_n_top: p_top,
                as_n_bot: p_bot,
                status: r.status,
            }
        })
        .collect()
}

/// Средневзвешенный диаметр сетки на фасетке.
fn eff_ds(ax: f64, dx: f64, ay: f64, dy: f64, c2: f64, s2: f64) -> f64 {
    let w = ax * c2 + ay * s2;
    if w < 1e-14 {
        (dx + dy) / 2.0
    } else {
        (ax * c2 * dx + ay * s2 * dy) / w
    }
}

/// Проверка ПС2 через ограничения; d_s ряда — средневзвешенный по сетке фасетки.
pub fn check_sls(ctx: &SlsContext, f: &PlateForces, p: &ProvidedMeshes, step_deg: f64) -> Vec<ConstraintFacet> {
    facets(f, step_deg)
        .into_iter()
        .map(|(theta, n, m)| {
            let (p_top, p_bot, c2, s2) = p.at(theta);
            let tension_top = m >= 0.0;
            let p_other = if tension_top { p_bot } else { p_top };
            let ds = if tension_top {
                eff_ds(p.asx_top, p.ds_x, p.asy_top, p.ds_y, c2, s2)
            } else {
                eff_ds(p.asx_bot, p.ds_x, p.asy_bot, p.ds_y, c2, s2)
            };
            let r = strip_design::required_sls(n, m, &ctx.strip(ds), p_other, ctx.acrc_lim_mm);
            let other = if r.as_compression > p_other * (1.0 + 1e-9) + 1e-12 { r.as_compression } else { 0.0 };
            ConstraintFacet {
                theta_deg: theta,
                n,
                m,
                r_top: if tension_top { r.as_tension } else { other },
                r_bot: if tension_top { other } else { r.as_tension },
                as_n_top: p_top,
                as_n_bot: p_bot,
                status: r.status,
            }
        })
        .collect()
}

/// Параллельное отображение по строкам усилий: порядок результатов совпадает с порядком входа,
/// результат каждой строки не зависит от числа потоков.
pub fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(items.len().max(1));
    if threads <= 1 {
        return items.iter().map(&f).collect();
    }
    let chunk = items.len().div_ceil(threads);
    std::thread::scope(|scope| {
        let handles: Vec<_> = items
            .chunks(chunk)
            .map(|part| {
                let f = &f;
                scope.spawn(move || part.iter().map(f).collect::<Vec<R>>())
            })
            .collect();
        handles.into_iter().flat_map(|h| h.join().expect("поток подбора завершился с ошибкой")).collect()
    })
}

/// Подбор по ПС1 для массива строк усилий на всех ядрах.
pub fn design_uls_batch(ctx: &UlsContext, forces: &[PlateForces], step_deg: f64) -> Vec<Design> {
    par_map(forces, |f| design_uls(ctx, f, step_deg))
}

/// Подбор по ПС2 для массива строк усилий на всех ядрах.
pub fn design_sls_batch(ctx: &SlsContext, forces: &[PlateForces], step_deg: f64) -> Vec<Design> {
    par_map(forces, |f| design_sls(ctx, f, step_deg))
}
