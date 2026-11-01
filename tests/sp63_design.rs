//! Подбор полосы, трещины, ЛП сетки и Капра–Мори.
//! Перенос `StripNdmDesignTests`, `StripNdmSlsTests`, `OrthogonalMeshLpTests`,
//! `CapraMauryDesignTests`, `CapraMauryDesignSlsTests` (OpenCS, коммит `387e157`).
//!
//! Проверка подобранных сеток в C# идёт через `ShellSimplSolver` (проверка по фасеткам); здесь —
//! та же НДМ полосы на каждой фасетке при обеспечиваемых площадях.

mod common;
use common::*;
use femgeo::sp63::capra_maury::{self, Design, PlateForces, ProvidedMeshes, SlsContext, SlsOptions, UlsContext};
use femgeo::sp63::materials::{CalcType, ConcreteDiagramType};
use femgeo::sp63::mesh_lp;
use femgeo::sp63::strip::{self, SigmaSCrcMethod};
use femgeo::sp63::strip_design::{self, DesignResult, DesignStatus, SlsStrip};

fn uls(n: f64, m: f64) -> DesignResult {
    let c = concrete_diagram(CalcType::Uls, ConcreteDiagramType::L3);
    let r = rebar_diagram(CalcType::Uls);
    strip_design::required_uls(n, m, &strip_design::uls_table(H, Z, -Z, &c, &r, &limits()))
}

fn eta_at(r: &DesignResult, n: f64, m: f64, at: f64, ac: f64) -> f64 {
    let (t, b) = if r.tension_top { (at, ac) } else { (ac, at) };
    strip::ultimate_factor(&strip_as(CalcType::Uls, t, b), n, m, &limits()).eta
}

// ── подбор полосы ПС1 ──

#[test]
fn tension_only_found_area_is_minimal() {
    for (n, m) in [(0.0, 30.0), (0.0, -20.0), (-200.0, 40.0), (80.0, 25.0)] {
        let r = uls(n, m);
        assert_eq!(r.status, DesignStatus::TensionOnly, "N={n}, M={m}");
        assert_eq!(m >= 0.0, r.tension_top);
        assert!((0.99..=1.0).contains(&r.value), "η = {}", r.value);
        assert!(eta_at(&r, n, m, 0.99 * r.as_tension, 0.0) > 1.0);
    }
}

#[test]
fn no_reinforcement_when_concrete_carries_compression() {
    let r = uls(-1000.0, 5.0);
    assert_eq!(r.status, DesignStatus::NotRequired);
    assert_eq!(r.as_tension, 0.0);
    assert!(r.value <= 1.0);
}

/// Сквозное растяжение: точный минимум As_t + As_c = N/Rs.
#[test]
fn through_tension_both_rows_at_yield() {
    let (n, m) = (300.0, 5.0);
    let r = uls(n, m);
    let ft = (n + m / Z) / 2.0;
    let fc = n - ft;
    assert_eq!(r.status, DesignStatus::WithCompression);
    assert!(close(n / RS, r.as_tension + r.as_compression, 0.01 * n / RS));
    assert!(close(ft / RS, r.as_tension, 0.03 * ft / RS));
    assert!(close(fc / RS, r.as_compression, 0.03 * fc / RS));
    assert!(r.value <= 1.0);
}

#[test]
fn large_moment_needs_compression_reinforcement() {
    let r = uls(0.0, 200.0);
    assert_eq!(r.status, DesignStatus::WithCompression);
    assert!(r.as_compression > 0.0);
    assert!(r.value <= 1.0);
    assert!(eta_at(&r, 0.0, 200.0, r.as_tension, 0.0) > 1.0);
}

#[test]
fn impossible_when_thickness_insufficient() {
    let r = uls(-20_000.0, 0.0);
    assert_eq!(r.status, DesignStatus::Impossible);
    assert!(!r.reason.is_empty());
}

/// Ручная проверка порядка величины: изгиб M = 30 кН·м/м, h0 = 0,17 м, B25/A500.
/// Прямоугольная эпюра: x = h0 − √(h0² − 2M/Rb), As = Rb·x/Rs. НДМ даёт близкое значение.
#[test]
fn bending_area_close_to_rectangular_block() {
    let (m, h0, rb) = (30.0, H / 2.0 + Z, 14_500.0);
    let x = h0 - (h0 * h0 - 2.0 * m / rb).sqrt();
    let as_block = rb * x / RS;
    let r = uls(0.0, m);
    let dev = (r.as_tension - as_block).abs() / as_block;
    assert!(dev < 0.03, "НДМ {:.3} см²/м, эпюра {:.3} см²/м", r.as_tension * 1e4, as_block * 1e4);
}

// ── трещины ──

fn sls_strip(ct: CalcType) -> (femgeo::sp63::diagram::Diagram, femgeo::sp63::diagram::Diagram) {
    (concrete_diagram(ct, ConcreteDiagramType::L3), rebar_diagram(ct))
}

fn sls(n: f64, m: f64, as_c: f64, lim: f64) -> DesignResult {
    let (c, r) = sls_strip(CalcType::Sls);
    let p = SlsStrip { h: H, z_top: Z, z_bot: -Z, concrete: &c, rebar: &r, ds_top: DS, ds_bot: DS, crack: crack(SigmaSCrcMethod::ReleasedConcrete8137) };
    strip_design::required_sls(n, m, &p, as_c, lim)
}

fn acrc(r: &DesignResult, n: f64, m: f64, at: f64, ac: f64) -> f64 {
    let (t, b) = if r.tension_top { (at, ac) } else { (ac, at) };
    strip::crack_width(&strip_as(CalcType::Sls, t, b), n, m, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137)).acrc_mm
}

/// Упругий Wred приведённого сечения для растянутой верхней грани (п. 8.2.12) и ex = Wred/Ared.
fn elastic() -> (f64, f64) {
    let a_s = as_per_meter();
    let (h0, a, alpha) = (H / 2.0 + Z, H / 2.0 - Z, ES / EB);
    let ared = H + 2.0 * alpha * a_s;
    let sred = H * H / 2.0 + alpha * a_s * h0 + alpha * a_s * a;
    let y = sred / ared;
    let ired = H * H * H / 12.0 + H * (H / 2.0 - y).powi(2) + alpha * a_s * (h0 - y).powi(2) + alpha * a_s * (a - y).powi(2);
    let w = ired / (H - y);
    (w, w / ared)
}

/// Эталон OpenCS (18.10.2026): γ_экв = 1,724 при N = 0, убывает с растяжением.
#[test]
fn mcrc_equivalent_gamma_matches_opencs_reference() {
    let (w, ex) = elastic();
    let mut prev = f64::MAX;
    for n in [0.0, 100.0, 200.0, 350.0] {
        let r = strip::crack_width(&strip(CalcType::Sls), n, 40.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
        let gamma = (r.mcrc + n * ex) / (1_550.0 * w);
        if n == 0.0 {
            assert!((1.70..=1.90).contains(&gamma), "γ = {gamma}");
            assert!(close(gamma, 1.724, 0.002), "γ = {gamma} (OpenCS 1,724)");
        }
        assert!(gamma < prev, "N = {n}: γ = {gamma}");
        prev = gamma;
    }
}

#[test]
fn below_mcrc_no_crack() {
    let r = strip::crack_width(&strip(CalcType::Sls), 0.0, 5.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
    assert!(!r.cracked);
    assert_eq!(r.acrc_mm, 0.0);
    assert!(r.mcrc > 5.0);
}

#[test]
fn cracking_moment_variant_8138_psi_from_moment_ratio() {
    let r = strip::crack_width(&strip(CalcType::Sls), 0.0, 30.0, &[DS, DS], &crack(SigmaSCrcMethod::CrackingMoment8138));
    assert!(r.cracked);
    assert!(close(1.0 - 0.8 * r.mcrc / 30.0, r.rows[0].psi_s, 1e-9));
}

#[test]
fn psi_one_variant_is_conservative() {
    let base = strip::crack_width(&strip(CalcType::Sls), 0.0, 30.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
    let one = strip::crack_width(&strip(CalcType::Sls), 0.0, 30.0, &[DS, DS], &crack(SigmaSCrcMethod::PsiOne));
    assert_eq!(one.rows[0].psi_s, 1.0);
    assert!(one.acrc_mm >= base.acrc_mm);
}

#[test]
fn through_tension_acrc_for_both_rows() {
    let r = strip::crack_width(&strip(CalcType::Sls), 300.0, 5.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
    assert!(r.cracked);
    assert_eq!(r.state.as_ref().unwrap().compression_height(H), 0.0);
    assert!(r.rows.iter().all(|x| x.acrc_mm > 0.0));
    assert!(r.rows[0].acrc_mm > r.rows[1].acrc_mm);
    assert_eq!(r.rows.iter().map(|x| x.acrc_mm).fold(0.0, f64::max), r.acrc_mm);
    assert!(r.rows.iter().all(|x| x.phi3 == 1.2));
}

#[test]
fn negative_moment_cracks_bottom_face() {
    let r = strip::crack_width(&strip(CalcType::Sls), 0.0, -30.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
    assert!(!r.crack_face_top);
    assert!(r.rows[1].acrc_mm > 0.0);
    assert_eq!(r.rows[0].acrc_mm, 0.0);
}

/// Несходимость сечения с трещиной не выдаётся за успех.
#[test]
fn unconverged_cracked_section_is_not_success() {
    let r = strip::crack_width(&strip(CalcType::Sls), 0.0, 500.0, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
    assert!(r.cracked);
    assert!(!r.converged);
    assert!(!r.reason.is_empty());
}

#[test]
fn sls_found_area_is_minimal() {
    for (n, m, lim) in [(0.0, 30.0, 0.3), (0.0, 30.0, 0.1), (-100.0, -32.0, 0.2)] {
        let r = sls(n, m, as_per_meter(), lim);
        assert_eq!(r.status, DesignStatus::TensionOnly, "N={n}, M={m}, lim={lim}");
        assert!(r.value <= lim);
        assert!(acrc(&r, n, m, 0.99 * r.as_tension, r.as_compression) > lim);
    }
}

#[test]
fn sls_no_crack_not_required() {
    let r = sls(0.0, 5.0, 0.0, 0.3);
    assert_eq!(r.status, DesignStatus::NotRequired);
    assert_eq!(r.as_tension, 0.0);
}

#[test]
fn sls_through_tension_grows_both_rows() {
    let r = sls(300.0, 5.0, 0.0, 0.2);
    assert_eq!(r.status, DesignStatus::WithCompression);
    assert!(r.as_compression > 0.0 && r.as_compression < r.as_tension);
    assert!(r.value <= 0.2);
}

// ── ЛП ортогональной сетки ──

fn facet_demands(r: impl Fn(f64) -> f64, step: f64) -> Vec<(f64, f64)> {
    let n = (180.0 / step).round() as usize;
    (0..n).map(|i| (i as f64 * step, r(i as f64 * step))).collect()
}

fn provided(asx: f64, asy: f64, theta: f64) -> f64 {
    let a = theta.to_radians();
    asx * a.cos() * a.cos() + asy * a.sin() * a.sin()
}

fn grid_minimum(d: &[(f64, f64)], x_max: f64, n: usize) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..=n {
        let x = x_max * i as f64 / n as f64;
        let mut y: f64 = 0.0;
        let mut ok = true;
        for &(theta, r) in d {
            let a = theta.to_radians();
            let (c2, s2) = (a.cos() * a.cos(), a.sin() * a.sin());
            let rest = r - c2 * x;
            if rest <= 0.0 {
                continue;
            }
            if s2 < 1e-12 {
                ok = false;
                break;
            }
            y = y.max(rest / s2);
        }
        if ok {
            best = best.min(x + y);
        }
    }
    best
}

#[test]
fn lp_special_cases() {
    let s = mesh_lp::solve(&facet_demands(|t| if t == 0.0 { 5.0 } else { 0.0 }, 10.0));
    assert!(close(s.asx, 5.0, 1e-9) && close(s.asy, 0.0, 1e-9));
    assert_eq!(s.active_thetas, vec![0.0]);

    let s = mesh_lp::solve(&facet_demands(|_| 3.0, 10.0));
    assert!(close(s.asx, 3.0, 1e-9) && close(s.asy, 3.0, 1e-9));
    assert!(s.active_thetas.contains(&0.0) && s.active_thetas.contains(&90.0));

    let s = mesh_lp::solve(&facet_demands(|_| -1.0, 10.0));
    assert_eq!((s.asx, s.asy), (0.0, 0.0));
    assert!(s.active_thetas.is_empty());

    let s = mesh_lp::solve(&[(45.0, 2.0)]);
    assert!(close(s.asx + s.asy, 4.0, 1e-9));
    assert!(close(s.asx, s.asy, 1e-9));
}

#[test]
fn lp_random_demands_match_grid_search_and_are_minimal() {
    for seed in 1..=8u64 {
        let mut rnd = Lcg::new(seed);
        let a = rnd.next_f64() * 5.0;
        let b = rnd.next_f64() * 6.0 - 3.0;
        let c = rnd.next_f64() * 6.0 - 3.0;
        let step: f64 = if rnd.next_f64() < 0.5 { 10.0 } else { 5.0 };
        let n = (180.0 / step).round() as usize;
        let d: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let t = i as f64 * step;
                let r = t.to_radians();
                (t, a + b * (2.0 * r).cos() + c * (2.0 * r).sin() + rnd.next_f64() - 0.5)
            })
            .collect();
        let s = mesh_lp::solve(&d);
        for &(theta, r) in &d {
            assert!(provided(s.asx, s.asy, theta) >= r - 1e-9, "seed {seed}, θ = {theta}");
        }
        let grid = grid_minimum(&d, 20.0, 200_000);
        assert!(s.asx + s.asy <= grid + 1e-9, "seed {seed}: ЛП хуже перебора");
        assert!(s.asx + s.asy >= grid - 1e-3, "seed {seed}: ЛП лучше перебора");
        if s.asx > 1e-9 {
            assert!(d.iter().any(|f| provided(0.99 * s.asx, s.asy, f.0) < f.1 - 1e-12));
        }
        if s.asy > 1e-9 {
            assert!(d.iter().any(|f| provided(s.asx, 0.99 * s.asy, f.0) < f.1 - 1e-12));
        }
        assert!(!s.active_thetas.is_empty());
    }
}

// ── Капра–Мори ПС1 ──

fn uls_ctx() -> UlsContext {
    UlsContext::new(&plate(), &concrete(CalcType::Uls), &rebar(CalcType::Uls))
}

fn forces(nx: f64, ny: f64, nxy: f64, mx: f64, my: f64, mxy: f64) -> PlateForces {
    PlateForces { nx, ny, nxy, mx, my, mxy }
}

const CASES_ULS: [[f64; 6]; 6] = [
    [0.0, 0.0, 0.0, 30.0, 10.0, 8.0],
    [0.0, 0.0, 0.0, -25.0, 15.0, -12.0],
    [-300.0, -200.0, 50.0, 40.0, 20.0, 15.0],
    [150.0, 80.0, 40.0, 20.0, -10.0, 6.0],
    [400.0, 300.0, 100.0, 5.0, 3.0, 2.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 25.0],
];

fn assert_minimal(d: &Design) {
    for top in [true, false] {
        let mesh = if top { &d.top } else { &d.bottom };
        let req = |f: &capra_maury::DesignFacet| if top { f.as_top } else { f.as_bot };
        if mesh.asx > 1e-9 {
            assert!(d.facets.iter().any(|f| provided(0.99 * mesh.asx, mesh.asy, f.theta_deg) < req(f) - 1e-12));
        }
        if mesh.asy > 1e-9 {
            assert!(d.facets.iter().any(|f| provided(mesh.asx, 0.99 * mesh.asy, f.theta_deg) < req(f) - 1e-12));
        }
    }
}

#[test]
fn capra_maury_uls_meshes_pass_facet_check_and_are_minimal() {
    let ctx = uls_ctx();
    let lim = limits();
    for c in CASES_ULS {
        let f = forces(c[0], c[1], c[2], c[3], c[4], c[5]);
        let d = capra_maury::design_uls(&ctx, &f, 10.0);
        assert!(d.feasible() && d.diagnostics.is_empty(), "{c:?}");
        let mut worst: f64 = 0.0;
        for fc in &d.facets {
            assert!(fc.as_n_top >= fc.as_top * (1.0 - 1e-9) - 1e-12, "{c:?} θ = {}: верх", fc.theta_deg);
            assert!(fc.as_n_bot >= fc.as_bot * (1.0 - 1e-9) - 1e-12, "{c:?} θ = {}: низ", fc.theta_deg);
            let eta = strip::ultimate_factor(&strip_as(CalcType::Uls, fc.as_n_top, fc.as_n_bot), fc.n, fc.m, &lim).eta;
            worst = worst.max(eta);
        }
        assert!(worst <= 1.005, "{c:?}: η = {worst}");
        assert_minimal(&d);

        // Проверка через ограничения при подобранных сетках — Кисп ≤ 1 с точностью подбора.
        let p = ProvidedMeshes { asx_top: d.top.asx, asy_top: d.top.asy, asx_bot: d.bottom.asx, asy_bot: d.bottom.asy, ds_x: DS, ds_y: DS };
        let chk = capra_maury::check_uls(&ctx, &f, &p, 10.0);
        let eta_c = chk.iter().map(|x| x.eta()).fold(0.0, f64::max);
        assert!(eta_c <= 1.0 + 1e-6, "{c:?}: проверка η = {eta_c}");
    }
}

#[test]
fn capra_maury_symmetry_sign_of_twist_and_shear() {
    let ctx = uls_ctx();
    let a = capra_maury::design_uls(&ctx, &forces(-100.0, 50.0, 40.0, 20.0, -8.0, 12.0), 10.0);
    let b = capra_maury::design_uls(&ctx, &forces(-100.0, 50.0, -40.0, 20.0, -8.0, -12.0), 10.0);
    for (x, y) in [(a.top.asx, b.top.asx), (a.top.asy, b.top.asy), (a.bottom.asx, b.bottom.asx), (a.bottom.asy, b.bottom.asy)] {
        assert!(close(x, y, 1e-9));
    }
}

/// Поворот осей на 90° (x ↔ y) меняет местами Asx и Asy.
#[test]
fn capra_maury_axis_swap_invariance() {
    let ctx = uls_ctx();
    let a = capra_maury::design_uls(&ctx, &forces(-120.0, 60.0, 25.0, 22.0, -9.0, 11.0), 10.0);
    let b = capra_maury::design_uls(&ctx, &forces(60.0, -120.0, -25.0, -9.0, 22.0, -11.0), 10.0);
    assert!(close(a.top.asx, b.top.asy, 1e-9) && close(a.top.asy, b.top.asx, 1e-9));
    assert!(close(a.bottom.asx, b.bottom.asy, 1e-9) && close(a.bottom.asy, b.bottom.asx, 1e-9));
}

#[test]
fn capra_maury_convergence_with_step() {
    let ctx = uls_ctx();
    let f = forces(-50.0, 20.0, 30.0, 25.0, -5.0, 14.0);
    let s10 = capra_maury::design_uls(&ctx, &f, 10.0).total();
    let s5 = capra_maury::design_uls(&ctx, &f, 5.0).total();
    let s25 = capra_maury::design_uls(&ctx, &f, 2.5).total();
    assert!(s5 >= s10 * (1.0 - 1e-9));
    assert!(s25 >= s5 * (1.0 - 1e-9));
    assert!((s25 - s5).abs() <= (s5 - s10).abs().max(0.002 * s25) + 1e-12);
    assert!((s25 - s10).abs() <= 0.03 * s25);
}

#[test]
fn capra_maury_no_twist_matches_single_strips() {
    let ctx = uls_ctx();
    for (nx, ny, mx, my) in [(0.0, 0.0, 30.0, 12.0), (0.0, 0.0, -20.0, 35.0), (-200.0, -100.0, 25.0, -15.0)] {
        let d = capra_maury::design_uls(&ctx, &forces(nx, ny, 0.0, mx, my, 0.0), 10.0);
        let (rx, ry) = (uls(nx, mx), uls(ny, my));
        assert!(close(rx.top(), d.top.asx, 1e-8));
        assert!(close(ry.top(), d.top.asy, 1e-8));
        assert!(close(rx.bot(), d.bottom.asx, 1e-8));
        assert!(close(ry.bot(), d.bottom.asy, 1e-8));
    }
}

#[test]
fn capra_maury_top_and_bottom_follow_moment_sign() {
    let ctx = uls_ctx();
    let up = capra_maury::design_uls(&ctx, &forces(0.0, 0.0, 0.0, 30.0, 20.0, 0.0), 10.0);
    assert!(up.top.asx > 0.0 && up.top.asy > 0.0);
    assert_eq!((up.bottom.asx, up.bottom.asy), (0.0, 0.0));
    let down = capra_maury::design_uls(&ctx, &forces(0.0, 0.0, 0.0, -30.0, -20.0, 0.0), 10.0);
    assert!(close(up.top.asx, down.bottom.asx, 1e-9) && close(up.top.asy, down.bottom.asy, 1e-9));
    assert_eq!(down.top.asx, 0.0);
}

#[test]
fn capra_maury_impossible_and_zero() {
    let ctx = uls_ctx();
    let d = capra_maury::design_uls(&ctx, &forces(-20_000.0, -20_000.0, 0.0, 0.0, 0.0, 0.0), 10.0);
    assert!(!d.feasible());
    assert!(!d.diagnostics.is_empty());
    let z = capra_maury::design_uls(&ctx, &PlateForces::default(), 10.0);
    assert_eq!(z.total(), 0.0);
    assert!(z.feasible());
}

// ── Капра–Мори ПС2 ──

fn sls_ctx(lim: f64, with_uls: bool) -> SlsContext {
    let (cu, ru) = (concrete(CalcType::Uls), rebar(CalcType::Uls));
    SlsContext::new(
        &plate(),
        &concrete(CalcType::Sls),
        &rebar(CalcType::Sls),
        if with_uls { Some((&cu, &ru)) } else { None },
        SlsOptions { acrc_lim_mm: lim, ..Default::default() },
    )
}

const CASES_SLS: [[f64; 6]; 6] = [
    [0.0, 0.0, 0.0, 30.0, 10.0, 8.0],
    [0.0, 0.0, 0.0, -25.0, 15.0, -12.0],
    [-100.0, -50.0, 20.0, 40.0, 20.0, 15.0],
    [150.0, 80.0, 40.0, 20.0, -10.0, 6.0],
    [400.0, 300.0, 100.0, 5.0, 3.0, 2.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 25.0],
];

#[test]
fn capra_maury_sls_crack_width_within_limit_and_minimal() {
    let ctx = sls_ctx(0.3, true);
    for c in CASES_SLS {
        let f = forces(c[0], c[1], c[2], c[3], c[4], c[5]);
        let d = capra_maury::design_sls(&ctx, &f, 10.0);
        assert!(d.feasible() && d.diagnostics.is_empty(), "{c:?}");
        assert!(d.total() > 0.0);
        assert_eq!(d.acrc_lim_mm, 0.3);
        for fc in &d.facets {
            assert!(fc.as_n_top >= fc.as_top * (1.0 - 1e-9) - 1e-12);
            assert!(fc.as_n_bot >= fc.as_bot * (1.0 - 1e-9) - 1e-12);
            assert!(fc.acrc_mm <= 0.3 + 1e-9);
            let s = strip::crack_width(&strip_as(CalcType::Sls, fc.as_n_top, fc.as_n_bot), fc.n, fc.m, &[DS, DS], &crack(SigmaSCrcMethod::ReleasedConcrete8137));
            assert!(s.converged || !s.cracked, "{c:?} θ = {}: {}", fc.theta_deg, s.reason);
            assert!(s.acrc_mm <= 0.3 * 1.005, "{c:?} θ = {}: a_crc = {}", fc.theta_deg, s.acrc_mm);
        }
        assert_minimal(&d);
    }
}

#[test]
fn capra_maury_sls_compression_row_starts_from_uls() {
    let with = capra_maury::design_sls(&sls_ctx(0.3, true), &forces(400.0, 300.0, 100.0, 5.0, 3.0, 2.0), 10.0);
    for f in &with.facets {
        let opposite = if f.m >= 0.0 { f.as_bot } else { f.as_top };
        assert!(opposite >= f.as_c_uls * (1.0 - 1e-9) - 1e-12);
    }
    assert!(with.facets.iter().any(|f| f.as_c_uls > 1e-6));
    let ctx = sls_ctx(0.3, false);
    assert!(!ctx.with_uls());
    let without = capra_maury::design_sls(&ctx, &forces(400.0, 300.0, 100.0, 5.0, 3.0, 2.0), 10.0);
    assert!(without.facets.iter().all(|f| f.as_c_uls == 0.0));
}

#[test]
fn capra_maury_sls_stricter_limit_needs_more_steel() {
    let f = forces(0.0, 0.0, 0.0, 30.0, 10.0, 8.0);
    let d03 = capra_maury::design_sls(&sls_ctx(0.3, true), &f, 10.0);
    let d01 = capra_maury::design_sls(&sls_ctx(0.1, true), &f, 10.0);
    assert!(d01.total() > d03.total());
}

#[test]
fn capra_maury_sls_no_twist_matches_single_strips() {
    let d = capra_maury::design_sls(&sls_ctx(0.3, false), &forces(0.0, 0.0, 0.0, 30.0, -20.0, 0.0), 10.0);
    let rx = sls(0.0, 30.0, 0.0, 0.3);
    let ry = sls(0.0, -20.0, 0.0, 0.3);
    assert!(close(rx.as_tension, d.top.asx, 1e-8));
    assert!(close(ry.as_tension, d.bottom.asy, 1e-8));
    assert!(close(0.0, d.top.asy, 1e-12));
    assert!(close(0.0, d.bottom.asx, 1e-12));
}

#[test]
fn capra_maury_sls_zero_and_impossible() {
    let z = capra_maury::design_sls(&sls_ctx(0.3, true), &PlateForces::default(), 10.0);
    assert_eq!(z.total(), 0.0);
    assert!(z.feasible());
    assert!(z.facets.iter().all(|f| f.status == DesignStatus::NotRequired));
    let d = capra_maury::design_sls(&sls_ctx(0.001, true), &forces(0.0, 0.0, 0.0, 120.0, 0.0, 0.0), 10.0);
    assert!(!d.feasible());
    assert!(!d.diagnostics.is_empty());
}

/// Проверка ПС2 через ограничения при подобранных сетках проходит.
#[test]
fn capra_maury_sls_check_of_designed_meshes() {
    let ctx = sls_ctx(0.3, true);
    let f = forces(0.0, 0.0, 0.0, 30.0, 10.0, 8.0);
    let d = capra_maury::design_sls(&ctx, &f, 10.0);
    let p = ProvidedMeshes { asx_top: d.top.asx, asy_top: d.top.asy, asx_bot: d.bottom.asx, asy_bot: d.bottom.asy, ds_x: DS, ds_y: DS };
    let eta = capra_maury::check_sls(&ctx, &f, &p, 10.0).iter().map(|x| x.eta()).fold(0.0, f64::max);
    assert!(eta <= 1.0 + 1e-6, "η = {eta}");
}

/// Замер, не утверждение: `cargo test --release -- --ignored --nocapture timing`.
#[test]
#[ignore]
fn timing_1000_rows_uls_and_100_rows_sls() {
    let ctx = uls_ctx();
    let t = std::time::Instant::now();
    for i in 0..1000 {
        let i = i as f64;
        capra_maury::design_uls(&ctx, &forces(-100.0 + 0.2 * i, 50.0 - 0.1 * i, 20.0, 10.0 + 0.03 * i, -5.0 + 0.02 * i, 8.0), 10.0);
    }
    println!("ПС1: 1000 строк × 18 фасеток: {:?}", t.elapsed());
    let ctx = sls_ctx(0.3, true);
    let t = std::time::Instant::now();
    for i in 0..100 {
        let i = i as f64;
        capra_maury::design_sls(&ctx, &forces(-100.0 + 2.0 * i, 50.0 - i, 20.0, 10.0 + 0.3 * i, -5.0 + 0.2 * i, 8.0), 10.0);
    }
    println!("ПС2: 100 строк × 18 фасеток: {:?}", t.elapsed());
}
