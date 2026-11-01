//! НДМ полосы: интегрирование, равновесие, ПС1 и таблица контура.
//! Перенос `StripNdmSectionTests`, `StripNdmSolveTests`, `StripNdmUlsTests`,
//! `StripNdmUlsTableTests` (OpenCS, коммит `387e157`).

mod common;
use common::*;
use femgeo::sp63::materials::{CalcType, ConcreteDiagramType};
use femgeo::sp63::strip::{self, Pivot};

// ── интегрирование ──

#[test]
fn elastic_range_matches_analytical_stiffness() {
    let s = strip(CalcType::Uls);
    let (e0, k) = (-5e-5, 2e-4);
    let f = s.integrate(e0, k, false);
    let a = as_per_meter();
    let mut n = EB * H * e0;
    let mut m = EB * H * H * H / 12.0 * k;
    for z in [Z, -Z] {
        n += ES * a * (e0 + k * z);
        m += ES * a * (e0 + k * z) * z;
    }
    assert!(close(n, f.n, n.abs() * 1e-9));
    assert!(close(m, f.m, m.abs() * 1e-9));
    assert!(close(EB * H + 2.0 * ES * a, f.k00, 1e-6 * f.k00));
}

#[test]
fn tangent_matches_finite_differences() {
    for (e0, k, ten) in [(0.0004, 0.011, false), (-0.0011, 0.004, false), (-0.00002, 0.0009, true)] {
        let s = strip(CalcType::Uls);
        let f = s.integrate(e0, k, ten);
        let (de, dk) = (1e-9, 1e-8);
        let (fe, fe_) = (s.integrate(e0 + de, k, ten), s.integrate(e0 - de, k, ten));
        let (fk, fk_) = (s.integrate(e0, k + dk, ten), s.integrate(e0, k - dk, ten));
        let k00 = (fe.n - fe_.n) / (2.0 * de);
        let k10 = (fe.m - fe_.m) / (2.0 * de);
        let k01 = (fk.n - fk_.n) / (2.0 * dk);
        let k11 = (fk.m - fk_.m) / (2.0 * dk);
        assert!(close(k00, f.k00, k00.abs() * 1e-3 + 1.0), "k00 {e0} {k}");
        assert!(close(k10, f.k01, f.k00.abs() * H * 1e-3), "k10 {e0} {k}");
        assert!(close(k01, f.k01, f.k00.abs() * H * 1e-3), "k01 {e0} {k}");
        assert!(close(k11, f.k11, k11.abs() * 1e-3 + 1e-3), "k11 {e0} {k}");
    }
}

/// Аналитика вместо сравнения со стержневым НДМ OpenCS: при чистом растяжении бетон выключен,
/// ряды работают по закону Гука, пока не потекли.
#[test]
fn rebar_only_response_in_tension() {
    let s = strip(CalcType::Uls);
    let e0 = 0.001;
    let f = s.integrate(e0, 0.0, false);
    assert!(close(f.n, 2.0 * ES * as_per_meter() * e0, 1e-6));
    assert!(close(f.m, 0.0, 1e-9));
}

// ── равновесие ──

#[test]
fn solve_reproduces_target_forces() {
    for (n, m) in [(0.0, 30.0), (0.0, -25.0), (-800.0, 20.0), (-2500.0, 5.0), (100.0, 20.0), (300.0, 5.0)] {
        let s = strip(CalcType::Uls);
        let st = strip::solve(&s, n, m, false);
        assert!(st.converged, "N={n}, M={m}: {}", st.reason);
        let f = s.integrate(st.eps0, st.kappa, false);
        assert!(close(n, f.n, 1e-3), "N={n}, M={m}");
        assert!(close(m, f.m, 1e-4), "N={n}, M={m}");
    }
}

#[test]
fn solve_through_tension_both_rows_tensile() {
    let s = strip(CalcType::Uls);
    let st = strip::solve(&s, 300.0, 5.0, false);
    assert!(st.converged, "{}", st.reason);
    assert!(st.eps_at(Z) > 0.0 && st.eps_at(-Z) > 0.0);
    assert_eq!(st.compression_height(H), 0.0);
    let ft = (300.0 + 5.0 / Z) / 2.0;
    let fc = 300.0 - ft;
    assert!(close(ft / as_per_meter() / ES, st.eps_at(Z), 1e-7));
    assert!(close(fc / as_per_meter() / ES, st.eps_at(-Z), 1e-7));
}

#[test]
fn solve_bending_compression_height_on_compressed_face() {
    let st = strip::solve(&strip(CalcType::Uls), 0.0, 30.0, false);
    assert!(st.eps_at(H / 2.0) > 0.0 && st.eps_at(-H / 2.0) < 0.0);
    let x = st.compression_height(H);
    assert!((0.005..=0.05).contains(&x));
    assert!(close(-H / 2.0 + x, st.neutral_axis(), 1e-9));
}

#[test]
fn solve_beyond_capacity_not_converged_with_reason() {
    for (n, m) in [(0.0, 500.0), (800.0, 0.0), (-6000.0, 0.0)] {
        let st = strip::solve(&strip(CalcType::Uls), n, m, false);
        assert!(!st.converged, "N={n}, M={m}");
        assert!(!st.reason.trim().is_empty());
    }
}

// ── ПС1 ──

/// Эталон OpenCS (18.10.2026): η НДМ L2 = 0,4968, L3 = 0,4960; шарнир B у нижней грани.
#[test]
fn pure_bending_eta_matches_opencs_reference() {
    let lim = limits();
    for (t, lo, hi) in [(ConcreteDiagramType::L2, 0.4963, 0.4973), (ConcreteDiagramType::L3, 0.4955, 0.4965)] {
        let s = strip_t(CalcType::Uls, t, as_per_meter(), as_per_meter());
        let u = strip::ultimate_factor(&s, 0.0, 20.0, &lim);
        assert_eq!(u.pivot, Pivot::B, "{t:?}");
        assert!(!u.top_face);
        assert!((lo..=hi).contains(&u.eta), "{t:?}: η = {}", u.eta);
        let es = u.ultimate.unwrap().eps_at(Z);
        assert!((0.01..=0.025).contains(&es));
    }
}

#[test]
fn ultimate_state_lies_on_boundary_and_on_ray() {
    let lim = limits();
    for (n, m) in [(0.0, 20.0), (0.0, -15.0), (-800.0, 20.0), (-1500.0, -10.0), (100.0, 15.0), (250.0, 3.0)] {
        let s = strip(CalcType::Uls);
        let u = strip::ultimate_factor(&s, n, m, &lim);
        let ult = u.ultimate.clone().expect("предельное состояние");
        assert!(close(1.0, strip::strain_utilization(&s, &ult, &lim), 1e-4), "N={n}, M={m}");

        let f = s.integrate(ult.eps0, ult.kappa, false);
        let scale = n.abs() * H / 6.0 + m.abs();
        assert!(close(u.lambda * n * H / 6.0, f.n * H / 6.0, 1e-6 * scale * u.lambda));
        assert!(close(u.lambda * m, f.m, 1e-6 * scale * u.lambda));

        let below = strip::solve(&s, 0.99 * u.lambda * n, 0.99 * u.lambda * m, false);
        assert!(below.converged, "{}", below.reason);
        assert!(strip::strain_utilization(&s, &below, &lim) <= 1.0);
        let above = strip::solve(&s, 1.02 * u.lambda * n, 1.02 * u.lambda * m, false);
        assert!(!above.converged || strip::strain_utilization(&s, &above, &lim) > 1.0);
    }
}

#[test]
fn pivots_for_compression_and_tension() {
    let s = strip(CalcType::Uls);
    let lim = limits();
    let t = strip::ultimate_factor(&s, 100.0, 15.0, &lim);
    assert_eq!(t.pivot, Pivot::A);
    assert!(t.top_face);
    let c = strip::ultimate_factor(&s, -1500.0, 10.0, &lim);
    assert!(matches!(c.pivot, Pivot::B | Pivot::C));
    assert!(!c.top_face);
}

#[test]
fn central_compression_lambda_is_squash_load_with_eps_b0() {
    let u = strip::ultimate_factor(&strip(CalcType::Uls), -2000.0, 0.0, &limits());
    let nult = 14_500.0 * H + 2.0 * 400_000.0 * as_per_meter();
    assert_eq!(u.pivot, Pivot::C);
    assert!(close(nult / 2000.0, u.lambda, 1e-6));
}

#[test]
fn zero_load_eta_is_zero() {
    let u = strip::ultimate_factor(&strip(CalcType::Uls), 0.0, 0.0, &limits());
    assert_eq!(u.eta, 0.0);
    assert_eq!(u.pivot, Pivot::None);
}

#[test]
fn eta_is_continuous_across_through_tension_transition() {
    let s = strip(CalcType::Uls);
    let lim = limits();
    let mut prev = f64::NAN;
    let mut max_jump: f64 = 0.0;
    let mut m = 0.0;
    while m <= 20.0001 {
        let eta = strip::ultimate_factor(&s, 200.0, m, &lim).eta;
        if !prev.is_nan() {
            max_jump = max_jump.max((eta - prev).abs());
        }
        prev = eta;
        m += 0.25;
    }
    assert!(max_jump < 0.03, "скачок {max_jump}");
}

// ── таблица контура ──

#[test]
fn table_eta_matches_direct_ray() {
    let lim = limits();
    let a = as_per_meter();
    let table = strip::tabulate(&strip_as(CalcType::Uls, 0.0, 0.0), &lim);
    for (n, m, kt, kb) in [
        (0.0, 30.0, 1.0, 1.0),
        (0.0, -45.0, 1.0, 1.0),
        (-800.0, 40.0, 1.0, 0.5),
        (-3000.0, 10.0, 2.0, 2.0),
        (150.0, 10.0, 1.0, 1.0),
        (250.0, -3.0, 0.7, 1.4),
        (60.0, 25.0, 3.0, 0.0),
        (-200.0, 120.0, 4.0, 2.0),
        (0.0, 15.0, 0.0, 1.0),
    ] {
        let (at, ab) = (kt * a, kb * a);
        let direct = strip::ultimate_factor(&strip_as(CalcType::Uls, at, ab), n, m, &lim).eta;
        let fast = table.eta(&[at, ab], n, m);
        assert!(close(direct, fast, (1e-3 * direct).max(1e-9)), "N={n}, M={m}: {direct} vs {fast}");
    }
}

#[test]
fn table_zero_load_eta_is_zero() {
    let table = strip::tabulate(&strip(CalcType::Uls), &limits());
    assert_eq!(table.eta(&[as_per_meter(), as_per_meter()], 0.0, 0.0), 0.0);
}

// ── инвариантность ──

/// Зеркальная полоса (верх ↔ низ, знак M) даёт тот же Кисп.
#[test]
fn mirror_symmetry_of_strip() {
    let lim = limits();
    let a = as_per_meter();
    for (n, m, kt, kb) in [(-300.0, 25.0, 1.0, 2.0), (80.0, 10.0, 3.0, 1.0), (0.0, 40.0, 2.0, 0.5)] {
        let s1 = strip_as(CalcType::Uls, kt * a, kb * a);
        let s2 = strip_as(CalcType::Uls, kb * a, kt * a);
        let e1 = strip::ultimate_factor(&s1, n, m, &lim).eta;
        let e2 = strip::ultimate_factor(&s2, n, -m, &lim).eta;
        assert!(close(e1, e2, 1e-6 * e1), "{e1} vs {e2}");
    }
}
