//! Поперечная сила в пластине. Перенос `PlateShearCheckTests` (OpenCS, коммит `387e157`).
//! h = 0,2 м, h0 = 0,17 м, B25 (вид C): Rb = 14 500 кПа, Rbt = 1 050 кПа.

mod common;
use common::*;
use femgeo::sp63::capra_maury::PlateForces;
use femgeo::sp63::materials::{CalcType, RebarClass};
use femgeo::sp63::plate_shear::{self, PlateShear};

const H0: f64 = 0.17;
const RB: f64 = 14_500.0;
const RBT: f64 = 1_050.0;

fn chk(qx: f64, qy: f64, f: PlateForces) -> Option<PlateShear> {
    plate_shear::check(qx, qy, &f, &plate(), &concrete(CalcType::Uls), &rebar(CalcType::Uls), 300_000.0, 2.0 * as_per_meter())
}

fn nf(nx: f64) -> PlateForces {
    PlateForces { nx, ..Default::default() }
}

#[test]
fn zero_shear_gives_none() {
    assert!(chk(0.0, 0.0, nf(0.0)).is_none());
    assert!(chk(1e-12, -1e-12, nf(0.0)).is_none());
}

#[test]
fn direction_and_facet_forces() {
    let f = PlateForces { nx: 100.0, ny: 200.0, nxy: 50.0, mx: 10.0, my: 20.0, mxy: 5.0 };
    let r = chk(30.0, 40.0, f).unwrap();
    let (c, s) = (0.6, 0.8);
    assert!(close(r.q0, 50.0, 1e-9));
    assert!(close(r.theta_deg, 53.130102, 1e-5));
    assert!(close(r.n, 100.0 * c * c + 200.0 * s * s + 100.0 * c * s, 1e-9));
    assert!(close(r.m, 10.0 * c * c + 20.0 * s * s + 10.0 * c * s, 1e-9));
}

#[test]
fn no_axial_force_concrete_carries_shear() {
    let r = chk(60.0, 0.0, nf(0.0)).unwrap();
    assert_eq!(r.phi_n, 1.0);
    assert!(close(r.h0, H0, 1e-9));
    assert!(close(r.q_strip, 0.3 * RB * H0, 1e-6)); // 739,5 кН/м
    assert!(close(r.qb1, 0.5 * RBT * H0, 1e-6)); // 89,25 кН/м
    assert!(close(r.eta(), 60.0 / 89.25, 1e-9));
    assert!(!r.stirrups_required());
    assert_eq!((r.qsw_required, r.asw_per_sw), (0.0, 0.0));
}

#[test]
fn stirrups_by_8_60_and_minimum_threshold() {
    let r = chk(150.0, 0.0, nf(0.0)).unwrap();
    let qsw = (150.0 - 89.25) / H0;
    assert!(r.stirrups_required());
    assert!(close(r.qsw_min, 0.25 * RBT, 1e-9));
    assert!(close(r.qsw_required, qsw, 1e-6));
    assert!(close(r.asw_per_sw, qsw / 300_000.0, 1e-12));
    let r = chk(100.0, 0.0, nf(0.0)).unwrap();
    assert!(close(r.qsw_required, 0.25 * RBT, 1e-9));
    assert!(close(RebarClass::A400.rsw(), 280_000.0, 1e-9));
}

#[test]
fn strip_failure_reported() {
    let r = chk(800.0, 0.0, nf(0.0)).unwrap();
    assert!(close(r.eta_strip, 800.0 / 739.5, 1e-9));
}

#[test]
fn compression_raises_phi_n_and_applies_to_strip() {
    let r = chk(100.0, 0.0, nf(-1000.0)).unwrap();
    assert!(close(r.phi_n, 1.25, 1e-12));
    assert!(r.phi_n_to_strip);
    assert!(close(r.q_strip, 1.25 * 0.3 * RB * H0, 1e-6));
    assert!(close(r.qb1, 1.25 * 0.5 * RBT * H0, 1e-6));
}

#[test]
fn tension_lowers_phi_n_not_strip() {
    let r = chk(100.0, 0.0, nf(200.0)).unwrap();
    assert!(r.phi_n < 1.0 && !r.phi_n_to_strip);
    assert!(close(r.q_strip, 0.3 * RB * H0, 1e-6));
    assert!(close(r.qb1, r.phi_n * 0.5 * RBT * H0, 1e-6));
}
