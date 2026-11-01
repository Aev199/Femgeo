//! Прочность нормального сечения круглой сваи. Эталон — отчёты АРБАТ v23.1 (СП 63.13330.2018),
//! предоставленные пользователем и сохранённые в `cscore-pile-poc`
//! (`examples/verification/arbat_matched/arbat_reference_v2_nm8.json`). Исходные данные — как в
//! PoC (сверка CScore ↔ АРБАТ < 1 %): B30 с γb1 = 0,9, A500, a1 = 70 мм до грани стержня,
//! l0 = 3 м, случайный эксцентриситет при сжатии.

use femgeo::sp63::materials::{ConcreteChars, ConcreteDiagramType, RebarChars};
use femgeo::sp63::pile::{self, LimitOptions, Status};
use femgeo::sp63::section::Section;

fn pile_section(d: f64, bars: usize, bar_d: f64, rotation: f64) -> Section {
    let c = ConcreteChars::short_term(32_500_000.0, 15_300.0, 1_035.0);
    let r = RebarChars::physical_yield(2e8, 435_000.0, 435_000.0);
    let bar_radius = d / 2.0 - 0.07 - bar_d / 2.0;
    Section::circular(d / 2.0, 96, 20, bars, bar_d, bar_radius, rotation, c.diagram(ConcreteDiagramType::L2), r.diagram(), c.eps_b2)
}

/// (N, M, Кисп АРБАТ).
const P800: [(f64, f64, f64); 8] = [
    (0.0, 200.0, 0.28), (0.0, 600.0, 0.839), (-1000.0, 200.0, 0.229), (-1000.0, 600.0, 0.655),
    (-2000.0, 200.0, 0.317), (-2000.0, 600.0, 0.627), (500.0, 200.0, 0.446), (500.0, 600.0, 0.998),
];
const P1000: [(f64, f64, f64); 8] = [
    (0.0, 200.0, 0.131), (0.0, 600.0, 0.393), (-1000.0, 200.0, 0.122), (-1000.0, 600.0, 0.3),
    (-2000.0, 200.0, 0.182), (-2000.0, 600.0, 0.316), (500.0, 200.0, 0.234), (500.0, 600.0, 0.484),
];

#[test]
fn matches_arbat_v2_reports() {
    let o = LimitOptions::default();
    let mut worst: f64 = 0.0;
    for (d, bars, bd, cases) in [(0.8, 12, 0.025, P800), (1.0, 16, 0.028, P1000)] {
        let s = pile_section(d, bars, bd, 0.0);
        for (n, m, arbat) in cases {
            let (n, mx, my) = pile::with_accidental_eccentricity(n, m, 0.0, 3.0, d, None).unwrap();
            let r = pile::all_factor(&s, n, mx, my, &o);
            assert_eq!(r.status, Status::Ok, "D={d} N={n} M={m}");
            let dev = (r.utilization - arbat) / arbat;
            println!("D={d} N={n} M={m}: Кисп {:.4}, АРБАТ {arbat}, {:+.2} %", r.utilization, 100.0 * dev);
            worst = worst.max(dev.abs());
        }
    }
    assert!(worst < 0.02, "наибольшее отклонение от АРБАТ {:.2} %", 100.0 * worst);
}

/// Центральное сжатие: в OpenCS допустима εb = εcu = −0,0035 по всему сечению — бетон на Rb,
/// стержни на Rsc (разностная диаграмма вычитает бетон под стержнем).
#[test]
fn central_compression_squash_load() {
    let s = pile_section(0.8, 12, 0.025, 0.0);
    let ac: f64 = s.concrete.iter().map(|f| f.a).sum();
    let as_: f64 = s.bars.iter().map(|f| f.a).sum();
    let nult = 15_300.0 * ac + (435_000.0 - 15_300.0) * as_;
    let r = pile::all_factor(&s, -1000.0, 0.0, 0.0, &LimitOptions::default());
    assert!((r.factor * 1000.0 - nult).abs() < 2e-3 * nult, "{} vs {nult}", r.factor * 1000.0);
}

#[test]
fn central_tension_yields_all_bars() {
    let s = pile_section(0.8, 12, 0.025, 0.0);
    let nult = 435_000.0 * s.bars.iter().map(|f| f.a).sum::<f64>();
    let r = pile::all_factor(&s, 1000.0, 0.0, 0.0, &LimitOptions::default());
    assert!((r.factor * 1000.0 - nult).abs() < 2e-3 * nult);
}

/// Поворот сечения вместе с моментом не меняет результат: 90°, а также косой изгиб под 45° ≡
/// прямой изгиб сечения, повёрнутого на 45° (для 12 стержней — 15° по модулю шага 30°).
#[test]
fn rotation_invariance() {
    let o = LimitOptions::default();
    let a = pile::all_factor(&pile_section(0.8, 12, 0.025, 0.0), -1000.0, 500.0, 0.0, &o);
    let b = pile::all_factor(&pile_section(0.8, 12, 0.025, 90.0), -1000.0, 0.0, 500.0, &o);
    assert!((a.utilization - b.utilization).abs() < 1e-3 * a.utilization);
    let m = 500.0 / 2f64.sqrt();
    let c = pile::all_factor(&pile_section(0.8, 12, 0.025, 0.0), -1000.0, m, m, &o);
    let e = pile::all_factor(&pile_section(0.8, 12, 0.025, 15.0), -1000.0, 500.0, 0.0, &o);
    assert!((c.utilization - e.utilization).abs() < 1e-3 * c.utilization, "{} vs {}", c.utilization, e.utilization);
}

#[test]
fn moment_factor_keeps_axial_force() {
    let s = pile_section(0.8, 12, 0.025, 0.0);
    let r = pile::moment_factor(&s, -1000.0, 300.0, 0.0, &LimitOptions::default());
    assert_eq!(r.status, Status::Ok);
    assert_eq!(r.n_limit, -1000.0);
    assert!((r.mx_limit - 300.0 * r.factor).abs() < 1e-9);
}

#[test]
fn zero_load_unbounded() {
    let r = pile::all_factor(&pile_section(0.8, 12, 0.025, 0.0), 0.0, 0.0, 0.0, &LimitOptions::default());
    assert_eq!(r.status, Status::Unbounded);
}

#[test]
fn accidental_eccentricity() {
    // ea = max(3/600, 0,8/30, 0,01) = 26,67 мм.
    let (_, mx, my) = pile::with_accidental_eccentricity(-1500.0, 30.0, 40.0, 3.0, 0.8, None).unwrap();
    let dm = 1500.0 * 0.8 / 30.0;
    assert!((mx - (30.0 + dm * 0.6)).abs() < 1e-9 && (my - (40.0 + dm * 0.8)).abs() < 1e-9);
    assert_eq!(pile::with_accidental_eccentricity(500.0, 30.0, 0.0, 3.0, 0.8, None).unwrap(), (500.0, 30.0, 0.0));
    assert!(pile::with_accidental_eccentricity(-500.0, 0.0, 0.0, 3.0, 0.8, None).is_err());
}
