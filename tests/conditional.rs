//! Условный фундамент: α по Буссинеску, расширение контура, осадка.

use femgeo::piles::conditional::{alpha, offset, settlement};
use femgeo::piles::{Layer, SoilColumn};

/// Угол прямоугольника B × L (Ньюмарк), σz/p на глубине z.
fn corner(b: f64, l: f64, z: f64) -> f64 {
    let (m, n) = (b / z, l / z);
    let r = (1.0 + m * m + n * n).sqrt();
    (m * n / r * (1.0 / (1.0 + m * m) + 1.0 / (1.0 + n * n)) + (m * n / r).atan()) / (2.0 * std::f64::consts::PI)
}

#[test]
fn alpha_matches_newmark() {
    let rect = [[0.0, 0.0], [12.0, 0.0], [12.0, 30.0], [0.0, 30.0]];
    for z in [0.5, 3.0, 10.0, 25.0] {
        let exact = 4.0 * corner(6.0, 15.0, z);
        let a = alpha(&rect, [6.0, 15.0], z);
        assert!((a - exact).abs() < 2e-3, "z = {z}: {a} {exact}");
    }
    // Точка вне контура (L-образный): сумма прямоугольников.
    let l = [[0.0, 0.0], [20.0, 0.0], [20.0, 10.0], [10.0, 10.0], [10.0, 20.0], [0.0, 20.0]];
    let z = 6.0;
    let exact = 4.0 * corner(5.0, 5.0, z) + 2.0 * (corner(15.0, 5.0, z) - corner(5.0, 5.0, z)) + 2.0 * (corner(5.0, 15.0, z) - corner(5.0, 5.0, z));
    assert!((alpha(&l, [5.0, 5.0], z) - exact).abs() < 2e-3);
}

#[test]
fn offset_square() {
    let o = offset(&[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]], 1.5);
    assert!(o.iter().any(|p| (p[0] + 1.5).abs() < 1e-12 && (p[1] + 1.5).abs() < 1e-12));
    assert!(o.iter().any(|p| (p[0] - 11.5).abs() < 1e-12 && (p[1] - 11.5).abs() < 1e-12));
}

/// Однородный слой: осадка по формуле послойного суммирования, независимо посчитанная здесь.
#[test]
fn settlement_uniform() {
    let col = SoilColumn { top: 0.0, layers: vec![Layer { bottom: -200.0, e: 20_000.0, nu: 0.3, gamma: 20.0, phi: 20.0 }] };
    let outline = [[0.0, 0.0], [20.0, 0.0], [20.0, 20.0], [0.0, 20.0]];
    let (head, l, n) = (-3.0, 15.0, 80_000.0);
    let c = settlement(&col, &outline, head, l, n).unwrap();
    let t = l * (20.0f64 / 4.0).to_radians().tan();
    let b = 20.0 + 2.0 * t;
    assert!((c.area - b * b).abs() < 1e-6);
    let p = (n + 20.0 * l * b * b) / (b * b);
    assert!((c.p - p).abs() < 1e-6);
    let s_tip = 20.0 * 18.0;
    let (mut s, mut z) = (0.0, 0.0);
    loop {
        let a = 4.0 * corner(b / 2.0, b / 2.0, z + 0.1);
        if a * p <= 0.5 * 20.0 * (18.0 + z + 0.1) {
            break;
        }
        s += 0.8 * ((a * p - a * s_tip) * 0.2 / 20_000.0 + a * s_tip * 0.2 / 100_000.0);
        z += 0.2;
    }
    assert!((c.s / s - 1.0).abs() < 5e-3, "{} {s}", c.s);
    assert!((c.hc - z).abs() < 0.41);
}
