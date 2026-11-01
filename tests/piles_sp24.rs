//! Свайное поле по СП 24, п. 7.4: свойства формул и связь с МКЭ.

use femgeo::fem::{LoadCase, Model, Shell, ShellMaterial, SpringGroup};
use femgeo::piles::{interaction, Field, Layer, Pile, SinglePile, SoilColumn};

fn soil() -> SoilColumn {
    SoilColumn { top: 0.0, layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3 }, Layer { bottom: -40.0, e: 40_000.0, nu: 0.3 }] }
}

fn pile(x: f64, y: f64) -> Pile {
    let d: f64 = 0.6;
    Pile { x, y, head: -1.0, length: 15.0, d, ea: 3.0e7 * std::f64::consts::PI * d * d / 4.0 }
}

#[test]
fn soil_averaging() {
    // По длине сваи −16…−1: 7 м слоя 1 и 8 м слоя 2.
    let (g, nu) = soil().average(-16.0, -1.0);
    let g1 = 15_000.0 / 2.6;
    let g2 = 40_000.0 / 2.6;
    assert!((g - (7.0 * g1 + 8.0 * g2) / 15.0).abs() < 1e-9);
    assert!((nu - 0.3).abs() < 1e-12);
}

/// Формула п. 7.4.3 по шагам, вручную.
#[test]
fn single_pile_formula() {
    let p = pile(0.0, 0.0);
    let sp = SinglePile::new(&p, &soil());
    let (g1, g2): (f64, f64) = ((7.0 * 15_000.0 + 8.0 * 40_000.0) / 15.0 / 2.6, 40_000.0 / 2.6);
    let kv: f64 = 2.82 - 3.78 * 0.3 + 2.18 * 0.09;
    let br = 0.17 * (kv * g1 * 15.0 / (g2 * 0.6)).ln();
    let ar = 0.17 * (kv * 15.0 / 0.6).ln();
    let chi = p.ea / (g1 * 225.0);
    let l1 = 2.12 * chi.powf(0.75) / (1.0 + 2.12 * chi.powf(0.75));
    let beta = br / l1 + (1.0 - br / ar) / chi;
    assert!((sp.beta - beta).abs() < 1e-12);
    assert!((sp.compliance - beta / (g1 * 15.0)).abs() < 1e-15);
    assert!(sp.warnings.is_empty());
    // Порядок величины: осадка от 1000 кН — миллиметры.
    let s = 1000.0 * sp.compliance;
    assert!((0.002..0.02).contains(&s), "s = {s}");
}

#[test]
fn stiffer_soil_less_settlement_and_interaction_decays() {
    let p = pile(0.0, 0.0);
    let soft = SinglePile::new(&p, &soil());
    let mut hard = soil();
    for l in &mut hard.layers {
        l.e *= 2.0;
    }
    assert!(SinglePile::new(&p, &hard).compliance < soft.compliance);
    let (d1, d2) = (interaction(&soft, 15.0, 1.8), interaction(&soft, 15.0, 3.6));
    assert!(d1 > d2 && d2 > 0.0);
    assert_eq!(interaction(&soft, 15.0, 1e3), 0.0);
}

#[test]
fn field_matrix_symmetric_positive_and_invariant() {
    let ps: Vec<Pile> = (0..9).map(|i| pile(1.8 * (i % 3) as f64, 1.8 * (i / 3) as f64)).collect();
    let f = Field::new(&ps, &|_, _| soil());
    let k = f.stiffness().unwrap();
    let n = ps.len();
    for i in 0..n {
        for j in 0..n {
            assert!((k[i * n + j] - k[j * n + i]).abs() < 1e-6 * k[i * n + i]);
        }
    }
    // Центральная свая в группе — податливее одиночной (взаимовлияние): её доля нагрузки
    // при равной осадке меньше, чем у угловой.
    let row = |i: usize| (0..n).map(|j| k[i * n + j]).sum::<f64>();
    assert!(row(4) < row(0));
    // Перестановка и перенос свай не меняют результат.
    let mut ps2: Vec<Pile> = ps.iter().rev().copied().collect();
    for p in &mut ps2 {
        p.x += 100.0;
        p.y -= 50.0;
    }
    let k2 = Field::new(&ps2, &|_, _| soil()).stiffness().unwrap();
    let r = |kk: &[f64], i: usize| (0..n).map(|j| kk[i * n + j]).sum::<f64>();
    assert!((r(&k, 4) - r(&k2, 4)).abs() < 1e-9 * r(&k, 4));
    assert!((r(&k, 0) - r(&k2, 8)).abs() < 1e-9 * r(&k, 0));
}

/// Жёсткий ростверк на 4 сваях: осадка совпадает с решением по матрице податливости.
#[test]
fn rigid_raft_on_pile_group_in_fem() {
    let ps: Vec<Pile> = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].iter().map(|&(x, y)| pile(x, y)).collect();
    let field = Field::new(&ps, &|_, _| soil());
    let k = field.stiffness().unwrap();
    let mut m = Model::default();
    m.nodes = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]];
    m.shells.push(Shell { nodes: vec![0, 1, 2, 3], mat: ShellMaterial { e: 3e10, nu: 0.2, h: 1.5 } });
    m.spring_groups.push(SpringGroup { nodes: vec![0, 1, 2, 3], k });
    m.fixed = vec![(0, [true, true, false, false, false, true]), (1, [false, true, false, false, false, true]), (2, [false, false, false, false, false, true]), (3, [false, false, false, false, false, true])];
    let lc = LoadCase { nodal: (0..4).map(|i| (i, [0.0, 0.0, -1000.0, 0.0, 0.0, 0.0])).collect(), ..Default::default() };
    let u = m.factorize().unwrap().solve(&m, &lc).unwrap();
    // Симметрия: все сваи по 1000 кН, s = Σ_j F_0j·1000.
    let s: f64 = (0..4).map(|j| field.flexibility[j] * 1000.0).sum();
    for i in 0..4 {
        assert!((-u[6 * i + 2] - s).abs() < 1e-3 * s, "{} vs {s}", -u[6 * i + 2]);
    }
}
