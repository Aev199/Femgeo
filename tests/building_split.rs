//! Здание целиком → плита + сконденсированная надземная часть.

use femgeo::building::{split, Mesh};
use femgeo::fem::{Beam, Shell, ShellMaterial};

fn grid(m: &mut Mesh, z: f64, h: f64, hole: bool) -> Vec<usize> {
    let base = m.nodes.len();
    for i in 0..=10 {
        for j in 0..=10 {
            m.nodes.push([j as f64, i as f64, z]);
        }
    }
    let mut ids = vec![];
    for i in 0..10 {
        for j in 0..10 {
            if hole && (4..6).contains(&i) && (4..6).contains(&j) {
                continue;
            }
            let n = |a: usize, b: usize| base + a * 11 + b;
            ids.push(m.shells.len());
            m.shells.push(Shell { nodes: vec![n(i, j), n(i, j + 1), n(i + 1, j + 1), n(i + 1, j)], mat: ShellMaterial { e: 3e7, nu: 0.2, h } });
        }
    }
    ids
}

#[test]
fn raft_with_columns_and_slab() {
    let mut m = Mesh::default();
    m.raft = grid(&mut m, 0.0, 1.0, true);
    grid(&mut m, 3.0, 0.2, false);
    for (x, y) in [(2, 2), (8, 2), (2, 8), (8, 8)] {
        let (a, b) = (y * 11 + x, 121 + y * 11 + x);
        m.beams.push(Beam { n1: a, n2: b, e: 3e7, g: 1.25e7, a: 0.16, iy: 2.1e-3, iz: 2.1e-3, j: 3.6e-3, ref_dir: None });
    }
    // Загружение: 1 кН вниз в каждом узле перекрытия и 2 кН в узле плиты.
    m.cases.push((121..242).map(|n| (n, [0.0, 0.0, -1.0, 0.0, 0.0, 0.0])).chain([(0, [0.0, 0.0, -2.0, 0.0, 0.0, 0.0])]).collect());
    m.piles.push((24, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0]));
    let s = split(&m).unwrap();
    assert_eq!(s.raft.shells.len(), 96);
    assert_eq!(s.region.holes.len(), 1);
    let area = |p: &[[f64; 2]]| (0..p.len()).map(|i| p[i][0] * p[(i + 1) % p.len()][1] - p[(i + 1) % p.len()][0] * p[i][1]).sum::<f64>().abs() / 2.0;
    assert!((area(&s.region.outer) - 100.0).abs() < 1e-9 && (area(&s.region.holes[0]) - 4.0).abs() < 1e-9);
    assert_eq!((s.z_mid, s.h), (0.0, 1.0));
    assert_eq!(s.piles.len(), 1);
    assert_eq!(s.raft_nodes[s.piles[0].0], 24);
    let b = s.building.unwrap();
    assert_eq!(b.block.dofs.len(), 4 * 3);
    // Перекрытие на колоннах: вся нагрузка (121 кН) приходит в 4 узла сопряжения.
    let fz: f64 = b.loads[0].iter().map(|l| l.1[2]).sum();
    assert!((fz + 121.0).abs() < 1e-6, "{fz}");
    assert_eq!(s.raft_cases[0].nodal.iter().map(|l| l.1[2]).sum::<f64>(), -2.0);
}
