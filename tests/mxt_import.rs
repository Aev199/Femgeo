//! Чтение MIDAS .mxt из topo_reconstruct_rs: «Тестовая плита» ЛИРА (10 × 10 м, h = 0,2 м,
//! 10 стержней 20 × 20 см, загружение 1 — 11,3188 т/м² на 5 пластинах 1 × 1 м).

use femgeo::fem::Spring;
use femgeo::mxt;

#[test]
fn test_plate() {
    let m = mxt::parse(include_str!("data/test_plate.mxt")).unwrap();
    assert_eq!((m.model.nodes.len(), m.model.shells.len(), m.model.beams.len()), (226, 195, 20));
    assert!(m.model.shells.iter().all(|s| s.mat.h == 0.2 && s.mat.nu == 0.17));
    assert_eq!(m.cases.len(), 2);
    let fz = |k: usize| m.cases[k].1.nodal.iter().map(|x| x.1[2]).sum::<f64>();
    // 5 м² × 11,3188 т/м² × 9,80665 кН/т.
    let q = 5.0 * 11.3188 * 9.80665;
    assert!((fz(0) + q).abs() < 0.01, "{}", fz(0));
    // Собственный вес: плита 100 м² × 0,2 м × γ + 10 стержней 1 м × 0,04 м² × γ (γ = 2,5 в файле).
    let sw = 2.5 * (100.0 * 0.2 + 10.0 * 0.04);
    assert!((fz(1) + sw).abs() < 1e-6, "{}", fz(1));
    // Равновесие на пружинах во всех узлах.
    let mut model = m.model.clone();
    let k = 1e6;
    for n in 0..model.nodes.len() {
        model.springs.push(Spring { node: n, k: [k; 6] });
    }
    let u = model.factorize().unwrap().solve(&model, &m.cases[0].1).unwrap();
    let r: f64 = (0..model.nodes.len()).map(|n| k * u[6 * n + 2]).sum();
    assert!((r - fz(0)).abs() < 1e-6 * q, "{r}");
}
