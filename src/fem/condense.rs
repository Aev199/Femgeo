//! Статическая конденсация модели на заданные степени свободы (суперэлемент):
//! K_c = K_bb − K_bi·K_ii⁻¹·K_ib, f_c = f_b − K_bi·K_ii⁻¹·f_i.
//!
//! Применение: надземная часть здания конденсируется на узлы сопряжения с плитой один раз до
//! подбора свай; в каждом варианте к плите добавляется [`super::DenseBlock`] и приведённые
//! нагрузки.

use std::collections::HashMap;

use faer::linalg::solvers::Solve;
use faer::sparse::{SparseColMat, Triplet};
use faer::{Mat, Side};

use super::{LoadCase, Model};

/// Результат конденсации: матрица по сохранённым степеням свободы (в порядке `dofs`,
/// построчно) и приведённые нагрузки для каждого загружения.
#[derive(Clone, Debug)]
pub struct Condensed {
    pub dofs: Vec<usize>,
    pub k: Vec<f64>,
    pub loads: Vec<Vec<f64>>,
}

/// Конденсация `m` на глобальные степени свободы `retained` (не закреплённые). Внутренние —
/// все остальные незакреплённые; при них схема должна быть геометрически неизменяемой.
pub fn condense(m: &Model, retained: &[usize], cases: &[LoadCase]) -> Result<Condensed, String> {
    let (eq_all, _) = m.equations();
    let nr = retained.len();
    let mut ret: HashMap<usize, usize> = HashMap::with_capacity(nr);
    for (r, &d) in retained.iter().enumerate() {
        if eq_all.get(d).copied().flatten().is_none() {
            return Err(format!("сохраняемая степень свободы {d} закреплена или вне модели"));
        }
        if ret.insert(d, r).is_some() {
            return Err(format!("степень свободы {d} указана дважды"));
        }
    }
    let mut int = vec![None; eq_all.len()];
    let mut ni = 0;
    for (d, e) in eq_all.iter().enumerate() {
        if e.is_some() && !ret.contains_key(&d) {
            int[d] = Some(ni);
            ni += 1;
        }
    }
    let mut trip = Vec::new();
    let mut kib: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nr];
    let mut kbb = vec![0.0; nr * nr];
    m.entries(&mut |a, b, v| match (int[a], ret.get(&b), int[b], ret.get(&a)) {
        (Some(ia), _, Some(ib), _) if ia >= ib => trip.push(Triplet::new(ia, ib, v)),
        (Some(ia), Some(&rb), None, _) => kib[rb].push((ia, v)),
        (None, _, None, Some(&ra)) => {
            if let Some(&rb) = ret.get(&b) {
                kbb[ra * nr + rb] += v;
            }
        }
        _ => {}
    })?;

    let mut kc = kbb;
    let mut loads = Vec::with_capacity(cases.len());
    if ni == 0 {
        for c in cases {
            let f = m.load_vector(c)?;
            loads.push(retained.iter().map(|&d| f[d]).collect());
        }
        return Ok(Condensed { dofs: retained.to_vec(), k: kc, loads });
    }
    let a = SparseColMat::<usize, f64>::try_new_from_triplets(ni, ni, &trip).map_err(|e| format!("сборка: {e:?}"))?;
    let llt = a
        .sp_cholesky(Side::Lower)
        .map_err(|_| "внутренняя часть изменяема: проверьте закрепления узлов сопряжения".to_string())?;

    const BLOCK: usize = 64;
    for start in (0..nr).step_by(BLOCK) {
        let cols = (start..(start + BLOCK).min(nr)).collect::<Vec<_>>();
        let mut rhs = Mat::<f64>::zeros(ni, cols.len());
        for (c, &j) in cols.iter().enumerate() {
            for &(i, v) in &kib[j] {
                rhs[(i, c)] -= v;
            }
        }
        let x = llt.solve(&rhs);
        for r in 0..nr {
            for (c, &j) in cols.iter().enumerate() {
                let s: f64 = kib[r].iter().map(|&(i, v)| v * x[(i, c)]).sum();
                kc[r * nr + j] += s;
            }
        }
    }
    for r in 0..nr {
        for j in r + 1..nr {
            let v = 0.5 * (kc[r * nr + j] + kc[j * nr + r]);
            kc[r * nr + j] = v;
            kc[j * nr + r] = v;
        }
    }
    for c in cases {
        let f = m.load_vector(c)?;
        let mut fi = Mat::<f64>::zeros(ni, 1);
        for (d, e) in int.iter().enumerate() {
            if let Some(i) = e {
                fi[(*i, 0)] = f[d];
            }
        }
        let u = llt.solve(&fi);
        loads.push((0..nr).map(|r| f[retained[r]] - kib[r].iter().map(|&(i, v)| v * u[(i, 0)]).sum::<f64>()).collect());
    }
    Ok(Condensed { dofs: retained.to_vec(), k: kc, loads })
}

/// Надземная часть, приведённая к плите.
#[derive(Clone, Debug)]
pub struct Building {
    /// Суперэлемент в степенях свободы плиты (uz, θx, θy узлов сопряжения).
    pub block: super::DenseBlock,
    /// Приведённые нагрузки для каждого загружения: (узел плиты, силы по 6 степеням свободы).
    pub loads: Vec<Vec<(usize, [f64; 6])>>,
}

/// Конденсация здания `m` на узлы сопряжения с плитой: `interface` — пары (узел здания, узел
/// плиты). В узлах сопряжения горизонтальные смещения и поворот вокруг вертикали удержаны
/// (плита в своей плоскости жёсткая) — в модели плиты те же степени свободы должны быть
/// закреплены или жёстки; сохраняются uz, θx, θy.
pub fn building(m: &Model, interface: &[(usize, usize)], cases: &[LoadCase]) -> Result<Building, String> {
    let mut s = m.clone();
    for &(n, _) in interface {
        s.fixed.push((n, [true, true, false, false, false, true]));
    }
    let retained: Vec<usize> = interface.iter().flat_map(|&(n, _)| [6 * n + 2, 6 * n + 3, 6 * n + 4]).collect();
    let c = condense(&s, &retained, cases)?;
    let dofs = interface.iter().flat_map(|&(_, r)| [6 * r + 2, 6 * r + 3, 6 * r + 4]).collect();
    let loads = c
        .loads
        .iter()
        .map(|f| interface.iter().enumerate().map(|(i, &(_, r))| (r, [0.0, 0.0, f[3 * i], f[3 * i + 1], f[3 * i + 2], 0.0])).collect())
        .collect();
    Ok(Building { block: super::DenseBlock { dofs, k: c.k }, loads })
}
