//! Здание целиком (сетка препроцессора topo_reconstruct_rs) → плита Femgeo и надземная часть,
//! сконденсированная на плиту.
//!
//! Плита — заданные пластины, горизонтальные (разброс Z срединной поверхности ≤ 1 мм).
//! Здание — остальные пластины и стержни; узлы сопряжения — общие с плитой. Нагрузки узлов
//! плиты (в том числе сопряжения) остаются плите, остальные приводятся конденсацией
//! (`fem::condense::building`: в узлах сопряжения удержаны горизонтальные смещения и поворот
//! вокруг вертикали). Существующие сваи — одноузловые КЭ (ЛИРА КЭ 56) в узлах плиты.

use crate::fem::condense::{self, Building};
use crate::fem::{Beam, LoadCase, Model, Shell, V3};
use crate::piles::layout::Region;
use std::collections::HashMap;

/// Сетка здания: силы — кН, длины — м, глобальные оси, Z вверх.
#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub nodes: Vec<V3>,
    pub shells: Vec<Shell>,
    pub beams: Vec<Beam>,
    /// Пластины фундаментной плиты (индексы `shells`).
    pub raft: Vec<usize>,
    /// Узловые нагрузки по загружениям: (узел, силы и моменты по 6 степеням свободы).
    pub cases: Vec<Vec<(usize, [f64; 6])>>,
    /// Одноузловые упругие связи (сваи, КЭ 56): точка и жёсткости по 6 степеням свободы.
    pub piles: Vec<([f64; 3], [f64; 6])>,
}

pub struct Split {
    /// Плита: свои узлы (с отметкой срединной поверхности) и пластины.
    pub raft: Model,
    /// Узел сетки для каждого узла плиты.
    pub raft_nodes: Vec<usize>,
    /// Нагрузки узлов плиты по загружениям.
    pub raft_cases: Vec<LoadCase>,
    /// Здание, приведённое к плите (None — надземной части нет).
    pub building: Option<Building>,
    pub region: Region,
    /// Отметка срединной поверхности и толщина плиты (наибольшая из пластин).
    pub z_mid: f64,
    pub h: f64,
    /// Существующие сваи под плитой (в контуре, по высоте — в пределах плиты ± 1 м): (X, Y), жёсткости.
    pub piles: Vec<([f64; 2], [f64; 6])>,
    /// Что исключено или приближено.
    pub notes: Vec<String>,
}

/// Поворотные пружины-стабилизаторы узлов здания, кН·м/рад: гасят местные механизмы
/// (вращение стержня вокруг оси, поворот из плоскости у пластин без «сверлящей» жёсткости).
/// Перемещения не стабилизируются — силы передаются на плиту полностью; поступательно
/// свободные части исключаются обходом связности.
pub const STABILIZER: f64 = 1.0;

fn signed_area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    (0..n).map(|i| p[i][0] * p[(i + 1) % n][1] - p[(i + 1) % n][0] * p[i][1]).sum::<f64>() / 2.0
}

/// Контуры границы пластин: рёбра, принадлежащие одной пластине, сцепленные в циклы.
fn boundary(nodes: &[V3], shells: &[Shell]) -> Result<Region, String> {
    let mut count: HashMap<(usize, usize), usize> = HashMap::new();
    for s in shells {
        for k in 0..s.nodes.len() {
            let (a, b) = (s.nodes[k], s.nodes[(k + 1) % s.nodes.len()]);
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let mut next: HashMap<usize, Vec<usize>> = HashMap::new();
    for (&(a, b), &c) in &count {
        if c == 1 {
            next.entry(a).or_default().push(b);
            next.entry(b).or_default().push(a);
        }
    }
    if next.values().any(|v| v.len() != 2) {
        return Err("плита: граница не распадается на простые контуры (касание в узле или неконформная сетка)".into());
    }
    let mut loops: Vec<Vec<[f64; 2]>> = vec![];
    let mut keys: Vec<usize> = next.keys().copied().collect();
    keys.sort_unstable();
    let mut seen = std::collections::HashSet::new();
    for &start in &keys {
        if !seen.insert(start) {
            continue;
        }
        let (mut prev, mut cur) = (start, next[&start][0]);
        let mut lp = vec![[nodes[start][0], nodes[start][1]]];
        while cur != start {
            seen.insert(cur);
            lp.push([nodes[cur][0], nodes[cur][1]]);
            let n = &next[&cur];
            let to = if n[0] != prev { n[0] } else { n[1] };
            (prev, cur) = (cur, to);
        }
        loops.push(lp);
    }
    let outer = (0..loops.len()).max_by(|&a, &b| signed_area(&loops[a]).abs().total_cmp(&signed_area(&loops[b]).abs())).ok_or("плита: нет пластин")?;
    let outer_loop = loops.remove(outer);
    Ok(Region { outer: outer_loop, holes: loops })
}

pub fn split(m: &Mesh) -> Result<Split, String> {
    if m.raft.is_empty() {
        return Err("плита: не выбраны пластины".into());
    }
    let mut is_raft = vec![false; m.shells.len()];
    for &s in &m.raft {
        *is_raft.get_mut(s).ok_or(format!("плита: нет пластины {s}"))? = true;
    }
    // Узлы плиты.
    let mut raft_index: HashMap<usize, usize> = HashMap::new();
    let mut raft_nodes = vec![];
    let mut raft = Model::default();
    let mut h: f64 = 0.0;
    for &s in &m.raft {
        let sh = &m.shells[s];
        h = h.max(sh.mat.h);
        let nodes = sh
            .nodes
            .iter()
            .map(|&n| {
                *raft_index.entry(n).or_insert_with(|| {
                    raft_nodes.push(n);
                    raft.nodes.push(m.nodes[n]);
                    raft_nodes.len() - 1
                })
            })
            .collect();
        raft.shells.push(Shell { nodes, mat: sh.mat });
    }
    let zs: Vec<f64> = raft.nodes.iter().map(|p| p[2]).collect();
    let (z0, z1) = zs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &z| (a.0.min(z), a.1.max(z)));
    if z1 - z0 > 1e-3 {
        return Err(format!("плита не горизонтальна: Z от {z0:.3} до {z1:.3} м"));
    }
    let z_mid = zs.iter().sum::<f64>() / zs.len() as f64;
    let region = boundary(&raft.nodes, &raft.shells)?;

    // Здание: элементы, связанные с плитой (обход от узлов плиты); остальные — исключить.
    let mut notes = vec![];
    let (ns, nb) = (m.shells.len(), m.beams.len());
    let elem_nodes = |e: usize| -> Vec<usize> { if e < ns { m.shells[e].nodes.clone() } else { vec![m.beams[e - ns].n1, m.beams[e - ns].n2] } };
    let mut by_node: HashMap<usize, Vec<usize>> = HashMap::new();
    for e in 0..ns + nb {
        if e < ns && is_raft[e] {
            continue;
        }
        for n in elem_nodes(e) {
            by_node.entry(n).or_default().push(e);
        }
    }
    let mut keep = vec![false; ns + nb];
    let mut stack: Vec<usize> = raft_index.keys().copied().collect();
    let mut seen: std::collections::HashSet<usize> = stack.iter().copied().collect();
    while let Some(n) = stack.pop() {
        for &e in by_node.get(&n).map_or(&[][..], |v| &v[..]) {
            if !keep[e] {
                keep[e] = true;
                for q in elem_nodes(e) {
                    if seen.insert(q) {
                        stack.push(q);
                    }
                }
            }
        }
    }
    let dropped = (0..ns + nb).filter(|&e| !(e < ns && is_raft[e]) && !keep[e]).count();
    if dropped > 0 {
        notes.push(format!("элементы, не связанные с плитой, исключены: {dropped}"));
    }

    // Здание: остальные элементы на своих узлах.
    let mut b_index: HashMap<usize, usize> = HashMap::new();
    let mut bm = Model::default();
    let mut node_of = |n: usize, bm: &mut Model| -> usize {
        *b_index.entry(n).or_insert_with(|| {
            bm.nodes.push(m.nodes[n]);
            bm.nodes.len() - 1
        })
    };
    for (s, sh) in m.shells.iter().enumerate() {
        if !is_raft[s] && keep[s] {
            let nodes = sh.nodes.iter().map(|&n| node_of(n, &mut bm)).collect();
            bm.shells.push(Shell { nodes, mat: sh.mat });
        }
    }
    for (i, b) in m.beams.iter().enumerate() {
        if !keep[ns + i] {
            continue;
        }
        let (n1, n2) = (node_of(b.n1, &mut bm), node_of(b.n2, &mut bm));
        bm.beams.push(Beam { n1, n2, ..b.clone() });
    }
    let mut interface: Vec<(usize, usize)> = b_index.iter().filter_map(|(n, &bn)| raft_index.get(n).map(|&r| (bn, r))).collect();
    interface.sort_unstable();
    if !b_index.is_empty() && interface.is_empty() {
        return Err("здание не связано с плитой (нет общих узлов)".into());
    }

    let mut raft_cases = vec![];
    let mut b_cases = vec![];
    for (c, loads) in m.cases.iter().enumerate() {
        let (mut rl, mut bl) = (LoadCase::default(), LoadCase::default());
        let mut lost = 0.0;
        for &(n, f) in loads {
            if let Some(&r) = raft_index.get(&n) {
                rl.nodal.push((r, f));
            } else if let Some(&b) = b_index.get(&n) {
                bl.nodal.push((b, f));
            } else {
                lost += f[2];
            }
        }
        if lost.abs() > 1e-6 {
            notes.push(format!("загружение {}: нагрузка вне здания, связанного с плитой, не учтена: Fz = {lost:.1} кН", c + 1));
        }
        raft_cases.push(rl);
        b_cases.push(bl);
    }
    for n in 0..bm.nodes.len() {
        bm.springs.push(crate::fem::Spring { node: n, k: [0.0, 0.0, 0.0, STABILIZER, STABILIZER, STABILIZER] });
    }
    let building = if bm.nodes.is_empty() { None } else { Some(condense::building(&bm, &interface, &b_cases)?) };

    let piles: Vec<([f64; 2], [f64; 6])> =
        m.piles.iter().filter(|(p, _)| (p[2] - z_mid).abs() <= h / 2.0 + 1.0 && region.contains([p[0], p[1]])).map(|&(p, k)| ([p[0], p[1]], k)).collect();
    if piles.len() < m.piles.len() {
        notes.push(format!("сваи КЭ 56 вне плиты: {} из {}", m.piles.len() - piles.len(), m.piles.len()));
    }
    Ok(Split { raft, raft_nodes, raft_cases, building, region, z_mid, h, piles, notes })
}
