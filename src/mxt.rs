//! Модель здания из MIDAS Civil `.mxt` в виде, который пишет topo_reconstruct_rs
//! (`midas.rs`, `midas_stiffness.rs`): кН, м. Читаются *NODE, *ELEMENT (BEAM, PLATE 3/4),
//! *MATERIAL (USER: E кН/м², ν, γ кН/м³), *SECTION (DBUSER: SB, SR, P), *THICKNESS (VALUE),
//! *STLDCASE и блоки *USE-STLD с *CONLOAD, *BEAMLOAD (UNILOAD, глобальные GX/GY/GZ, доли
//! длины), *PRESSURE (PLATE FACE, глобальные GX/GY/GZ на весь элемент).
//!
//! Распределённые нагрузки приводятся к узловым силам статически эквивалентно, без узловых
//! моментов (допущение: на сетке здания погрешность локальна). Собственный вес topo в `.mxt`
//! не передаёт — загружение «Собственный вес» строится здесь по γ материалов (узлам поровну
//! для пластин, по концам для стержней). Закрепления в `.mxt` не пишутся.

use crate::fem::{Beam, LoadCase, Model, Shell, ShellMaterial, V3};
use std::collections::HashMap;

pub struct Mxt {
    pub model: Model,
    /// Загружения файла по порядку *STLDCASE; последним — собственный вес.
    pub cases: Vec<(String, LoadCase)>,
    /// Номер узла MIDAS для каждого узла модели.
    pub node_ids: Vec<usize>,
    pub warnings: Vec<String>,
}

pub const SELF_WEIGHT: &str = "Собственный вес";

struct Mat {
    e: f64,
    nu: f64,
    gamma: f64,
}

struct Sect {
    a: f64,
    iy: f64,
    iz: f64,
    j: f64,
}

fn fields(line: &str) -> Vec<&str> {
    line.split(';').next().unwrap_or("").split(',').map(str::trim).collect()
}

fn num(f: &[&str], i: usize, what: &str, ln: usize) -> Result<f64, String> {
    f.get(i).and_then(|s| s.parse::<f64>().ok()).ok_or_else(|| format!("mxt, строка {ln}: {what}: не число в поле {}", i + 1))
}

fn int(f: &[&str], i: usize, what: &str, ln: usize) -> Result<usize, String> {
    num(f, i, what, ln).map(|x| x as usize)
}

/// Прямоугольник b × h (h — по локальной z): A, Iy (изгиб в xz), Iz, Iк.
fn rect(b: f64, h: f64) -> Sect {
    let (s, l) = (b.min(h), b.max(h));
    let j = s.powi(3) * l * (1.0 / 3.0 - 0.21 * s / l * (1.0 - (s / l).powi(4) / 12.0));
    Sect { a: b * h, iy: b * h.powi(3) / 12.0, iz: h * b.powi(3) / 12.0, j }
}

fn pipe(d: f64, t: f64) -> Sect {
    let di = (d - 2.0 * t).max(0.0);
    let i = std::f64::consts::PI * (d.powi(4) - di.powi(4)) / 64.0;
    Sect { a: std::f64::consts::PI * (d * d - di * di) / 4.0, iy: i, iz: i, j: 2.0 * i }
}

fn area(p: &[V3]) -> f64 {
    let c = |a: V3, b: V3| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let d = |a: V3, b: V3| [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let mut s = [0.0; 3];
    for k in 1..p.len() - 1 {
        let x = c(d(p[0], p[k]), d(p[0], p[k + 1]));
        (0..3).for_each(|i| s[i] += x[i]);
    }
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() / 2.0
}

fn add(lc: &mut HashMap<usize, [f64; 6]>, n: usize, f: [f64; 3]) {
    let e = lc.entry(n).or_insert([0.0; 6]);
    (0..3).for_each(|i| e[i] += f[i]);
}

fn dir(s: &str, ln: usize) -> Result<usize, String> {
    match s {
        "GX" => Ok(0),
        "GY" => Ok(1),
        "GZ" => Ok(2),
        _ => Err(format!("mxt, строка {ln}: направление «{s}» не поддерживается (только GX/GY/GZ)")),
    }
}

pub fn parse(text: &str) -> Result<Mxt, String> {
    let mut nodes: Vec<(usize, V3)> = vec![];
    // (id, тип, материал, сечение/толщина, узлы MIDAS)
    let mut elems: Vec<(usize, bool, usize, usize, Vec<usize>)> = vec![];
    let mut mats: HashMap<usize, Mat> = HashMap::new();
    let mut sects: HashMap<usize, Sect> = HashMap::new();
    let mut thick: HashMap<usize, f64> = HashMap::new();
    let mut case_names: Vec<String> = vec![];
    // (загружение, строка, блок, поля)
    let mut loads: Vec<(usize, usize, String, Vec<String>)> = vec![];
    let mut warnings = vec![];
    let (mut block, mut case) = (String::new(), None::<usize>);
    for (k, raw) in text.lines().enumerate() {
        let ln = k + 1;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('*') {
            let head = rest.split(';').next().unwrap().trim();
            if let Some(name) = head.strip_prefix("USE-STLD,") {
                let name = name.trim();
                case = Some(case_names.iter().position(|c| c == name).ok_or(format!("mxt, строка {ln}: загружение «{name}» не объявлено в *STLDCASE"))?);
            }
            block = head.split(',').next().unwrap().trim().to_string();
            continue;
        }
        let f = fields(line);
        match block.as_str() {
            "NODE" => nodes.push((int(&f, 0, "узел", ln)?, [num(&f, 1, "X", ln)?, num(&f, 2, "Y", ln)?, num(&f, 3, "Z", ln)?])),
            "ELEMENT" => {
                let id = int(&f, 0, "элемент", ln)?;
                let (mat, prop) = (int(&f, 2, "материал", ln)?, int(&f, 3, "сечение", ln)?);
                match f.get(1).copied() {
                    Some("BEAM") => elems.push((id, false, mat, prop, vec![int(&f, 4, "узел", ln)?, int(&f, 5, "узел", ln)?])),
                    Some("PLATE") => {
                        let mut n: Vec<usize> = (4..8).map(|i| int(&f, i, "узел", ln)).collect::<Result<_, _>>()?;
                        if n[3] == 0 {
                            n.pop();
                        }
                        elems.push((id, true, mat, prop, n));
                    }
                    Some(t) => return Err(format!("mxt, строка {ln}: тип элемента «{t}» не поддерживается")),
                    None => return Err(format!("mxt, строка {ln}: нет типа элемента")),
                }
            }
            "MATERIAL" => {
                if f.get(1) != Some(&"USER") {
                    return Err(format!("mxt, строка {ln}: материал не USER"));
                }
                mats.insert(int(&f, 0, "материал", ln)?, Mat { e: num(&f, 10, "E", ln)?, nu: num(&f, 11, "ν", ln)?, gamma: num(&f, 13, "γ", ln)? });
            }
            "SECTION" => {
                let id = int(&f, 0, "сечение", ln)?;
                if f.get(1) != Some(&"DBUSER") {
                    return Err(format!("mxt, строка {ln}: сечение «{}» задано не размерами (DBUSER) — не поддерживается", f.get(2).unwrap_or(&"")));
                }
                let s = match f.get(11).copied() {
                    // SB, 2, h, b
                    Some("SB") => rect(num(&f, 14, "b", ln)?, num(&f, 13, "h", ln)?),
                    Some("SR") => pipe(num(&f, 13, "D", ln)?, num(&f, 13, "D", ln)? / 2.0),
                    Some("P") => pipe(num(&f, 13, "D", ln)?, num(&f, 14, "t", ln)?),
                    Some(t) => return Err(format!("mxt, строка {ln}: форма сечения «{t}» не поддерживается")),
                    None => return Err(format!("mxt, строка {ln}: нет формы сечения")),
                };
                sects.insert(id, s);
            }
            "THICKNESS" => {
                thick.insert(int(&f, 0, "толщина", ln)?, num(&f, 3, "толщина", ln)?);
            }
            "STLDCASE" => case_names.push(f[0].to_string()),
            "CONLOAD" | "BEAMLOAD" | "PRESSURE" => {
                let c = case.ok_or(format!("mxt, строка {ln}: нагрузка вне *USE-STLD"))?;
                loads.push((c, ln, block.clone(), f.iter().map(|s| s.to_string()).collect()));
            }
            "VERSION" | "UNIT" | "PROJINFO" | "USE-STLD" => {
                if block == "UNIT" && !(f.first() == Some(&"KN") && f.get(1) == Some(&"M")) {
                    return Err(format!("mxt, строка {ln}: единицы «{line}», ожидаются KN, M"));
                }
            }
            b => {
                if !warnings.iter().any(|w: &String| w.contains(&format!("*{b}"))) {
                    warnings.push(format!("блок *{b} не читается"));
                }
            }
        }
    }
    if nodes.is_empty() {
        return Err("mxt: нет узлов".into());
    }
    let index: HashMap<usize, usize> = nodes.iter().enumerate().map(|(i, n)| (n.0, i)).collect();
    let node = |id: usize| index.get(&id).copied().ok_or(format!("mxt: нет узла {id}"));
    let mut model = Model { nodes: nodes.iter().map(|n| n.1).collect(), ..Default::default() };
    let mut sw: HashMap<usize, [f64; 6]> = HashMap::new();
    // Элемент MIDAS → (пластина?, номер в модели).
    let mut elem_index: HashMap<usize, (bool, usize)> = HashMap::new();
    for (id, plate, mat, prop, ns) in &elems {
        let m = mats.get(mat).ok_or(format!("mxt: элемент {id}: нет материала {mat}"))?;
        let ns: Vec<usize> = ns.iter().map(|&n| node(n)).collect::<Result<_, _>>()?;
        if *plate {
            let h = *thick.get(prop).ok_or(format!("mxt: элемент {id}: нет толщины {prop}"))?;
            let w = m.gamma * h * area(&ns.iter().map(|&n| model.nodes[n]).collect::<Vec<_>>()) / ns.len() as f64;
            ns.iter().for_each(|&n| add(&mut sw, n, [0.0, 0.0, -w]));
            elem_index.insert(*id, (true, model.shells.len()));
            model.shells.push(Shell { nodes: ns, mat: ShellMaterial { e: m.e, nu: m.nu, h } });
        } else {
            let s = sects.get(prop).ok_or(format!("mxt: элемент {id}: нет сечения {prop}"))?;
            let (a, b) = (model.nodes[ns[0]], model.nodes[ns[1]]);
            let l = ((0..3).map(|i| (b[i] - a[i]).powi(2)).sum::<f64>()).sqrt();
            ns.iter().for_each(|&n| add(&mut sw, n, [0.0, 0.0, -m.gamma * s.a * l / 2.0]));
            elem_index.insert(*id, (false, model.beams.len()));
            model.beams.push(Beam { n1: ns[0], n2: ns[1], e: m.e, g: m.e / (2.0 * (1.0 + m.nu)), a: s.a, iy: s.iy, iz: s.iz, j: s.j, ref_dir: None });
        }
    }
    let mut acc: Vec<HashMap<usize, [f64; 6]>> = vec![HashMap::new(); case_names.len()];
    for (c, ln, block, f) in &loads {
        let f: Vec<&str> = f.iter().map(String::as_str).collect();
        let lc = &mut acc[*c];
        match block.as_str() {
            "CONLOAD" => {
                let n = node(int(&f, 0, "узел", *ln)?)?;
                let e = lc.entry(n).or_insert([0.0; 6]);
                for i in 0..6 {
                    e[i] += num(&f, i + 1, "сила", *ln)?;
                }
            }
            "BEAMLOAD" => {
                if f.get(2) != Some(&"UNILOAD") {
                    return Err(format!("mxt, строка {ln}: нагрузка на стержень «{}» не поддерживается", f.get(2).unwrap_or(&"")));
                }
                let id = int(&f, 0, "элемент", *ln)?;
                let Some(&(false, b)) = elem_index.get(&id) else { return Err(format!("mxt, строка {ln}: нет стержня {id}")) };
                let k = dir(f[3], *ln)?;
                let (t1, q1, t2, q2) = (num(&f, 10, "x1", *ln)?, num(&f, 11, "q1", *ln)?, num(&f, 12, "x2", *ln)?, num(&f, 13, "q2", *ln)?);
                let bm = &model.beams[b];
                let (a, e) = (model.nodes[bm.n1], model.nodes[bm.n2]);
                let l = ((0..3).map(|i| (e[i] - a[i]).powi(2)).sum::<f64>()).sqrt();
                // Трапеция на [t1, t2] (доли длины): равнодействующая и её положение.
                let r = (q1 + q2) / 2.0 * (t2 - t1) * l;
                let tc = if (q1 + q2).abs() > 1e-12 { t1 + (t2 - t1) * (q1 + 2.0 * q2) / (3.0 * (q1 + q2)) } else { (t1 + t2) / 2.0 };
                let mut v = [0.0; 3];
                v[k] = r * (1.0 - tc);
                add(lc, bm.n1, v);
                v[k] = r * tc;
                add(lc, bm.n2, v);
            }
            _ => {
                if f.get(1) != Some(&"PRES") || f.get(2) != Some(&"PLATE") || f.get(3) != Some(&"FACE") {
                    return Err(format!("mxt, строка {ln}: давление этого вида не поддерживается"));
                }
                let id = int(&f, 0, "элемент", *ln)?;
                let Some(&(true, s)) = elem_index.get(&id) else { return Err(format!("mxt, строка {ln}: нет пластины {id}")) };
                let k = dir(f[4], *ln)?;
                let p = num(&f, 9, "давление", *ln)?;
                let ns = &model.shells[s].nodes;
                let w = p * area(&ns.iter().map(|&n| model.nodes[n]).collect::<Vec<_>>()) / ns.len() as f64;
                let mut v = [0.0; 3];
                v[k] = w;
                for &n in ns.clone().iter() {
                    add(lc, n, v);
                }
            }
        }
    }
    let to_case = |m: HashMap<usize, [f64; 6]>| {
        let mut nodal: Vec<(usize, [f64; 6])> = m.into_iter().collect();
        nodal.sort_by_key(|x| x.0);
        LoadCase { nodal, pressure: vec![] }
    };
    let mut cases: Vec<(String, LoadCase)> = case_names.into_iter().zip(acc).map(|(n, m)| (n, to_case(m))).collect();
    cases.push((SELF_WEIGHT.to_string(), to_case(sw)));
    Ok(Mxt { model, cases, node_ids: nodes.iter().map(|n| n.0).collect(), warnings })
}
