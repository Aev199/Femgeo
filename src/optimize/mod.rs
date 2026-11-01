//! Подбор расстановки свай под сплошной плитой (`docs/SPEC.md`, §3): фронт вариантов
//! «число свай — расход армирования плиты, кг/м³» при ограничениях по несущей способности
//! свай по грунту и материалу, контуру, правилам СП 24 и осадкам.
//!
//! Первая версия: сваи ставятся в узлы сетки плиты; каждый вариант считается заново
//! (свайное поле по СП 24 → плита МКЭ → нагрузки на сваи, осадки → подбор Капра–Мори).
//! Алгоритм детерминирован и не зависит от нумерации кандидатов.

use crate::fem::condense::Building;
use crate::fem::{self, LoadCase, Model, SpringGroup, V3};
use crate::piles::layout::{self, Region, Rules};
use crate::piles::{Field, Pile, SoilColumn};
use crate::sp63::capra_maury::{self, PlateForces, UlsContext};

/// Тип сваи (одинаковый для всех свай варианта).
#[derive(Clone, Copy, Debug)]
pub struct PileType {
    pub length: f64,
    pub d: f64,
    pub ea: f64,
    /// Несущая способность по грунту с учётом γk (задаёт пользователь), кН.
    pub fd_ground: f64,
    /// Кисп по материалу при N = 1 кН (со случайным эксцентриситетом); 0 — не проверять.
    pub material_util_per_kn: f64,
}

/// Предельные осадки.
#[derive(Clone, Copy, Debug)]
pub struct SettlementLimits {
    /// Наибольшая осадка, м.
    pub s_max: f64,
    /// Наибольшая относительная разность осадок свай |Δs|/L.
    pub ds_rel_max: f64,
}

/// Постановка.
pub struct Problem<'a> {
    /// Плита: оболочки, узлы, закрепления в плоскости. Сваи добавляются в узлы.
    pub raft: Model,
    /// Отметка голов свай (низ плиты), м.
    pub pile_head: f64,
    pub pile: PileType,
    pub region: Region,
    pub rules: Rules,
    pub soil: &'a (dyn Fn(f64, f64) -> SoilColumn + Sync),
    /// Расчётные сочетания (ПС1) — нагрузки на плиту.
    pub combinations: Vec<LoadCase>,
    pub limits: SettlementLimits,
    /// Подбор армирования плиты (ПС1) и ось x армирования.
    pub design: &'a UlsContext,
    pub rebar_axis: V3,
    /// Толщина плиты, м (для кг/м³).
    pub raft_h: f64,
    /// Надземная часть, сконденсированная на плиту (`fem::condense::building`); приведённые
    /// нагрузки — по одному набору на сочетание. None — здание не учитывается.
    pub building: Option<&'a Building>,
}

/// Оценка варианта.
#[derive(Clone, Debug)]
pub struct Evaluation {
    /// Узлы плиты со сваями.
    pub piles: Vec<usize>,
    /// Наибольшее сжатие сваи по сочетаниям, кН (для каждой сваи).
    pub pile_n: Vec<f64>,
    pub util_ground: f64,
    pub util_material: f64,
    pub s_max: f64,
    pub ds_rel: f64,
    /// Теоретический расход арматуры плиты, кг/м³ (огибающая по сочетаниям).
    pub kg_m3: f64,
    /// Огибающая требуемой арматуры по элементам плиты, м²/м: [верх x, верх y, низ x, низ y].
    pub rebar: Vec<[f64; 4]>,
    /// Осадка узлов плиты в сочетании с наибольшей осадкой, м (вниз — положительная).
    pub w_nodes: Vec<f64>,
    /// Подбор невозможен хотя бы в одном элементе.
    pub design_impossible: bool,
    pub feasible: bool,
}

impl Evaluation {
    /// Наибольший Кисп из ограничений (для сравнения недопустимых вариантов).
    pub fn worst(&self, lim: &SettlementLimits) -> f64 {
        self.util_ground.max(self.util_material).max(self.s_max / lim.s_max).max(self.ds_rel / lim.ds_rel_max)
    }
}

const STEEL_DENSITY: f64 = 7850.0;

impl Problem<'_> {
    /// Узлы плиты, допустимые для свай.
    pub fn candidates(&self) -> Vec<usize> {
        let mut c: Vec<usize> = (0..self.raft.nodes.len())
            .filter(|&i| self.raft.shells.iter().any(|s| s.nodes.contains(&i)))
            .filter(|&i| self.region.admits([self.raft.nodes[i][0], self.raft.nodes[i][1]], &self.rules))
            .collect();
        // Порядок — по координатам, а не по номерам узлов.
        c.sort_by(|&a, &b| {
            let (p, q) = (self.raft.nodes[a], self.raft.nodes[b]);
            p[1].total_cmp(&q[1]).then(p[0].total_cmp(&q[0]))
        });
        c
    }

    fn xy(&self, n: usize) -> [f64; 2] {
        [self.raft.nodes[n][0], self.raft.nodes[n][1]]
    }

    /// Расчёт варианта со сваями в узлах `piles`.
    pub fn evaluate(&self, piles: &[usize]) -> Result<Evaluation, String> {
        let ps: Vec<Pile> = piles
            .iter()
            .map(|&n| {
                let [x, y] = self.xy(n);
                Pile { x, y, head: self.pile_head, length: self.pile.length, d: self.pile.d, ea: self.pile.ea }
            })
            .collect();
        let field = Field::new(&ps, self.soil);
        let k = field.stiffness()?;
        let mut m = self.raft.clone();
        m.spring_groups.push(SpringGroup { nodes: piles.to_vec(), k: k.clone() });
        let mut combinations = self.combinations.clone();
        if let Some(b) = self.building {
            if b.loads.len() != combinations.len() {
                return Err("здание: число наборов приведённых нагрузок не равно числу сочетаний".into());
            }
            m.blocks.push(b.block.clone());
            for &d in b.block.dofs.iter().step_by(3) {
                m.fixed.push((d / 6, [true, true, false, false, false, true]));
            }
            for (lc, extra) in combinations.iter_mut().zip(&b.loads) {
                lc.nodal.extend(extra.iter().copied());
            }
        }
        let lin = m.factorize()?;
        let np = piles.len();
        let mut pile_n = vec![0.0f64; np];
        let mut s_max: f64 = 0.0;
        let mut ds_rel: f64 = 0.0;
        let ne = m.shells.len();
        let mut env = vec![[0.0f64; 4]; ne];
        let mut w_nodes = vec![0.0f64; m.nodes.len()];
        let mut impossible = false;
        for lc in &combinations {
            let u = lin.solve(&m, lc)?;
            let w: Vec<f64> = piles.iter().map(|&n| u[6 * n + 2]).collect();
            let s_lc = w.iter().map(|x| -x).fold(0.0, f64::max);
            if s_lc >= s_max {
                w_nodes = (0..m.nodes.len()).map(|n| -u[6 * n + 2]).collect();
            }
            for i in 0..np {
                // Реакция сваи R = −K·w (сжатие — положительное).
                let r: f64 = -(0..np).map(|j| k[i * np + j] * w[j]).sum::<f64>();
                pile_n[i] = pile_n[i].max(r);
                s_max = s_max.max(-w[i]);
                for j in i + 1..np {
                    let [a, b] = [self.xy(piles[i]), self.xy(piles[j])];
                    let l = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt();
                    ds_rel = ds_rel.max((w[i] - w[j]).abs() / l);
                }
            }
            let forces: Vec<PlateForces> = m
                .shells
                .iter()
                .map(|s| {
                    let c: Vec<V3> = s.nodes.iter().map(|&n| m.nodes[n]).collect();
                    let ue: Vec<f64> = s.nodes.iter().flat_map(|&n| u[6 * n..6 * n + 6].to_vec()).collect();
                    let (r, f) = fem::shell::forces(&c, &s.mat, &ue).expect("элемент проверен при сборке");
                    let f = fem::shell::rotate_forces(&r, &f, self.rebar_axis);
                    PlateForces { nx: f[0], ny: f[1], nxy: f[2], mx: f[3], my: f[4], mxy: f[5] }
                })
                .collect();
            let designs = capra_maury::design_uls_batch(self.design, &forces, 10.0);
            for (e, d) in designs.iter().enumerate() {
                impossible |= !d.feasible();
                let a = [d.top.asx, d.top.asy, d.bottom.asx, d.bottom.asy];
                for k in 0..4 {
                    env[e][k] = env[e][k].max(a[k]);
                }
            }
        }
        let (mut steel, mut area) = (0.0, 0.0);
        for (e, s) in m.shells.iter().enumerate() {
            let c: Vec<V3> = s.nodes.iter().map(|&n| m.nodes[n]).collect();
            let a = shell_area(&c);
            steel += env[e].iter().sum::<f64>() * a;
            area += a;
        }
        let kg_m3 = steel * STEEL_DENSITY / (area * self.raft_h);
        let n_max = pile_n.iter().copied().fold(0.0, f64::max);
        let util_ground = n_max / self.pile.fd_ground;
        let util_material = n_max * self.pile.material_util_per_kn;
        let mut ev = Evaluation {
            piles: piles.to_vec(),
            pile_n,
            util_ground,
            util_material,
            s_max,
            ds_rel,
            kg_m3,
            rebar: env,
            w_nodes,
            design_impossible: impossible,
            feasible: false,
        };
        ev.feasible = !impossible && ev.worst(&self.limits) <= 1.0;
        Ok(ev)
    }

    fn spacing_ok(&self, piles: &[usize], skip: usize, cand: usize) -> bool {
        let p = self.xy(cand);
        piles.iter().enumerate().all(|(i, &q)| {
            if i == skip {
                return true;
            }
            let q = self.xy(q);
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt() >= self.rules.min_spacing - 1e-9
        })
    }

    /// Начальная расстановка из ≈ `n` свай: жадный отбор кандидатов (по координатам) с
    /// шагом s, убывающим от √(A/n) до `min_spacing`, пока свай не станет не меньше `n`.
    fn initial(&self, cand: &[usize], n: usize) -> Vec<usize> {
        let area: f64 = self.raft.shells.iter().map(|s| shell_area(&s.nodes.iter().map(|&k| self.raft.nodes[k]).collect::<Vec<_>>())).sum();
        let mut step = (area / n.max(1) as f64).sqrt().max(self.rules.min_spacing);
        loop {
            let mut out: Vec<usize> = Vec::new();
            for &c in cand {
                let p = self.xy(c);
                if out.iter().all(|&q| dist(p, self.xy(q)) >= step - 1e-9) {
                    out.push(c);
                }
            }
            if out.len() >= n || step <= self.rules.min_spacing * 1.0001 {
                return out;
            }
            step = (step * 0.97).max(self.rules.min_spacing);
        }
    }

    /// Локальное улучшение при неизменном числе свай: самую недогруженную сваю переносим к
    /// самой нагруженной, если вариант лучше (сначала по допустимости, затем по кг/м³).
    fn improve(&self, cand: &[usize], mut cur: Evaluation, moves: usize) -> Result<Evaluation, String> {
        let better = |a: &Evaluation, b: &Evaluation| match (a.feasible, b.feasible) {
            (true, false) => true,
            (false, true) => false,
            (true, true) => a.kg_m3 < b.kg_m3 * (1.0 - 1e-6),
            (false, false) => a.worst(&self.limits) < b.worst(&self.limits) * (1.0 - 1e-6),
        };
        for _ in 0..moves {
            let np = cur.piles.len();
            let mut order: Vec<usize> = (0..np).collect();
            order.sort_by(|&a, &b| cur.pile_n[a].total_cmp(&cur.pile_n[b]).then(self.xy(cur.piles[a])[1].total_cmp(&self.xy(cur.piles[b])[1])));
            let (lo, hi) = (order[0], order[np - 1]);
            let target = self.xy(cur.piles[hi]);
            // Кандидаты у перегруженной сваи, ближайшие сначала.
            let mut near: Vec<usize> = cand.iter().copied().filter(|&c| self.spacing_ok(&cur.piles, lo, c) && !cur.piles.contains(&c)).collect();
            near.sort_by(|&a, &b| dist(self.xy(a), target).total_cmp(&dist(self.xy(b), target)));
            let mut moved = false;
            for &c in near.iter().take(3) {
                let mut p = cur.piles.clone();
                p[lo] = c;
                let ev = self.evaluate(&p)?;
                if better(&ev, &cur) {
                    cur = ev;
                    moved = true;
                    break;
                }
            }
            if !moved {
                break;
            }
        }
        Ok(cur)
    }

    /// Фронт вариантов: от начальной расстановки из ≈ `n_start` свай (больше, если она
    /// недопустима) снимаем по одной самую
    /// недогруженную сваю и улучшаем, пока вариант допустим. Возвращает лучшие варианты по
    /// числу свай (недопустимые — последним элементом, если допустимых нет).
    pub fn frontier(&self, n_start: usize, moves: usize, progress: &mut dyn FnMut(&Evaluation)) -> Result<Vec<Evaluation>, String> {
        let cand = self.candidates();
        if cand.is_empty() {
            return Err("нет допустимых позиций свай".into());
        }
        // Недопустимо — больше свай (на 10 %), пока есть место.
        let mut n = n_start;
        let mut cur = self.improve(&cand, self.evaluate(&self.initial(&cand, n))?, moves)?;
        while !cur.feasible {
            let next = self.initial(&cand, (n as f64 * 1.1).ceil() as usize);
            if next.len() <= cur.piles.len() {
                break;
            }
            n = next.len();
            cur = self.improve(&cand, self.evaluate(&next)?, moves)?;
        }
        let mut out = vec![];
        loop {
            if !cur.feasible {
                if out.is_empty() {
                    out.push(cur);
                }
                return Ok(out);
            }
            progress(&cur);
            out.push(cur.clone());
            if cur.piles.len() <= 1 {
                return Ok(out);
            }
            // Снимаем ≈ 3 % самых недогруженных свай (не меньше одной); у границы допустимости,
            // если шаг больше одной сваи и вариант стал недопустим, — повтор по одной.
            let step = ((cur.piles.len() as f64 * 0.03).round() as usize).max(1);
            let next = self.remove(&cand, &cur, step, moves)?;
            cur = if !next.feasible && step > 1 { self.remove(&cand, &cur, 1, moves)? } else { next };
        }
    }

    fn remove(&self, cand: &[usize], cur: &Evaluation, k: usize, moves: usize) -> Result<Evaluation, String> {
        let mut order: Vec<usize> = (0..cur.piles.len()).collect();
        order.sort_by(|&a, &b| cur.pile_n[a].total_cmp(&cur.pile_n[b]));
        let drop: std::collections::HashSet<usize> = order.into_iter().take(k.min(cur.piles.len() - 1)).collect();
        let p: Vec<usize> = cur.piles.iter().enumerate().filter(|(i, _)| !drop.contains(i)).map(|(_, &n)| n).collect();
        self.improve(cand, self.evaluate(&p)?, moves)
    }

    /// Нарушения правил расстановки варианта (проверка результата).
    pub fn layout_violations(&self, piles: &[usize]) -> Vec<layout::Violation> {
        let pts: Vec<[f64; 2]> = piles.iter().map(|&n| self.xy(n)).collect();
        layout::violations(&self.region, &self.rules, &pts)
    }
}

fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

fn shell_area(c: &[V3]) -> f64 {
    let sub = |a: V3, b: V3| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let d = if c.len() == 4 { fem::shell::cross(sub(c[2], c[0]), sub(c[3], c[1])) } else { fem::shell::cross(sub(c[1], c[0]), sub(c[2], c[0])) };
    0.5 * fem::shell::norm(d)
}
