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
use crate::sp63::capra_maury::{self, PlateForces, SlsContext, UlsContext};

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
    /// Загружения плиты (нормативные значения) и сочетания — коэффициенты к загружениям
    /// (`loads::combinations`): ПС1 — нагрузки на сваи и армирование; ПС2 полные — осадки
    /// и раскрытие трещин непродолжительное; ПС2 длительное — длительное раскрытие. Пустые
    /// ПС2 — осадки по ПС1.
    pub cases: Vec<LoadCase>,
    pub uls: Vec<Vec<f64>>,
    pub sls: Vec<Vec<f64>>,
    pub sls_long: Option<Vec<f64>>,
    pub limits: SettlementLimits,
    /// Подбор армирования плиты (ПС1) и ось x армирования.
    pub design: &'a UlsContext,
    pub rebar_axis: V3,
    /// Толщина плиты, м (для кг/м³).
    pub raft_h: f64,
    /// Надземная часть, сконденсированная на плиту (`fem::condense::building`); приведённые
    /// нагрузки — по одному набору на загружение. None — здание не учитывается.
    pub building: Option<&'a Building>,
    /// Трещины: непродолжительное раскрытие (полные сочетания ПС2) и продолжительное
    /// (длительное сочетание); контексты задают φ1 и предел.
    pub crack_short: Option<&'a SlsContext>,
    pub crack_long: Option<&'a SlsContext>,
    /// Продавливание плиты сваями (пп. 8.1.46–8.1.47 СП 63, без поперечной арматуры; у края — доля замкнутого контура внутри плиты (упрощение двух вариантов контура п. 8.1.47)).
    pub punching: Option<Punching>,
    /// Существующие сваи (ориентир, не правило): старт свободного режима и предпочтение
    /// позиций рядом с ними.
    pub existing: Vec<[f64; 2]>,
}

/// Данные проверки продавливания.
#[derive(Clone, Debug)]
pub struct Punching {
    /// Rbt, кПа.
    pub rbt: f64,
    /// Рабочая высота плиты, м.
    pub h0: f64,
    /// Контур плиты (часть расчётного контура вне плиты не учитывается).
    pub raft: Region,
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
    /// Наибольший Кисп продавливания плиты сваей (0 — не проверялось).
    pub util_punch: f64,
    /// Среднее расстояние от свай до ближайшей существующей сваи, м (NaN — ориентира нет).
    pub shift: f64,
    /// Подбор невозможен хотя бы в одном элементе.
    pub design_impossible: bool,
    pub feasible: bool,
}

impl Evaluation {
    /// Наибольший Кисп из ограничений (для сравнения недопустимых вариантов).
    pub fn worst(&self, lim: &SettlementLimits) -> f64 {
        self.util_ground.max(self.util_material).max(self.util_punch).max(self.s_max / lim.s_max).max(self.ds_rel / lim.ds_rel_max)
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
        self.evaluate_with(piles, true)
    }

    /// Расчёт без подбора армирования: нагрузки на сваи и осадки (`kg_m3` = NaN).
    pub fn evaluate_fast(&self, piles: &[usize]) -> Result<Evaluation, String> {
        self.evaluate_with(piles, false)
    }

    fn evaluate_with(&self, piles: &[usize], design: bool) -> Result<Evaluation, String> {
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
        let mut cases = self.cases.clone();
        if let Some(b) = self.building {
            if b.loads.len() != cases.len() {
                return Err("здание: число наборов приведённых нагрузок не равно числу загружений".into());
            }
            m.blocks.push(b.block.clone());
            for &d in b.block.dofs.iter().step_by(3) {
                m.fixed.push((d / 6, [true, true, false, false, false, true]));
            }
            for (lc, extra) in cases.iter_mut().zip(&b.loads) {
                lc.nodal.extend(extra.iter().copied());
            }
        }
        let lin = m.factorize()?;
        let u_cases: Vec<Vec<f64>> = cases.iter().map(|lc| lin.solve(&m, lc)).collect::<Result<_, _>>()?;
        let combine = |c: &[f64]| -> Vec<f64> {
            let mut u = vec![0.0; m.ndof()];
            for (f, uc) in c.iter().zip(&u_cases) {
                if *f != 0.0 {
                    for (a, b) in u.iter_mut().zip(uc) {
                        *a += f * b;
                    }
                }
            }
            u
        };
        let np = piles.len();
        let mut pile_n = vec![0.0f64; np];
        let mut s_max: f64 = 0.0;
        let mut ds_rel: f64 = 0.0;
        let ne = m.shells.len();
        let mut env = vec![[0.0f64; 4]; ne];
        let mut w_nodes = vec![0.0f64; m.nodes.len()];
        let mut impossible = false;
        let forces_of = |u: &[f64]| -> Vec<PlateForces> {
            m.shells
                .iter()
                .map(|s| {
                    let c: Vec<V3> = s.nodes.iter().map(|&n| m.nodes[n]).collect();
                    let ue: Vec<f64> = s.nodes.iter().flat_map(|&n| u[6 * n..6 * n + 6].to_vec()).collect();
                    let (r, f) = fem::shell::forces(&c, &s.mat, &ue).expect("элемент проверен при сборке");
                    let f = fem::shell::rotate_forces(&r, &f, self.rebar_axis);
                    PlateForces { nx: f[0], ny: f[1], nxy: f[2], mx: f[3], my: f[4], mxy: f[5] }
                })
                .collect()
        };
        let mut add = |designs: &[capra_maury::Design], env: &mut Vec<[f64; 4]>| {
            for (e, d) in designs.iter().enumerate() {
                impossible |= !d.feasible();
                let a = [d.top.asx, d.top.asy, d.bottom.asx, d.bottom.asy];
                for k in 0..4 {
                    env[e][k] = env[e][k].max(a[k]);
                }
            }
        };
        // ПС1: нагрузки на сваи (R = −K·w, сжатие — положительное) и армирование.
        for c in &self.uls {
            let u = combine(c);
            let w: Vec<f64> = piles.iter().map(|&n| u[6 * n + 2]).collect();
            for i in 0..np {
                let r: f64 = -(0..np).map(|j| k[i * np + j] * w[j]).sum::<f64>();
                pile_n[i] = pile_n[i].max(r);
            }
            if design {
                add(&capra_maury::design_uls_batch(self.design, &forces_of(&u), 10.0), &mut env);
            }
        }
        // ПС2: осадки и разность осадок; непродолжительное раскрытие трещин.
        let settle = if self.sls.is_empty() { &self.uls } else { &self.sls };
        for c in settle {
            let u = combine(c);
            let w: Vec<f64> = piles.iter().map(|&n| u[6 * n + 2]).collect();
            let s_lc = w.iter().map(|x| -x).fold(0.0, f64::max);
            if s_lc >= s_max {
                w_nodes = (0..m.nodes.len()).map(|n| -u[6 * n + 2]).collect();
            }
            for i in 0..np {
                s_max = s_max.max(-w[i]);
                for j in i + 1..np {
                    let l = dist(self.xy(piles[i]), self.xy(piles[j]));
                    ds_rel = ds_rel.max((w[i] - w[j]).abs() / l);
                }
            }
            if design {
                if let Some(ctx) = self.crack_short {
                    add(&capra_maury::design_sls_batch(ctx, &forces_of(&u), 10.0), &mut env);
                }
            }
        }
        if design {
            if let (Some(ctx), Some(c)) = (self.crack_long, &self.sls_long) {
                add(&capra_maury::design_sls_batch(ctx, &forces_of(&combine(c)), 10.0), &mut env);
            }
        }
        let (mut steel, mut area) = (0.0, 0.0);
        for (e, s) in m.shells.iter().enumerate() {
            let c: Vec<V3> = s.nodes.iter().map(|&n| m.nodes[n]).collect();
            let a = shell_area(&c);
            steel += env[e].iter().sum::<f64>() * a;
            area += a;
        }
        let kg_m3 = if design { steel * STEEL_DENSITY / (area * self.raft_h) } else { f64::NAN };
        let n_max = pile_n.iter().copied().fold(0.0, f64::max);
        let util_ground = n_max / self.pile.fd_ground;
        let util_material = n_max * self.pile.material_util_per_kn;
        let util_punch = match &self.punching {
            None => 0.0,
            Some(p) => {
                let r = 0.5 * (self.pile.d + p.h0);
                let u = std::f64::consts::PI * 2.0 * r;
                piles
                    .iter()
                    .zip(&pile_n)
                    .map(|(&n, &f)| {
                        let [x, y] = self.xy(n);
                        let inside = (0..72)
                            .filter(|&k| {
                                let a = k as f64 * std::f64::consts::PI / 36.0;
                                p.raft.contains([x + r * a.cos(), y + r * a.sin()])
                            })
                            .count() as f64
                            / 72.0;
                        f / (p.rbt * u * inside.max(1e-9) * p.h0)
                    })
                    .fold(0.0, f64::max)
            }
        };
        let shift = if self.existing.is_empty() {
            f64::NAN
        } else {
            piles.iter().map(|&n| self.existing_distance(self.xy(n))).sum::<f64>() / piles.len().max(1) as f64
        };
        let mut ev = Evaluation {
            util_punch,
            shift,
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

    fn existing_distance(&self, p: [f64; 2]) -> f64 {
        self.existing.iter().map(|&q| dist(p, q)).fold(f64::INFINITY, f64::min)
    }

    /// Существующие сваи, приведённые к допустимым позициям: каждая — к ближайшему
    /// кандидату не дальше max(min_spacing, 2·шаг сетки), с соблюдением расстояний.
    pub fn snap_existing(&self) -> Vec<usize> {
        let cand = self.candidates();
        let reach = self.rules.min_spacing.max(2.0 * self.mesh_step());
        let mut pts = self.existing.clone();
        pts.sort_by(|a, b| a[1].total_cmp(&b[1]).then(a[0].total_cmp(&b[0])));
        let mut out: Vec<usize> = Vec::new();
        for p in pts {
            let best = cand
                .iter()
                .copied()
                .filter(|&c| dist(self.xy(c), p) <= reach && !out.contains(&c) && self.spacing_ok(&out, usize::MAX, c))
                .min_by(|&a, &b| dist(self.xy(a), p).total_cmp(&dist(self.xy(b), p)));
            if let Some(c) = best {
                out.push(c);
            }
        }
        out
    }

    /// Расчёт существующей расстановки (для сравнения); None — ориентира нет.
    pub fn reference(&self) -> Result<Option<Evaluation>, String> {
        if self.existing.is_empty() {
            return Ok(None);
        }
        let p = self.snap_existing();
        if p.is_empty() {
            return Ok(None);
        }
        self.evaluate(&p).map(Some)
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

    /// Локальное улучшение при неизменном числе свай по быстрой оценке: самую недогруженную
    /// сваю переносим к самой нагруженной, если вариант лучше (сначала по допустимости, затем
    /// по наибольшему коэффициенту использования — выравнивание нагрузок). Возвращает полный
    /// расчёт (с подбором армирования) итоговой расстановки.
    fn improve(&self, cand: &[usize], piles: Vec<usize>, moves: usize) -> Result<Evaluation, String> {
        let better = |a: &Evaluation, b: &Evaluation| match (a.feasible, b.feasible) {
            (true, false) => true,
            (false, true) => false,
            _ => a.worst(&self.limits) < b.worst(&self.limits) * (1.0 - 1e-6),
        };
        let mut cur = self.evaluate_fast(&piles)?;
        for _ in 0..moves {
            let np = cur.piles.len();
            let mut order: Vec<usize> = (0..np).collect();
            order.sort_by(|&a, &b| cur.pile_n[a].total_cmp(&cur.pile_n[b]).then(self.xy(cur.piles[a])[1].total_cmp(&self.xy(cur.piles[b])[1])));
            let (lo, hi) = (order[0], order[np - 1]);
            let target = self.xy(cur.piles[hi]);
            let mut near: Vec<usize> = cand.iter().copied().filter(|&c| self.spacing_ok(&cur.piles, lo, c) && !cur.piles.contains(&c)).collect();
            // Ориентир: при наличии существующих свай предпочтение позициям рядом с ними.
            let key = |c: usize| dist(self.xy(c), target) + if self.existing.is_empty() { 0.0 } else { 0.5 * self.existing_distance(self.xy(c)) };
            near.sort_by(|&a, &b| key(a).total_cmp(&key(b)));
            let mut moved = false;
            for &c in near.iter().take(3) {
                let mut p = cur.piles.clone();
                p[lo] = c;
                let ev = self.evaluate_fast(&p)?;
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
        self.evaluate(&cur.piles)
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
        // Старт — существующие сваи (если заданы), иначе жадная расстановка.
        let snapped = self.snap_existing();
        let mut n = if snapped.is_empty() { n_start } else { snapped.len() };
        let mut cur = self.improve(&cand, if snapped.is_empty() { self.initial(&cand, n) } else { snapped }, moves)?;
        while !cur.feasible {
            let next = self.initial(&cand, (n as f64 * 1.1).ceil() as usize);
            if next.len() <= cur.piles.len() {
                break;
            }
            n = next.len();
            cur = self.improve(&cand, next, moves)?;
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
        self.improve(cand, p, moves)
    }

    /// Расстановка рядами: начальная сетка — ряды вдоль X и Y с наименьшим шагом, кратным
    /// шагу сетки плиты и не меньшим `min_spacing`, по центру допустимых позиций; далее по
    /// одному снимается целый ряд или столбец — из всех вариантов снятия (быстрая оценка)
    /// выбирается допустимый с наименьшим коэффициентом использования. Каждый шаг — вариант
    /// фронта (с подбором армирования); в выдаче — недоминируемые по (число свай, кг/м³), по
    /// убыванию числа свай. Если допустима не даже самая плотная сетка — она одна.
    pub fn regular(&self, progress: &mut dyn FnMut(&Evaluation)) -> Result<Vec<Evaluation>, String> {
        let cand = self.candidates();
        if cand.is_empty() {
            return Err("нет допустимых позиций свай".into());
        }
        let step = self.mesh_step();
        let (x0, x1, y0, y1) = cand.iter().fold((f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY), |a, &c| {
            let [x, y] = self.xy(c);
            (a.0.min(x), a.1.max(x), a.2.min(y), a.3.max(y))
        });
        let pitch = (self.rules.min_spacing / step - 1e-9).ceil().max(1.0) * step;
        let lines = |lo: f64, hi: f64| -> Vec<f64> {
            let o = lo + (((hi - lo) % pitch) / 2.0 / step).floor() * step;
            (0..).map(|i| o + i as f64 * pitch).take_while(|&v| v <= hi + 1e-9).collect()
        };
        let tol = 0.25 * step;
        let layout = |xs: &[f64], ys: &[f64]| -> Vec<usize> {
            cand.iter()
                .copied()
                .filter(|&c| {
                    let [x, y] = self.xy(c);
                    xs.iter().any(|&v| (x - v).abs() < tol) && ys.iter().any(|&v| (y - v).abs() < tol)
                })
                .collect()
        };
        let (mut xs, mut ys) = (lines(x0, x1), lines(y0, y1));
        let mut cur = self.evaluate_fast(&layout(&xs, &ys))?;
        if cur.piles.is_empty() {
            return Err("ни один ряд не помещается в допустимую область".into());
        }
        let mut found: Vec<Evaluation> = Vec::new();
        while cur.feasible {
            let full = self.evaluate(&cur.piles)?;
            progress(&full);
            found.push(full);
            // Варианты снятия одного ряда (X) или столбца (Y).
            let mut best: Option<(Evaluation, Vec<f64>, Vec<f64>)> = None;
            for axis in 0..2 {
                let n = if axis == 0 { xs.len() } else { ys.len() };
                if n <= 2 {
                    continue;
                }
                for i in 0..n {
                    let (mut nx, mut ny) = (xs.clone(), ys.clone());
                    if axis == 0 { nx.remove(i); } else { ny.remove(i); }
                    let piles = layout(&nx, &ny);
                    if piles.len() < 3 {
                        continue;
                    }
                    let ev = self.evaluate_fast(&piles)?;
                    let better = match &best {
                        None => true,
                        Some((b, _, _)) => match (ev.feasible, b.feasible) {
                            (true, false) => true,
                            (false, true) => false,
                            _ => ev.worst(&self.limits) < b.worst(&self.limits),
                        },
                    };
                    if better {
                        best = Some((ev, nx, ny));
                    }
                }
            }
            match best {
                Some((ev, nx, ny)) if ev.feasible => {
                    cur = ev;
                    xs = nx;
                    ys = ny;
                }
                _ => break,
            }
        }
        if found.is_empty() {
            return Ok(vec![self.evaluate(&cur.piles)?]);
        }
        // От меньшего числа свай к большему, при равном числе — меньший расход первым:
        // вариант остаётся, только если его расход меньше всех оставленных (с меньшим числом).
        found.sort_by(|a, b| a.piles.len().cmp(&b.piles.len()).then(a.kg_m3.total_cmp(&b.kg_m3)));
        let mut out: Vec<Evaluation> = Vec::new();
        for e in found {
            if out.last().is_none_or(|o| o.piles.len() != e.piles.len() && e.kg_m3 < o.kg_m3 * (1.0 - 1e-9)) {
                out.push(e);
            }
        }
        out.reverse();
        Ok(out)
    }

    /// Шаг сетки плиты: наименьшая ненулевая длина ребра элементов.
    fn mesh_step(&self) -> f64 {
        self.raft
            .shells
            .iter()
            .flat_map(|s| (0..s.nodes.len()).map(move |i| (s.nodes[i], s.nodes[(i + 1) % s.nodes.len()])))
            .map(|(a, b)| dist(self.xy(a), self.xy(b)))
            .filter(|&d| d > 1e-9)
            .fold(f64::INFINITY, f64::min)
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
