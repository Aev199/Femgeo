//! Подбор расстановки свай под сплошной плитой (`docs/SPEC.md`, §3): фронт вариантов
//! «число свай — расход армирования плиты, кг/м³» при ограничениях по несущей способности
//! свай по грунту и материалу, контуру, правилам СП 24 и осадкам.
//!
//! Первая версия: сваи ставятся в узлы сетки плиты; каждый вариант считается заново
//! (свайное поле по СП 24 → плита МКЭ → нагрузки на сваи, осадки → подбор Капра–Мори).
//! Алгоритм детерминирован и не зависит от нумерации кандидатов.

use crate::fem::condense::Building;
use std::cell::OnceCell;
use std::collections::HashMap;

use crate::fem::{self, LoadCase, Model, Spring, V3};
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
    /// Предел сваи Fu, кН: в каждом сочетании сваи с N > Fu держат постоянную силу Fu,
    /// избыток перераспределяется (активное множество). None — без предела.
    pub pile_fu: Option<f64>,
    /// Жёсткость свай по осадке условного фундамента (`piles::conditional`): одинаковая
    /// k = N / (n·s), N — наибольшая суммарная нагрузка длительного сочетания ПС2 (иначе
    /// полных ПС2). None — взаимовлияние по СП 24, п. 7.4.
    pub conditional: Option<ConditionalModel>,
    /// Продавливание плиты вертикальными конструкциями (нужны `punching` и `building`).
    pub structure_punch: Vec<PunchArea>,
    /// Подготовка расчёта (заполняется при первой оценке): `Default::default()`.
    pub cache: OnceCell<Prepared>,
}

/// Один раз для постановки: разложение модели без свай K0 (плита, здание, регуляризующие
/// пружины 1e-6 от жёсткости сваи), перемещения по загружениям u0 и столбцы G = K0⁻¹·e_uz
/// для всех кандидатных узлов — на степенях свободы узлов плиты. Сваи учитываются по
/// Вудбери: (F + G_SS)·t = w0, u = u0 − G·t, F — податливость свайного поля (СП 24).
pub struct Prepared {
    index: HashMap<usize, usize>,
    /// Узел плиты → номер в `nodes`.
    pos: HashMap<usize, usize>,
    nodes: Vec<usize>,
    /// По загружениям: перемещения узлов плиты (6 на узел).
    u0: Vec<Vec<f64>>,
    /// По кандидатам: перемещения узлов плиты от единичной силы +z в кандидате.
    g: Vec<Vec<f64>>,
    /// Пары узлов соседних вертикальных конструкций (узлы сопряжения здания, разные
    /// конструкции, не дальше `STRUCTURE_SPAN`) и расстояние между ними — для Δs/L.
    structure_pairs: Vec<(usize, usize, f64)>,
}

/// Наибольшая доля свай на пределе Fu в одном сочетании для допустимого варианта
/// (решение пользователя 2026-11-02).
pub const CAPPED_SHARE_MAX: f64 = 0.10;

/// Наибольшее расстояние между узлами соседних конструкций для Δs/L, м (допущение: пролёт).
pub const STRUCTURE_SPAN: f64 = 12.0;

/// Узлы сопряжения → конструкции (связные цепочки узлов не дальше 1,5 шага сетки) → пары
/// узлов разных конструкций на расстоянии до `STRUCTURE_SPAN`.
fn structure_pairs(xy: &[[f64; 2]], step: f64) -> Vec<(usize, usize, f64)> {
    let n = xy.len();
    let mut comp: Vec<usize> = (0..n).collect();
    fn root(c: &mut Vec<usize>, i: usize) -> usize {
        let mut r = i;
        while c[r] != r {
            r = c[r];
        }
        c[i] = r;
        r
    }
    for i in 0..n {
        for j in i + 1..n {
            if dist(xy[i], xy[j]) <= 1.5 * step {
                let (a, b) = (root(&mut comp, i), root(&mut comp, j));
                comp[a] = b;
            }
        }
    }
    let c: Vec<usize> = (0..n).map(|i| root(&mut comp, i)).collect();
    let mut out = vec![];
    for i in 0..n {
        for j in i + 1..n {
            let l = dist(xy[i], xy[j]);
            if c[i] != c[j] && l <= STRUCTURE_SPAN {
                out.push((i, j, l));
            }
        }
    }
    out
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

/// Продавливание плиты вертикальной конструкцией (колонна, торец или угол стены).
#[derive(Clone, Debug)]
pub struct PunchArea {
    /// "column" | "wall_end" | "wall_corner" — для отчёта.
    pub kind: String,
    pub at: [f64; 2],
    /// Узлы плиты, сила здания в которых — нагрузка на загруженную площадь.
    pub nodes: Vec<usize>,
    /// Длина расчётного контура (м) без учёта края плиты; точки контура — доля внутри плиты.
    pub u: f64,
    pub samples: Vec<[f64; 2]>,
    /// Область внутри контура: сваи с центром в ней разгружают продавливание.
    pub inside: Region,
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
    /// Наибольший Кисп продавливания плиты колоннами, торцами и углами стен (0 — не проверялось).
    pub util_punch_struct: f64,
    /// Среднее расстояние от свай до ближайшей существующей сваи, м (NaN — ориентира нет).
    pub shift: f64,
    /// Подбор невозможен хотя бы в одном элементе.
    pub design_impossible: bool,
    /// Хотя бы в одном сочетании сваи ограничены пределом Fu; наибольшая доля таких свай.
    pub capped: bool,
    pub capped_share: f64,
    /// Условный фундамент (если жёсткость свай по нему).
    pub conditional: Option<crate::piles::conditional::Conditional>,
    pub feasible: bool,
}

/// Исходные данные условного фундамента: контур по наружным граням свай и грунт в его центре.
#[derive(Clone, Debug)]
pub struct ConditionalModel {
    pub outline: Vec<[f64; 2]>,
    pub column: crate::piles::SoilColumn,
}

impl Evaluation {
    /// Наибольший Кисп из ограничений (для сравнения недопустимых вариантов).
    pub fn worst(&self, lim: &SettlementLimits) -> f64 {
        self.util_ground.max(self.util_material).max(self.util_punch).max(self.util_punch_struct).max(self.s_max / lim.s_max).max(self.ds_rel / lim.ds_rel_max)
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
        let field = if self.conditional.is_some() { None } else { Some(Field::new(&ps, self.soil)) };
        let p = self.prepared()?;
        let m = &self.raft;
        let np = piles.len();
        let idx: Vec<usize> = piles.iter().map(|n| p.index.get(n).copied().ok_or("свая не в допустимой позиции")).collect::<Result<_, _>>()?;
        let at = |v: &[f64], node: usize| v[6 * p.pos[&node] + 2];
        // A = F + G_SS (м/кН), симметричная положительно определённая.
        let mut gss = vec![0.0; np * np];
        for i in 0..np {
            for j in 0..np {
                gss[i * np + j] = 0.5 * (at(&p.g[idx[j]], piles[i]) + at(&p.g[idx[i]], piles[j]));
            }
        }
        let w0: Vec<Vec<f64>> = p.u0.iter().map(|u0| piles.iter().map(|&n| at(u0, n)).collect()).collect();
        let mut conditional = None;
        let a: Vec<f64> = match (&field, &self.conditional) {
            (Some(f), _) => f.flexibility.iter().zip(&gss).map(|(x, y)| x + y).collect(),
            (None, Some(cm)) => {
                // Суммарные нагрузки загружений — из равновесия при любой одинаковой жёсткости.
                let probe = 1e-5;
                let mut a = gss.clone();
                (0..np).for_each(|i| a[i * np + i] += probe);
                let ch = cholesky(&a, np)?;
                let totals: Vec<f64> = w0.iter().map(|w| -chol_solve(&ch, np, w).iter().sum::<f64>()).collect();
                let combos: &[Vec<f64>] = match &self.sls_long {
                    Some(c) => std::slice::from_ref(c),
                    None => &self.sls,
                };
                let n_sls = combos.iter().map(|c| c.iter().zip(&totals).map(|(f, t)| f * t).sum::<f64>()).fold(0.0, f64::max);
                let cf = crate::piles::conditional::settlement(&cm.column, &cm.outline, self.pile_head, self.pile.length, n_sls)?;
                let flex = cf.s * np as f64 / n_sls.max(1e-9);
                conditional = Some(cf);
                let mut a = gss.clone();
                (0..np).for_each(|i| a[i * np + i] += flex);
                a
            }
            (None, None) => unreachable!(),
        };
        let chol = cholesky(&a, np)?;
        let mut u_cases = Vec::with_capacity(self.cases.len());
        let mut t_cases = Vec::with_capacity(self.cases.len());
        for u0 in &p.u0 {
            let w0: Vec<f64> = piles.iter().map(|&n| at(u0, n)).collect();
            let t = chol_solve(&chol, np, &w0);
            let mut u = u0.clone();
            for (j, &tj) in t.iter().enumerate() {
                for (x, gx) in u.iter_mut().zip(&p.g[idx[j]]) {
                    *x -= tj * gx;
                }
            }
            u_cases.push(u);
            t_cases.push(t);
        }
        // Перемещения в глобальной нумерации узлов плиты.
        let global = |u: &[f64]| -> Vec<f64> {
            let mut g = vec![0.0; m.ndof()];
            for (k, &n) in p.nodes.iter().enumerate() {
                g[6 * n..6 * n + 6].copy_from_slice(&u[6 * k..6 * k + 6]);
            }
            g
        };
        let combine = |c: &[f64]| -> Vec<f64> {
            let mut u = vec![0.0; 6 * p.nodes.len()];
            for (f, uc) in c.iter().zip(&u_cases) {
                if *f != 0.0 {
                    for (a, b) in u.iter_mut().zip(uc) {
                        *a += f * b;
                    }
                }
            }
            global(&u)
        };
        // Сочетание: перемещения плиты и нагрузки на сваи (сжатие > 0); при пределе Fu сваи с
        // N > Fu держат Fu, остальные — из совместности (активное множество, только добавление).
        let mut capped_any = false;
        let mut overloaded = false;
        let mut capped_share: f64 = 0.0;
        let mut solve = |c: &[f64]| -> Result<(Vec<f64>, Vec<f64>), String> {
            let t_lin: Vec<f64> = (0..np).map(|i| c.iter().zip(&t_cases).map(|(f, t)| f * t[i]).sum::<f64>()).collect();
            let mut u = combine(c);
            let Some(fu) = self.pile_fu else { return Ok((u, t_lin.iter().map(|t| -t).collect())) };
            let mut t = t_lin.clone();
            let mut capped = vec![false; np];
            loop {
                let over: Vec<usize> = (0..np).filter(|&i| !capped[i] && -t[i] > fu * (1.0 + 1e-9)).collect();
                if over.is_empty() {
                    break;
                }
                over.iter().for_each(|&i| capped[i] = true);
                let free: Vec<usize> = (0..np).filter(|&i| !capped[i]).collect();
                if free.is_empty() {
                    // Сваи не держат нагрузку даже на пределе: вариант недопустим, нагрузки — линейно.
                    overloaded = true;
                    return Ok((u, t_lin.iter().map(|t| -t).collect()));
                }
                let w: Vec<f64> = (0..np).map(|i| c.iter().zip(&w0).map(|(f, w)| f * w[i]).sum()).collect();
                let nf = free.len();
                let mut aff = vec![0.0; nf * nf];
                let mut rhs = vec![0.0; nf];
                for (x, &i) in free.iter().enumerate() {
                    rhs[x] = w[i] + (0..np).filter(|&j| capped[j]).map(|j| a[i * np + j] * fu).sum::<f64>();
                    for (y, &j) in free.iter().enumerate() {
                        aff[x * nf + y] = a[i * np + j];
                    }
                }
                let tf = chol_solve(&cholesky(&aff, nf)?, nf, &rhs);
                (0..np).filter(|&j| capped[j]).for_each(|j| t[j] = -fu);
                free.iter().zip(&tf).for_each(|(&i, &v)| t[i] = v);
            }
            let nc = capped.iter().filter(|&&x| x).count();
            capped_share = capped_share.max(nc as f64 / np as f64);
            if nc > 0 {
                capped_any = true;
                let mut du = vec![0.0; 6 * p.nodes.len()];
                for j in 0..np {
                    let dt = t[j] - t_lin[j];
                    if dt != 0.0 {
                        for (x, gx) in du.iter_mut().zip(&p.g[idx[j]]) {
                            *x -= dt * gx;
                        }
                    }
                }
                let dg = global(&du);
                u.iter_mut().zip(&dg).for_each(|(x, d)| *x += d);
            }
            Ok((u, t.iter().map(|x| -x).collect()))
        };
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
        // Вертикальная сила здания на узлы плиты в сочетании (вниз > 0): r = f − K·u суперэлемента.
        let building_down = |c: &[f64], u: &[f64]| -> HashMap<usize, f64> {
            let mut out: HashMap<usize, f64> = HashMap::new();
            let Some(b) = self.building else { return out };
            let nd = b.block.dofs.len();
            for (i, &di) in b.block.dofs.iter().enumerate() {
                if di % 6 != 2 {
                    continue;
                }
                let ku: f64 = (0..nd).map(|j| b.block.k[i * nd + j] * u[b.block.dofs[j]]).sum();
                *out.entry(di / 6).or_default() += ku;
            }
            for (f, loads) in c.iter().zip(&b.loads) {
                for &(n, l) in loads {
                    *out.entry(n).or_default() -= f * l[2];
                }
            }
            out
        };
        let mut util_punch_struct: f64 = 0.0;
        for c in &self.uls {
            let (u, n) = solve(c)?;
            if let (Some(pp), false) = (&self.punching, self.structure_punch.is_empty()) {
                let down = building_down(c, &u);
                for a in &self.structure_punch {
                    let f: f64 = a.nodes.iter().map(|k| down.get(k).copied().unwrap_or(0.0)).sum();
                    let relief: f64 = piles.iter().zip(&n).filter(|(&q, _)| a.inside.contains(self.xy(q))).map(|(_, &x)| x).sum();
                    let frac = a.samples.iter().filter(|&&q| pp.raft.contains(q)).count() as f64 / a.samples.len().max(1) as f64;
                    let util = (f - relief).max(0.0) / (pp.rbt * a.u * frac.max(1e-9) * pp.h0);
                    util_punch_struct = util_punch_struct.max(util);
                }
            }
            for (i, r) in n.into_iter().enumerate() {
                pile_n[i] = pile_n[i].max(r);
            }
            if design {
                add(&capra_maury::design_uls_batch(self.design, &forces_of(&u), 10.0), &mut env);
            }
        }
        // ПС2: осадки и разность осадок; непродолжительное раскрытие трещин.
        let settle = if self.sls.is_empty() { &self.uls } else { &self.sls };
        for c in settle {
            let (u, _) = solve(c)?;
            let w: Vec<f64> = piles.iter().map(|&n| u[6 * n + 2]).collect();
            let s_lc = w.iter().map(|x| -x).fold(0.0, f64::max);
            if s_lc >= s_max {
                w_nodes = (0..m.nodes.len()).map(|n| -u[6 * n + 2]).collect();
            }
            for i in 0..np {
                s_max = s_max.max(-w[i]);
            }
            // Δs/L: между соседними вертикальными конструкциями здания; без здания — между сваями.
            if p.structure_pairs.is_empty() {
                for i in 0..np {
                    for j in i + 1..np {
                        let l = dist(self.xy(piles[i]), self.xy(piles[j]));
                        ds_rel = ds_rel.max((w[i] - w[j]).abs() / l);
                    }
                }
            } else {
                for &(a, b, l) in &p.structure_pairs {
                    ds_rel = ds_rel.max((u[6 * a + 2] - u[6 * b + 2]).abs() / l);
                }
            }
            if design {
                if let Some(ctx) = self.crack_short {
                    let d = crack_designs(ctx, &forces_of(&u), &env);
                    add_some(&d, &mut env, &mut impossible);
                }
            }
        }
        if design {
            if let (Some(ctx), Some(c)) = (self.crack_long, &self.sls_long) {
                let d = crack_designs(ctx, &forces_of(&solve(c)?.0), &env);
                add_some(&d, &mut env, &mut impossible);
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
            util_punch_struct,
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
            capped: capped_any,
            capped_share,
            conditional,
            feasible: false,
        };
        ev.feasible = !impossible && !overloaded && capped_share <= CAPPED_SHARE_MAX + 1e-12 && ev.worst(&self.limits) <= 1.0;
        Ok(ev)
    }

    /// Подготовка (один раз): K0, u0, G — см. [`Prepared`].
    fn prepared(&self) -> Result<&Prepared, String> {
        if let Some(p) = self.cache.get() {
            return Ok(p);
        }
        let mut m = self.raft.clone();
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
        let mut nodes: Vec<usize> = self.raft.shells.iter().flat_map(|s| s.nodes.iter().copied()).collect();
        nodes.sort_unstable();
        nodes.dedup();
        // Регуляризация: без свай плита не закреплена по вертикали.
        let c = nodes.iter().fold([0.0, 0.0], |a, &n| [a[0] + self.xy(n)[0], a[1] + self.xy(n)[1]]);
        let probe = Pile { x: c[0] / nodes.len() as f64, y: c[1] / nodes.len() as f64, head: self.pile_head, length: self.pile.length, d: self.pile.d, ea: self.pile.ea };
        let k_pile = 1.0 / crate::piles::SinglePile::new(&probe, &(self.soil)(probe.x, probe.y)).compliance;
        for &n in &nodes {
            m.springs.push(Spring { node: n, k: [0.0, 0.0, 1e-6 * k_pile, 0.0, 0.0, 0.0] });
        }
        let lin = m.factorize()?;
        let pick = |u: &[f64]| -> Vec<f64> { nodes.iter().flat_map(|&n| u[6 * n..6 * n + 6].to_vec()).collect() };
        let fs: Vec<Vec<f64>> = cases.iter().map(|lc| m.load_vector(lc)).collect::<Result<_, _>>()?;
        let u0 = lin.solve_many(&fs).iter().map(|u| pick(u)).collect();
        let cand = self.candidates();
        let units: Vec<Vec<f64>> = cand
            .iter()
            .map(|&n| {
                let mut f = vec![0.0; m.ndof()];
                f[6 * n + 2] = 1.0;
                f
            })
            .collect();
        let g = lin.solve_many(&units).iter().map(|u| pick(u)).collect();
        let index = cand.iter().enumerate().map(|(i, &n)| (n, i)).collect();
        let pos = nodes.iter().enumerate().map(|(i, &n)| (n, i)).collect();
        let structure_pairs = match self.building {
            Some(b) => {
                let mut ns: Vec<usize> = b.block.dofs.iter().step_by(3).map(|d| d / 6).collect();
                ns.sort_unstable();
                ns.dedup();
                let xy: Vec<[f64; 2]> = ns.iter().map(|&n| self.xy(n)).collect();
                structure_pairs(&xy, self.mesh_step()).into_iter().map(|(i, j, l)| (ns[i], ns[j], l)).collect()
            }
            None => vec![],
        };
        let _ = self.cache.set(Prepared { index, pos, nodes, u0, g, structure_pairs });
        Ok(self.cache.get().expect("только что заполнено"))
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
            // Варианты снятия одного ряда (X) или столбца (Y): допустимые, с наименьшим Кисп;
            // среди близких (в пределах 2 %) — наиболее симметричный набор рядов.
            let mut options: Vec<(Evaluation, Vec<f64>, Vec<f64>)> = Vec::new();
            for axis in 0..2 {
                let n = if axis == 0 { xs.len() } else { ys.len() };
                if n <= 2 {
                    continue;
                }
                let (v, lo, hi) = if axis == 0 { (&xs, x0, x1) } else { (&ys, y0, y1) };
                // Одиночные ряды и пары зеркальных рядов (относительно середины области).
                let mut sets: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
                for i in 0..n {
                    for j in i + 1..n {
                        if (v[i] + v[j] - lo - hi).abs() < tol && n - 2 >= 2 {
                            sets.push(vec![i, j]);
                        }
                    }
                }
                for set in sets {
                    let keep = |w: &[f64]| -> Vec<f64> { w.iter().enumerate().filter(|(k, _)| !set.contains(k)).map(|(_, &x)| x).collect() };
                    let (nx, ny) = if axis == 0 { (keep(&xs), ys.clone()) } else { (xs.clone(), keep(&ys)) };
                    let piles = layout(&nx, &ny);
                    if piles.len() < 3 {
                        continue;
                    }
                    options.push((self.evaluate_fast(&piles)?, nx, ny));
                }
            }
            let asym = |v: &[f64], lo: f64, hi: f64| -> f64 {
                v.iter().map(|&a| v.iter().map(|&b| (lo + hi - a - b).abs()).fold(f64::INFINITY, f64::min)).sum()
            };
            let feasible: Vec<&(Evaluation, Vec<f64>, Vec<f64>)> = options.iter().filter(|o| o.0.feasible).collect();
            let symmetric: Vec<&&(Evaluation, Vec<f64>, Vec<f64>)> =
                feasible.iter().filter(|o| asym(&o.1, x0, x1) + asym(&o.2, y0, y1) < tol).collect();
            let best = if feasible.is_empty() {
                None
            } else if !symmetric.is_empty() {
                symmetric.into_iter().min_by(|a, b| a.0.worst(&self.limits).total_cmp(&b.0.worst(&self.limits))).map(|o| (*o).clone())
            } else {
                let w_min = feasible.iter().map(|o| o.0.worst(&self.limits)).fold(f64::INFINITY, f64::min);
                feasible
                    .into_iter()
                    .filter(|o| o.0.worst(&self.limits) <= w_min * 1.02)
                    .min_by(|a, b| {
                        let sa = asym(&a.1, x0, x1) + asym(&a.2, y0, y1);
                        let sb = asym(&b.1, x0, x1) + asym(&b.2, y0, y1);
                        sa.total_cmp(&sb).then(a.0.worst(&self.limits).total_cmp(&b.0.worst(&self.limits)))
                    })
                    .cloned()
            };
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

/// Плотный Холецкий n×n (построчно); L в нижнем треугольнике.
fn cholesky(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut l = a.to_vec();
    for j in 0..n {
        let mut d = l[j * n + j];
        for k in 0..j {
            d -= l[j * n + k] * l[j * n + k];
        }
        if !(d > 0.0) {
            return Err("система свайного поля не положительно определена".into());
        }
        let d = d.sqrt();
        l[j * n + j] = d;
        for i in j + 1..n {
            let mut s = l[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = s / d;
        }
    }
    Ok(l)
}

fn chol_solve(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut y = b.to_vec();
    for i in 0..n {
        for k in 0..i {
            y[i] -= l[i * n + k] * y[k];
        }
        y[i] /= l[i * n + i];
    }
    for i in (0..n).rev() {
        for k in i + 1..n {
            y[i] -= l[k * n + i] * y[k];
        }
        y[i] /= l[i * n + i];
    }
    y
}

/// Подбор по трещинам только там, где сетки ПС1 (`env`) не проходят проверку раскрытия.
fn crack_designs(ctx: &SlsContext, forces: &[PlateForces], env: &[[f64; 4]]) -> Vec<(usize, capra_maury::Design)> {
    let fail: Vec<usize> = capra_maury::par_map(&(0..forces.len()).collect::<Vec<_>>(), |&e| {
        let a = env[e];
        let p = capra_maury::ProvidedMeshes { asx_top: a[0], asy_top: a[1], asx_bot: a[2], asy_bot: a[3], ds_x: ctx.ds, ds_y: ctx.ds };
        (!capra_maury::cracks_ok(ctx, &forces[e], &p, 10.0)).then_some(e)
    })
    .into_iter()
    .flatten()
    .collect();
    let f: Vec<PlateForces> = fail.iter().map(|&e| forces[e]).collect();
    fail.into_iter().zip(capra_maury::design_sls_batch(ctx, &f, 10.0)).collect()
}

fn add_some(d: &[(usize, capra_maury::Design)], env: &mut [[f64; 4]], impossible: &mut bool) {
    for (e, des) in d {
        *impossible |= !des.feasible();
        let a = [des.top.asx, des.top.asy, des.bottom.asx, des.bottom.asy];
        for k in 0..4 {
            env[*e][k] = env[*e][k].max(a[k]);
        }
    }
}
