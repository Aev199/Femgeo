//! Команды приложения: одна точка входа `Service::dispatch` для окна Tauri и отладочного
//! моста (`examples/ui_server.rs`). Аргументы и ответы — JSON.

use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::fem::{LoadCase, Model, Shell, ShellMaterial};
use crate::geobase;
use crate::optimize::{Evaluation, PileType, Problem, Punching, SettlementLimits};
use crate::piles::layout::{PileKind, Region, Rules};
use crate::piles::{Layer, SoilColumn};
use crate::sp63::capra_maury::{PlateSection, SlsContext, SlsOptions, UlsContext};
use crate::sp63::materials::{CalcType, ConcreteClass, ConcreteDiagramType, RebarClass};
use crate::sp63::{pile, section::Section};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RaftInput {
    pub outline: Vec<[f64; 2]>,
    #[serde(default)]
    pub holes: Vec<Vec<[f64; 2]>>,
    pub h: f64,
    /// Отметка низа плиты (голов свай), м.
    pub z_bottom: f64,
    /// Размер элемента сетки, м.
    pub mesh: f64,
    pub concrete: String,
    pub rebar: String,
    /// Расстояние от грани до центра арматуры, м.
    pub cover: f64,
    pub bar_d: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PointLoad {
    pub x: f64,
    pub y: f64,
    /// Расчётная вертикальная сила, кН (вниз — положительная).
    pub n: f64,
}

/// Загружение плиты (нормативные значения).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CaseInput {
    pub name: String,
    /// "permanent" | "long" | "short".
    pub kind: String,
    pub gamma_f: f64,
    /// Длительная доля кратковременной нагрузки (СП 20); для остальных не используется.
    #[serde(default)]
    pub psi_long: f64,
    /// Равномерная нагрузка, кПа (вниз — положительная).
    pub q: f64,
    /// Значения заданы расчётными (нормативные = расчётные / γf).
    #[serde(default)]
    pub design: bool,
    #[serde(default)]
    pub points: Vec<PointLoad>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PileInput {
    pub d: f64,
    pub length: f64,
    /// Несущая способность по грунту с учётом γk, кН.
    pub fd: f64,
    /// "bored" | "driven" | "end_bearing".
    pub kind: String,
    pub concrete: String,
    pub rebar: String,
    pub bars: usize,
    pub bar_d: f64,
    /// Расстояние от грани до центра стержней, м.
    pub cover: f64,
    /// Расчётная длина для случайного эксцентриситета, м.
    pub l0: f64,
    /// Наименьшее расстояние между осями и до края, м; None — по СП 24.
    #[serde(default)]
    pub min_spacing: Option<f64>,
    #[serde(default)]
    pub min_edge: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "source", rename_all = "snake_case")]
pub enum SoilInput {
    Geobase { path: String, model: Option<String> },
    Layers { top: f64, layers: Vec<Layer> },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Input {
    pub raft: RaftInput,
    pub cases: Vec<CaseInput>,
    pub pile: PileInput,
    /// Допустимая область свай; None — контур плиты.
    #[serde(default)]
    pub region: Option<Region>,
    pub soil: SoilInput,
    pub s_max: f64,
    pub ds_rel: f64,
    /// "rows" — рядами (по умолчанию), "free" — свободная.
    #[serde(default = "rows")]
    pub layout: String,
    /// Предельная ширина раскрытия трещин, мм (п. 8.2.6 СП 63): продолжительное и
    /// непродолжительное; 0 — без проверки.
    #[serde(default = "acrc_long")]
    pub acrc_long: f64,
    #[serde(default = "acrc_short")]
    pub acrc_short: f64,
    /// Проверять продавливание плиты сваями.
    #[serde(default = "yes")]
    pub punching: bool,
    /// Существующие сваи (X, Y) — ориентир.
    #[serde(default)]
    pub existing: Vec<[f64; 2]>,
}

fn acrc_long() -> f64 {
    0.3
}
fn acrc_short() -> f64 {
    0.4
}
fn yes() -> bool {
    true
}

fn rows() -> String {
    "rows".into()
}

impl serde::Serialize for Region {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        json!({ "outer": self.outer, "holes": self.holes }).serialize(s)
    }
}

impl<'de> serde::Deserialize<'de> for Region {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct R {
            outer: Vec<[f64; 2]>,
            #[serde(default)]
            holes: Vec<Vec<[f64; 2]>>,
        }
        let r = R::deserialize(d)?;
        Ok(Region { outer: r.outer, holes: r.holes })
    }
}

fn concrete(name: &str) -> Result<ConcreteClass, String> {
    use ConcreteClass::*;
    Ok(match name.trim() {
        "B15" => B15,
        "B20" => B20,
        "B25" => B25,
        "B30" => B30,
        "B35" => B35,
        "B40" => B40,
        "B45" => B45,
        "B50" => B50,
        "B55" => B55,
        "B60" => B60,
        n => return Err(format!("неизвестный класс бетона «{n}»")),
    })
}

fn rebar(name: &str) -> Result<RebarClass, String> {
    Ok(match name.trim() {
        "A240" => RebarClass::A240,
        "A400" => RebarClass::A400,
        "A500" => RebarClass::A500,
        n => return Err(format!("неизвестный класс арматуры «{n}»")),
    })
}

impl Serialize for Layer {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        json!({ "bottom": self.bottom, "e": self.e, "nu": self.nu }).serialize(s)
    }
}

impl<'de> Deserialize<'de> for Layer {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct L {
            bottom: f64,
            e: f64,
            nu: f64,
        }
        let l = L::deserialize(d)?;
        Ok(Layer { bottom: l.bottom, e: l.e, nu: l.nu })
    }
}

/// Сетка плиты: квадратные элементы по габариту, в плиту входят элементы с центром внутри
/// контура и вне отверстий.
fn raft_model(r: &RaftInput, e: f64) -> Result<Model, String> {
    if r.outline.len() < 3 {
        return Err("контур плиты: нужно не меньше трёх точек".into());
    }
    if !(r.mesh > 0.0 && r.h > 0.0) {
        return Err("толщина плиты и размер элемента должны быть положительными".into());
    }
    let region = Region { outer: r.outline.clone(), holes: r.holes.clone() };
    let (x0, x1) = r.outline.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, p| (a.0.min(p[0]), a.1.max(p[0])));
    let (y0, y1) = r.outline.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, p| (a.0.min(p[1]), a.1.max(p[1])));
    let nx = ((x1 - x0) / r.mesh).round().max(1.0) as usize;
    let ny = ((y1 - y0) / r.mesh).round().max(1.0) as usize;
    if (nx + 1) * (ny + 1) > 200_000 {
        return Err("слишком мелкая сетка плиты".into());
    }
    let (dx, dy) = ((x1 - x0) / nx as f64, (y1 - y0) / ny as f64);
    let mut m = Model::default();
    let mut node = vec![usize::MAX; (nx + 1) * (ny + 1)];
    let mat = ShellMaterial { e, nu: 0.2, h: r.h };
    for i in 0..ny {
        for j in 0..nx {
            let c = [x0 + (j as f64 + 0.5) * dx, y0 + (i as f64 + 0.5) * dy];
            if !region.contains(c) {
                continue;
            }
            let mut ids = vec![];
            for (a, b) in [(i, j), (i, j + 1), (i + 1, j + 1), (i + 1, j)] {
                let k = a * (nx + 1) + b;
                if node[k] == usize::MAX {
                    node[k] = m.nodes.len();
                    m.nodes.push([x0 + b as f64 * dx, y0 + a as f64 * dy, 0.0]);
                }
                ids.push(node[k]);
            }
            m.shells.push(Shell { nodes: ids, mat });
        }
    }
    if m.shells.is_empty() {
        return Err("в контур плиты не попал ни один элемент: уменьшите размер элемента".into());
    }
    plane_fix(&mut m);
    Ok(m)
}

/// Плита в своей плоскости: поворот вокруг вертикали во всех узлах, сдвиг и поворот в плане
/// — в двух удалённых узлах (сваи держат только вертикально).
fn plane_fix(m: &mut Model) {
    for k in 0..m.nodes.len() {
        m.fixed.push((k, [false, false, false, false, false, true]));
    }
    // Плоскость плиты: закрепить сдвиг и поворот в плане (сваи держат только вертикально).
    let a = 0;
    let b = (1..m.nodes.len()).max_by(|&p, &q| dist2(m.nodes[a], m.nodes[p]).total_cmp(&dist2(m.nodes[a], m.nodes[q]))).unwrap_or(0);
    m.fixed.push((a, [true, true, false, false, false, false]));
    let ab = [m.nodes[b][0] - m.nodes[a][0], m.nodes[b][1] - m.nodes[a][1]];
    m.fixed.push((b, if ab[0].abs() >= ab[1].abs() { [false, true, false, false, false, false] } else { [true, false, false, false, false, false] }));
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)
}

#[derive(Default)]
pub struct Service {
    geobase: Option<(String, geobase::Object)>,
    last: Option<(Input, Model, Vec<Evaluation>)>,
    /// Первый вариант последнего расчёта — существующая расстановка.
    last_reference: bool,
}

impl Service {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn dispatch(&mut self, command: &str, args: Value, progress: &mut dyn FnMut(&str)) -> Result<Value, String> {
        match command {
            "defaults" => Ok(serde_json::to_value(default_input()).unwrap()),
            "open_geobase" => self.open_geobase(&args),
            "run" => self.run(args, progress),
            "export_csv" => self.export(&args, false),
            "export_dxf" => self.export(&args, true),
            "export_rebar_dxf" => self.export_rebar(&args),
            "export_report" => self.export_report(&args),
            "save_project" => {
                let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
                let input = args.get("input").ok_or("нет исходных данных")?;
                let text = serde_json::to_string_pretty(&json!({ "femgeo": 1, "input": input })).unwrap();
                std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
                Ok(json!({ "path": path }))
            }
            "load_project" => {
                let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
                let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
                let v: Value = serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
                let input = v.get("input").cloned().ok_or("файл не является проектом Femgeo")?;
                // Проверка формата: те же правила, что при расчёте.
                serde_json::from_value::<Input>(input.clone()).map_err(|e| format!("проект: {e}"))?;
                Ok(input)
            }
            _ => Err(format!("неизвестная команда «{command}»")),
        }
    }

    fn open_geobase(&mut self, args: &Value) -> Result<Value, String> {
        let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к объекту geobase")?;
        let model = args.get("model").and_then(Value::as_str);
        let o = geobase::Object::open(Path::new(path), model)?;
        let out = json!({
            "path": path,
            "model": o.model,
            "soils": o.soils.len(),
            "boreholes": o.boreholes.iter().map(|b| json!({ "name": b.name, "x": b.x, "y": b.y, "z": b.z })).collect::<Vec<_>>(),
        });
        self.geobase = Some((path.to_string(), o));
        Ok(out)
    }

    fn soil_fn(&mut self, s: &SoilInput) -> Result<Box<dyn Fn(f64, f64) -> SoilColumn + Sync>, String> {
        match s {
            SoilInput::Layers { top, layers } => {
                if layers.is_empty() {
                    return Err("грунт: нет слоёв".into());
                }
                let c = SoilColumn { top: *top, layers: layers.clone() };
                Ok(Box::new(move |_, _| c.clone()))
            }
            SoilInput::Geobase { path, model } => {
                let reopen = !matches!(&self.geobase, Some((p, o)) if p == path && model.as_deref().is_none_or(|m| m == o.model));
                if reopen {
                    let o = geobase::Object::open(Path::new(path), model.as_deref())?;
                    self.geobase = Some((path.clone(), o));
                }
                let o = self.geobase.as_ref().unwrap().1.clone();
                // Колонки скважин — заранее, ошибки ИГЭ — сразу.
                let cols: Vec<([f64; 2], SoilColumn)> =
                    o.boreholes.iter().map(|b| o.column(b).map(|c| ([b.x, b.y], c))).collect::<Result<_, _>>()?;
                if cols.is_empty() {
                    return Err("geobase: в объекте нет скважин".into());
                }
                Ok(Box::new(move |x, y| {
                    cols.iter().min_by(|a, b| ((a.0[0] - x).powi(2) + (a.0[1] - y).powi(2)).total_cmp(&((b.0[0] - x).powi(2) + (b.0[1] - y).powi(2)))).unwrap().1.clone()
                }))
            }
        }
    }

    fn run(&mut self, args: Value, progress: &mut dyn FnMut(&str)) -> Result<Value, String> {
        let input: Input = serde_json::from_value(args).map_err(|e| format!("исходные данные: {e}"))?;
        self.run_with(input, None, progress)
    }

    /// Расчёт по сетке здания (препроцессор topo): плита, её контур, толщина и отметка, нагрузки
    /// и существующие сваи — из `mesh`; из `input` — материалы, сваи, грунт, ограничения и вид
    /// загружений (`input.cases` по порядку `mesh.cases`; q и сосредоточенные силы не
    /// используются, признак «расчётные» относится к нагрузкам сетки).
    pub fn run_mesh(&mut self, input: Input, mesh: &crate::building::Mesh, progress: &mut dyn FnMut(&str)) -> Result<Value, String> {
        self.run_with(input, Some(mesh), progress)
    }

    fn run_with(&mut self, mut input: Input, mesh: Option<&crate::building::Mesh>, progress: &mut dyn FnMut(&str)) -> Result<Value, String> {
        let split = match mesh {
            Some(m) => {
                if m.cases.len() != input.cases.len() {
                    return Err(format!("загружений в сетке {}, в исходных данных {}", m.cases.len(), input.cases.len()));
                }
                let s = crate::building::split(m)?;
                input.raft.outline = s.region.outer.clone();
                input.raft.holes = s.region.holes.clone();
                input.raft.h = s.h;
                input.raft.z_bottom = s.z_mid - s.h / 2.0;
                input.existing.extend(s.piles.iter().map(|p| p.0));
                Some(s)
            }
            None => None,
        };
        let r = &input.raft;
        let cc = concrete(&r.concrete)?;
        let (c_uls, r_uls) = (cc.chars(CalcType::Uls), rebar(&r.rebar)?.chars(CalcType::Uls));
        // Плита — всегда регулярная сетка Femgeo (сваи рядами в её узлах). Для сетки здания узлы
        // её плиты (сопряжение, нагрузки) переносятся в ближайшие узлы регулярной сетки.
        let raft = raft_model(r, c_uls.eb)?;
        let to_grid: Vec<usize> = split.as_ref().map_or(vec![], |s| {
            s.raft.nodes.iter().map(|p| (0..raft.nodes.len()).min_by(|&a, &b| dist2(raft.nodes[a], *p).total_cmp(&dist2(raft.nodes[b], *p))).unwrap()).collect()
        });
        progress("сетка плиты");

        // Загружения: давление на все элементы, сосредоточенные силы — в ближайший узел.
        if input.cases.is_empty() {
            return Err("нет загружений".into());
        }
        let mut cases = Vec::new();
        let mut kinds = Vec::new();
        let mut norms = Vec::new();
        for (ci, c) in input.cases.iter().enumerate() {
            let kind = match c.kind.as_str() {
                "permanent" => crate::loads::Kind::Permanent,
                "long" => crate::loads::Kind::Long,
                "short" => crate::loads::Kind::Short,
                k => return Err(format!("загружение «{}»: неизвестный вид «{k}»", c.name)),
            };
            if !(c.gamma_f >= 1.0) {
                return Err(format!("загружение «{}»: γf должен быть не меньше 1", c.name));
            }
            kinds.push(crate::loads::Case { kind, gamma_f: c.gamma_f, psi_long: c.psi_long.clamp(0.0, 1.0) });
            // Загружения хранятся нормативными; расчётные значения приводятся делением на γf.
            let to_norm = if c.design { 1.0 / c.gamma_f } else { 1.0 };
            norms.push(to_norm);
            if let Some(s) = &split {
                let nodal = s.raft_cases[ci].nodal.iter().map(|&(n, f)| (to_grid[n], f.map(|x| x * to_norm))).collect();
                cases.push(LoadCase { nodal, pressure: vec![] });
                continue;
            }
            let mut lc = LoadCase { pressure: (0..raft.shells.len()).map(|s| (s, -c.q * to_norm)).collect(), ..Default::default() };
            for p in &c.points {
                let k = (0..raft.nodes.len()).min_by(|&a, &b| dist2(raft.nodes[a], [p.x, p.y, 0.0]).total_cmp(&dist2(raft.nodes[b], [p.x, p.y, 0.0]))).unwrap();
                lc.nodal.push((k, [0.0, 0.0, -p.n * to_norm, 0.0, 0.0, 0.0]));
            }
            cases.push(lc);
        }
        let combos = crate::loads::combinations(&kinds);
        let building = split.as_ref().and_then(|s| s.building.as_ref()).map(|b| crate::fem::condense::Building {
            block: crate::fem::DenseBlock { dofs: b.block.dofs.iter().map(|&d| 6 * to_grid[d / 6] + d % 6).collect(), k: b.block.k.clone() },
            loads: b.loads.iter().zip(&norms).map(|(l, &k)| l.iter().map(|&(n, f)| (to_grid[n], f.map(|x| x * k))).collect()).collect(),
        });
        let area: f64 = polygon_area(&r.outline) - r.holes.iter().map(|h| polygon_area(h)).sum::<f64>();
        // Нормативная вертикальная нагрузка загружений (вниз — положительная).
        let case_total: Vec<f64> = match &building {
            None if split.is_none() => input.cases.iter().map(|k| (k.q * area + k.points.iter().map(|p| p.n).sum::<f64>()) / if k.design { k.gamma_f } else { 1.0 }).collect(),
            _ => (0..cases.len())
                .map(|i| -cases[i].nodal.iter().chain(building.as_ref().map_or(&[][..], |b| &b.loads[i][..])).map(|l| l.1[2]).sum::<f64>())
                .collect(),
        };

        let pi = &input.pile;
        let kind = match pi.kind.as_str() {
            "driven" => PileKind::DrivenFriction,
            "end_bearing" => PileKind::EndBearing,
            _ => PileKind::Bored,
        };
        let mut rules = Rules::sp24(kind, pi.d);
        if let Some(s) = pi.min_spacing {
            rules.min_spacing = s;
        }
        if let Some(e) = pi.min_edge {
            rules.min_edge = e;
        }
        let pc = concrete(&pi.concrete)?.chars(CalcType::Uls);
        let pr = rebar(&pi.rebar)?.chars(CalcType::Uls);
        let sec = Section::circular(pi.d / 2.0, 96, 20, pi.bars, pi.bar_d, pi.d / 2.0 - pi.cover, 0.0, pc.diagram(ConcreteDiagramType::L2), pr.diagram(), pc.eps_b2);
        let (n, mx, my) = pile::with_accidental_eccentricity(-1000.0, 0.0, 0.0, pi.l0, pi.d, Some((1.0, 0.0)))?;
        let lim = pile::all_factor(&sec, n, mx, my, &Default::default());
        if lim.status != pile::Status::Ok {
            return Err("сечение сваи: предельное усилие не найдено".into());
        }
        progress("сечение сваи");

        let psec = PlateSection { h: r.h, a_top: r.cover, a_bot: r.cover, ds_x: r.bar_d, ds_y: r.bar_d, concrete_diagram: ConcreteDiagramType::L3 };
        let design = UlsContext::new(&psec, &c_uls, &r_uls);
        // Трещины (п. 8.2.6 СП 63), φ1 = 1,4: продолжительное — длительное сочетание;
        // непродолжительное — полные сочетания (верхняя оценка acrc1 + acrc2 − acrc3).
        let crack = |lim: f64| {
            (lim > 0.0).then(|| {
                SlsContext::new(
                    &psec,
                    &cc.chars(CalcType::Sls),
                    &rebar(&r.rebar).expect("проверено").chars(CalcType::Sls),
                    Some((&c_uls, &r_uls)),
                    SlsOptions { acrc_lim_mm: lim, phi1: 1.4, ..Default::default() },
                )
            })
        };
        let (crack_long, crack_short) = (crack(input.acrc_long), crack(input.acrc_short));
        let soil = self.soil_fn(&input.soil)?;
        let region = input.region.clone().unwrap_or(Region { outer: r.outline.clone(), holes: r.holes.clone() });
        let problem = Problem {
            raft: raft.clone(),
            pile_head: r.z_bottom,
            pile: PileType {
                length: pi.length,
                d: pi.d,
                ea: pc.eb * std::f64::consts::PI * pi.d * pi.d / 4.0,
                fd_ground: pi.fd,
                material_util_per_kn: lim.utilization / 1000.0,
            },
            region,
            rules,
            soil: soil.as_ref(),
            cases,
            uls: combos.uls.clone(),
            sls: combos.sls.clone(),
            sls_long: Some(combos.sls_long.clone()),
            limits: SettlementLimits { s_max: input.s_max, ds_rel_max: input.ds_rel },
            design: &design,
            rebar_axis: [1.0, 0.0, 0.0],
            raft_h: r.h,
            building: building.as_ref(),
            crack_short: crack_short.as_ref(),
            crack_long: crack_long.as_ref(),
            punching: input.punching.then(|| Punching { rbt: c_uls.rbt, h0: r.h - r.cover, raft: Region { outer: r.outline.clone(), holes: r.holes.clone() } }),
            existing: input.existing.clone(),
            cache: Default::default(),
        };
        // Начальное число свай — по суммарной нагрузке и Fd с запасом 1,5.
        // Наибольшая суммарная расчётная нагрузка по сочетаниям ПС1.
        let total = combos
            .uls
            .iter()
            .map(|c| c.iter().zip(&case_total).map(|(f, t)| f * t).sum::<f64>())
            .fold(0.0, f64::max);
        let n_start = ((1.5 * total / pi.fd).ceil() as usize).max(4);
        progress("подбор расстановки");
        let mut report = |e: &crate::optimize::Evaluation| progress(&format!("вариант: {} свай, {:.1} кг/м³", e.piles.len(), e.kg_m3));
        let mut variants = if input.layout == "free" { problem.frontier(n_start, 8, &mut report)? } else { problem.regular(&mut report)? };
        // Существующая расстановка — первой строкой, для сравнения.
        let reference = problem.reference()?;
        let has_reference = reference.is_some();
        if let Some(r) = reference {
            variants.insert(0, r);
        }
        let mut out = result_json(&raft, &variants, &input, has_reference);
        if let Some(s) = &split {
            out["notes"] = json!(s.notes);
        }
        self.last = Some((input, raft, variants));
        self.last_reference = has_reference;
        Ok(out)
    }

    fn export(&self, args: &Value, dxf: bool) -> Result<Value, String> {
        let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
        let i = args.get("variant").and_then(Value::as_u64).ok_or("не выбран вариант")? as usize;
        let (input, raft, vs) = self.last.as_ref().ok_or("нет результата расчёта")?;
        let v = vs.get(i).ok_or("нет такого варианта")?;
        let text = if dxf { dxf_text(input, raft, v) } else { csv_text(raft, v) };
        std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
        Ok(json!({ "path": path }))
    }
}

impl Service {
    /// Поля требуемой арматуры плиты выбранного варианта: DXF, слои верх/низ × X/Y, см²/м.
    fn export_rebar(&self, args: &Value) -> Result<Value, String> {
        let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
        let i = args.get("variant").and_then(Value::as_u64).ok_or("не выбран вариант")? as usize;
        let (_, raft, vs) = self.last.as_ref().ok_or("нет результата расчёта")?;
        let v = vs.get(i).ok_or("нет такого варианта")?;
        let f = |k: usize| -> Vec<f64> { v.rebar.iter().map(|a| a[k] * 1e4).collect() };
        let (tx, ty, bx, by) = (f(0), f(1), f(2), f(3));
        let text = crate::reinforce::dxf_fields(raft, &[("AS_TOP_X", &tx), ("AS_TOP_Y", &ty), ("AS_BOT_X", &bx), ("AS_BOT_Y", &by)], "см²/м");
        std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
        Ok(json!({ "path": path }))
    }
}

impl Service {
    /// Отчёт по выбранному варианту: автономный HTML.
    fn export_report(&self, args: &Value) -> Result<Value, String> {
        let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
        let i = args.get("variant").and_then(Value::as_u64).ok_or("не выбран вариант")? as usize;
        let (input, raft, vs) = self.last.as_ref().ok_or("нет результата расчёта")?;
        if i >= vs.len() {
            return Err("нет такого варианта".into());
        }
        let text = crate::report::html(input, raft, vs, i, self.last_reference);
        std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
        Ok(json!({ "path": path }))
    }
}

fn polygon_area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    (0..n).map(|i| p[i][0] * p[(i + 1) % n][1] - p[(i + 1) % n][0] * p[i][1]).sum::<f64>().abs() / 2.0
}

fn result_json(raft: &Model, vs: &[Evaluation], input: &Input, has_reference: bool) -> Value {
    json!({
        "nodes": raft.nodes.iter().map(|p| [p[0], p[1]]).collect::<Vec<_>>(),
        "elements": raft.shells.iter().map(|s| s.nodes.clone()).collect::<Vec<_>>(),
        "outline": input.raft.outline,
        "holes": input.raft.holes,
        "fd": input.pile.fd,
        "d": input.pile.d,
        "existing": input.existing,
        "variants": vs.iter().enumerate().map(|(i, v)| json!({
            "reference": has_reference && i == 0,
            "util_punch": v.util_punch,
            "shift": if v.shift.is_finite() { json!(v.shift) } else { Value::Null },
            "count": v.piles.len(),
            "kg_m3": v.kg_m3,
            "util_ground": v.util_ground,
            "util_material": v.util_material,
            "s_max": v.s_max,
            "ds_rel": v.ds_rel,
            "feasible": v.feasible,
            "design_impossible": v.design_impossible,
            "piles": v.piles.iter().zip(&v.pile_n).map(|(&n, &f)| json!({ "x": raft.nodes[n][0], "y": raft.nodes[n][1], "n": f })).collect::<Vec<_>>(),
            // см²/м
            "rebar": v.rebar.iter().map(|a| a.map(|x| x * 1e4)).collect::<Vec<_>>(),
            "w": v.w_nodes,
        })).collect::<Vec<_>>(),
    })
}

fn csv_text(raft: &Model, v: &Evaluation) -> String {
    let mut s = String::from("№;X, м;Y, м;N, кН\n");
    for (k, (&n, &f)) in v.piles.iter().zip(&v.pile_n).enumerate() {
        s += &format!("{};{:.3};{:.3};{:.1}\n", k + 1, raft.nodes[n][0], raft.nodes[n][1], f).replace('.', ",");
    }
    s
}

/// DXF R12: контур плиты и отверстия (слой RAFT), сваи — окружности диаметра сваи (слой PILES).
fn dxf_text(input: &Input, raft: &Model, v: &Evaluation) -> String {
    let mut s = String::from("0\nSECTION\n2\nENTITIES\n");
    let poly = |pts: &[[f64; 2]], s: &mut String| {
        for i in 0..pts.len() {
            let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
            *s += &format!("0\nLINE\n8\nRAFT\n10\n{}\n20\n{}\n30\n0\n11\n{}\n21\n{}\n31\n0\n", a[0], a[1], b[0], b[1]);
        }
    };
    poly(&input.raft.outline, &mut s);
    for h in &input.raft.holes {
        poly(h, &mut s);
    }
    for &n in &v.piles {
        s += &format!("0\nCIRCLE\n8\nPILES\n10\n{}\n20\n{}\n30\n0\n40\n{}\n", raft.nodes[n][0], raft.nodes[n][1], input.pile.d / 2.0);
    }
    s + "0\nENDSEC\n0\nEOF\n"
}

/// Пример исходных данных для нового окна.
pub fn default_input() -> Input {
    Input {
        raft: RaftInput {
            outline: vec![[0.0, 0.0], [24.0, 0.0], [24.0, 18.0], [0.0, 18.0]],
            holes: vec![],
            h: 1.0,
            z_bottom: -3.0,
            mesh: 1.0,
            concrete: "B30".into(),
            rebar: "A500".into(),
            cover: 0.06,
            bar_d: 0.02,
        },
        cases: vec![
            CaseInput { name: "Собственный вес и постоянные".into(), kind: "permanent".into(), gamma_f: 1.1, psi_long: 0.0, q: 130.0, design: false, points: vec![] },
            CaseInput { name: "Полезная".into(), kind: "short".into(), gamma_f: 1.2, psi_long: 0.35, q: 30.0, design: false, points: vec![] },
        ],
        pile: PileInput {
            d: 0.6,
            length: 16.0,
            fd: 2500.0,
            kind: "bored".into(),
            concrete: "B25".into(),
            rebar: "A500".into(),
            bars: 8,
            bar_d: 0.016,
            cover: 0.06,
            l0: 3.0,
            min_spacing: None,
            min_edge: None,
        },
        region: None,
        soil: SoilInput::Layers {
            top: 0.0,
            layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3 }, Layer { bottom: -60.0, e: 40_000.0, nu: 0.3 }],
        },
        s_max: 0.15,
        ds_rel: 0.002,
        layout: rows(),
        acrc_long: acrc_long(),
        acrc_short: acrc_short(),
        punching: true,
        existing: vec![],
    }
}
