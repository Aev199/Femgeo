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
    /// Предел сваи Fu, кН: сваи с N > Fu держат Fu (None — без предела).
    #[serde(default)]
    pub pile_fu: Option<f64>,
    /// Жёсткость свай: "sp24" — взаимовлияние по СП 24, п. 7.4 (по умолчанию);
    /// "conditional" — одинаковая, по осадке условного фундамента (нужны γ и φ слоёв).
    #[serde(default = "sp24")]
    pub pile_model: String,
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

fn sp24() -> String {
    "sp24".into()
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
        json!({ "bottom": self.bottom, "e": self.e, "nu": self.nu, "gamma": self.gamma, "phi": self.phi }).serialize(s)
    }
}

impl<'de> Deserialize<'de> for Layer {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct L {
            bottom: f64,
            e: f64,
            nu: f64,
            #[serde(default)]
            gamma: f64,
            #[serde(default)]
            phi: f64,
        }
        let l = L::deserialize(d)?;
        Ok(Layer { bottom: l.bottom, e: l.e, nu: l.nu, gamma: l.gamma, phi: l.phi })
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
    /// Последний расчёт по сетке здания — для полной модели (армирование подземной части).
    last_mesh: Option<LastMesh>,
    /// Армирование подземной части последнего варианта: полная модель и поля, м²/м.
    last_underground: Option<(Model, Vec<[f64; 4]>)>,
}

struct LastMesh {
    split: crate::building::Split,
    to_grid: Vec<usize>,
    /// Нормативные нагрузки: плита (узлы сетки) и здание (узлы модели здания).
    raft_cases: Vec<LoadCase>,
    building_cases: Vec<LoadCase>,
    combos: crate::loads::Combinations,
    /// Отметка поверхности грунта в центре плиты: ниже — подземная часть.
    ground: f64,
    /// Здание, приведённое к плите (узлы сетки, нормативные нагрузки).
    building: Option<crate::fem::condense::Building>,
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
            "design_underground" => self.design_underground(&args, progress),
            "export_underground_dxf" => {
                let path = args.get("path").and_then(Value::as_str).ok_or("нет пути к файлу")?;
                let (m, f) = self.last_underground.as_ref().ok_or("нет армирования подземной части")?;
                let c = |k: usize| -> Vec<f64> { f.iter().map(|a| a[k] * 1e4).collect() };
                let (tx, ty, bx, by) = (c(0), c(1), c(2), c(3));
                let text = crate::reinforce::dxf_fields(m, &[("AS_TOP_X", &tx), ("AS_TOP_Y", &ty), ("AS_BOT_X", &bx), ("AS_BOT_Y", &by)], "см²/м");
                std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
                Ok(json!({ "path": path }))
            }
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
        // Жёсткость свай по условному фундаменту: контур — контур плиты (сваи у края; допущение),
        // грунт — в центре контура.
        let conditional = match input.pile_model.as_str() {
            "conditional" => {
                let n = r.outline.len() as f64;
                let c = r.outline.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
                Some(crate::optimize::ConditionalModel { outline: r.outline.clone(), column: soil(c[0], c[1]) })
            }
            "sp24" => None,
            m => return Err(format!("неизвестная модель свай «{m}»")),
        };
        // Продавливание вертикальными конструкциями (сетка здания): загруженные площади и контуры
        // на h0/2 (СП 63, п. 8.1.46); торцы и углы стен — СП 430.1325800.2018, п. 6.3.8, рис. 6.1–6.2.
        let h0 = r.h - r.cover;
        let structure_punch: Vec<crate::optimize::PunchArea> = match (&split, input.punching) {
            (Some(s), true) => s
                .structures
                .iter()
                .flat_map(|st| {
                    let p = s.raft.nodes[st.node];
                    let c = [p[0], p[1]];
                    let mut nodes: Vec<usize> = st.nodes.iter().map(|&k| to_grid[k]).collect();
                    nodes.sort_unstable();
                    nodes.dedup();
                    let rect = |x0: f64, x1: f64, y0: f64, y1: f64| vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
                    let along = |poly: Vec<[f64; 2]>, k: usize| -> Vec<[f64; 2]> {
                        (0..poly.len()).flat_map(|i| {
                            let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
                            (0..k).map(move |j| { let t = (j as f64 + 0.5) / k as f64; [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])] })
                        }).collect()
                    };
                    // Расчётные контуры: (вид, длина, область внутри, отрезки (начало, конец, наружная нормаль)
                    // для силы по поперечным силам плиты — торцы и углы стен).
                    type Seg = ([f64; 2], [f64; 2], [f64; 2]);
                    let checks: Vec<(&str, f64, Vec<[f64; 2]>, Vec<Seg>, Vec<[f64; 2]>)> = match st.kind {
                        crate::building::StructureKind::Column { b, h } => {
                            let (hx, hy) = ((b + h0) / 2.0, (h + h0) / 2.0);
                            vec![("column", 2.0 * (b + h) + 4.0 * h0, rect(c[0] - hx, c[0] + hx, c[1] - hy, c[1] + hy), vec![], vec![[1.0, 0.0], [0.0, 1.0]])]
                        }
                        crate::building::StructureKind::WallEnd { t, dir } => {
                            // СП 430, п. 6.3.8: участок у торца на h0/2 и два вдоль стены той же длины t + h0.
                            let l = t + h0;
                            let nrm = [-dir[1], dir[0]];
                            let w = l / 2.0;
                            let q = |s: f64, n: f64| [c[0] + s * dir[0] + n * nrm[0], c[1] + s * dir[1] + n * nrm[1]];
                            let (s0, s1) = (-h0 / 2.0, -h0 / 2.0 + l);
                            let segs = vec![(q(s0, -w), q(s0, w), [-dir[0], -dir[1]]), (q(s0, w), q(s1, w), nrm), (q(s0, -w), q(s1, -w), [-nrm[0], -nrm[1]])];
                            let mut v = vec![("wall_end", 3.0 * l, vec![q(s0, -w), q(s1, -w), q(s1, w), q(s0, w)], segs, vec![nrm])];
                            // Торец у края плиты (СП 430, п. 6.3.8): дополнительно — контур из двух участков
                            // вдоль граней стены на h0/2 от края плиты.
                            let region = Region { outer: r.outline.clone(), holes: r.holes.clone() };
                            let edge = (1..=40).map(|k| -(k as f64) * h0 / 40.0).find(|&sv| !region.contains(q(sv, 0.0)));
                            if let Some(se) = edge {
                                // Рис. 6.1, б (контур 1′): от края плиты вдоль граней стены до a1 = t за торцом
                                // (s отсчитывается от торца внутрь стены).
                                let se = se + h0 / 80.0;
                                let s2 = t;
                                let segs2 = vec![(q(se, w), q(s2, w), nrm), (q(se, -w), q(s2, -w), [-nrm[0], -nrm[1]])];
                                // Контур 1′ — для стены, подходящей торцом к краю; у стены вдоль края (участок
                                // грани вне плиты) он не применяется, остаётся контур 1, открытый к краю.
                                let inside = |a: [f64; 2], b: [f64; 2]| (0..10).filter(|&k| { let t = (k as f64 + 0.5) / 10.0; region.contains([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]) }).count() >= 5;
                                if segs2.iter().all(|sg| inside(sg.0, sg.1)) {
                                v.push(("wall_end_edge", 2.0 * (s2 - se), vec![q(se, -w), q(s2, -w), q(s2, w), q(se, w)], segs2, vec![nrm]));
                                }
                            }
                            v
                        }
                        crate::building::StructureKind::WallCorner { t, d1, d2 } => {
                            // СП 430, рис. 6.2 (a1 = a2 = t): наружный контур — два участка a1 + a2 на h0/2 от
                            // наружных граней; внутренний — диагональ во входящем углу.
                            let pt = |o: [f64; 2], k1: f64, k2: f64| [o[0] + k1 * d1[0] + k2 * d2[0], o[1] + k1 * d1[1] + k2 * d2[1]];
                            let e = t / 2.0 + h0 / 2.0;
                            let o = pt(c, -e, -e);
                            let outer = vec![(o, pt(o, 2.0 * t, 0.0), [-d2[0], -d2[1]]), (o, pt(o, 0.0, 2.0 * t), [-d1[0], -d1[1]])];
                            let i = pt(c, t / 2.0, t / 2.0);
                            let (p1, p2) = (pt(i, h0 / 2.0, h0 / 2.0 + t), pt(i, h0 / 2.0 + t, h0 / 2.0));
                            let dn = [d1[0] + d2[0], d1[1] + d2[1]];
                            let ln = (dn[0] * dn[0] + dn[1] * dn[1]).sqrt().max(1e-12);
                            let inner = vec![(p1, p2, [dn[0] / ln, dn[1] / ln])];
                            let area = vec![o, pt(o, 2.0 * t, 0.0), p2, p1, pt(o, 0.0, 2.0 * t)];
                            vec![("wall_corner_outer", 4.0 * t, area.clone(), outer, vec![]), ("wall_corner_inner", (2.0f64).sqrt() * t, area, inner, vec![])]
                        }
                    };
                    checks
                        .into_iter()
                        .map(|(kind, u, inside, segs, moment_axes)| {
                            let samples = along(inside.clone(), 18);
                            // Точки контура в элементах плиты (вне плиты — не учитываются).
                            let shear = (!segs.is_empty()).then(|| {
                                let k = 12;
                                segs.iter()
                                    .flat_map(|&(a, b, nrm)| {
                                        let ds = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt() / k as f64;
                                        (0..k).map(move |j| { let t = (j as f64 + 0.5) / k as f64; ([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])], nrm, ds) })
                                    })
                                    .filter_map(|(p, nrm, ds)| element_at(&raft, p).map(|e| (e, nrm, ds, p)))
                                    .collect()
                            });
                            crate::optimize::PunchArea { kind: kind.into(), at: c, nodes: nodes.clone(), u, samples, inside: Region { outer: inside, holes: vec![] }, shear, moment_axes }
                        })
                        .collect::<Vec<_>>()
                })
                .collect(),
            _ => vec![],
        };
        let raft_cases_norm = cases.clone();
        let ground = {
            let n = r.outline.len() as f64;
            let c = r.outline.iter().fold([0.0, 0.0], |a, p| [a[0] + p[0] / n, a[1] + p[1] / n]);
            soil(c[0], c[1]).top
        };
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
            punching: input.punching.then(|| Punching { rbt: c_uls.rbt, h0: r.h - r.cover, rsw: RSW_PUNCH, raft: Region { outer: r.outline.clone(), holes: r.holes.clone() } }),
            existing: input.existing.clone(),
            pile_fu: input.pile_fu,
            conditional,
            structure_punch,
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
            // Вертикальные конструкции: узлы сопряжения здания с плитой (X, Y).
            let mut st: Vec<[f64; 2]> = s.building.as_ref().map_or(vec![], |b| b.block.dofs.iter().step_by(3).map(|d| { let p = s.raft.nodes[d / 6]; [p[0], p[1]] }).collect());
            st.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
            st.dedup();
            out["structures"] = json!(st);
        }
        self.last_mesh = split.map(|s| {
            let building_cases = s.building_cases.iter().zip(&norms).map(|(lc, &k)| LoadCase { nodal: lc.nodal.iter().map(|&(n, f)| (n, f.map(|x| x * k))).collect(), pressure: vec![] }).collect();
            LastMesh { split: s, to_grid, raft_cases: raft_cases_norm, building_cases, combos, ground, building: building.clone() }
        });
        self.last_underground = None;
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

/// Элемент плиты, в плане содержащий точку (регулярная сетка: по габариту элемента).
fn element_at(m: &Model, p: [f64; 2]) -> Option<usize> {
    m.shells.iter().position(|s| {
        let (x0, x1) = s.nodes.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &n| (a.0.min(m.nodes[n][0]), a.1.max(m.nodes[n][0])));
        let (y0, y1) = s.nodes.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &n| (a.0.min(m.nodes[n][1]), a.1.max(m.nodes[n][1])));
        p[0] >= x0 && p[0] <= x1 && p[1] >= y0 && p[1] <= y1
    })
}

/// Rsw поперечной арматуры продавливания, кПа: А500 — 300 МПа (СП 63, табл. 6.15; по памяти,
/// сверить), класс арматуры не задаётся.
const RSW_PUNCH: f64 = 300_000.0;

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
            "util_punch_struct": v.util_punch_struct,
            "punch_kg_m3": v.punch_kg_m3,
            "punch_struct_worst": v.punch_struct_worst.as_ref().map(|w| json!({ "kind": w.0, "x": w.1[0], "y": w.1[1], "f": w.2, "cap": w.3 })),
            "shift": if v.shift.is_finite() { json!(v.shift) } else { Value::Null },
            "count": v.piles.len(),
            "kg_m3": v.kg_m3,
            "util_ground": v.util_ground,
            "util_material": v.util_material,
            "s_max": v.s_max,
            "ds_rel": v.ds_rel,
            "feasible": v.feasible,
            "design_impossible": v.design_impossible,
            "capped": v.capped,
            "capped_share": v.capped_share,
            "conditional": v.conditional.as_ref().map(|c| json!({ "s": c.s, "p": c.p, "area": c.area, "phi_mt": c.phi_mt, "hc": c.hc })),
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
            layers: vec![Layer { bottom: -8.0, e: 15_000.0, nu: 0.3, gamma: 0.0, phi: 0.0 }, Layer { bottom: -60.0, e: 40_000.0, nu: 0.3, gamma: 0.0, phi: 0.0 }],
        },
        s_max: 0.15,
        ds_rel: 0.002,
        layout: rows(),
        acrc_long: acrc_long(),
        acrc_short: acrc_short(),
        punching: true,
        existing: vec![],
        pile_fu: None,
        pile_model: sp24(),
    }
}

impl Service {
    /// Сверка с расчётом исходной модели (ЛИРА): плита на регулярной сетке со зданием
    /// суперэлементом и сваями-пружинами `mesh.piles` (жёсткость Rz КЭ 56) в ближайших узлах.
    /// `combos` — коэффициенты к загружениям сетки (как они заданы, без γf). Для каждой сваи —
    /// Rz (кН, сжатие < 0, как в ЛИРА) по сочетаниям и смещение до узла; для плиты — наибольшая
    /// и наименьшая вертикальная осадка (м, вниз < 0). Грунт и правила СП 24 не участвуют.
    pub fn verify_mesh(&mut self, mut input: Input, mesh: &crate::building::Mesh, combos: &[Vec<f64>]) -> Result<Value, String> {
        let s = crate::building::split(mesh)?;
        if combos.iter().any(|c| c.len() != mesh.cases.len()) {
            return Err(format!("сочетание: нужно {} коэффициентов", mesh.cases.len()));
        }
        input.raft.outline = s.region.outer.clone();
        input.raft.holes = s.region.holes.clone();
        input.raft.h = s.h;
        let eb = concrete(&input.raft.concrete)?.chars(CalcType::Uls).eb;
        let mut m = raft_model(&input.raft, eb)?;
        let nearest = |p: [f64; 3], m: &Model| (0..m.nodes.len()).min_by(|&a, &b| dist2(m.nodes[a], p).total_cmp(&dist2(m.nodes[b], p))).unwrap();
        let to_grid: Vec<usize> = s.raft.nodes.iter().map(|&p| nearest(p, &m)).collect();
        let mut loads: Vec<LoadCase> = s.raft_cases.iter().map(|lc| LoadCase { nodal: lc.nodal.iter().map(|&(n, f)| (to_grid[n], f)).collect(), pressure: vec![] }).collect();
        if let Some(b) = &s.building {
            m.blocks.push(crate::fem::DenseBlock { dofs: b.block.dofs.iter().map(|&d| 6 * to_grid[d / 6] + d % 6).collect(), k: b.block.k.clone() });
            for &d in b.block.dofs.iter().step_by(3) {
                m.fixed.push((to_grid[d / 6], [true, true, false, false, false, true]));
            }
            for (lc, l) in loads.iter_mut().zip(&b.loads) {
                lc.nodal.extend(l.iter().map(|&(n, f)| (to_grid[n], f)));
            }
        }
        let piles: Vec<(usize, f64, f64)> = mesh
            .piles
            .iter()
            .map(|(p, k)| {
                let n = nearest([p[0], p[1], 0.0], &m);
                (n, k[2], dist2(m.nodes[n], [p[0], p[1], 0.0]).sqrt())
            })
            .collect();
        for &(n, k, _) in &piles {
            m.springs.push(crate::fem::Spring { node: n, k: [0.0, 0.0, k, 0.0, 0.0, 0.0] });
        }
        let lin = m.factorize()?;
        let fs: Vec<Vec<f64>> = loads.iter().map(|lc| m.load_vector(lc)).collect::<Result<_, _>>()?;
        let us = lin.solve_many(&fs);
        let mut rz = vec![vec![0.0; combos.len()]; piles.len()];
        let mut w = vec![];
        for (c, f) in combos.iter().enumerate() {
            let u: Vec<f64> = (0..us[0].len()).map(|i| f.iter().zip(&us).map(|(k, u)| k * u[i]).sum()).collect();
            for (i, &(n, k, _)) in piles.iter().enumerate() {
                rz[i][c] = k * u[6 * n + 2];
            }
            let ws: Vec<f64> = (0..m.nodes.len()).map(|n| u[6 * n + 2]).collect();
            w.push(json!({ "min": ws.iter().cloned().fold(f64::INFINITY, f64::min), "max": ws.iter().cloned().fold(f64::NEG_INFINITY, f64::max) }));
        }
        Ok(json!({
            "piles": piles.iter().zip(&rz).map(|(p, r)| json!({ "node": p.0, "snap": p.2, "rz": r })).collect::<Vec<_>>(),
            "settlement": w,
            "notes": s.notes,
            "raft": { "area": polygon_area(&s.region.outer) - s.region.holes.iter().map(|h| polygon_area(h)).sum::<f64>(), "h": s.h, "z": s.z_mid, "holes": s.region.holes.len() },
        }))
    }
}

impl Service {
    /// Армирование подземной части (этап 2): модель этапа 1 варианта (плита, здание суперэлементом,
    /// сваи — пружины с податливостью варианта, без предела Fu) даёт перемещения узлов сопряжения;
    /// здание на своих узлах решается с этими перемещениями (пружины-штрафы uz, θx, θy; сдвиги и
    /// поворот вокруг вертикали закреплены) — без искажения геометрии слиянием узлов. Подбор (`reinforce::design`, СП 63) —
    /// пластины здания ниже поверхности грунта; сечение — толщина пластины, защитный слой и
    /// диаметр — как у плиты, бетон и арматура — как у плиты (допущение).
    fn design_underground(&mut self, args: &Value, progress: &mut dyn FnMut(&str)) -> Result<Value, String> {
        let i = args.get("variant").and_then(Value::as_u64).ok_or("не выбран вариант")? as usize;
        let (input, raft, vs) = self.last.as_ref().ok_or("нет результата расчёта")?;
        let lm = self.last_mesh.as_ref().ok_or("нет расчёта по сетке здания")?;
        let v = vs.get(i).ok_or("нет такого варианта")?;
        let bm = lm.split.building_model.as_ref().ok_or("в сетке нет здания")?;
        progress("плита варианта");
        // 1. Модель этапа 1: плита, здание суперэлементом, сваи — группа пружин; перемещения по загружениям.
        let np = v.piles.len();
        let mut rm = raft.clone();
        let blk = lm.building.as_ref().ok_or("в сетке нет здания")?;
        rm.blocks.push(blk.block.clone());
        for &d in blk.block.dofs.iter().step_by(3) {
            rm.fixed.push((d / 6, [true, true, false, false, false, true]));
        }
        rm.spring_groups.push(crate::fem::SpringGroup { nodes: v.piles.clone(), k: invert(&v.pile_flex, np)? });
        let lin = rm.factorize()?;
        let u_raft: Vec<Vec<f64>> = lm
            .raft_cases
            .iter()
            .zip(&blk.loads)
            .map(|(rc, bl)| lin.solve(&rm, &LoadCase { nodal: rc.nodal.iter().copied().chain(bl.iter().copied()).collect(), pressure: vec![] }))
            .collect::<Result<_, _>>()?;
        // 2. Здание на своих узлах: в узлах сопряжения — перемещения плиты (uz, θx, θy) пружинами-
        // штрафами, сдвиги и поворот вокруг вертикали закреплены (как при конденсации).
        progress("здание: перемещения узлов сопряжения");
        const PENALTY: f64 = 1e10;
        let mut m = bm.clone();
        let mut ifs: Vec<(usize, usize)> = lm.split.interface.iter().map(|&(b, r)| (b, lm.to_grid[r])).collect();
        ifs.sort_unstable();
        for &(b, _) in &ifs {
            m.fixed.push((b, [true, true, false, false, false, true]));
            m.springs.push(crate::fem::Spring { node: b, k: [0.0, 0.0, PENALTY, PENALTY, PENALTY, 0.0] });
        }
        // Перемещения плиты в точке узла сопряжения: билинейно по элементу сетки (иначе — ближайший узел).
        let weights: Vec<Vec<(usize, f64)>> = ifs
            .iter()
            .map(|&(b, g)| {
                let p = [bm.nodes[b][0], bm.nodes[b][1]];
                match element_at(raft, p) {
                    Some(e) if raft.shells[e].nodes.len() == 4 => {
                        let ns = &raft.shells[e].nodes;
                        let (x0, x1) = ns.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &n| (a.0.min(raft.nodes[n][0]), a.1.max(raft.nodes[n][0])));
                        let (y0, y1) = ns.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |a, &n| (a.0.min(raft.nodes[n][1]), a.1.max(raft.nodes[n][1])));
                        let (sx, sy) = ((p[0] - x0) / (x1 - x0), (p[1] - y0) / (y1 - y0));
                        ns.iter()
                            .map(|&n| {
                                let fx = if (raft.nodes[n][0] - x0).abs() < 1e-9 { 1.0 - sx } else { sx };
                                let fy = if (raft.nodes[n][1] - y0).abs() < 1e-9 { 1.0 - sy } else { sy };
                                (n, fx * fy)
                            })
                            .collect()
                    }
                    _ => vec![(g, 1.0)],
                }
            })
            .collect();
        let cases: Vec<LoadCase> = lm
            .building_cases
            .iter()
            .zip(&u_raft)
            .map(|(bc, u)| {
                let mut nodal = bc.nodal.clone();
                for (&(b, _), w) in ifs.iter().zip(&weights) {
                    let v = |k: usize| w.iter().map(|&(g, c)| c * u[6 * g + k]).sum::<f64>();
                    nodal.push((b, [0.0, 0.0, PENALTY * v(2), PENALTY * v(3), PENALTY * v(4), 0.0]));
                }
                LoadCase { nodal, pressure: vec![] }
            })
            .collect();
        let raft_shells = 0;
        // Подземные пластины здания и группы по толщине.
        let r = &input.raft;
        let cc = concrete(&r.concrete)?;
        let rc = rebar(&r.rebar)?;
        let (c_uls, r_uls) = (cc.chars(CalcType::Uls), rc.chars(CalcType::Uls));
        let (c_sls, r_sls) = (cc.chars(CalcType::Sls), rc.chars(CalcType::Sls));
        let under: Vec<bool> = (0..m.shells.len()).map(|e| e >= raft_shells && m.shells[e].nodes.iter().all(|&n| m.nodes[n][2] <= lm.ground + 1e-6)).collect();
        let mut hs: Vec<f64> = (0..m.shells.len()).filter(|&e| under[e]).map(|e| (m.shells[e].mat.h * 1000.0).round() / 1000.0).collect();
        hs.sort_by(|a, b| a.total_cmp(b));
        hs.dedup();
        if hs.is_empty() {
            return Err(format!("нет пластин здания ниже поверхности грунта (отметка {:.2})", lm.ground));
        }
        let secs: Vec<PlateSection> = hs.iter().map(|&h| PlateSection { h, a_top: r.cover, a_bot: r.cover, ds_x: r.bar_d, ds_y: r.bar_d, concrete_diagram: ConcreteDiagramType::L3 }).collect();
        let ulss: Vec<UlsContext> = secs.iter().map(|s| UlsContext::new(s, &c_uls, &r_uls)).collect();
        let crack = |lim: f64, s: &PlateSection| (lim > 0.0).then(|| SlsContext::new(s, &c_sls, &r_sls, Some((&c_uls, &r_uls)), SlsOptions { acrc_lim_mm: lim, phi1: 1.4, ..Default::default() }));
        let cl: Vec<Option<SlsContext>> = secs.iter().map(|s| crack(input.acrc_long, s)).collect();
        let cs: Vec<Option<SlsContext>> = secs.iter().map(|s| crack(input.acrc_short, s)).collect();
        let groups: Vec<crate::reinforce::Group> = (0..secs.len())
            .map(|g| crate::reinforce::Group {
                section: secs[g].clone(),
                concrete: c_uls.clone(),
                rebar: r_uls.clone(),
                rsw: rc.rsw(),
                uls: &ulss[g],
                crack_short: cs[g].as_ref(),
                crack_long: cl[g].as_ref(),
                background: [0.0; 4],
            })
            .collect();
        let group_of: Vec<Option<usize>> = (0..m.shells.len())
            .map(|e| under[e].then(|| hs.iter().position(|&h| (h - (m.shells[e].mat.h * 1000.0).round() / 1000.0).abs() < 1e-9).unwrap()))
            .collect();
        progress("армирование подземной части");
        let out = crate::reinforce::design(&m, &groups, &group_of, &cases, &lm.combos.uls, &lm.combos.sls, Some(&lm.combos.sls_long))?;
        // Масса: продольная — Σ площадей на метр × площадь элемента; поперечная — Asw/sw × площадь × h.
        let (mut kg, mut vol, mut n, mut impossible, mut strip, mut asw_max) = (0.0, 0.0, 0usize, 0usize, 0usize, 0.0f64);
        let mut fields = vec![[0.0; 4]; m.shells.len()];
        for e in 0..m.shells.len() {
            if !under[e] {
                continue;
            }
            let c: Vec<crate::fem::V3> = m.shells[e].nodes.iter().map(|&k| m.nodes[k]).collect();
            let a = shell_area3(&c);
            let h = m.shells[e].mat.h;
            let el = &out[e];
            kg += (el.required.iter().sum::<f64>() * a + el.asw * a * h) * 7850.0;
            vol += a * h;
            n += 1;
            impossible += el.impossible as usize;
            strip += el.strip_fails as usize;
            asw_max = asw_max.max(el.asw);
            fields[e] = el.required;
        }
        let mx: Vec<f64> = (0..4).map(|k| fields.iter().map(|f| f[k] * 1e4).fold(0.0, f64::max)).collect();
        let diag = {
            // Диагностика худшего элемента: перемещения узлов в первом сочетании ПС1, примыкающие
            // стержни и пластины, узловые нагрузки.
            let wn = m.shells.iter().enumerate().max_by(|a, b| fields[a.0].iter().cloned().fold(0.0, f64::max).total_cmp(&fields[b.0].iter().cloned().fold(0.0, f64::max))).unwrap().0;
            let ns = m.shells[wn].nodes.clone();
            let lin = m.factorize()?;
            let c0 = &lm.combos.uls[0];
            let us: Vec<Vec<f64>> = cases.iter().map(|lc| lin.solve(&m, lc)).collect::<Result<_, _>>()?;
            let u = |n: usize, k: usize| c0.iter().zip(&us).map(|(f, u)| f * u[6 * n + k]).sum::<f64>();
            json!(ns.iter().map(|&n| json!({
                "node": n, "xyz": m.nodes[n], "uz_mm": u(n, 2) * 1e3, "rx": u(n, 3), "ry": u(n, 4),
                "beams": m.beams.iter().filter(|b| b.n1 == n || b.n2 == n).count(),
                "shells": m.shells.iter().filter(|s| s.nodes.contains(&n)).count(),
                "fz": cases.iter().zip(c0).map(|(lc, f)| f * lc.nodal.iter().filter(|x| x.0 == n).map(|x| x.1[2]).sum::<f64>()).sum::<f64>(),
            })).collect::<Vec<_>>())
        };
        // Элемент с наибольшим армированием: центр, толщина, требуемое (см²/м).
        let worst = (0..m.shells.len()).max_by(|&a, &b| fields[a].iter().cloned().fold(0.0, f64::max).total_cmp(&fields[b].iter().cloned().fold(0.0, f64::max))).unwrap();
        let wc = m.shells[worst].nodes.iter().fold([0.0; 3], |a, &k| { let p = m.nodes[k]; let n = m.shells[worst].nodes.len() as f64; [a[0] + p[0] / n, a[1] + p[1] / n, a[2] + p[2] / n] });
        let summary = json!({
            "variant": i, "elements": n, "ground": lm.ground, "kg": kg, "kg_m3": kg / vol.max(1e-9), "volume": vol,
            "impossible": impossible, "strip_fails": strip, "asw_max": asw_max, "thicknesses": hs,
            "max": mx,
            "worst": { "x": wc[0], "y": wc[1], "z": wc[2], "h": m.shells[worst].mat.h, "as": fields[worst].map(|x| x * 1e4), "diag": diag },
        });
        self.last_underground = Some((m, fields));
        Ok(summary)
    }
}

fn shell_area3(p: &[crate::fem::V3]) -> f64 {
    let mut s = [0.0; 3];
    for k in 1..p.len() - 1 {
        let (u, v) = ([p[k][0] - p[0][0], p[k][1] - p[0][1], p[k][2] - p[0][2]], [p[k + 1][0] - p[0][0], p[k + 1][1] - p[0][1], p[k + 1][2] - p[0][2]]);
        let c = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        (0..3).for_each(|i| s[i] += c[i]);
    }
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() / 2.0
}

/// Обращение симметричной положительно определённой матрицы n × n (Гаусс — Жордан).
fn invert(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; n * n];
    (0..n).for_each(|i| inv[i * n + i] = 1.0);
    for c in 0..n {
        let p = (c..n).max_by(|&x, &y| m[x * n + c].abs().total_cmp(&m[y * n + c].abs())).unwrap();
        if m[p * n + c].abs() < 1e-300 {
            return Err("матрица податливости свай вырождена".into());
        }
        for k in 0..n {
            m.swap(c * n + k, p * n + k);
            inv.swap(c * n + k, p * n + k);
        }
        let d = m[c * n + c];
        for k in 0..n {
            m[c * n + k] /= d;
            inv[c * n + k] /= d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r * n + c];
                if f != 0.0 {
                    for k in 0..n {
                        m[r * n + k] -= f * m[c * n + k];
                        inv[r * n + k] -= f * inv[c * n + k];
                    }
                }
            }
        }
    }
    Ok(inv)
}
