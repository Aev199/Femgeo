//! Чтение объекта geobase (`*.db`, SQLite) только на чтение: скважины, слои, ИГЭ и их
//! параметры по моделям грунта. Схема — geobase `create_db.py`.
//!
//! Глубина слоя (`borehole_layers.depth`) — подошва слоя от устья скважины, как в экспорте
//! geobase в PLAXIS. Модуль деформации — свойство «Е» (кириллица) или «E», кПа; ν — «ν» или
//! «nu». Грунтовая колонка в точке — по ближайшей скважине (поверхности кригинга geobase —
//! следующий шаг).

use std::collections::HashMap;
use std::path::Path;

use rusqlite::{Connection, OpenFlags};

use crate::piles::{Layer, SoilColumn};

#[derive(Clone, Debug, PartialEq)]
pub struct Soil {
    pub id: i64,
    pub name: String,
    /// Модуль деформации, кПа.
    pub e: Option<f64>,
    pub nu: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Borehole {
    pub name: String,
    pub x: f64,
    pub y: f64,
    /// Отметка устья, м.
    pub z: f64,
    /// (ИГЭ, глубина подошвы от устья, м) сверху вниз.
    pub layers: Vec<(i64, f64)>,
}

#[derive(Clone, Debug)]
pub struct Object {
    pub soils: HashMap<i64, Soil>,
    pub boreholes: Vec<Borehole>,
    /// Модель грунта, из которой взяты E и ν.
    pub model: String,
}

const E_KEYS: [&str; 2] = ["Е", "E"];
const NU_KEYS: [&str; 3] = ["ν", "nu", "v"];

impl Object {
    /// Открыть объект. `model` — имя модели грунта (`ConstitutiveModels.name`); None — первая
    /// по алфавиту модель, в которой у ИГЭ задан модуль деформации.
    pub fn open(path: &Path, model: Option<&str>) -> Result<Self, String> {
        let c = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let q = |e: rusqlite::Error| format!("geobase: {e}");

        let mut props: HashMap<(String, i64), (Option<f64>, Option<f64>)> = HashMap::new();
        let mut st = c
            .prepare("SELECT m.name, p.soil_id, p.property_name, p.property_value FROM SoilProperties p JOIN ConstitutiveModels m ON m.id = p.model_id")
            .map_err(q)?;
        let rows = st
            .query_map([], |r| Ok((r.get::<_, String>(0)?, r.get::<_, i64>(1)?, r.get::<_, String>(2)?, r.get::<_, Option<f64>>(3)?)))
            .map_err(q)?;
        for row in rows {
            let (m, soil, name, value) = row.map_err(q)?;
            let entry = props.entry((m, soil)).or_default();
            let name = name.trim();
            if E_KEYS.contains(&name) {
                entry.0 = value.or(entry.0);
            } else if NU_KEYS.contains(&name) {
                entry.1 = value.or(entry.1);
            }
        }
        let model = match model {
            Some(m) => m.to_string(),
            None => {
                let mut names: Vec<&String> = props.iter().filter(|(_, v)| v.0.is_some()).map(|((m, _), _)| m).collect();
                names.sort();
                names.first().map(|s| s.to_string()).ok_or("geobase: ни у одного ИГЭ не задан модуль деформации")?
            }
        };

        let mut soils = HashMap::new();
        let mut st = c.prepare("SELECT id, name FROM Soils").map_err(q)?;
        for row in st.query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, String>(1)?))).map_err(q)? {
            let (id, name) = row.map_err(q)?;
            let (e, nu) = props.get(&(model.clone(), id)).copied().unwrap_or((None, None));
            soils.insert(id, Soil { id, name, e, nu });
        }

        let mut boreholes = Vec::new();
        let mut st = c.prepare("SELECT id, name, x_coord, y_coord, z_coord FROM boreholes ORDER BY id").map_err(q)?;
        let bh: Vec<(i64, String, f64, f64, f64)> = st
            .query_map([], |r| Ok((r.get(0)?, r.get(1)?, r.get::<_, Option<f64>>(2)?.unwrap_or(f64::NAN), r.get::<_, Option<f64>>(3)?.unwrap_or(f64::NAN), r.get::<_, Option<f64>>(4)?.unwrap_or(f64::NAN))))
            .map_err(q)?
            .collect::<Result<_, _>>()
            .map_err(q)?;
        let mut st = c
            .prepare("SELECT soil_id, depth FROM borehole_layers WHERE borehole_id = ?1 ORDER BY depth, layer_number")
            .map_err(q)?;
        for (id, name, x, y, z) in bh {
            if !(x.is_finite() && y.is_finite() && z.is_finite()) {
                return Err(format!("geobase: у скважины «{name}» не заданы координаты"));
            }
            let layers: Vec<(i64, f64)> = st
                .query_map([id], |r| Ok((r.get::<_, Option<i64>>(0)?.unwrap_or(-1), r.get::<_, Option<f64>>(1)?.unwrap_or(f64::NAN))))
                .map_err(q)?
                .collect::<Result<_, _>>()
                .map_err(q)?;
            boreholes.push(Borehole { name, x, y, z, layers });
        }
        Ok(Self { soils, boreholes, model })
    }

    /// Грунтовая колонка скважины: слои нулевой толщины пропускаются; ИГЭ без E или ν — ошибка.
    pub fn column(&self, b: &Borehole) -> Result<SoilColumn, String> {
        let mut layers = Vec::new();
        let mut top_depth = 0.0;
        for &(soil, depth) in &b.layers {
            if !(depth > top_depth) {
                continue;
            }
            let s = self.soils.get(&soil).ok_or_else(|| format!("скважина «{}»: ИГЭ {soil} нет в справочнике", b.name))?;
            let (Some(e), Some(nu)) = (s.e, s.nu) else {
                return Err(format!("ИГЭ «{}»: в модели «{}» не заданы E и ν", s.name, self.model));
            };
            layers.push(Layer { bottom: b.z - depth, e, nu });
            top_depth = depth;
        }
        if layers.is_empty() {
            return Err(format!("скважина «{}»: нет слоёв", b.name));
        }
        Ok(SoilColumn { top: b.z, layers })
    }

    /// Колонка по ближайшей скважине в плане.
    pub fn column_at(&self, x: f64, y: f64) -> Result<SoilColumn, String> {
        let b = self
            .boreholes
            .iter()
            .min_by(|a, b| ((a.x - x).powi(2) + (a.y - y).powi(2)).total_cmp(&((b.x - x).powi(2) + (b.y - y).powi(2))))
            .ok_or("geobase: в объекте нет скважин")?;
        self.column(b)
    }
}
