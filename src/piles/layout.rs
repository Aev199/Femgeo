//! Ограничения расстановки свай: допустимый контур с отверстиями, минимальные расстояния
//! между осями и от осей до границы контура.
//!
//! Расстояния между осями — СП 24.13330.2021, п. 8.13 (сверено 2026-11-01): висячие забивные
//! и вдавливаемые — 3d, сваи-стойки — 1,5d, буровые, набивные — в свету не менее 1,0 м.
//! Отступ от края ростверка в СП 24 числом не задан (п. 8.3 — с учётом допускаемых отклонений
//! свай): 0,5d + 0,1 м — допущение, параметр пользователя.

/// Многоугольник: вершины по порядку (замыкание подразумевается).
pub type Polygon = Vec<[f64; 2]>;

/// Допустимая область: внешний контур и отверстия.
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub outer: Polygon,
    pub holes: Vec<Polygon>,
}

/// Тип сваи для правил расстояний.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PileKind {
    /// Висячая забивная.
    DrivenFriction,
    /// Свая-стойка.
    EndBearing,
    /// Буровая, набивная.
    Bored,
}

/// Правила расстановки, м.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rules {
    /// Наименьшее расстояние между осями свай.
    pub min_spacing: f64,
    /// Наименьшее расстояние от оси сваи до границы допустимой области.
    pub min_edge: f64,
}

impl Rules {
    /// Шаг по п. 8.13 СП 24; до границы — 0,5d + 0,1 м (допущение, см. модуль).
    pub fn sp24(kind: PileKind, d: f64) -> Self {
        let min_spacing = match kind {
            PileKind::DrivenFriction => 3.0 * d,
            PileKind::EndBearing => 1.5 * d,
            PileKind::Bored => d + 1.0,
        };
        Self { min_spacing, min_edge: 0.5 * d + 0.1 }
    }
}

fn inside(poly: &Polygon, p: [f64; 2]) -> bool {
    let mut c = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
            if p[0] < x {
                c = !c;
            }
        }
    }
    c
}

fn dist_to_boundary(poly: &Polygon, p: [f64; 2]) -> f64 {
    let n = poly.len();
    (0..n)
        .map(|i| {
            let (a, b) = (poly[i], poly[(i + 1) % n]);
            let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
            let l2 = dx * dx + dy * dy;
            let t = if l2 > 0.0 { (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
            ((p[0] - a[0] - t * dx).powi(2) + (p[1] - a[1] - t * dy).powi(2)).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

impl Region {
    /// Точка внутри области (вне отверстий).
    pub fn contains(&self, p: [f64; 2]) -> bool {
        inside(&self.outer, p) && !self.holes.iter().any(|h| inside(h, p))
    }

    /// Расстояние до ближайшей границы (внешней или отверстия).
    pub fn edge_distance(&self, p: [f64; 2]) -> f64 {
        self.holes.iter().map(|h| dist_to_boundary(h, p)).fold(dist_to_boundary(&self.outer, p), f64::min)
    }

    /// Допустима ли свая в точке по контуру и отступу.
    pub fn admits(&self, p: [f64; 2], rules: &Rules) -> bool {
        self.contains(p) && self.edge_distance(p) >= rules.min_edge - 1e-9
    }
}

/// Нарушение правил расстановки.
#[derive(Clone, Debug, PartialEq)]
pub enum Violation {
    /// Свая вне области или ближе `min_edge` к границе: (свая, расстояние до границы; < 0 — вне).
    Edge(usize, f64),
    /// Сваи ближе `min_spacing`: (i, j, расстояние).
    Spacing(usize, usize, f64),
}

/// Все нарушения для набора свай. Пары по близости ищутся через сетку ячеек размером
/// `min_spacing` — O(n) в среднем.
pub fn violations(region: &Region, rules: &Rules, piles: &[[f64; 2]]) -> Vec<Violation> {
    let mut out = Vec::new();
    for (i, &p) in piles.iter().enumerate() {
        if !region.contains(p) {
            out.push(Violation::Edge(i, -region.edge_distance(p)));
        } else {
            let d = region.edge_distance(p);
            if d < rules.min_edge - 1e-9 {
                out.push(Violation::Edge(i, d));
            }
        }
    }
    let s = rules.min_spacing.max(1e-9);
    let mut cells: std::collections::HashMap<(i64, i64), Vec<usize>> = std::collections::HashMap::new();
    for (i, p) in piles.iter().enumerate() {
        cells.entry(((p[0] / s).floor() as i64, (p[1] / s).floor() as i64)).or_default().push(i);
    }
    for (i, p) in piles.iter().enumerate() {
        let (cx, cy) = ((p[0] / s).floor() as i64, (p[1] / s).floor() as i64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for &j in cells.get(&(cx + dx, cy + dy)).map_or(&[][..], |v| v) {
                    if j > i {
                        let d = ((p[0] - piles[j][0]).powi(2) + (p[1] - piles[j][1]).powi(2)).sqrt();
                        if d < rules.min_spacing - 1e-9 {
                            out.push(Violation::Spacing(i, j, d));
                        }
                    }
                }
            }
        }
    }
    out
}
