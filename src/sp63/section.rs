//! Фибровое сечение стержня (НДМ, пп. 8.1.20–8.1.30 СП 63): бетонные волокна, стержни
//! арматуры, плоскость деформаций ε(x, y) = ε₀ + ky·y + kz·x.
//!
//! Источник: OpenCS `CScore/CrossSection.cs` (интегрирование), `CScore/Diagramm.cs`
//! (разностная диаграмма стержня в бетоне), `CScore/LimitSectionStrainSolver.cs` (Ньютон 3×3),
//! `CScore/CrossSectionLimitAdapter.cs` (εcu, εsu), коммит `387e157`. Отличия: касательная
//! матрица аналитическая (в OpenCS — конечные разности) и шаг Ньютона с демпфированием, как в
//! НДМ полосы; круглое сечение разбивается на кольца и секторы точно по площади многоугольника
//! (в OpenCS — сетка `SliceXY`). Сходимость и допуски — как в OpenCS.

use std::sync::Arc;

use super::diagram::Diagram;

/// Волокно: центр (x, y), м, и площадь, м².
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fiber {
    pub x: f64,
    pub y: f64,
    pub a: f64,
}

/// Плоскость деформаций.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Plane {
    pub e0: f64,
    /// Кривизна при y (момент Mx = Σσ·A·y).
    pub ky: f64,
    /// Кривизна при x (момент My = Σσ·A·x).
    pub kz: f64,
}

impl Plane {
    pub fn eps(&self, x: f64, y: f64) -> f64 {
        self.e0 + self.ky * y + self.kz * x
    }
}

/// Усилия и касательная матрица ∂(N, Mx, My)/∂(ε₀, ky, kz).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Response {
    pub n: f64,
    pub mx: f64,
    pub my: f64,
    pub k: [[f64; 3]; 3],
}

/// Фибровое сечение.
#[derive(Clone, Debug)]
pub struct Section {
    pub concrete: Vec<Fiber>,
    pub bars: Vec<Fiber>,
    pub concrete_diagram: Arc<Diagram>,
    pub rebar_diagram: Arc<Diagram>,
    /// Вершины контура бетона — для проверки εcu.
    pub contour: Vec<(f64, f64)>,
    /// Стержни работают по разностной диаграмме σs − σb (бетон не вычитается из волокон).
    pub differential: bool,
    /// Предельная деформация сжатия бетона (< 0): min(εb2, первый узел ветви сжатия).
    pub eps_cu: f64,
    /// Предельная деформация растяжения арматуры: последний узел ветви растяжения.
    pub eps_su: f64,
}

impl Section {
    pub fn new(
        concrete: Vec<Fiber>,
        bars: Vec<Fiber>,
        contour: Vec<(f64, f64)>,
        concrete_diagram: Diagram,
        rebar_diagram: Diagram,
        eps_b2: f64,
    ) -> Self {
        let eps_cu = (-eps_b2.abs()).min(concrete_diagram.compression.xs()[0]);
        let eps_su = rebar_diagram.tension.last_x();
        Self {
            concrete,
            bars,
            concrete_diagram: Arc::new(concrete_diagram),
            rebar_diagram: Arc::new(rebar_diagram),
            contour,
            differential: true,
            eps_cu,
            eps_su,
        }
    }

    /// Круглое сечение: правильный многоугольник с `sides` сторонами радиуса `radius`,
    /// `rings` колец; `bar_count` стержней диаметром `bar_d` на радиусе `bar_radius`,
    /// первый — под углом `rotation_deg` от оси x.
    #[allow(clippy::too_many_arguments)]
    pub fn circular(
        radius: f64,
        sides: usize,
        rings: usize,
        bar_count: usize,
        bar_d: f64,
        bar_radius: f64,
        rotation_deg: f64,
        concrete_diagram: Diagram,
        rebar_diagram: Diagram,
        eps_b2: f64,
    ) -> Self {
        use std::f64::consts::PI;
        let v: Vec<(f64, f64)> = (0..sides)
            .map(|i| {
                let a = 2.0 * PI * i as f64 / sides as f64;
                (radius * a.cos(), radius * a.sin())
            })
            .collect();
        let mut fibers = Vec::with_capacity(sides * rings);
        for i in 0..sides {
            let (p, q) = (v[i], v[(i + 1) % sides]);
            let tri = 0.5 * (p.0 * q.1 - p.1 * q.0);
            let (gx, gy) = ((p.0 + q.0) / 3.0, (p.1 + q.1) / 3.0);
            for j in 0..rings {
                let (t1, t2) = (j as f64 / rings as f64, (j + 1) as f64 / rings as f64);
                let (s2, s3) = (t2 * t2 - t1 * t1, t2 * t2 * t2 - t1 * t1 * t1);
                let f = s3 / s2;
                fibers.push(Fiber { x: gx * f, y: gy * f, a: tri * s2 });
            }
        }
        let area = PI * bar_d * bar_d / 4.0;
        let rot = rotation_deg.to_radians();
        let bars = (0..bar_count)
            .map(|i| {
                let a = rot + 2.0 * PI * i as f64 / bar_count as f64;
                Fiber { x: bar_radius * a.cos(), y: bar_radius * a.sin(), a: area }
            })
            .collect();
        Self::new(fibers, bars, v, concrete_diagram, rebar_diagram, eps_b2)
    }

    /// Усилия и касательная матрица; `ten` — учитывать растянутый бетон.
    pub fn integrate(&self, p: &Plane, ten: bool) -> Response {
        let mut r = Response::default();
        let mut add = |x: f64, y: f64, sig: f64, e2: f64, a: f64| {
            let f = sig * a;
            r.n += f;
            r.mx += f * y;
            r.my += f * x;
            let w = e2 * a;
            let g = [1.0, y, x];
            for i in 0..3 {
                for j in 0..3 {
                    r.k[i][j] += w * g[i] * g[j];
                }
            }
        };
        for c in &self.concrete {
            let (s, e) = self.concrete_diagram.sig(p.eps(c.x, c.y), ten, true);
            add(c.x, c.y, s, e, c.a);
        }
        for b in &self.bars {
            let eps = p.eps(b.x, b.y);
            let (mut s, mut e) = self.rebar_diagram.sig(eps, true, true);
            if self.differential {
                let (sc, ec) = self.concrete_diagram.sig(eps, ten, true);
                s -= sc;
                e -= ec;
            }
            add(b.x, b.y, s, e, b.a);
        }
        r
    }

    fn elastic(&self) -> [[f64; 3]; 3] {
        let (_, eb) = self.concrete_diagram.sig(-1e-9, false, true);
        let (_, es) = self.rebar_diagram.sig(1e-9, true, true);
        let mut k = [[0.0; 3]; 3];
        let mut add = |x: f64, y: f64, w: f64| {
            let g = [1.0, y, x];
            for i in 0..3 {
                for j in 0..3 {
                    k[i][j] += w * g[i] * g[j];
                }
            }
        };
        for c in &self.concrete {
            add(c.x, c.y, eb * c.a);
        }
        for b in &self.bars {
            add(b.x, b.y, (es - if self.differential { eb } else { 0.0 }) * b.a);
        }
        k
    }

    /// Минимальная деформация по вершинам контура.
    pub fn eps_contour_min(&self, p: &Plane) -> f64 {
        self.contour.iter().map(|&(x, y)| p.eps(x, y)).fold(f64::INFINITY, f64::min)
    }

    /// Максимальная деформация стержней.
    pub fn eps_bars_max(&self, p: &Plane) -> f64 {
        self.bars.iter().map(|b| p.eps(b.x, b.y)).fold(f64::NEG_INFINITY, f64::max)
    }
}

/// Решение 3×3 по Гауссу с выбором ведущего элемента; None — вырождено.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let (mut m, mut v) = (a, b);
    for col in 0..3 {
        let piv = (col..3).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs())).unwrap();
        if !m[piv][col].is_finite() || m[piv][col].abs() < 1e-15 {
            return None;
        }
        m.swap(col, piv);
        v.swap(col, piv);
        for row in col + 1..3 {
            let f = m[row][col] / m[col][col];
            for k in col..3 {
                m[row][k] -= f * m[col][k];
            }
            v[row] -= f * v[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        let mut s = v[row];
        for k in row + 1..3 {
            s -= m[row][k] * x[k];
        }
        x[row] = s / m[row][row];
    }
    Some(x)
}

/// Результат поиска плоскости деформаций.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Equilibrium {
    pub plane: Plane,
    pub converged: bool,
    /// √(ΔN² + ΔMx² + ΔMy²), кН и кН·м.
    pub residual: f64,
    pub iterations: usize,
}

/// Плоскость деформаций при (N, Mx, My): Ньютон с аналитической касательной, демпфированием и
/// упругой матрицей при вырождении. Сходимость — невязка < `tol` (в OpenCS 0,5).
pub fn solve(s: &Section, n: f64, mx: f64, my: f64, ten: bool, tol: f64, max_iter: usize) -> Equilibrium {
    let el = s.elastic();
    let x0 = solve3(el, [n, mx, my]).unwrap_or([0.0; 3]);
    let mut p = Plane { e0: x0[0], ky: x0[1], kz: x0[2] };
    let res = |f: &Response| [n - f.n, mx - f.mx, my - f.my];
    let norm = |r: [f64; 3]| (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let mut f = s.integrate(&p, ten);
    let mut r = res(&f);
    let mut nr = norm(r);
    for it in 0..max_iter {
        if nr < tol {
            return Equilibrium { plane: p, converged: true, residual: nr, iterations: it };
        }
        let d = solve3(f.k, r).or_else(|| solve3(el, r));
        let Some(d) = d else { break };
        let mut alpha = 1.0;
        let mut accepted = false;
        for _ in 0..30 {
            let q = Plane { e0: p.e0 + alpha * d[0], ky: p.ky + alpha * d[1], kz: p.kz + alpha * d[2] };
            let fq = s.integrate(&q, ten);
            let rq = res(&fq);
            let nq = norm(rq);
            if nq < nr {
                p = q;
                f = fq;
                r = rq;
                nr = nq;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            // Касательная не даёт спуска (сечение потекло) — шаг по упругой матрице.
            let Some(d) = solve3(el, r) else { break };
            let mut alpha = 1.0;
            for _ in 0..30 {
                let q = Plane { e0: p.e0 + alpha * d[0], ky: p.ky + alpha * d[1], kz: p.kz + alpha * d[2] };
                let fq = s.integrate(&q, ten);
                let nq = norm(res(&fq));
                if nq < nr {
                    p = q;
                    f = fq;
                    r = res(&f);
                    nr = nq;
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted {
                break;
            }
        }
    }
    Equilibrium { plane: p, converged: nr < tol, residual: nr, iterations: max_iter }
}
