//! Линейный МКЭ: оболочки, стержни, пружины, связи; разреженное разложение Холецкого (`faer`).
//!
//! Глобальные степени свободы узла: ux, uy, uz, θx, θy, θz.

pub mod condense;
pub mod shell;

use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::Llt;
use faer::sparse::{SparseColMat, Triplet};
use faer::{Mat, Side};

pub use shell::{ShellMaterial, V3};
use shell::{cross, norm};

/// Оболочка: 3 или 4 узла (обход против часовой стрелки задаёт нормаль).
#[derive(Clone, Debug)]
pub struct Shell {
    pub nodes: Vec<usize>,
    pub mat: ShellMaterial,
}

/// Стержень Эйлера–Бернулли. Локальная ось x — от `n1` к `n2`; ось y ⟂ `ref_dir` и x
/// (`ref_dir` по умолчанию — глобальная Z, для вертикальных стержней — глобальная X).
#[derive(Clone, Debug)]
pub struct Beam {
    pub n1: usize,
    pub n2: usize,
    pub e: f64,
    pub g: f64,
    pub a: f64,
    /// Момент инерции относительно локальной y (изгиб в плоскости xz).
    pub iy: f64,
    /// Момент инерции относительно локальной z (изгиб в плоскости xy).
    pub iz: f64,
    /// Момент инерции при кручении.
    pub j: f64,
    pub ref_dir: Option<V3>,
}

/// Упругая связь узла с землёй по глобальным осям (кН/м, кН·м/рад).
#[derive(Clone, Copy, Debug)]
pub struct Spring {
    pub node: usize,
    pub k: [f64; 6],
}

/// Связанная группа вертикальных пружин: полная матрица жёсткости n×n по uz узлов
/// (свайное поле с взаимовлиянием свай), построчно.
#[derive(Clone, Debug)]
pub struct SpringGroup {
    pub nodes: Vec<usize>,
    pub k: Vec<f64>,
}

/// Суперэлемент: полная матрица жёсткости по заданным глобальным степеням свободы, построчно
/// (например, сконденсированная надземная часть здания).
#[derive(Clone, Debug)]
pub struct DenseBlock {
    pub dofs: Vec<usize>,
    pub k: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct Model {
    pub nodes: Vec<V3>,
    pub shells: Vec<Shell>,
    pub beams: Vec<Beam>,
    pub springs: Vec<Spring>,
    pub spring_groups: Vec<SpringGroup>,
    pub blocks: Vec<DenseBlock>,
    /// Закреплённые степени свободы.
    pub fixed: Vec<(usize, [bool; 6])>,
}

/// Загружение: узловые силы (глобальные оси) и давление на оболочки (по нормали элемента).
#[derive(Clone, Debug, Default)]
pub struct LoadCase {
    pub nodal: Vec<(usize, [f64; 6])>,
    pub pressure: Vec<(usize, f64)>,
}

/// Локальные оси стержня (строки).
pub fn beam_frame(nodes: &[V3], b: &Beam) -> Result<([V3; 3], f64), String> {
    let (p, q) = (nodes[b.n1], nodes[b.n2]);
    let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
    let l = norm(d);
    if l < 1e-12 {
        return Err("стержень нулевой длины".into());
    }
    let ex = d.map(|x| x / l);
    let r = b.ref_dir.unwrap_or(if ex[2].abs() > 0.999 { [1.0, 0.0, 0.0] } else { [0.0, 0.0, 1.0] });
    let ey = cross(r, ex);
    let n = norm(ey);
    if n < 1e-9 {
        return Err("опорное направление стержня параллельно его оси".into());
    }
    let ey = ey.map(|x| x / n);
    Ok(([ex, ey, cross(ex, ey)], l))
}

fn beam_local(b: &Beam, l: f64) -> [[f64; 12]; 12] {
    let mut k = [[0.0; 12]; 12];
    let mut set = |i: usize, j: usize, v: f64| {
        k[i][j] += v;
        if i != j {
            k[j][i] += v;
        }
    };
    let (ea, gj) = (b.e * b.a / l, b.g * b.j / l);
    set(0, 0, ea);
    set(6, 6, ea);
    set(0, 6, -ea);
    set(3, 3, gj);
    set(9, 9, gj);
    set(3, 9, -gj);
    // Плоскость xy (v, θz), Iz.
    let (a, c, e, f) = (12.0 * b.e * b.iz / l.powi(3), 6.0 * b.e * b.iz / l.powi(2), 4.0 * b.e * b.iz / l, 2.0 * b.e * b.iz / l);
    set(1, 1, a);
    set(7, 7, a);
    set(1, 7, -a);
    set(1, 5, c);
    set(1, 11, c);
    set(5, 7, -c);
    set(7, 11, -c);
    set(5, 5, e);
    set(11, 11, e);
    set(5, 11, f);
    // Плоскость xz (w, θy), Iy.
    let (a, c, e, f) = (12.0 * b.e * b.iy / l.powi(3), 6.0 * b.e * b.iy / l.powi(2), 4.0 * b.e * b.iy / l, 2.0 * b.e * b.iy / l);
    set(2, 2, a);
    set(8, 8, a);
    set(2, 8, -a);
    set(2, 4, -c);
    set(2, 10, -c);
    set(4, 8, c);
    set(8, 10, c);
    set(4, 4, e);
    set(10, 10, e);
    set(4, 10, f);
    k
}

/// Глобальная матрица стержня 12×12.
pub fn beam_stiffness(nodes: &[V3], b: &Beam) -> Result<[[f64; 12]; 12], String> {
    let (r, l) = beam_frame(nodes, b)?;
    let kl = beam_local(b, l);
    // K = Tᵀ·Kl·T, T — блоки R (строки — локальные оси).
    let t = |i: usize, j: usize| if i / 3 == j / 3 { r[i % 3][j % 3] } else { 0.0 };
    let mut kt = [[0.0; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            let blk = j / 3 * 3;
            kt[i][j] = (blk..blk + 3).map(|p| kl[i][p] * t(p, j)).sum();
        }
    }
    let mut kg = [[0.0; 12]; 12];
    for i in 0..12 {
        let blk = i / 3 * 3;
        for j in 0..12 {
            kg[i][j] = (blk..blk + 3).map(|p| t(p, i) * kt[p][j]).sum();
        }
    }
    Ok(kg)
}

/// Внутренние усилия стержня в начале и в конце в локальных осях:
/// [N, Vy, Vz, T, My, Mz]; «+» N — растяжение.
pub fn beam_forces(nodes: &[V3], b: &Beam, u: &[f64]) -> Result<([f64; 6], [f64; 6]), String> {
    let (r, l) = beam_frame(nodes, b)?;
    let kl = beam_local(b, l);
    let mut ul = [0.0; 12];
    for (i, uli) in ul.iter_mut().enumerate() {
        let node = if i < 6 { b.n1 } else { b.n2 };
        let blk = i / 3 * 3;
        let g = &u[6 * node + (blk % 6)..6 * node + (blk % 6) + 3];
        *uli = (0..3).map(|k| r[i % 3][k] * g[k]).sum();
    }
    let f: Vec<f64> = (0..12).map(|i| (0..12).map(|j| kl[i][j] * ul[j]).sum()).collect();
    let start = [-f[0], -f[1], -f[2], -f[3], -f[4], -f[5]];
    let end = [f[6], f[7], f[8], f[9], f[10], f[11]];
    Ok((start, end))
}

/// Разложенная система: решение для любого числа загружений.
pub struct Linear {
    llt: Llt<usize, f64>,
    /// Номер уравнения для каждой глобальной степени свободы; None — закреплена.
    eq: Vec<Option<usize>>,
    neq: usize,
}

impl Model {
    pub fn ndof(&self) -> usize {
        6 * self.nodes.len()
    }

    pub fn equations(&self) -> (Vec<Option<usize>>, usize) {
        let mut fixed = vec![false; self.ndof()];
        for &(node, f) in &self.fixed {
            for k in 0..6 {
                fixed[6 * node + k] |= f[k];
            }
        }
        let mut eq = vec![None; fixed.len()];
        let mut n = 0;
        for (i, &f) in fixed.iter().enumerate() {
            if !f {
                eq[i] = Some(n);
                n += 1;
            }
        }
        (eq, n)
    }

    /// Все вклады в глобальную матрицу жёсткости: (степень свободы, степень свободы, значение),
    /// обе половины симметричной матрицы.
    pub fn entries(&self, put: &mut dyn FnMut(usize, usize, f64)) -> Result<(), String> {
        let block = |dofs: &[usize], k: &dyn Fn(usize, usize) -> f64, put: &mut dyn FnMut(usize, usize, f64)| {
            for (a, &da) in dofs.iter().enumerate() {
                for (b, &db) in dofs.iter().enumerate() {
                    let v = k(a, b);
                    if v != 0.0 {
                        put(da, db, v);
                    }
                }
            }
        };
        let dofs_of = |nodes: &[usize]| -> Vec<usize> { nodes.iter().flat_map(|&n| (0..6).map(move |k| 6 * n + k)).collect() };
        for (i, s) in self.shells.iter().enumerate() {
            let c: Vec<V3> = s.nodes.iter().map(|&n| self.nodes[n]).collect();
            let k = shell::stiffness(&c, &s.mat).map_err(|e| format!("оболочка {i}: {e}"))?;
            let nd = 6 * s.nodes.len();
            block(&dofs_of(&s.nodes), &|a, b| k[a * nd + b], put);
        }
        for (i, b) in self.beams.iter().enumerate() {
            let k = beam_stiffness(&self.nodes, b).map_err(|e| format!("стержень {i}: {e}"))?;
            block(&dofs_of(&[b.n1, b.n2]), &|x, y| k[x][y], put);
        }
        for s in &self.springs {
            for k in 0..6 {
                if s.k[k] != 0.0 {
                    put(6 * s.node + k, 6 * s.node + k, s.k[k]);
                }
            }
        }
        for g in &self.spring_groups {
            let n = g.nodes.len();
            if g.k.len() != n * n {
                return Err("группа пружин: размер матрицы не совпадает с числом узлов".into());
            }
            let dofs: Vec<usize> = g.nodes.iter().map(|&n| 6 * n + 2).collect();
            block(&dofs, &|a, b| g.k[a * n + b], put);
        }
        for (i, d) in self.blocks.iter().enumerate() {
            let n = d.dofs.len();
            if d.k.len() != n * n {
                return Err(format!("суперэлемент {i}: размер матрицы не совпадает с числом степеней свободы"));
            }
            block(&d.dofs, &|a, b| d.k[a * n + b], put);
        }
        Ok(())
    }

    /// Сборка и разложение. Ошибка — вырожденная схема (механизм) или плохой элемент.
    pub fn factorize(&self) -> Result<Linear, String> {
        let (eq, neq) = self.equations();
        let mut trip: Vec<Triplet<usize, usize, f64>> = Vec::new();
        self.entries(&mut |a, b, v| {
            if let (Some(ia), Some(ib)) = (eq[a], eq[b]) {
                if ia >= ib {
                    trip.push(Triplet::new(ia, ib, v));
                }
            }
        })?;
        let a = SparseColMat::<usize, f64>::try_new_from_triplets(neq, neq, &trip).map_err(|e| format!("сборка: {e:?}"))?;
        let llt = a.sp_cholesky(Side::Lower).map_err(|_| "матрица не положительно определена: схема изменяема или есть незакреплённые степени свободы".to_string())?;
        Ok(Linear { llt, eq, neq })
    }

    /// Узловой вектор нагрузки (глобальные оси).
    pub fn load_vector(&self, lc: &LoadCase) -> Result<Vec<f64>, String> {
        let mut f = vec![0.0; self.ndof()];
        for &(node, v) in &lc.nodal {
            for k in 0..6 {
                f[6 * node + k] += v[k];
            }
        }
        for &(si, q) in &lc.pressure {
            let s = &self.shells[si];
            let c: Vec<V3> = s.nodes.iter().map(|&n| self.nodes[n]).collect();
            let r = shell::local_frame(&c)?;
            let d = if c.len() == 4 {
                cross([c[2][0] - c[0][0], c[2][1] - c[0][1], c[2][2] - c[0][2]], [c[3][0] - c[1][0], c[3][1] - c[1][1], c[3][2] - c[1][2]])
            } else {
                cross([c[1][0] - c[0][0], c[1][1] - c[0][1], c[1][2] - c[0][2]], [c[2][0] - c[0][0], c[2][1] - c[0][1], c[2][2] - c[0][2]])
            };
            let area = 0.5 * norm(d);
            let per = q * area / c.len() as f64;
            for &n in &s.nodes {
                for k in 0..3 {
                    f[6 * n + k] += per * r[2][k];
                }
            }
        }
        Ok(f)
    }
}

impl Linear {
    pub fn equations(&self) -> usize {
        self.neq
    }

    /// Перемещения по глобальному вектору нагрузки (закреплённые — 0).
    pub fn solve_vector(&self, f: &[f64]) -> Vec<f64> {
        let mut rhs = Mat::<f64>::zeros(self.neq, 1);
        for (g, e) in self.eq.iter().enumerate() {
            if let Some(i) = e {
                rhs[(*i, 0)] = f[g];
            }
        }
        let x = self.llt.solve(&rhs);
        self.eq.iter().map(|e| e.map_or(0.0, |i| x[(i, 0)])).collect()
    }

    /// Решение для нескольких векторов нагрузки сразу (блоками по 64).
    pub fn solve_many(&self, fs: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(fs.len());
        for chunk in fs.chunks(64) {
            let mut rhs = Mat::<f64>::zeros(self.neq, chunk.len());
            for (c, f) in chunk.iter().enumerate() {
                for (g, e) in self.eq.iter().enumerate() {
                    if let Some(i) = e {
                        rhs[(*i, c)] = f[g];
                    }
                }
            }
            let x = self.llt.solve(&rhs);
            for c in 0..chunk.len() {
                out.push(self.eq.iter().map(|e| e.map_or(0.0, |i| x[(i, c)])).collect());
            }
        }
        out
    }

    pub fn solve(&self, model: &Model, lc: &LoadCase) -> Result<Vec<f64>, String> {
        Ok(self.solve_vector(&model.load_vector(lc)?))
    }
}
