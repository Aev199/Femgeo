//! Плоские оболочки Q4 (мембрана Q4 + пластина MITC4) и T3 (CST + MITC3-подобный сдвиг),
//! изотропный материал. Локальные степени свободы узла: u, v, w, θx, θy; в глобальных — 6
//! с малой жёсткостью на вращение вокруг нормали (drilling).
//!
//! Источник: OpenCS `CSfea.Core/Elements/Shell4.cs`, `Shell3.cs`, `ShellGeometry.cs`,
//! `ShellElementMatrices.cs` (KLinear), коммит `387e157`.
//!
//! Знаки: κx = θy,x, κy = −θx,y, κxy = θy,y − θx,x; «+» M растягивает грань +z локальной
//! системы (как в [`crate::sp63`]).

pub type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Изотропная оболочка толщиной `h`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellMaterial {
    pub e: f64,
    pub nu: f64,
    pub h: f64,
}

impl ShellMaterial {
    /// A (мембрана), D (изгиб), As (сдвиг, κ = 5/6) — матрицы 3×3, 3×3, 2×2.
    fn abd(&self) -> ([[f64; 3]; 3], [[f64; 3]; 3], [[f64; 2]; 2]) {
        let (e, nu, h) = (self.e, self.nu, self.h);
        let q = e / (1.0 - nu * nu);
        let c = [[q, q * nu, 0.0], [q * nu, q, 0.0], [0.0, 0.0, q * (1.0 - nu) / 2.0]];
        let a = c.map(|r| r.map(|x| x * h));
        let d = c.map(|r| r.map(|x| x * h * h * h / 12.0));
        let g = e / (2.0 * (1.0 + nu));
        let s = 5.0 / 6.0 * g * h;
        (a, d, [[s, 0.0], [0.0, s]])
    }
}

/// Локальная система элемента (строки — оси x, y, n в глобальных координатах).
pub fn local_frame(c: &[V3]) -> Result<[V3; 3], String> {
    let (v1, v2) = match c.len() {
        3 => (sub(c[1], c[0]), sub(c[2], c[0])),
        4 => (sub(c[2], c[0]), sub(c[3], c[1])),
        _ => return Err("поддерживаются только 3 или 4 узла".into()),
    };
    let vn = cross(v1, v2);
    let nn = norm(vn);
    if nn < 1e-14 {
        return Err("вырожденный элемент (нулевая площадь)".into());
    }
    let vn = scale(vn, 1.0 / nn);
    let vx = sub(c[1], c[0]);
    let vx = sub(vx, scale(vn, dot(vx, vn)));
    let vx = scale(vx, 1.0 / norm(vx));
    Ok([vx, cross(vn, vx), vn])
}

fn inv2(j: [[f64; 2]; 2]) -> ([[f64; 2]; 2], f64) {
    let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    ([[j[1][1] / det, -j[0][1] / det], [-j[1][0] / det, j[0][0] / det]], det)
}

/// Матрицы B: мембрана, изгиб (3×5n), сдвиг (2×5n), плюс вес точки (detJ·w).
struct BSet {
    bm: Vec<[f64; 3]>,
    bb: Vec<[f64; 3]>,
    bs: Vec<[f64; 2]>,
    w: f64,
}

fn q4_shape(xi: f64, eta: f64) -> ([f64; 4], [[f64; 2]; 4]) {
    let n = [
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    ];
    let dn = [
        [-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)],
        [0.25 * (1.0 - eta), -0.25 * (1.0 + xi)],
        [0.25 * (1.0 + eta), 0.25 * (1.0 + xi)],
        [-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)],
    ];
    (n, dn)
}

/// (dN/dx, detJ, J, invJ); J[a][b] = Σ dN_i/dξ_a · x_i,b.
fn q4_jac(xy: &[[f64; 2]], dn: &[[f64; 2]; 4]) -> Result<([[f64; 2]; 4], f64, [[f64; 2]; 2], [[f64; 2]; 2]), String> {
    let mut j = [[0.0; 2]; 2];
    for a in 0..2 {
        for b in 0..2 {
            j[a][b] = (0..4).map(|i| dn[i][a] * xy[i][b]).sum();
        }
    }
    let (inv, det) = inv2(j);
    if det <= 0.0 {
        return Err("отрицательный якобиан (плохая форма элемента)".into());
    }
    let mut d = [[0.0; 2]; 4];
    for i in 0..4 {
        for k in 0..2 {
            d[i][k] = dn[i][0] * inv[k][0] + dn[i][1] * inv[k][1];
        }
    }
    Ok((d, det, j, inv))
}

fn bend_membrane(d: &[[f64; 2]], nd: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut bm = vec![[0.0; 3]; 5 * nd];
    let mut bb = vec![[0.0; 3]; 5 * nd];
    for i in 0..nd {
        let (nx, ny, c) = (d[i][0], d[i][1], 5 * i);
        bm[c][0] = nx;
        bm[c + 1][1] = ny;
        bm[c][2] = ny;
        bm[c + 1][2] = nx;
        bb[c + 4][0] = nx;
        bb[c + 3][1] = -ny;
        bb[c + 3][2] = -nx;
        bb[c + 4][2] = ny;
    }
    (bm, bb)
}

fn q4_bset(xy: &[[f64; 2]], xi: f64, eta: f64, w: f64) -> Result<BSet, String> {
    let (_, dn) = q4_shape(xi, eta);
    let (d, det, _, inv) = q4_jac(xy, &dn)?;
    let (bm, bb) = bend_membrane(&d, 4);
    // MITC4: естественные сдвиговые деформации в тай-точках A(0,−1), C(0,1), B(1,0), D(−1,0).
    let natural = |xs: f64, es: f64| -> Result<Vec<[f64; 2]>, String> {
        let (n, dns) = q4_shape(xs, es);
        let (ds, _, js, _) = q4_jac(xy, &dns)?;
        let mut phys = vec![[0.0; 2]; 20];
        for i in 0..4 {
            let c = 5 * i;
            phys[c + 2] = [ds[i][0], ds[i][1]];
            phys[c + 4][0] = n[i];
            phys[c + 3][1] = -n[i];
        }
        // γ_nat = J · γ_phys (J[a][b] = ∂x_b/∂ξ_a).
        Ok(phys.iter().map(|g| [js[0][0] * g[0] + js[0][1] * g[1], js[1][0] * g[0] + js[1][1] * g[1]]).collect())
    };
    let (ga, gc, gb, gd) = (natural(0.0, -1.0)?, natural(0.0, 1.0)?, natural(1.0, 0.0)?, natural(-1.0, 0.0)?);
    let mut bs = vec![[0.0; 2]; 20];
    for k in 0..20 {
        let nat = [
            0.5 * (1.0 - eta) * ga[k][0] + 0.5 * (1.0 + eta) * gc[k][0],
            0.5 * (1.0 - xi) * gd[k][1] + 0.5 * (1.0 + xi) * gb[k][1],
        ];
        // γ_phys = J⁻¹·γ_nat. В OpenCS здесь и при γ_nat = J·γ_phys матрица транспонирована
        // (как в исправленном автором `Jacobian`): на прямоугольниках вдоль осей разницы нет,
        // на скошенных элементах связываются не те компоненты. Здесь — ковариантные компоненты.
        bs[k] = [inv[0][0] * nat[0] + inv[0][1] * nat[1], inv[1][0] * nat[0] + inv[1][1] * nat[1]];
    }
    Ok(BSet { bm, bb, bs, w: w * det })
}

fn t3_bset(xy: &[[f64; 2]]) -> Result<BSet, String> {
    let dn = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
    let mut j = [[0.0; 2]; 2];
    for a in 0..2 {
        for b in 0..2 {
            j[a][b] = (0..3).map(|i| dn[i][a] * xy[i][b]).sum();
        }
    }
    let (inv, det) = inv2(j);
    if det <= 0.0 {
        return Err("отрицательный якобиан в T3".into());
    }
    let mut d = [[0.0; 2]; 3];
    for i in 0..3 {
        for k in 0..2 {
            d[i][k] = dn[i][0] * inv[k][0] + dn[i][1] * inv[k][1];
        }
    }
    let (bm, bb) = bend_membrane(&d, 3);
    // Сдвиг по касательным рёбер: γ_t = (w_j − w_i)/l + ½(θ_i + θ_j)·(−t_y, t_x), затем
    // псевдообратная матрица касательных (3×2 → 2×3).
    let edges = [(0usize, 1usize), (1, 2), (2, 0)];
    let mut t = [[0.0; 2]; 3];
    let mut rows = vec![[0.0; 3]; 15];
    for (ke, &(i, jn)) in edges.iter().enumerate() {
        let dv = [xy[jn][0] - xy[i][0], xy[jn][1] - xy[i][1]];
        let l = (dv[0] * dv[0] + dv[1] * dv[1]).sqrt();
        let (tx, ty) = (dv[0] / l, dv[1] / l);
        t[ke] = [tx, ty];
        let (ci, cj) = (5 * i, 5 * jn);
        rows[ci + 2][ke] += -1.0 / l;
        rows[cj + 2][ke] += 1.0 / l;
        for c in [ci, cj] {
            rows[c + 3][ke] += 0.5 * -ty;
            rows[c + 4][ke] += 0.5 * tx;
        }
    }
    // pinv(T) = (TᵀT)⁻¹Tᵀ.
    let mut tt = [[0.0; 2]; 2];
    for r in &t {
        for a in 0..2 {
            for b in 0..2 {
                tt[a][b] += r[a] * r[b];
            }
        }
    }
    let (ti, _) = inv2(tt);
    let mut bs = vec![[0.0; 2]; 15];
    for k in 0..15 {
        for a in 0..2 {
            bs[k][a] = (0..3).map(|ke| (ti[a][0] * t[ke][0] + ti[a][1] * t[ke][1]) * rows[k][ke]).sum();
        }
    }
    Ok(BSet { bm, bb, bs, w: 0.5 * det })
}

fn bsets(xy: &[[f64; 2]]) -> Result<Vec<BSet>, String> {
    match xy.len() {
        4 => {
            let g = 1.0 / 3f64.sqrt();
            [(-g, -g), (g, -g), (g, g), (-g, g)].iter().map(|&(a, b)| q4_bset(xy, a, b, 1.0)).collect()
        }
        3 => Ok(vec![t3_bset(xy)?]),
        _ => Err("поддерживаются только 3 или 4 узла".into()),
    }
}

/// Проекция узлов на плоскость элемента.
pub fn project(c: &[V3]) -> Result<([V3; 3], Vec<[f64; 2]>), String> {
    let r = local_frame(c)?;
    let xy = c.iter().map(|&p| {
        let d = sub(p, c[0]);
        [dot(d, r[0]), dot(d, r[1])]
    });
    Ok((r, xy.collect()))
}

/// Глобальная матрица жёсткости (6n × 6n, построчно).
pub fn stiffness(c: &[V3], m: &ShellMaterial) -> Result<Vec<f64>, String> {
    let (r, xy) = project(c)?;
    let n = c.len();
    let nl = 5 * n;
    let (a, d, s) = m.abd();
    let mut kl = vec![0.0; nl * nl];
    for b in bsets(&xy)? {
        for i in 0..nl {
            let (mi, bi, si) = (b.bm[i], b.bb[i], b.bs[i]);
            let am = [0, 1, 2].map(|p| (0..3).map(|q| a[p][q] * mi[q]).sum::<f64>());
            let db = [0, 1, 2].map(|p| (0..3).map(|q| d[p][q] * bi[q]).sum::<f64>());
            let ss = [0, 1].map(|p| (0..2).map(|q| s[p][q] * si[q]).sum::<f64>());
            for j in 0..nl {
                let (mj, bj, sj) = (b.bm[j], b.bb[j], b.bs[j]);
                let v = am[0] * mj[0] + am[1] * mj[1] + am[2] * mj[2]
                    + db[0] * bj[0] + db[1] * bj[1] + db[2] * bj[2]
                    + ss[0] * sj[0] + ss[1] * sj[1];
                kl[i * nl + j] += v * b.w;
            }
        }
    }
    // T (6n × 5n): перемещения — Rᵀ, повороты θx, θy — первые два столбца Rᵀ.
    let ng = 6 * n;
    let tm = |gi: usize, lj: usize| -> f64 {
        let (node_g, i) = (gi / 6, gi % 6);
        let (node_l, j) = (lj / 5, lj % 5);
        if node_g != node_l {
            return 0.0;
        }
        match (i < 3, j < 3) {
            (true, true) => r[j][i],
            (false, false) => r[j - 3][i - 3],
            _ => 0.0,
        }
    };
    let mut tk = vec![0.0; ng * nl];
    for gi in 0..ng {
        let node = gi / 6;
        for lj in 0..nl {
            tk[gi * nl + lj] = (5 * node..5 * node + 5).map(|p| tm(gi, p) * kl[p * nl + lj]).sum();
        }
    }
    let mut kg = vec![0.0; ng * ng];
    for gi in 0..ng {
        for gj in 0..ng {
            let node = gj / 6;
            kg[gi * ng + gj] = (5 * node..5 * node + 5).map(|p| tk[gi * nl + p] * tm(gj, p)).sum();
        }
    }
    let kref = (0..nl).map(|i| kl[i * nl + i].abs()).fold(0.0, f64::max) * 1e-6;
    let ez = r[2];
    for node in 0..n {
        let b = 6 * node + 3;
        for i in 0..3 {
            for j in 0..3 {
                kg[(b + i) * ng + b + j] += kref * ez[i] * ez[j];
            }
        }
    }
    Ok(kg)
}

/// Усилия оболочки в центре элемента в локальной системе: Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy.
pub fn forces(c: &[V3], m: &ShellMaterial, u_global: &[f64]) -> Result<([V3; 3], [f64; 8]), String> {
    let (r, xy) = project(c)?;
    let n = c.len();
    let mut ul = vec![0.0; 5 * n];
    for node in 0..n {
        let g = &u_global[6 * node..6 * node + 6];
        for k in 0..3 {
            ul[5 * node + k] = dot(r[k], [g[0], g[1], g[2]]);
        }
        for k in 0..2 {
            ul[5 * node + 3 + k] = dot(r[k], [g[3], g[4], g[5]]);
        }
    }
    let b = match n {
        4 => q4_bset(&xy, 0.0, 0.0, 1.0)?,
        _ => t3_bset(&xy)?,
    };
    let (a, d, s) = m.abd();
    let mut e = [0.0; 3];
    let mut k = [0.0; 3];
    let mut g = [0.0; 2];
    for i in 0..5 * n {
        for p in 0..3 {
            e[p] += b.bm[i][p] * ul[i];
            k[p] += b.bb[i][p] * ul[i];
        }
        for p in 0..2 {
            g[p] += b.bs[i][p] * ul[i];
        }
    }
    let mul3 = |mat: [[f64; 3]; 3], v: [f64; 3]| [0, 1, 2].map(|p| (0..3).map(|q| mat[p][q] * v[q]).sum::<f64>());
    let nn = mul3(a, e);
    let mm = mul3(d, k);
    Ok((r, [nn[0], nn[1], nn[2], mm[0], mm[1], mm[2], s[0][0] * g[0], s[1][1] * g[1]]))
}

/// Перевод усилий [Nx, Ny, Nxy, Mx, My, Mxy, Qx, Qy] из системы `from` в систему, где ось x —
/// проекция `axis` на плоскость элемента (оси армирования).
pub fn rotate_forces(from: &[V3; 3], f: &[f64; 8], axis: V3) -> [f64; 8] {
    let n = from[2];
    let ax = sub(axis, scale(n, dot(axis, n)));
    let ax = scale(ax, 1.0 / norm(ax));
    let (c, s) = (dot(ax, from[0]), dot(ax, from[1]));
    let t = |xx: f64, yy: f64, xy: f64| {
        (
            xx * c * c + yy * s * s + 2.0 * xy * c * s,
            xx * s * s + yy * c * c - 2.0 * xy * c * s,
            (yy - xx) * c * s + xy * (c * c - s * s),
        )
    };
    let (nx, ny, nxy) = t(f[0], f[1], f[2]);
    let (mx, my, mxy) = t(f[3], f[4], f[5]);
    [nx, ny, nxy, mx, my, mxy, f[6] * c + f[7] * s, -f[6] * s + f[7] * c]
}
