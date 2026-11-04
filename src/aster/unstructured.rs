//! Неструктурированная сетка основания (Gmsh, тетраэдры второго порядка TETRA10): размер
//! элемента у ствола сваи — `h_pile·d`, в зоне свай — `h`, к границам модели растёт до `h_far`
//! (поле Box с переходом на ширине боковой зоны). Сваи — призмы квадратного сечения равной площади (E = EA/A по сетке), над
//! отметкой голов — оголовок высотой d; слои — секущими плоскостями по колонке в центре поля
//! (изменчивость грунта в плане не учитывается — допущение). Сетка строится скриптом Python с
//! модулем gmsh (`Options::python`, переменная `FEMGEO_PYTHON`).

use super::{layer_at, nu_geostatic, Material, Model3, Options};
use crate::piles::{Pile, SoilColumn};
use serde_json::json;

const SCRIPT: &str = r#"
import sys, json, math, gmsh
p = json.load(open(sys.argv[1]))
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.NumThreads", p["threads"])
occ = gmsh.model.occ
X0, X1, Y0, Y1, ZB = p["box"]
H, ST = p["head"], p["stub"]
box = occ.addBox(X0, Y0, ZB, X1 - X0, Y1 - Y0, H - ZB)
# Сваи — призмы квадратного сечения равной площади (b = d·√π/2): без кривых поверхностей.
B = [d * math.sqrt(math.pi) / 2 for x, y, d, L in p["piles"]]
cyl = [occ.addBox(x - b / 2, y - b / 2, H - L, b, b, L + ST) for (x, y, d, L), b in zip(p["piles"], B)]
rects = [occ.addRectangle(X0, Y0, z, X1 - X0, Y1 - Y0) for z in p["planes"]]
occ.fragment([(3, box)], [(3, c) for c in cyl] + [(2, r) for r in rects])
occ.synchronize()
planes = p["planes"]
def interval(z):
    return sum(1 for q in planes if q > z)
groups = {}
piece = {}
pile_vols = []
heads = [[] for _ in p["piles"]]
for dim, v in gmsh.model.getEntities(3):
    xa, ya, za, xb, yb, zb = gmsh.model.getBoundingBox(3, v)
    cx, cy, cz = (xa + xb) / 2, (ya + yb) / 2, (za + zb) / 2
    m = None
    for k, (x, y, d, L) in enumerate(p["piles"]):
        if xb - xa < 1.01 * B[k] and abs(cx - x) < 0.01 * d and abs(cy - y) < 0.01 * d:
            pile_vols.append(v)
            if cz > H:
                m = "PS%d" % k
                heads[k].append(v)
            else:
                m = "PL%d_%d" % (k, interval(cz))
                piece[v] = k
            break
    if m is None:
        m = "MAT%d" % p["mat_soil"][interval(cz)]
    groups.setdefault(m, []).append(v)
# Размеры: у свай — Threshold от расстояния до поверхностей свай, зона — Box с переходом.
surf = sorted({abs(t) for v in pile_vols for _, t in gmsh.model.getBoundary([(3, v)], oriented=False)})
f = gmsh.model.mesh.field
fd = f.add("Distance"); f.setNumbers(fd, "SurfacesList", surf); f.setNumber(fd, "Sampling", 12)
ft = f.add("Threshold"); f.setNumber(ft, "InField", fd)
f.setNumber(ft, "SizeMin", p["h_pile"]); f.setNumber(ft, "SizeMax", p["h_far"])
# Рост размера от ствола: `grow` на метр.
f.setNumber(ft, "DistMin", 0.0); f.setNumber(ft, "DistMax", (p["h_far"] - p["h_pile"]) / p["grow"])
fb = f.add("Box"); zx0, zx1, zy0, zy1, ztip = p["zone"]
f.setNumber(fb, "VIn", p["h"]); f.setNumber(fb, "VOut", p["h_far"])
f.setNumber(fb, "XMin", zx0); f.setNumber(fb, "XMax", zx1); f.setNumber(fb, "YMin", zy0); f.setNumber(fb, "YMax", zy1)
f.setNumber(fb, "ZMin", ztip); f.setNumber(fb, "ZMax", H + ST); f.setNumber(fb, "Thickness", p["transition"])
fm = f.add("Min"); f.setNumbers(fm, "FieldsList", [ft, fb]); f.setAsBackgroundMesh(fm)
for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
             ("Mesh.Algorithm", 6), ("Mesh.Algorithm3D", 10), ("General.NumThreads", 1), ("Mesh.ElementOrder", 2), ("Mesh.SecondOrderLinear", 1)):
    gmsh.option.setNumber(k, v)
gmsh.model.mesh.generate(3)
tags, xyz, _ = gmsh.model.mesh.getNodes()
out = ["COOR_3D"]
out += ["N%d %.6f %.6f %.6f" % (t, xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]) for i, t in enumerate(tags)]
out += ["FINSF", "TETRA10"]
perm = [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
ne = 0
xyz_of = {t: (xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]) for i, t in enumerate(tags)}
vol = [0.0] * len(p["piles"])
def tetvol(q):
    a, b, c, d = (xyz_of[q[0]], xyz_of[q[1]], xyz_of[q[2]], xyz_of[q[3]])
    u = [b[i] - a[i] for i in range(3)]; w = [c[i] - a[i] for i in range(3)]; z = [d[i] - a[i] for i in range(3)]
    return abs(u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) + u[2] * (w[0] * z[1] - w[1] * z[0])) / 6.0
gm = {}
for m, vols in sorted(groups.items()):
    for v in vols:
        et, etags, enodes = gmsh.model.mesh.getElements(3, v)
        for t, tg, nd in zip(et, etags, enodes):
            if t != 11:
                raise SystemExit("unexpected element type %d" % t)
            nd = list(nd)
            for k, e in enumerate(tg):
                q = nd[10 * k:10 * k + 10]
                # Формат ASTER: строки не длиннее 80 символов — узлы на двух строках.
                out.append("M%d " % e + " ".join("N%d" % q[j] for j in perm[:5]))
                out.append(" ".join("N%d" % q[j] for j in perm[5:]))
                gm.setdefault(m, []).append(e)
                if v in piece:
                    vol[piece[v]] += tetvol(q)
                ne += 1
out.append("FINSF")
for m, es in sorted(gm.items()):
    out.append("GROUP_MA\n%s" % m)
    out += [" ".join("M%d" % e for e in es[i:i + 8]) for i in range(0, len(es), 8)]
    out.append("FINSF")
for k, vols in enumerate(heads):
    top = []
    for v in vols:
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            s = abs(s)
            if abs(gmsh.model.getBoundingBox(2, s)[2] - (H + ST)) < 1e-6 and abs(gmsh.model.getBoundingBox(2, s)[5] - (H + ST)) < 1e-6:
                top += list(gmsh.model.mesh.getNodes(2, s, includeBoundary=True)[0])
    if len(top) < 3:
        raise SystemExit("pile %d: head nodes not found" % k)
    top = sorted(set(top))
    out.append("GROUP_NO\nH%d" % k)
    out += [" ".join("N%d" % n for n in top[i:i + 8]) for i in range(0, len(top), 8)]
    out.append("FINSF")
out.append("FIN")
open(p["out"], "w").write("\n".join(out) + "\n")
area = [vol[k] / L for k, (x, y, d, L) in enumerate(p["piles"])]
json.dump({"nodes": len(tags), "elements": ne, "area": area}, open(p["out"] + ".json", "w"))
gmsh.finalize()
"#;

pub fn model(piles: &[Pile], soil: &dyn Fn(f64, f64) -> SoilColumn, o: &Options) -> Result<Model3, String> {
    if piles.is_empty() {
        return Err("нет свай".into());
    }
    let head = piles[0].head;
    if piles.iter().any(|p| (p.head - head).abs() > 1e-6) {
        return Err("3D-модель: головы свай на разных отметках не поддерживаются".into());
    }
    let (mut x0, mut x1, mut y0, mut y1) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for p in piles {
        (x0, x1, y0, y1) = (x0.min(p.x - p.d / 2.0), x1.max(p.x + p.d / 2.0), y0.min(p.y - p.d / 2.0), y1.max(p.y + p.d / 2.0));
    }
    let tip = piles.iter().map(|p| p.head - p.length).fold(f64::MAX, f64::min);
    let bottom = tip - o.below;
    let dmax = piles.iter().map(|p| p.d).fold(0.0, f64::max);
    let col = soil(0.5 * (x0 + x1), 0.5 * (y0 + y1));
    // Плоскости слоёв внутри модели (сверху вниз) и интервалы между ними.
    let mut planes: Vec<f64> = col.layers.iter().map(|l| l.bottom).filter(|&z| z < head - 1e-6 && z > bottom + 1e-6).collect();
    planes.sort_by(|a, b| b.total_cmp(a));
    planes.dedup_by(|a, b| (*a - *b).abs() < 0.05);
    let mut bounds = vec![head];
    bounds.extend(&planes);
    bounds.push(bottom);
    let mut materials: Vec<Material> = vec![];
    let mut key = |m: Material| -> usize {
        if let Some(i) = materials.iter().position(|q| *q == m) {
            return i;
        }
        materials.push(m);
        materials.len() - 1
    };
    let mut layers = vec![];
    for w in bounds.windows(2) {
        let l = layer_at(&col, 0.5 * (w[0] + w[1])).ok_or("грунт: нет слоёв")?;
        if !(l.e > 0.0) || !(l.nu >= 0.0 && l.nu < 0.5) {
            return Err(format!("грунт на отметке {:.2}: E {}, ν {}", 0.5 * (w[0] + w[1]), l.e, l.nu));
        }
        layers.push(l);
    }
    let mat_soil: Vec<usize> = layers.iter().map(|l| key(Material { e: l.e, nu: l.nu, gamma: l.gamma, phi: l.phi, c: l.c, pile: false, e0: l.e, nu0: nu_geostatic(l) })).collect();
    std::fs::create_dir_all(&o.workdir).map_err(|e| format!("{}: {e}", o.workdir.display()))?;
    let wd = o.workdir.canonicalize().map_err(|e| e.to_string())?;
    let out = wd.join("model.mail");
    let _ = std::fs::remove_file(&out);
    let box_ = [x0 - o.lateral, x1 + o.lateral, y0 - o.lateral, y1 + o.lateral, bottom];
    let args = json!({
        "piles": piles.iter().map(|p| [p.x, p.y, p.d, p.length]).collect::<Vec<_>>(),
        "head": head, "stub": dmax, "box": box_, "planes": planes,
        "zone": [x0 - dmax, x1 + dmax, y0 - dmax, y1 + dmax, tip - dmax],
        "h_pile": o.h_pile * piles.iter().map(|p| super::square_side(p.d)).fold(f64::MAX, f64::min), "h": o.h_zone, "h_far": o.h_far, "grow": o.grow, "transition": o.lateral,
        "mat_soil": mat_soil,
        "out": out, "threads": o.ncpus,
    });
    std::fs::write(wd.join("mesh.json"), args.to_string()).map_err(|e| e.to_string())?;
    std::fs::write(wd.join("mesh.py"), SCRIPT).map_err(|e| e.to_string())?;
    let r = std::process::Command::new(&o.python)
        .arg("-I")
        .arg(wd.join("mesh.py"))
        .arg(wd.join("mesh.json"))
        .current_dir(&wd)
        .output()
        .map_err(|e| format!("{}: {e} (нужен Python с модулем gmsh)", o.python.display()))?;
    if !r.status.success() {
        let err = String::from_utf8_lossy(&r.stderr);
        return Err(format!("Gmsh: {}", err.lines().rev().find(|l| !l.trim().is_empty()).unwrap_or("ошибка")));
    }
    let stats: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(out.with_extension("mail.json")).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let mail = std::fs::read_to_string(&out).map_err(|e| e.to_string())?;
    // Сваи: E = EA / A сетки (прямые рёбра TETRA10 — многоугольник, площадь меньше круга).
    let areas: Vec<f64> = stats["area"].as_array().map(|a| a.iter().filter_map(|v| v.as_f64()).collect()).unwrap_or_default();
    if areas.len() != piles.len() || areas.iter().any(|a| !(*a > 0.0)) {
        return Err("Gmsh: не определена площадь сечения свай".into());
    }
    let mut rename: Vec<(String, usize)> = vec![];
    for (k, p) in piles.iter().enumerate() {
        let e = p.ea / areas[k];
        for (i, l) in layers.iter().enumerate() {
            rename.push((format!("PL{k}_{i}"), key(Material { e, nu: 0.2, gamma: l.gamma, phi: 0.0, c: 0.0, pile: true, e0: l.e, nu0: nu_geostatic(l) })));
        }
        rename.push((format!("PS{k}"), key(Material { e, nu: 0.2, gamma: 0.0, phi: 0.0, c: 0.0, pile: true, e0: e, nu0: 0.2 })));
    }
    let mut groups: Vec<Vec<String>> = vec![vec![]; materials.len()];
    for (m, g) in groups.iter_mut().enumerate() {
        if mail.contains(&format!("\nGROUP_MA\nMAT{m}\n")) {
            g.push(format!("MAT{m}"));
        }
    }
    for (name, m) in rename {
        if mail.contains(&format!("\nGROUP_MA\n{name}\n")) {
            groups[m].push(name);
        }
    }
    Ok(Model3 {
        mail,
        materials,
        groups,
        bounds: box_,
        heads: piles.len(),
        nodes: stats["nodes"].as_u64().unwrap_or(0) as usize,
        elements: stats["elements"].as_u64().unwrap_or(0) as usize,
        quadratic: true,
    })
}
