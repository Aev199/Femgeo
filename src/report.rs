//! Отчёт по выбранному варианту свайного поля: автономный HTML (исходные данные, методики
//! со ссылками на пункты норм, результаты, план свай SVG). Печать в PDF — средствами браузера.

use crate::fem::Model;
use crate::optimize::Evaluation;
use crate::service::{Input, SoilInput};

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

/// Число с запятой.
fn f(x: f64, d: usize) -> String {
    format!("{x:.d$}").replace('.', ",")
}

fn row(s: &mut String, k: &str, v: &str) {
    *s += &format!("<tr><th>{}</th><td>{}</td></tr>\n", esc(k), v);
}

fn ok(util: f64) -> &'static str {
    if util <= 1.0 { "" } else { " class=\"bad\"" }
}

fn kind_ru(k: &str) -> &str {
    match k {
        "permanent" => "постоянное",
        "long" => "длительное",
        "short" => "кратковременное",
        "bored" => "буровая",
        "driven" => "забивная",
        "end_bearing" => "свая-стойка",
        _ => k,
    }
}

/// План: контур, отверстия, сваи (цвет — доля от Fd), существующие сваи — крестики.
fn plan_svg(input: &Input, raft: &Model, v: &Evaluation) -> String {
    let pts = input.raft.outline.iter().chain(input.raft.holes.iter().flatten());
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for p in pts {
        x0 = x0.min(p[0]);
        y0 = y0.min(p[1]);
        x1 = x1.max(p[0]);
        y1 = y1.max(p[1]);
    }
    let m = 0.05 * (x1 - x0).max(y1 - y0).max(1.0);
    let (w, h) = (x1 - x0 + 2.0 * m, y1 - y0 + 2.0 * m);
    // Ось Y вверх: y → y0 + y1 − y.
    let path = |p: &[[f64; 2]]| {
        let mut s = String::new();
        for (i, q) in p.iter().enumerate() {
            s += &format!("{}{:.3},{:.3}", if i == 0 { "M" } else { "L" }, q[0], y0 + y1 - q[1]);
        }
        s + "Z"
    };
    let mut d = path(&input.raft.outline);
    for hole in &input.raft.holes {
        d += &path(hole);
    }
    let lw = w.max(h) / 600.0;
    let mut s = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{:.3} {:.3} {:.3} {:.3}\" width=\"100%\" style=\"max-height:70vh\">\n\
         <path d=\"{d}\" fill=\"#f4f4f0\" fill-rule=\"evenodd\" stroke=\"#333\" stroke-width=\"{:.4}\"/>\n",
        x0 - m, y0 - m, w, h, 2.0 * lw
    );
    let r = input.pile.d / 2.0;
    for e in &input.existing {
        let (x, y) = (e[0], y0 + y1 - e[1]);
        s += &format!(
            "<path d=\"M{:.3},{:.3}L{:.3},{:.3}M{:.3},{:.3}L{:.3},{:.3}\" stroke=\"#888\" stroke-width=\"{lw:.4}\"/>\n",
            x - r, y - r, x + r, y + r, x - r, y + r, x + r, y - r
        );
    }
    for (k, (&n, &nf)) in v.piles.iter().zip(&v.pile_n).enumerate() {
        let u = nf / input.pile.fd;
        let c = if u > 1.0 { "#c62828" } else if u > 0.9 { "#ef8f00" } else { "#2e7d32" };
        let (x, y) = (raft.nodes[n][0], y0 + y1 - raft.nodes[n][1]);
        s += &format!(
            "<circle cx=\"{x:.3}\" cy=\"{y:.3}\" r=\"{r:.3}\" fill=\"{c}\"><title>№{} N = {} кН</title></circle>\n",
            k + 1,
            f(nf, 0)
        );
    }
    s + "</svg>"
}

pub fn html(input: &Input, raft: &Model, vs: &[Evaluation], i: usize, reference: bool) -> String {
    let v = &vs[i];
    let p = &input.pile;
    let r = &input.raft;
    let mut s = String::from(
        "<!DOCTYPE html><html lang=\"ru\"><head><meta charset=\"utf-8\"><title>Femgeo — отчёт</title><style>\
         body{font:14px/1.45 sans-serif;max-width:980px;margin:24px auto;padding:0 16px;color:#222}\
         h1{font-size:22px}h2{font-size:17px;border-bottom:1px solid #ccc;padding-bottom:2px;margin-top:28px}\
         table{border-collapse:collapse;margin:8px 0}th,td{border:1px solid #ccc;padding:3px 8px;text-align:left;vertical-align:top}\
         th{background:#f3f3f3;font-weight:normal}.bad{color:#c62828;font-weight:bold}td.n{text-align:right}\
         @media print{h2{break-after:avoid}svg{break-inside:avoid}}</style></head><body>\n",
    );
    let name = if reference && i == 0 { "существующая расстановка".to_string() } else { format!("вариант {}", if reference { i } else { i + 1 }) };
    s += &format!("<h1>Свайное поле под плитой: {}</h1>\n", esc(&name));
    s += &format!(
        "<p>Итог: {} свай, расход арматуры плиты {} кг/м³ — {}.</p>\n",
        v.piles.len(),
        f(v.kg_m3, 1),
        if v.feasible { "все проверки выполнены" } else { "<span class=\"bad\">проверки не выполнены</span>" }
    );

    s += "<h2>1. Исходные данные</h2>\n<table>\n";
    let area = |q: &[[f64; 2]]| {
        let n = q.len();
        (0..n).map(|k| q[k][0] * q[(k + 1) % n][1] - q[(k + 1) % n][0] * q[k][1]).sum::<f64>().abs() / 2.0
    };
    let a = area(&r.outline) - r.holes.iter().map(|h| area(h)).sum::<f64>();
    row(&mut s, "Плита: площадь / толщина", &format!("{} м² / {} м, отверстий: {}", f(a, 1), f(r.h, 2), r.holes.len()));
    row(&mut s, "Низ плиты (головы свай)", &format!("{} м", f(r.z_bottom, 2)));
    row(&mut s, "Бетон / арматура плиты", &format!("{} / {}, ⌀{} мм, a = {} мм", esc(&r.concrete), esc(&r.rebar), f(r.bar_d * 1e3, 0), f(r.cover * 1e3, 0)));
    row(&mut s, "Сетка МКЭ", &format!("{} м, узлов {}, элементов {}", f(r.mesh, 2), raft.nodes.len(), raft.shells.len()));
    row(
        &mut s,
        "Свая",
        &format!(
            "{}, d = {} м, L = {} м, Fd = {} кН; {} / {}, {}⌀{} мм, l0 = {} м",
            kind_ru(&p.kind),
            f(p.d, 2),
            f(p.length, 1),
            f(p.fd, 0),
            esc(&p.concrete),
            esc(&p.rebar),
            p.bars,
            f(p.bar_d * 1e3, 0),
            f(p.l0, 1)
        ),
    );
    let spacing = |o: Option<f64>| o.map(|x| format!("{} м", f(x, 2))).unwrap_or_else(|| "по СП 24".into());
    row(&mut s, "Шаг свай / отступ от края", &format!("{} / {}", spacing(p.min_spacing), spacing(p.min_edge)));
    row(&mut s, "Расстановка", if input.layout == "free" { "свободная" } else { "рядами" });
    row(&mut s, "Предельная осадка / относительная разность", &format!("{} м / {}", f(input.s_max, 3), f(input.ds_rel, 4)));
    let crack = |x: f64| if x > 0.0 { format!("{} мм", f(x, 2)) } else { "не проверяется".into() };
    row(&mut s, "acrc,ult продолж. / непродолж.", &format!("{} / {}", crack(input.acrc_long), crack(input.acrc_short)));
    row(&mut s, "Продавливание", if input.punching { "проверяется" } else { "не проверяется" });
    let soil = match &input.soil {
        SoilInput::Geobase { path, model } => {
            let file = std::path::Path::new(path).file_name().map(|x| x.to_string_lossy().into_owned()).unwrap_or_default();
            format!("geobase «{}», модель {} (колонка ближайшей скважины)", esc(&file), esc(model.as_deref().unwrap_or("по умолчанию")))
        }
        SoilInput::Layers { top, layers } => format!("слои вручную, верх {} м, слоёв {}", f(*top, 2), layers.len()),
    };
    row(&mut s, "Грунт", &soil);
    row(&mut s, "Существующие сваи (ориентир)", &input.existing.len().to_string());
    s += "</table>\n<p>Загружения:</p>\n<table><tr><th>Имя</th><th>Вид</th><th>γf</th><th>ψдл</th><th>q, кПа</th><th>Сосредоточенные, кН</th><th>Значения</th></tr>\n";
    for c in &input.cases {
        s += &format!(
            "<tr><td>{}</td><td>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td>{}</td></tr>\n",
            esc(&c.name),
            kind_ru(&c.kind),
            f(c.gamma_f, 2),
            if c.kind == "short" { f(c.psi_long, 2) } else { "—".into() },
            f(c.q, 2),
            if c.points.is_empty() { "—".into() } else { format!("{} шт., Σ {}", c.points.len(), f(c.points.iter().map(|q| q.n).sum::<f64>(), 0)) },
            if c.design { "расчётные" } else { "нормативные" }
        );
    }
    s += "</table>\n";
    if v.soil_3d {
        s += "<p><b>Вариант пересчитан на 3D-основании</b> (этап 3): податливость свай с взаимовлиянием — упругий слоистый массив \
              в Code_Aster (гексаэдры HEXA20, сваи — объёмные стойки равной площади, полное сцепление с грунтом; контакт плиты с \
              грунтом не учитывается). Методика жёсткости свай ниже к этому варианту не относится. Несущая способность сваи по грунту — \
              заданное Fd (не по 3D-модели); перегрузка краевых и угловых свай закономерна: сваи держат Fu, избыток перераспределяется, доля таких свай \
              не ограничивается (ограничение 10 % — только для расчёта по СП 24 / условному фундаменту).</p>\n";
        if !v.nl_factor.is_empty() {
            let (lo, hi) = v.nl_factor.iter().fold((f64::MAX, 0.0f64), |a, &d| (a.0.min(d), a.1.max(d)));
            s += &format!("<p>Основание упругопластическое (Кулон–Мор, ψ = 0; природные напряжения K0 = 1 − sin φ), гибридный метод: \
                  нелинейность — по кривой одиночной сваи (Code_Aster), взаимовлияние — упругое; отношение осадок к упругим {}–{}.</p>\n", f(lo, 2), f(hi, 2));
        }
    }

    s += "<h2>2. Методики</h2>\n<ul>\n\
          <li>Сочетания нагрузок — СП 20.13330.2016, п. 6 (основные), коэффициенты п. 6.4; ПС1 — расчётные, \
          ПС2 — нормативные значения, длительное ПС2 — с длительными частями кратковременных. Особые сочетания не формируются.</li>\n\
          <li>Плита — оболочка МКЭ (MITC4/T3) на сваях-пружинах; жёсткость свай и осадка — СП 24.13330.2021, п. 7.4 \
          (одиночная свая — п. 7.4.2, взаимовлияние — п. 7.4.4).</li>\n\
          <li>Несущая способность сваи по грунту — заданное Fd (с учётом γk); по материалу — СП 63.13330, внецентренное сжатие \
          со случайным эксцентриситетом (l0).</li>\n\
          <li>Расстояния между сваями — СП 24.13330.2021, п. 8.13; отступ от края 0,5d + 0,1 м — допущение (п. 8.3 числа не задаёт), если не задан.</li>\n\
          <li>Армирование плиты — подбор по СП 63.13330 (перенос OpenCS) по огибающей усилий; трещины — п. 8.2.6 \
          (продолжительное — длительное сочетание, φ1 = 1,4; непродолжительное — полные сочетания ПС2).</li>\n\
          <li>Продавливание плиты сваями — СП 63.13330, пп. 8.1.46–8.1.47, без поперечной арматуры и без моментов; \
          у края — доля замкнутого контура внутри плиты (упрощение).</li>\n\
          </ul>\n";

    s += "<h2>3. Результаты</h2>\n<table>\n";
    let lim = |name: &str, val: String, util: f64| format!("<tr><th>{}</th><td{}>{}</td><td{} class=n>{}</td></tr>\n", name, ok(util), val, ok(util), f(util, 2));
    s += "<tr><th>Проверка</th><th>Значение</th><th>Кисп</th></tr>\n";
    s += &lim("Свая по грунту", format!("Nmax = {} кН", f(v.pile_n.iter().cloned().fold(0.0, f64::max), 0)), v.util_ground);
    s += &lim("Свая по материалу", "—".into(), v.util_material);
    s += &lim("Осадка", format!("{} м", f(v.s_max, 3)), v.s_max / input.s_max);
    s += &lim("Относительная разность осадок", f(v.ds_rel, 5), v.ds_rel / input.ds_rel);
    if input.punching {
        s += &lim("Продавливание", "—".into(), v.util_punch);
    }
    s += "</table>\n";
    if v.design_impossible {
        s += "<p class=\"bad\">Подбор арматуры плиты невозможен хотя бы в одном элементе.</p>\n";
    }
    if v.shift.is_finite() {
        s += &format!("<p>Среднее смещение от существующих свай: {} м.</p>\n", f(v.shift, 2));
    }
    let mx = |k: usize| v.rebar.iter().map(|a| a[k] * 1e4).fold(0.0, f64::max);
    s += &format!(
        "<p>Наибольшая требуемая арматура плиты, см²/м: верх X {}, верх Y {}, низ X {}, низ Y {}; расход {} кг/м³ (теоретический).</p>\n",
        f(mx(0), 2),
        f(mx(1), 2),
        f(mx(2), 2),
        f(mx(3), 2),
        f(v.kg_m3, 1)
    );

    s += "<h2>4. План свай</h2>\n<p>Цвет: зелёный — N ≤ 0,9Fd, оранжевый — до Fd, красный — больше Fd; крестики — существующие сваи.</p>\n";
    s += &plan_svg(input, raft, v);
    s += "\n<h2>5. Ведомость свай</h2>\n<table><tr><th>№</th><th>X, м</th><th>Y, м</th><th>N, кН</th><th>N/Fd</th></tr>\n";
    for (k, (&n, &nf)) in v.piles.iter().zip(&v.pile_n).enumerate() {
        let u = nf / p.fd;
        s += &format!(
            "<tr><td>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n{}>{}</td></tr>\n",
            k + 1,
            f(raft.nodes[n][0], 3),
            f(raft.nodes[n][1], 3),
            f(nf, 1),
            if u > 1.0 { " style=\"color:#c62828\"" } else { "" },
            f(u, 2)
        );
    }
    s += "</table>\n";
    if vs.len() > 1 {
        s += "<h2>6. Сравнение вариантов</h2>\n<table><tr><th></th><th>Свай</th><th>кг/м³</th><th>Кисп грунт</th><th>Кисп матер.</th><th>Осадка, м</th><th>Допустим</th></tr>\n";
        for (k, e) in vs.iter().enumerate() {
            let label = if reference && k == 0 { "существ.".into() } else { (if reference { k } else { k + 1 }).to_string() };
            s += &format!(
                "<tr{}><td>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td class=n>{}</td><td>{}</td></tr>\n",
                if k == i { " style=\"background:#eef5ff\"" } else { "" },
                label,
                e.piles.len(),
                f(e.kg_m3, 1),
                f(e.util_ground, 2),
                f(e.util_material, 2),
                f(e.s_max, 3),
                if e.feasible { "да" } else { "нет" }
            );
        }
        s += "</table>\n";
    }
    s + "<p style=\"color:#888;margin-top:32px\">Сформировано Femgeo.</p></body></html>\n"
}
