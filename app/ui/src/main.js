// Окно Femgeo: исходные данные → расчёт фронта вариантов → план и таблицы.
import { call, onProgress, pickFile } from './api.js';
import { Plan, color, GRADIENT } from './plan.js';

const $ = (id) => document.getElementById(id);
const plan = new Plan($('plan'));
const state = { result: null, current: -1, mode: 'piles', soil: 'layers', geobase: null, layout: 'rows' };

const CONCRETE = ['B15', 'B20', 'B25', 'B30', 'B35', 'B40', 'B45', 'B50', 'B55', 'B60'];
const REBAR = ['A240', 'A400', 'A500'];
for (const s of document.querySelectorAll('select.concrete')) s.innerHTML = CONCRETE.map((c) => `<option>${c}</option>`).join('');
for (const s of document.querySelectorAll('select.rebar')) s.innerHTML = REBAR.map((c) => `<option>${c}</option>`).join('');

const fmt = (x, d) => Number(x).toFixed(d).replace('.', ',');

function status(text, error = false) {
  $('status').textContent = text;
  $('status').parentElement.classList.toggle('error', error);
}

// ── Таблицы чисел в текстовых полях ──

const rows = (text) => text.split('\n').map((l) => l.trim()).filter((l) => l)
  .map((l) => l.replace(/;/g, ' ').replace(/,(?=\d)/g, '.').split(/\s+/).map(Number));
const writeRows = (rs) => rs.map((r) => r.join(' ')).join('\n');

function table(id, width, what) {
  const area = $(id);
  const rs = rows(area.value);
  const bad = rs.some((r) => r.length !== width || r.some((x) => !Number.isFinite(x)));
  area.classList.toggle('invalid', bad);
  if (bad) throw new Error(`${what}: в каждой строке ${width} числа`);
  return rs;
}

function rings(id) {
  const blocks = $(id).value.split(/\n\s*\n/).map((b) => b.trim()).filter((b) => b);
  return blocks.map((b) => {
    const rs = rows(b);
    if (rs.length < 3 || rs.some((r) => r.length !== 2 || r.some((x) => !Number.isFinite(x)))) {
      $(id).classList.add('invalid');
      throw new Error('Отверстие: не менее трёх точек X Y');
    }
    return rs;
  });
}

const num = (id) => {
  const v = $(id).valueAsNumber;
  if (!Number.isFinite(v)) throw new Error(`Не заполнено поле «${$(id).closest('label').firstChild.textContent.trim()}»`);
  return v;
};
const optional = (id) => (Number.isFinite($(id).valueAsNumber) ? $(id).valueAsNumber : null);

// ── Форма ↔ исходные данные ──

function fill(input) {
  const r = input.raft, p = input.pile;
  $('outline').value = writeRows(r.outline);
  $('holes').value = r.holes.map(writeRows).join('\n\n');
  $('raft-h').value = r.h; $('raft-z').value = r.z_bottom; $('raft-mesh').value = r.mesh;
  $('raft-cover').value = r.cover; $('raft-concrete').value = r.concrete; $('raft-rebar').value = r.rebar;
  $('raft-bar').value = Math.round(r.bar_d * 1000);
  const KIND = { permanent: 'П', long: 'Д', short: 'К' };
  $('cases').value = input.cases.map((c) => [KIND[c.kind], c.gamma_f, c.psi_long ?? 0, c.q, c.design ? 'р' : ''].join(' ').trim()).join('\n');
  $('points').value = writeRows(input.cases.flatMap((c, i) => c.points.map((x) => [i + 1, x.x, x.y, x.n])));
  $('pile-d').value = p.d; $('pile-l').value = p.length; $('pile-fd').value = p.fd; $('pile-kind').value = p.kind;
  $('pile-concrete').value = p.concrete; $('pile-rebar').value = p.rebar; $('pile-bars').value = p.bars;
  $('pile-bar').value = Math.round(p.bar_d * 1000); $('pile-cover').value = p.cover; $('pile-l0').value = p.l0;
  $('pile-spacing').value = p.min_spacing ?? ''; $('pile-edge').value = p.min_edge ?? '';
  if (input.soil.source === 'layers') {
    $('soil-top').value = input.soil.top;
    $('layers').value = writeRows(input.soil.layers.map((l) => [l.bottom, l.e, l.nu]));
  }
  $('s-max').value = Math.round(input.s_max * 100); $('ds-rel').value = input.ds_rel;
  $('acrc-long').value = input.acrc_long; $('acrc-short').value = input.acrc_short; $('punching').checked = input.punching;
  $('existing').value = writeRows(input.existing);
  state.layout = input.layout || 'rows';
  for (const x of document.querySelectorAll('#layout button')) x.classList.toggle('current', x.dataset.layout === state.layout);
}

function collectCases() {
  const KIND = { 'П': 'permanent', 'Д': 'long', 'К': 'short' };
  const lines = $('cases').value.split('\n').map((l) => l.trim()).filter((l) => l);
  const cases = lines.map((l, i) => {
    const [k, ...rest] = l.split(/\s+/);
    const design = rest.length === 4 && /^[рp]$/i.test(rest[3]);
    if (design) rest.pop();
    const v = rest.map((x) => Number(x.replace(',', '.')));
    if (!KIND[k.toUpperCase()] || v.length !== 3 || v.some((x) => !Number.isFinite(x))) {
      $('cases').classList.add('invalid');
      throw new Error(`Загружение ${i + 1}: «вид γf ψдл q [р]», вид — П, Д или К`);
    }
    return { name: `${i + 1}`, kind: KIND[k.toUpperCase()], gamma_f: v[0], psi_long: v[1], q: v[2], design, points: [] };
  });
  if (!cases.length) throw new Error('Задайте хотя бы одно загружение');
  for (const [n, x, y, f] of table('points', 4, 'Сосредоточенные нагрузки')) {
    if (!cases[n - 1]) throw new Error(`Сосредоточенная нагрузка: нет загружения ${n}`);
    cases[n - 1].points.push({ x, y, n: f });
  }
  return cases;
}

function collect() {
  for (const a of document.querySelectorAll('textarea')) a.classList.remove('invalid');
  const outline = table('outline', 2, 'Контур плиты');
  if (outline.length < 3) throw new Error('Контур плиты: не менее трёх точек');
  let soil;
  if (state.soil === 'geobase') {
    if (!state.geobase) throw new Error('Грунт: откройте объект GeoBase');
    soil = { source: 'geobase', path: state.geobase.path, model: state.geobase.model };
  } else {
    const ls = table('layers', 3, 'Слои грунта');
    if (!ls.length) throw new Error('Грунт: задайте хотя бы один слой');
    soil = { source: 'layers', top: num('soil-top'), layers: ls.map(([bottom, e, nu]) => ({ bottom, e, nu })) };
  }
  return {
    raft: {
      outline, holes: rings('holes'), h: num('raft-h'), z_bottom: num('raft-z'), mesh: num('raft-mesh'),
      concrete: $('raft-concrete').value, rebar: $('raft-rebar').value, cover: num('raft-cover'), bar_d: num('raft-bar') / 1000,
    },
    cases: collectCases(),
    pile: {
      d: num('pile-d'), length: num('pile-l'), fd: num('pile-fd'), kind: $('pile-kind').value,
      concrete: $('pile-concrete').value, rebar: $('pile-rebar').value, bars: num('pile-bars'), bar_d: num('pile-bar') / 1000,
      cover: num('pile-cover'), l0: num('pile-l0'), min_spacing: optional('pile-spacing'), min_edge: optional('pile-edge'),
    },
    soil,
    s_max: num('s-max') / 100,
    ds_rel: num('ds-rel'),
    layout: state.layout,
    acrc_long: num('acrc-long'),
    acrc_short: num('acrc-short'),
    punching: $('punching').checked,
    existing: table('existing', 2, 'Существующие сваи'),
  };
}

// ── Предпросмотр до расчёта: контур и скважины ──

function preview() {
  if (state.result) return;
  try {
    const outline = table('outline', 2, 'Контур плиты');
    if (outline.length < 3) return;
    $('empty').hidden = false;
    plan.draw({ outline, holes: rings('holes'), boreholes: state.soil === 'geobase' ? state.geobase?.boreholes : [] });
  } catch { /* неполный ввод — без предпросмотра */ }
}

// ── Результат ──

const designVariants = () => state.result.variants.filter((v) => !v.reference && !v.checked);

// Классы нагрузки свай (упругий спрос N/Fd): ≤ 0,8 / ≤ 1,0 / ≤ 1,2 / > 1,2.
const CLASSES = [[0.8, '#3aa655', '≤ 0,8'], [1.0, '#e3c22b', '0,8–1,0'], [1.2, '#ee8a26', '1,0–1,2'], [Infinity, '#d63b2f', '> 1,2']];
const classColor = (k) => CLASSES.find((c) => k <= c[0] + 1e-9)[1];

function variantTable() {
  const tbody = $('variants').tBodies[0];
  tbody.replaceChildren();
  state.result.variants.forEach((v, i) => {
    const tr = tbody.insertRow();
    tr.className = [i === state.current ? 'current' : '', v.feasible ? '' : 'infeasible', v.reference ? 'reference' : ''].join(' ');
    const util = Math.max(v.util_ground, v.util_material);
    const cells = [[v.reference ? `${v.count} сущ.` : v.checked ? `${v.count} пров.` : v.count], [fmt(v.kg_m3, 1)], [fmt(util, 2)], [v.util_punch ? fmt(v.util_punch, 2) : '—'], [fmt(v.s_max * 100, 1)], [fmt(v.ds_rel, 4)],
      [v.over_fd, v.over_fd > 0], [v.over_12fd, v.over_12fd > 0], [fmt(v.demand_max, 2), v.demand_max > 1.2], [v.util_punch_struct ? fmt(v.util_punch_struct, 2) : '—', v.util_punch_struct > 1]];
    for (const [text, warn] of cells) {
      const td = tr.insertCell(); td.textContent = text; td.className = warn ? 'num warn' : 'num';
    }
    tr.title = v.reference ? 'Существующая расстановка (сваи приведены к узлам сетки плиты)' : (v.shift != null ? `Среднее смещение от существующих свай ${fmt(v.shift, 2)} м` : '');
    tr.onclick = () => select(i);
  });
  for (const th of $('variants').tHead.rows[0].cells) th.classList.add('num');
  const vs = designVariants();
  const feasible = vs.filter((v) => v.feasible);
  let note = feasible.length
    ? `Допустимых вариантов: ${feasible.length}. Меньше ${feasible[feasible.length - 1].count} свай — нарушается ${limiting(vs)}. ${state.input.pile_fu ? 'Перегрузка свай по грунту (>Fd — упругий спрос до перераспределения по Fu) — справочно, не отсекает варианты.' : 'Предел Fu не задан: N ≤ Fd — условие допустимости.'}`
    : `Допустимого варианта нет: ${reasons(vs[0]).join(', ')}.`;
  const ref = state.result.variants.find((v) => v.reference);
  if (ref) note += ref.feasible ? ' Существующая расстановка допустима.' : ` Существующая расстановка: ${reasons(ref).join(', ')}.`;
  $('variants-note').textContent = note;
}

function reasons(v) {
  const r = [];
  if (v.util_ground > 1 && !state.input.pile_fu) r.push(`нагрузка на сваю ${fmt(v.util_ground, 2)}·Fd (предел Fu не задан)`);
  if (v.util_material > 1) r.push(`прочность сваи по материалу (Кисп ${fmt(v.util_material, 2)})`);
  if (v.util_punch > 2) r.push(`продавливание плиты сваей (Кисп ${fmt(v.util_punch, 2)} > 2)`);
  if (v.util_punch_struct > 2) r.push(`продавливание плиты конструкциями (Кисп ${fmt(v.util_punch_struct, 2)} > 2)`);
  if (v.s_max > state.input.s_max) r.push(`осадка ${fmt(v.s_max * 100, 1)} см`);
  if (v.ds_rel > state.input.ds_rel) r.push(`Δs/L = ${fmt(v.ds_rel, 4)}`);
  if (v.design_impossible) r.push('армирование плиты не подбирается (мала толщина)');
  return r.length ? r : ['ограничения'];
}

function limiting(vs) {
  const v = vs[vs.length - 1];
  const lim = state.input;
  const terms = [[lim.pile_fu ? 0 : v.util_ground, 'несущая способность свай по грунту'], [v.util_material, 'прочность свай по материалу'], [v.util_punch / 2, 'продавливание плиты сваей'], [v.util_punch_struct / 2, 'продавливание плиты конструкциями'],
    [v.s_max / lim.s_max, 'предел осадки'], [v.ds_rel / lim.ds_rel, 'предел Δs/L']];
  return terms.sort((a, b) => b[0] - a[0])[0][1];
}

function pileTable() {
  const tbody = $('piles').tBodies[0];
  tbody.replaceChildren();
  const v = state.result.variants[state.current];
  if (!v) return;
  const order = v.piles.map((p, i) => i).sort((a, b) => v.piles[b].demand - v.piles[a].demand);
  for (const i of order) {
    const p = v.piles[i];
    const tr = tbody.insertRow();
    const k = p.demand / state.result.fd;
    tr.title = p.capped ? 'На пределе Fu: держит Fu, избыток перераспределён' : '';
    for (const [text, cls] of [[i + 1, 'num'], [fmt(p.x, 2), 'num'], [fmt(p.y, 2), 'num'], [fmt(p.n, 0) + (p.capped ? ' Fu' : ''), 'num'], [fmt(k, 2), k > 1 ? 'num over' : 'num']]) {
      const td = tr.insertCell(); td.textContent = text; td.className = cls;
    }
  }
  for (const th of $('piles').tHead.rows[0].cells) th.classList.add('num');
}

function legend(title, lo, hi, unit, note = '') {
  $('legend').hidden = false;
  $('legend').innerHTML = `<div class="title">${title}</div><div class="scale"><span>${lo}</span><span class="bar" style="background:${GRADIENT}"></span><span>${hi} ${unit}</span></div>${note ? `<div class="note">${note}</div>` : ''}`;
}

function render(refit = false) {
  const R = state.result;
  const v = R?.variants[state.current];
  if (!v) return preview();
  $('empty').hidden = true;
  const scene = { outline: R.outline, holes: R.holes, nodes: R.nodes, elements: R.elements, existing: R.existing, boreholes: state.soil === 'geobase' ? state.geobase?.boreholes : [] };
  const r = R.d / 2;
  if (state.mode === 'piles') {
    const fd = R.fd;
    const lab = $('pile-label').value;
    scene.piles = v.piles.map((p, i) => {
      const k = p.demand / fd;
      return { x: p.x, y: p.y, r, fill: classColor(k), over: p.capped, label: lab === 'k' ? fmt(k, 2) : lab === 'n' ? fmt(p.n, 0) : '',
        title: `Свая ${i + 1}: N = ${fmt(p.n, 0)} кН${p.capped ? ' (на пределе Fu)' : ''}, упругий спрос N/Fd = ${fmt(k, 2)}` };
    });
    $('legend').hidden = false;
    $('legend').innerHTML = `<div class="title">Нагрузка на сваю N/Fd (упругий спрос)</div><div class="classes">${CLASSES.map((c) => `<span><i style="background:${c[1]}"></i>${c[2]}</span>`).join('')}</div>`
      + `<div class="note">Fd = ${fmt(fd, 0)} кН; красная обводка — на пределе Fu${R.existing.length ? '; пунктир — существующие сваи' : ''}</div>`;
  } else if (state.mode === 'settlement') {
    const cellW = R.elements.map((e) => e.reduce((s, n) => s + v.w[n], 0) / e.length);
    const lo = Math.min(...cellW), hi = Math.max(...cellW);
    scene.cell = (i) => color((cellW[i] - lo) / Math.max(hi - lo, 1e-12));
    scene.piles = v.piles.map((p) => ({ x: p.x, y: p.y, r, fill: 'none', title: '' }));
    legend('Осадка', fmt(lo * 100, 2), fmt(hi * 100, 2), 'см', 'сочетание с наибольшей осадкой');
  } else {
    const k = Number($('rebar-dir').value);
    const vals = v.rebar.map((a) => a[k]);
    const hi = Math.max(...vals, 1e-9);
    scene.cell = (i) => (vals[i] > 1e-9 ? color(vals[i] / hi) : null);
    scene.piles = v.piles.map((p) => ({ x: p.x, y: p.y, r, fill: 'none', title: '' }));
    legend(`Требуемая арматура, ${$('rebar-dir').selectedOptions[0].textContent}`, '0', fmt(hi, 1), 'см²/м', 'ПС1, Капра–Мори; без фоновой арматуры');
  }
  plan.draw(scene, refit);
}

function checksTable() {
  const tbody = $('checks').tBodies[0];
  tbody.replaceChildren();
  const v = state.result.variants[state.current];
  if (!v) return;
  const lim = state.input;
  const rows = [
    [`Нагрузка на сваю N/Fd${lim.pile_fu ? ' (справочно)' : ''}`, `${fmt(v.demand_max, 2)}; >Fd: ${v.over_fd}, >1,2Fd: ${v.over_12fd}${v.capped_count ? `, на Fu: ${v.capped_count}` : ''}`, lim.pile_fu ? '' : '1,00', lim.pile_fu ? true : v.demand_max <= 1],
    ['Свая по материалу, Кисп', fmt(v.util_material, 2), '1,00', v.util_material <= 1],
    ['Продавливание сваей, Кисп', v.util_punch ? fmt(v.util_punch, 2) : 'не проверялось', '2,00', v.util_punch <= 2],
    ['Продавливание конструкциями, Кисп', v.util_punch_struct ? fmt(v.util_punch_struct, 2) : 'не проверялось', '2,00', v.util_punch_struct <= 2],
    ['Осадка, см', fmt(v.s_max * 100, 1), fmt(lim.s_max * 100, 1), v.s_max <= lim.s_max],
    ['Разность осадок Δs/L', fmt(v.ds_rel, 4), fmt(lim.ds_rel, 4), v.ds_rel <= lim.ds_rel],
    ['Армирование плиты', v.design_impossible ? 'не подбирается' : `${fmt(v.kg_m3, 1)} кг/м³`, '', !v.design_impossible],
  ];
  for (const [name, value, limit, ok] of rows) {
    const tr = tbody.insertRow();
    tr.className = ok ? '' : 'fail';
    for (const [text, cls] of [[name, ''], [value, 'num'], [limit ? `≤ ${limit}` : '', 'num muted']]) {
      const td = tr.insertCell(); td.textContent = text; td.className = cls;
    }
  }
}

function select(i) {
  state.current = i;
  variantTable();
  checksTable();
  pileTable();
  render();
  $('export-csv').disabled = $('export-dxf').disabled = $('export-rebar').disabled = $('export-report').disabled = i < 0;
  const v = state.result.variants[i];
  $('subtitle').textContent = `${v.reference ? 'существующая: ' : v.checked ? 'проверка: ' : ''}${v.count} свай, ${fmt(v.kg_m3, 1)} кг/м³`;
}

// ── Действия ──

async function busy(text, job) {
  $('busy-text').textContent = text;
  $('busy').hidden = false;
  try { return await job(); } finally { $('busy').hidden = true; }
}

$('run').onclick = async () => {
  let input;
  try { input = collect(); } catch (e) { return status(e.message, true); }
  try {
    const t = performance.now();
    const result = await busy('Расчёт…', () => call('run', input));
    state.input = input;
    state.result = result;
    const feasible = result.variants.filter((v) => v.feasible && !v.reference);
    // По умолчанию — вариант с наименьшим числом свай среди допустимых.
    select(feasible.length ? result.variants.indexOf(feasible[feasible.length - 1]) : result.variants.findIndex((v) => !v.reference));
    render(true);
    status(`Рассчитано вариантов: ${result.variants.length}, ${fmt((performance.now() - t) / 1000, 1)} с`);
  } catch (e) {
    status(`Ошибка: ${e.message ?? e}`, true);
  }
};

// Проверка заданной расстановки: координаты (по умолчанию — выбранный вариант) → расчёт без
// подбора; на тех же исходных данных результат добавляется к вариантам.
$('check').onclick = () => {
  const v = state.result?.variants[state.current];
  $('check-text').value = v ? v.piles.map((p) => `${fmt(p.x, 3)}; ${fmt(p.y, 3)}`).join('\n') : '';
  $('check-dialog').hidden = false;
};
$('check-cancel').onclick = () => { $('check-dialog').hidden = true; };
$('check-run').onclick = async () => {
  let input;
  try { input = collect(); } catch (e) { return status(e.message, true); }
  const piles = [];
  for (const line of $('check-text').value.split('\n')) {
    if (!line.trim()) continue;
    const xy = line.trim().split(/[;\s]+/).map((t) => Number(t.replace(',', '.')));
    if (xy.length < 2 || xy.some((x) => !Number.isFinite(x))) return status(`Не разобрана строка: «${line}»`, true);
    piles.push([xy[0], xy[1]]);
  }
  if (!piles.length) return status('Не задано ни одной сваи', true);
  $('check-dialog').hidden = true;
  try {
    const t = performance.now();
    const result = await busy('Проверка…', () => call('check_piles', { input, piles }));
    state.input = input;
    state.result = result;
    select(result.variants.length - 1);
    render(true);
    const notes = result.check_notes || [];
    status(`Проверено: ${result.variants[result.variants.length - 1].count} свай из ${piles.length}, ${fmt((performance.now() - t) / 1000, 1)} с${notes.length ? '. ' + notes.join('; ') : ''}`, notes.length > 0);
  } catch (e) {
    status(`Ошибка: ${e.message ?? e}`, true);
  }
};

$('save-project').onclick = async () => {
  let input;
  try { input = collect(); } catch (e) { return status(e.message, true); }
  const path = await pickFile('project', true);
  if (!path) return;
  try { await call('save_project', { path, input }); status(`Проект сохранён: ${path}`); } catch (e) { status(`Ошибка: ${e.message ?? e}`, true); }
};

$('open-project').onclick = async () => {
  const path = await pickFile('project');
  if (!path) return;
  try {
    const input = await call('load_project', { path });
    fill(input);
    state.result = null; state.current = -1;
    for (const id of ['export-csv', 'export-dxf', 'export-rebar', 'export-report']) $(id).disabled = true;
    $('variants').tBodies[0].replaceChildren(); $('piles').tBodies[0].replaceChildren(); $('checks').tBodies[0].replaceChildren();
    $('variants-note').textContent = ''; $('legend').hidden = true;
    $('subtitle').textContent = path.split(/[\\/]/).pop();
    preview();
    status(`Проект открыт: ${path}`);
  } catch (e) { status(`Ошибка: ${e.message ?? e}`, true); }
};

onProgress((stage) => { $('busy-text').textContent = `Расчёт: ${stage}`; status(stage); });

for (const b of $('mode').querySelectorAll('button')) {
  b.onclick = () => {
    state.mode = b.dataset.mode;
    for (const x of $('mode').querySelectorAll('button')) x.classList.toggle('current', x === b);
    $('rebar-dir').hidden = state.mode !== 'rebar';
    render();
  };
}
$('rebar-dir').onchange = () => render();
$('pile-label').onchange = () => render();

for (const b of $('layout').querySelectorAll('button')) {
  b.onclick = () => {
    state.layout = b.dataset.layout;
    for (const x of $('layout').querySelectorAll('button')) x.classList.toggle('current', x === b);
  };
}

for (const b of $('soil-source').querySelectorAll('button')) {
  b.onclick = () => {
    state.soil = b.dataset.source;
    for (const x of $('soil-source').querySelectorAll('button')) x.classList.toggle('current', x === b);
    $('soil-geobase').hidden = state.soil !== 'geobase';
    $('soil-layers').hidden = state.soil !== 'layers';
    render();
  };
}

$('open-geobase').onclick = async () => {
  const path = await pickFile('geobase');
  if (!path) return;
  try {
    state.geobase = await call('open_geobase', { path });
    const g = state.geobase;
    $('geobase-info').textContent = `${path.split(/[\\/]/).pop()}: скважин ${g.boreholes.length}, ИГЭ ${g.soils}, модель «${g.model}»`;
    status('Объект GeoBase открыт');
    render();
  } catch (e) { status(`GeoBase: ${e.message ?? e}`, true); }
};

for (const [id, kind, cmd] of [['export-csv', 'csv', 'export_csv'], ['export-dxf', 'dxf', 'export_dxf'], ['export-rebar', 'dxf', 'export_rebar_dxf'], ['export-report', 'html', 'export_report']]) {
  $(id).onclick = async () => {
    const path = await pickFile(kind, true);
    if (!path) return;
    try { await call(cmd, { path, variant: state.current }); status(`Сохранено: ${path}`); } catch (e) { status(`Ошибка: ${e.message ?? e}`, true); }
  };
}

$('outline').addEventListener('input', preview);
$('holes').addEventListener('input', preview);
window.addEventListener('resize', () => render(true));

fill(await call('defaults'));
preview();
