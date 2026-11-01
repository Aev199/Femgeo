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
  $('q').value = input.loads.q;
  $('points').value = writeRows(input.loads.points.map((x) => [x.x, x.y, x.n]));
  $('pile-d').value = p.d; $('pile-l').value = p.length; $('pile-fd').value = p.fd; $('pile-kind').value = p.kind;
  $('pile-concrete').value = p.concrete; $('pile-rebar').value = p.rebar; $('pile-bars').value = p.bars;
  $('pile-bar').value = Math.round(p.bar_d * 1000); $('pile-cover').value = p.cover; $('pile-l0').value = p.l0;
  $('pile-spacing').value = p.min_spacing ?? ''; $('pile-edge').value = p.min_edge ?? '';
  if (input.soil.source === 'layers') {
    $('soil-top').value = input.soil.top;
    $('layers').value = writeRows(input.soil.layers.map((l) => [l.bottom, l.e, l.nu]));
  }
  $('s-max').value = Math.round(input.s_max * 100); $('ds-rel').value = input.ds_rel;
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
    loads: { q: num('q'), points: table('points', 3, 'Сосредоточенные нагрузки').map(([x, y, n]) => ({ x, y, n })) },
    pile: {
      d: num('pile-d'), length: num('pile-l'), fd: num('pile-fd'), kind: $('pile-kind').value,
      concrete: $('pile-concrete').value, rebar: $('pile-rebar').value, bars: num('pile-bars'), bar_d: num('pile-bar') / 1000,
      cover: num('pile-cover'), l0: num('pile-l0'), min_spacing: optional('pile-spacing'), min_edge: optional('pile-edge'),
    },
    soil,
    s_max: num('s-max') / 100,
    ds_rel: num('ds-rel'),
    layout: state.layout,
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

function variantTable() {
  const tbody = $('variants').tBodies[0];
  tbody.replaceChildren();
  const vs = state.result.variants;
  vs.forEach((v, i) => {
    const tr = tbody.insertRow();
    tr.className = (i === state.current ? 'current ' : '') + (v.feasible ? '' : 'infeasible');
    const util = Math.max(v.util_ground, v.util_material);
    for (const [text, cls] of [[v.count, 'num'], [fmt(v.kg_m3, 1), 'num'], [fmt(util, 2), 'num'], [fmt(v.s_max * 100, 1), 'num'], [fmt(v.ds_rel, 4), 'num']]) {
      const td = tr.insertCell(); td.textContent = text; td.className = cls;
    }
    tr.onclick = () => select(i);
  });
  for (const th of $('variants').tHead.rows[0].cells) th.classList.add('num');
  const feasible = vs.filter((v) => v.feasible);
  $('variants-note').textContent = feasible.length
    ? `Допустимых вариантов: ${feasible.length}. Меньше ${feasible[feasible.length - 1].count} свай — нарушается ${limiting(vs)}.`
    : `Допустимого варианта нет: ${reasons(vs[0]).join(', ')}.`;
}

function reasons(v) {
  const r = [];
  if (v.util_ground > 1) r.push(`нагрузка на сваю ${fmt(v.util_ground, 2)}·Fd`);
  if (v.util_material > 1) r.push(`прочность сваи по материалу (Кисп ${fmt(v.util_material, 2)})`);
  if (v.s_max > state.input.s_max) r.push(`осадка ${fmt(v.s_max * 100, 1)} см`);
  if (v.ds_rel > state.input.ds_rel) r.push(`Δs/L = ${fmt(v.ds_rel, 4)}`);
  if (v.design_impossible) r.push('армирование плиты не подбирается (мала толщина)');
  return r.length ? r : ['ограничения'];
}

function limiting(vs) {
  const v = vs[vs.length - 1];
  const lim = state.input;
  const terms = [[v.util_ground, 'несущая способность свай по грунту'], [v.util_material, 'прочность свай по материалу'],
    [v.s_max / lim.s_max, 'предел осадки'], [v.ds_rel / lim.ds_rel, 'предел Δs/L']];
  return terms.sort((a, b) => b[0] - a[0])[0][1];
}

function pileTable() {
  const tbody = $('piles').tBodies[0];
  tbody.replaceChildren();
  const v = state.result.variants[state.current];
  if (!v) return;
  const order = v.piles.map((p, i) => i).sort((a, b) => v.piles[b].n - v.piles[a].n);
  for (const i of order) {
    const p = v.piles[i];
    const tr = tbody.insertRow();
    const k = p.n / state.result.fd;
    for (const [text, cls] of [[i + 1, 'num'], [fmt(p.x, 2), 'num'], [fmt(p.y, 2), 'num'], [fmt(p.n, 0), 'num'], [fmt(k, 2), k > 1 ? 'num over' : 'num']]) {
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
  const scene = { outline: R.outline, holes: R.holes, nodes: R.nodes, elements: R.elements, boreholes: state.soil === 'geobase' ? state.geobase?.boreholes : [] };
  const r = R.d / 2;
  if (state.mode === 'piles') {
    const fd = R.fd;
    scene.piles = v.piles.map((p, i) => ({ x: p.x, y: p.y, r, fill: color(p.n / fd), over: p.n > fd, title: `Свая ${i + 1}: N = ${fmt(p.n, 0)} кН, N/Fd = ${fmt(p.n / fd, 2)}` }));
    legend('Нагрузка на сваю N/Fd', '0', '1', '', `Fd = ${fmt(fd, 0)} кН; обводка красным — больше Fd`);
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

function select(i) {
  state.current = i;
  variantTable();
  pileTable();
  render();
  $('export-csv').disabled = $('export-dxf').disabled = i < 0;
  const v = state.result.variants[i];
  $('subtitle').textContent = `${v.count} свай, ${fmt(v.kg_m3, 1)} кг/м³`;
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
    const feasible = result.variants.filter((v) => v.feasible);
    // По умолчанию — вариант с наименьшим числом свай среди допустимых.
    select(feasible.length ? result.variants.indexOf(feasible[feasible.length - 1]) : 0);
    render(true);
    status(`Рассчитано вариантов: ${result.variants.length}, ${fmt((performance.now() - t) / 1000, 1)} с`);
  } catch (e) {
    status(`Ошибка: ${e.message ?? e}`, true);
  }
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

for (const [id, kind, cmd] of [['export-csv', 'csv', 'export_csv'], ['export-dxf', 'dxf', 'export_dxf']]) {
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
