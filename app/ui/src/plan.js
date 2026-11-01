// План плиты в SVG: контур, поле по элементам (осадки, армирование), сваи, скважины.
// Ось Y вверх; масштаб колесом, сдвиг перетаскиванием.

const NS = 'http://www.w3.org/2000/svg';

// Последовательная шкала: светлый → тёмный синий (одна тональность, без радуги).
const STOPS = [[0xF5, 0xF7, 0xF8], [0xBB, 0xDE, 0xFB], [0x64, 0xB5, 0xF6], [0x19, 0x76, 0xD2], [0x0D, 0x47, 0xA1]];
export function color(t) {
  const x = Math.min(Math.max(t, 0), 1) * (STOPS.length - 1);
  const i = Math.min(Math.floor(x), STOPS.length - 2);
  const f = x - i;
  const c = STOPS[i].map((v, k) => Math.round(v + (STOPS[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}
export const GRADIENT = `linear-gradient(to right, ${STOPS.map((c) => `rgb(${c.join(',')})`).join(', ')})`;

function el(name, attrs = {}, parent) {
  const e = document.createElementNS(NS, name);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (parent) parent.appendChild(e);
  return e;
}

export class Plan {
  constructor(svg) {
    this.svg = svg;
    this.view = null;
    svg.addEventListener('wheel', (e) => {
      if (!this.view) return;
      e.preventDefault();
      const k = e.deltaY > 0 ? 1.15 : 1 / 1.15;
      const [px, py] = this.toModel(e);
      const v = this.view;
      v.x = px - (px - v.x) * k; v.y = py - (py - v.y) * k; v.w *= k; v.h *= k;
      this.apply();
    }, { passive: false });
    let drag = null;
    svg.addEventListener('pointerdown', (e) => {
      if (!this.view) return;
      drag = { x: e.clientX, y: e.clientY, v: { ...this.view } };
      svg.setPointerCapture(e.pointerId);
      svg.classList.add('dragging');
    });
    svg.addEventListener('pointermove', (e) => {
      if (!drag) return;
      const s = this.view.w / svg.clientWidth;
      this.view.x = drag.v.x - (e.clientX - drag.x) * s;
      this.view.y = drag.v.y - (e.clientY - drag.y) * s;
      this.apply();
    });
    const end = () => { drag = null; svg.classList.remove('dragging'); };
    svg.addEventListener('pointerup', end);
    svg.addEventListener('pointercancel', end);
    svg.addEventListener('dblclick', () => this.fit());
  }

  toModel(e) {
    const r = this.svg.getBoundingClientRect();
    const v = this.view;
    return [v.x + (e.clientX - r.left) / r.width * v.w, v.y + (e.clientY - r.top) / r.height * v.h];
  }

  apply() {
    const v = this.view;
    this.svg.setAttribute('viewBox', `${v.x} ${v.y} ${v.w} ${v.h}`);
  }

  fit() {
    if (!this.bounds) return;
    const [x0, y0, x1, y1] = this.bounds;
    const r = this.svg.getBoundingClientRect();
    const m = 0.08 * Math.max(x1 - x0, y1 - y0);
    let w = x1 - x0 + 2 * m, h = y1 - y0 + 2 * m;
    const aspect = r.width / Math.max(r.height, 1);
    if (w / h < aspect) w = h * aspect; else h = w / aspect;
    // Y модели вверх: в SVG рисуем −y.
    this.view = { x: (x0 + x1) / 2 - w / 2, y: -(y0 + y1) / 2 - h / 2, w, h };
    this.apply();
  }

  /// `scene`: { outline, holes, nodes, elements, cell(i) → цвет|null, piles: [{x,y,r,fill,over,title}], boreholes }
  draw(scene, refit) {
    const svg = this.svg;
    svg.replaceChildren();
    const xs = scene.outline.map((p) => p[0]), ys = scene.outline.map((p) => p[1]);
    for (const b of scene.boreholes || []) { xs.push(b.x); ys.push(b.y); }
    const bounds = [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
    if (refit || !this.view || String(bounds) !== String(this.bounds)) { this.bounds = bounds; this.fit(); }
    const g = el('g', {}, svg);
    const pt = (p) => `${p[0]},${-p[1]}`;
    if (scene.elements && scene.cell) {
      const cells = el('g', {}, g);
      scene.elements.forEach((e, i) => {
        const c = scene.cell(i);
        const poly = el('polygon', { points: e.map((n) => pt(scene.nodes[n])).join(' '), class: c ? 'cell' : 'mesh' }, cells);
        if (c) poly.setAttribute('fill', c);
      });
    }
    for (const ring of [scene.outline, ...(scene.holes || [])]) {
      el('polygon', { points: ring.map(pt).join(' '), class: 'raft' }, g);
    }
    const size = Math.max(bounds[2] - bounds[0], bounds[3] - bounds[1]);
    // Габаритные размеры плиты: под планом и слева, как на чертеже.
    const [ox0, oy0, ox1, oy1] = [Math.min(...scene.outline.map((p) => p[0])), Math.min(...scene.outline.map((p) => p[1])), Math.max(...scene.outline.map((p) => p[0])), Math.max(...scene.outline.map((p) => p[1]))];
    const off = size * 0.04, fs = size * 0.022;
    const dim = (x1, y1, x2, y2, text, vertical) => {
      el('line', { x1, y1: -y1, x2, y2: -y2, class: 'dim' }, g);
      for (const [x, y] of [[x1, y1], [x2, y2]]) {
        el('line', vertical ? { x1: x - off * 0.3, y1: -y, x2: x + off * 0.3, y2: -y, class: 'dim' } : { x1: x, y1: -y - off * 0.3, x2: x, y2: -y + off * 0.3, class: 'dim' }, g);
      }
      const t = el('text', vertical
        ? { x: x1 - fs * 0.5, y: -(y1 + y2) / 2, 'font-size': fs, 'text-anchor': 'middle', transform: `rotate(-90 ${x1 - fs * 0.5} ${-(y1 + y2) / 2})` }
        : { x: (x1 + x2) / 2, y: -y1 + fs * 1.3, 'font-size': fs, 'text-anchor': 'middle' }, g);
      t.textContent = text;
    };
    dim(ox0, oy0 - off, ox1, oy0 - off, (ox1 - ox0).toFixed(2).replace('.', ','), false);
    dim(ox0 - off, oy0, ox0 - off, oy1, (oy1 - oy0).toFixed(2).replace('.', ','), true);
    for (const b of scene.boreholes || []) {
      const s = size * 0.012;
      el('polygon', { points: `${b.x},${-b.y - s} ${b.x - s},${-b.y + s * 0.7} ${b.x + s},${-b.y + s * 0.7}`, class: 'borehole' }, g);
      const t = el('text', { x: b.x + s * 1.3, y: -b.y + s * 0.4, 'font-size': s * 1.6 }, g);
      t.textContent = b.name;
    }
    for (const e of scene.existing || []) {
      el('circle', { cx: e[0], cy: -e[1], r: (scene.piles?.[0]?.r ?? 0.3), class: 'existing' }, g);
    }
    for (const p of scene.piles || []) {
      const c = el('circle', { cx: p.x, cy: -p.y, r: p.r, fill: p.fill, class: p.over ? 'pile over' : 'pile' }, g);
      el('title', {}, c).textContent = p.title;
    }
  }
}
