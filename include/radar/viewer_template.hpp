#pragma once

#include <string_view>

namespace radar {

/// Self-contained range-Doppler viewer. `__RADAR_DATA__` is replaced with the JSON payload.
inline constexpr std::string_view kViewerTemplate = R"HTML(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Range-Doppler Viewer</title>
<style>
:root{--bg:#0f1419;--panel:#171e26;--text:#e6edf3;--muted:#8b98a5;--accent:#4fc3f7;--line:#2a3440}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
header{padding:16px 20px;border-bottom:1px solid var(--line)}
h1{margin:0;font-size:18px;font-weight:600}
.sub{color:var(--muted);font-size:13px;margin-top:4px}
main{display:grid;grid-template-columns:minmax(0,1fr) 340px;gap:16px;padding:16px 20px;max-width:1280px}
@media (max-width:900px){main{grid-template-columns:1fr}}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px}
canvas{display:block;width:100%;image-rendering:pixelated;cursor:crosshair}
.controls{display:flex;flex-wrap:wrap;gap:12px 18px;align-items:center;margin-top:10px;color:var(--muted)}
.controls input[type=range]{width:160px;vertical-align:middle}
button{background:var(--accent);color:#06222e;border:0;border-radius:4px;padding:4px 12px;font-weight:600;cursor:pointer}
#readout{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12px;color:var(--accent);min-height:1.5em;margin-top:8px}
table{width:100%;border-collapse:collapse;font-size:12px;font-variant-numeric:tabular-nums}
th,td{text-align:right;padding:3px 6px;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:500}
h2{font-size:13px;margin:4px 0 8px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
dl{display:grid;grid-template-columns:auto 1fr;gap:2px 12px;margin:0;font-size:12px}
dt{color:var(--muted)}dd{margin:0;text-align:right}
.side>.panel+.panel{margin-top:16px}
</style>
</head>
<body>
<header><h1 id="title">Range-Doppler map</h1><div class="sub" id="subtitle"></div></header>
<main>
<section class="panel">
<canvas id="map" width="900" height="640"></canvas>
<div id="readout">Hover the map for range, velocity and power.</div>
<div class="controls">
<label>CPI <input id="cpi" type="range" min="0" value="0"> <span id="cpiLabel"></span></label>
<button id="play">Play</button>
<label>Dynamic range <input id="dr" type="range" min="20" max="70" value="60"> <span id="drLabel"></span></label>
<label><input id="showDet" type="checkbox" checked> detections</label>
<label><input id="showTrk" type="checkbox" checked> tracks</label>
</div>
</section>
<aside class="side">
<div class="panel"><h2>Radar</h2><dl id="params"></dl></div>
<div class="panel"><h2>Detections (this CPI)</h2><table><thead><tr><th>range m</th><th>vel m/s</th><th>SNR dB</th></tr></thead><tbody id="dets"></tbody></table></div>
<div class="panel"><h2>Confirmed tracks</h2><table><thead><tr><th>id</th><th>range m</th><th>vel m/s</th><th>hits</th></tr></thead><tbody id="trks"></tbody></table></div>
</aside>
</main>
<script>
const D = __RADAR_DATA__;
const cv = document.getElementById('map'), ctx = cv.getContext('2d');
const M = {l: 64, r: 70, t: 14, b: 44};
const R = D.rangeBins, N = D.dopplerBins;
const frames = D.frames.map(b64 => Uint8Array.from(atob(b64), c => c.charCodeAt(0)));
const turbo = x => {
  const p = (c) => c.reduceRight((a, k) => a * x + k, 0);
  const f = v => Math.max(0, Math.min(255, Math.round(v * 255)));
  return [f(p([0.13572138,4.61539260,-42.66032258,132.13108234,-152.94239396,59.28637943])),
          f(p([0.09140261,2.19418839,4.84296658,-14.18503333,4.27729857,2.82956604])),
          f(p([0.10667330,12.64194608,-60.58204836,110.36276771,-89.90310912,27.34824973]))];
};
const lut = Array.from({length: 256}, (_, i) => turbo(i / 255));
const off = document.createElement('canvas'); off.width = N; off.height = R;
const offCtx = off.getContext('2d'); const img = offCtx.createImageData(N, R);
const vel = d => (d - Math.floor(N / 2)) * D.velocityBinMps;
const el = id => document.getElementById(id);
let cpi = 0, timer = null;
el('cpi').max = frames.length - 1;
el('title').textContent = 'Range-Doppler map: ' + D.name;
const P = D.params;
el('subtitle').textContent = `${P.carrierGHz} GHz, B = ${P.bandwidthMHz} MHz, PRF = ${P.prfKHz} kHz, ${P.pulses} pulses, ${P.cfar} (Pfa ${P.pfa}), MTI ${P.mti}`;
el('params').innerHTML = [['Carrier', P.carrierGHz + ' GHz'], ['Bandwidth', P.bandwidthMHz + ' MHz'], ['Pulse width', P.pulseWidthUs + ' us'],
  ['PRF', P.prfKHz + ' kHz'], ['Pulses / CPI', P.pulses], ['Range resolution', P.rangeResM.toFixed(1) + ' m'],
  ['Unambiguous velocity', '±' + P.vmaxMps.toFixed(1) + ' m/s'], ['Range window', P.rangeWindow], ['Doppler window', P.dopplerWindow],
  ['Detector', P.cfar + ', Pfa ' + P.pfa], ['MTI', P.mti]].map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
function draw() {
  const W = cv.width - M.l - M.r, H = cv.height - M.t - M.b;
  const q = frames[cpi], dr = +el('dr').value, span = D.dbMax - D.dbMin, lo = 255 * (1 - dr / span);
  for (let r = 0; r < R; r++) for (let d = 0; d < N; d++) {
    const v = q[r * N + d], x = Math.max(0, Math.min(255, Math.round((v - lo) / (255 - lo) * 255)));
    const c = lut[x], o = ((R - 1 - r) * N + d) * 4;
    img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
  }
  offCtx.putImageData(img, 0, 0);
  ctx.fillStyle = '#171e26'; ctx.fillRect(0, 0, cv.width, cv.height);
  ctx.imageSmoothingEnabled = false; ctx.drawImage(off, M.l, M.t, W, H);
  ctx.strokeStyle = '#2a3440'; ctx.strokeRect(M.l, M.t, W, H);
  ctx.fillStyle = '#8b98a5'; ctx.font = '12px system-ui'; ctx.textAlign = 'center';
  const vmax = vel(N - 1);
  const vstep = vmax > 60 ? 20 : vmax > 30 ? 10 : 5;
  for (let v = -Math.floor(-vel(0) / vstep) * vstep; v <= vmax; v += vstep) {
    const x = M.l + (v / D.velocityBinMps + Math.floor(N / 2) + 0.5) / N * W;
    ctx.fillText(v, x, M.t + H + 16);
  }
  ctx.fillText('radial velocity (m/s, + closing)', M.l + W / 2, M.t + H + 36);
  ctx.textAlign = 'right';
  const rmax = (R - 1) * D.rangeBinM, rstep = rmax > 8000 ? 2000 : rmax > 4000 ? 1000 : 500;
  for (let rg = 0; rg <= rmax; rg += rstep) ctx.fillText((rg / 1000).toFixed(1), M.l - 6, M.t + H - rg / D.rangeBinM / R * H + 4);
  ctx.save(); ctx.translate(16, M.t + H / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = 'center'; ctx.fillText('range (km)', 0, 0); ctx.restore();
  for (let i = 0; i < 256; i++) { const c = lut[i]; ctx.fillStyle = `rgb(${c})`; ctx.fillRect(cv.width - 46, M.t + H - (i + 1) / 256 * H, 14, H / 256 + 1); }
  ctx.fillStyle = '#8b98a5'; ctx.textAlign = 'left';
  ctx.fillText(D.dbMax.toFixed(0) + ' dB', cv.width - 56, M.t + 10);
  ctx.fillText((D.dbMax - dr).toFixed(0) + ' dB', cv.width - 56, M.t + H);
  const toXY = (rb, db) => [M.l + (db + 0.5) / N * W, M.t + H - (rb + 0.5) / R * H];
  const dets = D.detections.filter(d => d[0] === cpi);
  if (el('showDet').checked) {
    ctx.strokeStyle = '#ffffff'; ctx.lineWidth = 1.5;
    for (const d of dets) { const [x, y] = toXY(d[1], d[2]); ctx.strokeRect(x - 7, y - 7, 14, 14); }
  }
  if (el('showTrk').checked) {
    ctx.fillStyle = '#ff4081'; ctx.font = 'bold 12px system-ui'; ctx.textAlign = 'left';
    for (const t of D.tracks) {
      if (cpi < t[5] || cpi > t[6]) continue;
      const near = dets.reduce((b, d) => (!b || Math.abs(d[3] - t[1]) < Math.abs(b[3] - t[1])) && Math.abs(d[4] - t[2]) < 5 ? d : b, null);
      if (!near) continue;
      const [x, y] = toXY(near[1], near[2]); ctx.fillText('T' + t[0], x + 10, y - 8);
    }
  }
  el('cpiLabel').textContent = `${cpi + 1} / ${frames.length}`;
  el('drLabel').textContent = dr + ' dB';
  el('dets').innerHTML = dets.map(d => `<tr><td>${d[3].toFixed(1)}</td><td>${d[4].toFixed(2)}</td><td>${d[5].toFixed(1)}</td></tr>`).join('') || '<tr><td colspan="3">none</td></tr>';
  el('trks').innerHTML = D.tracks.map(t => `<tr><td>T${t[0]}</td><td>${t[1].toFixed(1)}</td><td>${t[2].toFixed(2)}</td><td>${t[4]}</td></tr>`).join('') || '<tr><td colspan="4">none</td></tr>';
}
cv.addEventListener('mousemove', e => {
  const b = cv.getBoundingClientRect(), sx = cv.width / b.width, sy = cv.height / b.height;
  const W = cv.width - M.l - M.r, H = cv.height - M.t - M.b;
  const x = (e.clientX - b.left) * sx - M.l, y = (e.clientY - b.top) * sy - M.t;
  if (x < 0 || y < 0 || x >= W || y >= H) return;
  const d = Math.floor(x / W * N), r = R - 1 - Math.floor(y / H * R);
  const db = D.dbMin + frames[cpi][r * N + d] / 255 * (D.dbMax - D.dbMin);
  el('readout').textContent = `range ${(r * D.rangeBinM).toFixed(1)} m | velocity ${vel(d).toFixed(2)} m/s | ${db.toFixed(1)} dB (bin ${r}, ${d})`;
});
el('cpi').oninput = e => { cpi = +e.target.value; draw(); };
el('dr').oninput = draw; el('showDet').onchange = draw; el('showTrk').onchange = draw;
el('play').onclick = () => {
  if (timer) { clearInterval(timer); timer = null; el('play').textContent = 'Play'; return; }
  el('play').textContent = 'Pause';
  timer = setInterval(() => { cpi = (cpi + 1) % frames.length; el('cpi').value = cpi; draw(); }, 500);
};
draw();
</script>
</body>
</html>
)HTML";

}  // namespace radar
