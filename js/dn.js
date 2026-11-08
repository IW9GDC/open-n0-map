// Open N0 — read the map of the refractivity lapse rate DN (data/dn.bin, tools/build_dn.py) in a browser or Node: the
// median of DN = N(surface) - N(surface + 1 km), N-units/km, computed every 6 hours in 1991-2020 from NCEP-DOE
// Reanalysis 2 (public domain): the quantity ITU-R P.452 and P.1812 call ΔN (map DN50). One layer of 121 x 241 points,
// 1.5 degrees, rows from +90 down, columns from 0 east; uint16 = DN x 10, little-endian. MIT license.
const ROWS = 121, COLS = 241, STEP = 1.5, LAYER = ROWS * COLS;
let grid = null;

export function loadDN(url) {
  if (!grid) grid = fetch(url).then(r => {
    if (!r.ok) throw new Error(`DN map: HTTP ${r.status}`);
    return r.arrayBuffer();
  }).then(b => {
    const a = new Uint16Array(b);
    if (a.length !== LAYER) throw new Error('DN map: unexpected size');
    return a;
  }).catch(e => { grid = null; throw e; });
  return grid;
}

// DN at a point, bilinear between the four nearest grid points (as ITU-R P.1144 describes), N-units/km
export async function dnAt(url, lat, lon) {
  const g = await loadDN(url);
  const y = Math.min(Math.max((90 - lat) / STEP, 0), ROWS - 1.000001);
  const x = (((lon % 360) + 360) % 360) / STEP;
  const i = Math.floor(y), j = Math.min(Math.floor(x), COLS - 2), fy = y - i, fx = x - j;
  const v = (r, c) => g[r * COLS + c] / 10;
  return (v(i, j) * (1 - fx) + v(i, j + 1) * fx) * (1 - fy) + (v(i + 1, j) * (1 - fx) + v(i + 1, j + 1) * fx) * fy;
}
