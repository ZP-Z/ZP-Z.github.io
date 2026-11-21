// ---- Paths of light around the hole (Schwarzschild, rs = 1) ----
// A ray is described in its own plane by u = 1/r as a function of the angle
// psi swept around the hole, counted from the direction of the camera:
// u'' = 1.5 u^2 - u.
export const B_CRIT = 1.5 * Math.sqrt(3);      // rays with a smaller impact parameter are captured
export const LUT_W = 768, LUT_H = 96;          // angle x radius
export const PSI_MAX = 3 * Math.PI;
export const R_MIN = 1.55, R_MAX = 24;
export const SKY_W = 512, SKY_L0 = Math.log(1e-5), SKY_L1 = Math.log(38);
// The ray table: impact parameter x angle. Its columns crowd around the
// critical impact parameter, where a ray is most sensitive.
export const RAY_W = 512, RAY_H = 384, RAY_PSI = 2.3 * Math.PI;
export const RAY_B = 24, RAY_EPS = 0.002, RAY_X0 = 0.3;
export const RAY_LA = Math.log(1 + (RAY_B - B_CRIT) / RAY_EPS), RAY_LB = Math.log(1 + B_CRIT / RAY_EPS);
export const GAS_W = 512, GAS_H = 512;

// Follows the ray with impact parameter `b` from the camera. Calls
// onCross(j, psi, outgoing) each time it crosses radius rj[j], and returns
// the angle of closest approach (or -1 if the hole captures it).
function shoot(distance, b, uj, onCross, psiStop) {
  let u = 1 / distance, psi = 0;
  let v = Math.sqrt(Math.max(1 / (b * b) - u * u * (1 - u), 0));
  const h = Math.min(0.003, 0.004 * b);
  const H = uj.length;
  let jIn = H - 1, jOut = H, outgoing = false, apex = -1;
  for (let n = 0; n < 400000; n++) {
    const k1v = 1.5 * u * u - u;
    const u2 = u + 0.5 * h * v, v2 = v + 0.5 * h * k1v;
    const k2v = 1.5 * u2 * u2 - u2;
    const u3 = u + 0.5 * h * v2, v3 = v + 0.5 * h * k2v;
    const k3v = 1.5 * u3 * u3 - u3;
    const u4 = u + h * v3, v4 = v + h * k3v;
    const k4v = 1.5 * u4 * u4 - u4;
    const un = u + h / 6 * (v + 2 * v2 + 2 * v3 + v4);
    const vn = v + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
    if (!outgoing) {
      while (jIn >= 0 && un >= uj[jIn]) { onCross(jIn, psi + h * (uj[jIn] - u) / (un - u), false); jIn--; }
      if (vn < 0) {
        outgoing = true;
        apex = psi + h * v / (v - vn);
        jOut = jIn + 1;
      }
    }
    if (outgoing) {
      while (jOut < H && un <= uj[jOut]) { onCross(jOut, psi + h * (uj[jOut] - u) / (un - u), true); jOut++; }
      if (jOut >= H) break;
    }
    u = un; v = vn; psi += h;
    if (u > 0.985 || psi > psiStop) break;
  }
  return apex;
}

// Interpolates y(x) through points with increasing x (cubic Hermite).
function hermite(xs, ys, n, x, hint) {
  let k = hint[0];
  while (k < n - 2 && xs[k + 1] < x) k++;
  while (k > 0 && xs[k] > x) k--;
  hint[0] = k;
  const x0 = xs[k], x1 = xs[k + 1], y0 = ys[k], y1 = ys[k + 1];
  const dx = x1 - x0;
  if (dx <= 0) return y0;
  const m0 = k > 0 ? (y1 - ys[k - 1]) / (x1 - xs[k - 1]) : (y1 - y0) / dx;
  const m1 = k < n - 2 ? (ys[k + 2] - y0) / (xs[k + 2] - x0) : (y1 - y0) / dx;
  const t = Math.min(Math.max((x - x0) / dx, 0), 1), t2 = t * t, t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * dx * m0 + (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * dx * m1;
}

// The lens table. For a point at radius r whose direction from the hole
// makes the angle psi with the direction to the camera, it gives the impact
// parameter b of the ray that reaches it, and how b changes with psi and r.
// psi up to pi is the direct image; from pi to 2 pi the light went round the
// far side of the hole.
function buildLens(distance) {
  const rj = new Float64Array(LUT_H), uj = new Float64Array(LUT_H);
  for (let j = 0; j < LUT_H; j++) { rj[j] = R_MIN * Math.pow(R_MAX / R_MIN, j / (LUT_H - 1)); uj[j] = 1 / rj[j]; }
  const inc = [], out = [];
  for (let j = 0; j < LUT_H; j++) { inc.push([]); out.push([]); }
  const stop = PSI_MAX + 0.4;
  const run = b => shoot(distance, b, uj, (j, psi, outgoing) => (outgoing ? out : inc)[j].push(b, psi), stop);
  const bTop = R_MAX / Math.sqrt(1 - 1 / R_MAX) * 1.002;
  for (let k = 1; k <= 900; k++) run(bTop * k / 900);
  for (let k = 0; k <= 112; k++) { const d = Math.pow(10, -k / 16); run(B_CRIT + d); run(B_CRIT - d); }
  // The ray that just touches each radius marks where the two branches meet.
  const apexB = new Float64Array(LUT_H), apexPsi = new Float64Array(LUT_H);
  for (let j = 0; j < LUT_H; j++) {
    apexB[j] = rj[j] / Math.sqrt(1 - 1 / rj[j]);
    apexPsi[j] = shoot(distance, apexB[j] * (1 - 1e-9), uj, () => {}, stop);
  }

  const B = new Float64Array(LUT_W * LUT_H);
  const xs = new Float64Array(2400), ys = new Float64Array(2400), hint = [0];
  for (let j = 0; j < LUT_H; j++) {
    const a = [], pin = inc[j], pout = out[j];
    for (let k = 0; k < pin.length; k += 2) if (pin[k] < apexB[j]) a.push([pin[k + 1], pin[k]]);
    a.push([apexPsi[j], apexB[j]]);
    for (let k = 0; k < pout.length; k += 2) if (pout[k] < apexB[j]) a.push([pout[k + 1], pout[k]]);
    a.sort((p, q) => p[0] - q[0]);
    let n = 0;
    xs[n] = 0; ys[n] = 0; n++;
    for (const p of a) if (p[0] > xs[n - 1] + 1e-7) { xs[n] = p[0]; ys[n] = p[1]; n++; }
    hint[0] = 0;
    for (let i = 0; i < LUT_W; i++) {
      const psi = PSI_MAX * i / (LUT_W - 1);
      B[j * LUT_W + i] = psi >= xs[n - 1] ? ys[n - 1] : hermite(xs, ys, n, psi, hint);
    }
  }
  const data = new Float32Array(LUT_W * LUT_H * 4);
  const dpsi = PSI_MAX / (LUT_W - 1);
  for (let j = 0; j < LUT_H; j++) {
    const j0 = Math.max(j - 1, 0), j1 = Math.min(j + 1, LUT_H - 1);
    for (let i = 0; i < LUT_W; i++) {
      const i0 = Math.max(i - 1, 0), i1 = Math.min(i + 1, LUT_W - 1);
      const o = (j * LUT_W + i) * 4;
      data[o] = B[j * LUT_W + i];
      data[o + 1] = (B[j * LUT_W + i1] - B[j * LUT_W + i0]) / ((i1 - i0) * dpsi);
      data[o + 2] = (B[j1 * LUT_W + i] - B[j0 * LUT_W + i]) / (rj[j1] - rj[j0]);
      data[o + 3] = 0;
    }
  }
  return data;
}

// The sky table: the direction a ray ends up travelling in, as an angle
// around the hole from the direction of the camera, against its impact
// parameter. Without a hole the angle would be pi minus the angle of view.
function buildSky(distance) {
  const data = new Float32Array(SKY_W * 4);
  const far = [1 / 400];
  for (let i = 0; i < SKY_W; i++) {
    const b = B_CRIT + Math.exp(SKY_L0 + (SKY_L1 - SKY_L0) * i / (SKY_W - 1));
    let psiFar = -1;
    shoot(distance, b, far, (j, psi, outgoing) => { if (outgoing) psiFar = psi; }, 40);
    data[i * 4] = psiFar < 0 ? 40 : psiFar + Math.asin(Math.min(b / 400 * Math.sqrt(1 - 1 / 400), 1));
  }
  return data;
}

// The ray table. For the line of sight with impact parameter b it gives
// u = 1/r at the angle psi it has swept around the hole, and how fast u
// changes with psi. The gas is drawn with it: it tells where a line of
// sight meets the disk, and at what slant.
function buildRays(distance) {
  const data = new Float32Array(RAY_W * RAY_H * 2);
  const sub = 8, rowStep = RAY_PSI / ((RAY_H - 1) * sub);
  for (let i = 0; i < RAY_W; i++) {
    const x = i / (RAY_W - 1);
    const b = Math.max(x >= RAY_X0
      ? B_CRIT + RAY_EPS * (Math.exp((x - RAY_X0) / (1 - RAY_X0) * RAY_LA) - 1)
      : B_CRIT - RAY_EPS * (Math.exp((RAY_X0 - x) / RAY_X0 * RAY_LB) - 1), 1e-3);
    // A nearly radial ray reaches the horizon after an angle proportional
    // to b. Give each column a comparable number of samples before capture.
    const scale = Math.sqrt(1 + b * b);
    const h = rowStep * b / scale;
    let u = 1 / distance, v = Math.sqrt(Math.max(1 / (b * b) - u * u * (1 - u), 0)), done = false;
    for (let j = 0; j < RAY_H; j++) {
      // Both channels vary smoothly between columns, including at the axis.
      const o = (j * RAY_W + i) * 2;
      data[o] = u * scale; data[o + 1] = v * b;
      for (let n = 0; n < sub && !done; n++) {
        const k1 = 1.5 * u * u - u;
        const u2 = u + 0.5 * h * v, v2 = v + 0.5 * h * k1;
        const k2 = 1.5 * u2 * u2 - u2;
        const u3 = u + 0.5 * h * v2, v3 = v + 0.5 * h * k2;
        const k3 = 1.5 * u3 * u3 - u3;
        const u4 = u + h * v3, v4 = v + h * k3;
        const k4 = 1.5 * u4 * u4 - u4;
        u += h / 6 * (v + 2 * v2 + 2 * v3 + v4);
        v += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        // Past the horizon, or gone back out to the far sky: the ray is over.
        if (u > 3) { u = 3; v = 0; done = true; }
        if (u < -0.05) { u = -0.05; v = 0; done = true; }
      }
    }
  }
  return data;
}

function mulberry(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// The gas texture. Red and green: long fibres, for the two sheets of the
// disk. Blue: short tufts. Alpha: large clouds. A fibre is long along the
// flow (across the texture) and thin across it, and it wanders a little.
function buildGas() {
  const rand = mulberry(40977);
  // Smooth noise on a grid of cx by cy cells that repeats at the edges.
  const grid = (cx, cy) => {
    const g = new Float32Array(cx * cy);
    for (let i = 0; i < g.length; i++) g[i] = rand();
    return (x, y) => {
      const fx = x * cx, fy = y * cy;
      const ix = Math.floor(fx), iy = Math.floor(fy);
      let tx = fx - ix, ty = fy - iy;
      tx = tx * tx * (3 - 2 * tx); ty = ty * ty * (3 - 2 * ty);
      const x0 = ((ix % cx) + cx) % cx, x1 = (x0 + 1) % cx, y0 = ((iy % cy) + cy) % cy, y1 = (y0 + 1) % cy;
      const a = g[y0 * cx + x0], b = g[y0 * cx + x1], c = g[y1 * cx + x0], d = g[y1 * cx + x1];
      return a + (b - a) * tx + (c - a) * ty + (a - b - c + d) * tx * ty;
    };
  };
  // Several grids, each twice as fine as the one before. `ridge` turns part
  // of each into thin bright lines.
  const stack = (cx, cy, octaves, gain, ridge) => {
    const layers = [];
    for (let k = 0; k < octaves; k++) layers.push(grid(cx << k, cy << k));
    return (x, y) => {
      let v = 0, a = 1, norm = 0;
      for (let k = 0; k < octaves; k++) {
        const n = layers[k](x, y);
        v += a * (n + ridge * (1 - Math.abs(2 * n - 1) - n));
        norm += a; a *= gain;
      }
      return v / norm;
    };
  };
  const fibreA = stack(8, 16, 5, 0.8, 0.4), fibreB = stack(8, 16, 5, 0.8, 0.4);
  const tufts = stack(32, 8, 3, 0.75, 0.3), clouds = stack(8, 4, 4, 0.6, 0);
  const bendA = grid(7, 3), bendB = grid(9, 4), ripple = grid(18, 9), kinkA = grid(48, 20), kinkB = grid(56, 24);
  const N = GAS_W * GAS_H;
  const ch = [new Float32Array(N), new Float32Array(N), new Float32Array(N), new Float32Array(N)];
  for (let j = 0; j < GAS_H; j++) {
    const y = (j + 0.5) / GAS_H;
    for (let i = 0; i < GAS_W; i++) {
      const x = (i + 0.5) / GAS_W, o = j * GAS_W + i;
      const rp = 0.012 * (ripple(x, y) - 0.5);
      const ya = y + 0.05 * (bendA(x, y) - 0.5) + rp + 0.006 * (kinkA(x, y) - 0.5);
      const yb = y + 0.05 * (bendB(x, y) - 0.5) - rp + 0.006 * (kinkB(x, y) - 0.5);
      ch[0][o] = fibreA(x, ya);
      ch[1][o] = fibreB(x, yb);
      ch[2][o] = tufts(x, 0.5 * (ya + yb));
      ch[3][o] = clouds(x, y + 0.5 * (ya - y));
    }
  }
  // Every channel is spread to the same mean and contrast.
  const data = new Uint8Array(N * 4);
  for (let c = 0; c < 4; c++) {
    let m = 0, q = 0;
    for (let o = 0; o < N; o++) m += ch[c][o];
    m /= N;
    for (let o = 0; o < N; o++) q += (ch[c][o] - m) * (ch[c][o] - m);
    const k = 0.17 / Math.sqrt(q / N);
    for (let o = 0; o < N; o++) data[o * 4 + c] = Math.max(0, Math.min(255, Math.round(255 * (0.5 + (ch[c][o] - m) * k))));
  }
  return data;
}

export function buildTables(distance) {
  return { lens: buildLens(distance), sky: buildSky(distance), rays: buildRays(distance), gas: buildGas() };
}
