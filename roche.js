/*
 * Roche-lobe overflow background — a black-hole low-mass X-ray binary in a
 * modeled on a V404 Cygni–like system.
 *
 * Inner binary: a 9.0 M☉ black hole (M2) accretes from a donor (M1) that starts
 * at its 2.0 M☉ progenitor mass and is stripped by Roche-lobe overflow down to
 * its present-day mass. The black-hole mass is held fixed (non-conservative
 * transfer: most of the gas is lost in outflows). As the mass ratio
 * q = M1/M2 falls, the Roche lobes and L1 are recomputed live.
 *
 * Gas leaves L1 at the local sound speed and follows the restricted
 * three-body equations of motion in the corotating frame (gravity of both
 * bodies + centrifugal + Coriolis). The stream is ballistic until it meets disk
 * gas, where it dissipates energy at a rate set by the local disk density
 * (stream–disk shock), relaxes onto Keplerian orbits, spreads by viscous
 * diffusion out to the tidal radius, drifts inward and is accreted; accreted
 * gas re-launches from L1. Colour = projected column density Σ (darkest red =
 * highest Σ). The donor's face toward the black hole is shaded as irradiated.
 *
 * Units for the dynamics: G = M1 + M2 = a = Ω = 1.
 */
(function () {
  'use strict';

  var canvas = document.getElementById('roche');
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var reduceMotion = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ---- System parameters ----
  var M_BH = 9.0;                      // M☉, accretor (held fixed)
  var M_DONOR_0 = 2.0;                 // M☉, donor progenitor mass at onset of RLOF
  var M_DONOR_NOW = 0.54;              // M☉, present-day donor (q ≈ 0.06)
  var T_END = 1400;                    // Myr from onset of RLOF to present day
  var TAU = 400;                       // Myr, e-folding time of the mass-transfer rate
  var MYR_PER_S = 10;                  // displayed Myr per real second
  var INC = 67 * Math.PI / 180;        // orbital inclination

  // ---- Gas parameters ----
  var R_ACC = 0.03;                    // inner disk edge: gas inside is accreted
  var CS = 0.035;                      // sound speed at L1 (launch speed)
  var INFLOW = 0.03;                   // mean viscous inflow, as a fraction of v_Kepler
  var DIFFUSE = 0.18;                  // turbulent radial velocity spread (viscous diffusion)
  var KMAX = 10, DENS0 = 3;            // stream–disk dissipation rate and its density scale
  var PERIOD_S = 36;                   // real seconds per orbit (display)
  var RATE = 2 * Math.PI / PERIOD_S;   // simulation time per real second
  var HMAX = 0.0015;                   // max integration step
  var NP = 3000;                       // gas particles

  var cosI = Math.cos(INC), sinI = Math.sin(INC);
  // Position angle on the sky: the whole system slowly turns on the screen
  var PA0 = -0.35, PA_PERIOD_S = 240;
  var PA = PA0, cosP = Math.cos(PA), sinP = Math.sin(PA);

  // ---- Evolving geometry (recomputed as q changes) ----
  var Q, MU, MA, MD, XA, XD, L1, PHI_L1, LOBE_D, LOBE_A, R_DONOR, R_TIDAL, R_EGG_D;
  var tMyr = 0, geomAt = -1;

  function donorMass(tm) {
    var f = (1 - Math.exp(-tm / TAU)) / (1 - Math.exp(-T_END / TAU));
    return M_DONOR_0 - (M_DONOR_0 - M_DONOR_NOW) * f;
  }

  function eggleton(q) { // Roche-lobe radius / a for the star with mass ratio q (its mass / companion)
    var q13 = Math.cbrt(q), q23 = q13 * q13;
    return 0.49 * q23 / (0.6 * q23 + Math.log(1 + q13));
  }

  function phi(x, y) {
    var d1 = x - XA, d2 = x - XD;
    return -MA / Math.sqrt(d1 * d1 + y * y) - MD / Math.sqrt(d2 * d2 + y * y) - 0.5 * (x * x + y * y);
  }
  function dphidx(x) {
    var d1 = x - XA, d2 = x - XD;
    return MA * d1 / Math.abs(d1 * d1 * d1) + MD * d2 / Math.abs(d2 * d2 * d2) - x;
  }

  // Lobe outline in (u along the axis toward the companion, v perpendicular)
  function lobe(xc, dir) {
    var pts = [], n = 120, dist = Math.abs(L1 - xc);
    for (var k = 0; k <= n; k++) {
      var th = k / n * 2 * Math.PI, c = Math.cos(th), s = Math.sin(th);
      var lo = 1e-3, hi = dist;
      if (phi(xc + hi * c * dir, hi * s) < PHI_L1) { pts.push(hi * c, hi * s); continue; }
      for (var i = 0; i < 36; i++) {
        var m = (lo + hi) / 2;
        if (phi(xc + m * c * dir, m * s) < PHI_L1) lo = m; else hi = m;
      }
      pts.push(lo * c, lo * s);
    }
    return pts;
  }

  function setGeometry(tm) {
    Q = donorMass(tm) / M_BH;
    MU = Q / (1 + Q); MA = 1 - MU; MD = MU;
    XA = -MU; XD = 1 - MU;
    var lo = XA + 1e-4, hi = XD - 1e-4;
    for (var i = 0; i < 70; i++) { var m = (lo + hi) / 2; if (dphidx(m) > 0) lo = m; else hi = m; }
    L1 = (lo + hi) / 2;
    PHI_L1 = phi(L1, 0);
    LOBE_D = lobe(XD, -1);
    LOBE_A = lobe(XA, 1);
    R_DONOR = Math.abs(L1 - XD);
    R_EGG_D = eggleton(Q);
    R_TIDAL = 0.8 * eggleton(1 / Q);
  }

  // ---- Particle state (corotating frame) ----
  var PX = new Float64Array(NP), PY = new Float64Array(NP), PZ = new Float64Array(NP);
  var VX = new Float64Array(NP), VY = new Float64Array(NP), VZ = new Float64Array(NP);
  var ACT = new Uint8Array(NP), DISK = new Uint8Array(NP), WAIT = new Float32Array(NP), ETA = new Float32Array(NP);
  var SXs = new Float32Array(NP), SYs = new Float32Array(NP), BK = new Uint8Array(NP);
  var A = new Float64Array(3);

  // Disk-density grid in the orbital plane, centred on the accretor
  var NG = 52, GSPAN = 1.3, GCELL = GSPAN / NG;
  var GRID = new Uint16Array(NG * NG);

  // ---- Screen state ----
  var W = 0, H = 0, dpr = 1, CX = 0, CY = 0, S = 1, fadeStart = 0, fadeEnd = 0, mobile = false;
  var CELLPX = 4, sgw = 0, sgh = 0, sgx0 = 0, sgy0 = 0, SGRID = new Uint16Array(1);
  var starsCanvas = null, accSprite = null;
  var dark = true, ink = '', RAMP = [], RAMP_A = [];
  var phase = 0.6, t = 0, last = 0, running = false, wantRunning = false, rafId = 0;
  var cph = 1, sph = 0;
  var timeEl = document.getElementById('bgTime'), qEl = document.getElementById('bgQ'), lastLabel = -1;

  function gauss() { return (Math.random() + Math.random() + Math.random() - 1.5) * 1.15; }

  // ---- Initial conditions ----
  function launch(i) {
    PX[i] = L1 - 0.004; PY[i] = gauss() * 0.008; PZ[i] = gauss() * 0.008;
    VX[i] = -CS * (0.8 + 0.4 * Math.random()); VY[i] = gauss() * 0.008; VZ[i] = gauss() * 0.008;
    DISK[i] = 0; ACT[i] = 1;
  }

  function seedDisk(i) {
    var r = 0.04 + (R_TIDAL * 0.75 - 0.04) * Math.pow(Math.random(), 0.8), a = Math.random() * 2 * Math.PI;
    var dx = r * Math.cos(a), dy = r * Math.sin(a);
    var x = XA + dx, y = dy, vk = Math.sqrt(MA / r);
    PX[i] = x; PY[i] = y; PZ[i] = gauss() * 0.004;
    VX[i] = -vk * dy / r + y;   // inertial circular velocity minus Ω × r
    VY[i] = vk * dx / r - x;
    VZ[i] = 0; DISK[i] = 1; ACT[i] = 1;
  }

  function init() {
    setGeometry(tMyr); geomAt = tMyr;
    for (var i = 0; i < NP; i++) {
      if (i < NP * 0.68) seedDisk(i);
      else { ACT[i] = 0; WAIT[i] = Math.random() * 4; } // stream builds up over the first seconds
    }
  }

  // ---- Dynamics ----
  // Acceleration in the corotating frame: gravity + centrifugal + Coriolis
  function acc(x, y, z, vx, vy) {
    var d1 = x - XA, d2 = x - XD;
    var r1s = d1 * d1 + y * y + z * z + 1e-6, r2s = d2 * d2 + y * y + z * z + 1e-6;
    var a1 = MA / (r1s * Math.sqrt(r1s)), a2 = MD / (r2s * Math.sqrt(r2s));
    A[0] = 2 * vy + x - a1 * d1 - a2 * d2;
    A[1] = -2 * vx + y - (a1 + a2) * y;
    A[2] = -(a1 + a2) * z;
  }

  function rk4(i, h) {
    var x = PX[i], y = PY[i], z = PZ[i], vx = VX[i], vy = VY[i], vz = VZ[i];
    acc(x, y, z, vx, vy);
    var ax1 = A[0], ay1 = A[1], az1 = A[2];
    var x2 = x + 0.5 * h * vx, y2 = y + 0.5 * h * vy, z2 = z + 0.5 * h * vz;
    var vx2 = vx + 0.5 * h * ax1, vy2 = vy + 0.5 * h * ay1, vz2 = vz + 0.5 * h * az1;
    acc(x2, y2, z2, vx2, vy2);
    var ax2 = A[0], ay2 = A[1], az2 = A[2];
    var x3 = x + 0.5 * h * vx2, y3 = y + 0.5 * h * vy2, z3 = z + 0.5 * h * vz2;
    var vx3 = vx + 0.5 * h * ax2, vy3 = vy + 0.5 * h * ay2, vz3 = vz + 0.5 * h * az2;
    acc(x3, y3, z3, vx3, vy3);
    var ax3 = A[0], ay3 = A[1], az3 = A[2];
    var x4 = x + h * vx3, y4 = y + h * vy3, z4 = z + h * vz3;
    var vx4 = vx + h * ax3, vy4 = vy + h * ay3, vz4 = vz + h * az3;
    acc(x4, y4, z4, vx4, vy4);
    var h6 = h / 6;
    PX[i] = x + h6 * (vx + 2 * vx2 + 2 * vx3 + vx4);
    PY[i] = y + h6 * (vy + 2 * vy2 + 2 * vy3 + vy4);
    PZ[i] = z + h6 * (vz + 2 * vz2 + 2 * vz3 + vz4);
    VX[i] = vx + h6 * (ax1 + 2 * ax2 + 2 * ax3 + A[0]);
    VY[i] = vy + h6 * (ay1 + 2 * ay2 + 2 * ay3 + A[1]);
    VZ[i] = vz + h6 * (az1 + 2 * az2 + 2 * az3 + A[2]);
  }

  function buildDiskGrid() {
    GRID.fill(0);
    for (var i = 0; i < NP; i++) {
      if (!ACT[i] || !DISK[i]) continue;
      var gx = ((PX[i] - XA + GSPAN / 2) / GCELL) | 0, gy = ((PY[i] + GSPAN / 2) / GCELL) | 0;
      if (gx >= 0 && gy >= 0 && gx < NG && gy < NG) GRID[gy * NG + gx]++;
    }
  }

  // Stream–disk shock + viscosity: relax toward a Keplerian orbit
  function dissipate(i, h) {
    var x = PX[i], y = PY[i];
    var dx = x - XA, r = Math.sqrt(dx * dx + y * y);
    var gx = ((dx + GSPAN / 2) / GCELL) | 0, gy = ((y + GSPAN / 2) / GCELL) | 0;
    var dens = (gx >= 0 && gy >= 0 && gx < NG && gy < NG) ? GRID[gy * NG + gx] : 0;
    var k = KMAX * Math.min(1, dens / DENS0);
    if (k <= 0) return;
    var vk = Math.sqrt(MA / r);
    // Viscosity acts as radial diffusion: a slow mean inflow plus a random
    // turbulent radial velocity that changes on roughly an orbital timescale.
    // Outside the tidal radius the companion's torques push gas back in.
    if (Math.random() < h * 2) ETA[i] = gauss();
    var vr = vk * (-INFLOW + DIFFUSE * ETA[i]);
    if (r > R_TIDAL && vr > 0) vr = -vr;
    var tx = vk * (-y / r) + vr * dx / r + y;   // target velocity, corotating frame
    var ty = vk * (dx / r) + vr * y / r - x;
    var f = 1 - Math.exp(-k * h);
    var ex = tx - VX[i], ey = ty - VY[i];
    VX[i] += ex * f; VY[i] += ey * f; VZ[i] -= VZ[i] * f;
    if (!DISK[i] && ex * ex + ey * ey < 0.09 * vk * vk) DISK[i] = 1;
  }

  function step(dtReal) {
    // Secular evolution: donor mass, q and Roche geometry
    if (tMyr < T_END) tMyr = Math.min(T_END, tMyr + dtReal * MYR_PER_S);
    if (Math.abs(tMyr - geomAt) >= 2 || (tMyr === T_END && geomAt !== T_END)) { setGeometry(tMyr); geomAt = tMyr; }

    var dts = dtReal * RATE;
    if (dts <= 0) return;
    var nsub = Math.ceil(dts / HMAX), h = dts / nsub;
    buildDiskGrid();
    for (var i = 0; i < NP; i++) {
      if (!ACT[i]) {
        WAIT[i] -= dtReal;
        if (WAIT[i] <= 0) launch(i);
        continue;
      }
      for (var s = 0; s < nsub; s++) {
        var dx = PX[i] - XA, r1 = Math.sqrt(dx * dx + PY[i] * PY[i]);
        if (r1 < 0.1) { rk4(i, h / 3); rk4(i, h / 3); rk4(i, h / 3); } else rk4(i, h);
        dissipate(i, h);
        dx = PX[i] - XA;
        r1 = Math.sqrt(dx * dx + PY[i] * PY[i] + PZ[i] * PZ[i]);
        var d2 = PX[i] - XD, r2 = Math.sqrt(d2 * d2 + PY[i] * PY[i] + PZ[i] * PZ[i]);
        var rc = PX[i] * PX[i] + PY[i] * PY[i];
        if (r1 < R_ACC || r2 < R_DONOR * 0.6 || rc > 4.8) { launch(i); break; } // accreted, fell back, or lost
      }
    }
    phase += dts;
    PA = PA0 + t * 2 * Math.PI / PA_PERIOD_S;
    cosP = Math.cos(PA); sinP = Math.sin(PA);
  }

  // ---- Labels ----
  function updateLabels(force) {
    if (!force && t - lastLabel < 0.1) return;
    lastLabel = t;
    if (timeEl) timeEl.textContent = 't = ' + (tMyr < 1000 ? Math.floor(tMyr) + ' Myr' : (tMyr / 1000).toFixed(2) + ' Gyr');
    if (qEl) qEl.textContent = 'q = ' + Q.toFixed(3);
  }

  // ---- Theme & layout ----
  function glow(stops) {
    var c = document.createElement('canvas'), sz = 64;
    c.width = c.height = sz;
    var g = c.getContext('2d'), gr = g.createRadialGradient(sz / 2, sz / 2, 0, sz / 2, sz / 2, sz / 2);
    for (var i = 0; i < stops.length; i++) gr.addColorStop(stops[i][0], stops[i][1]);
    g.fillStyle = gr; g.fillRect(0, 0, sz, sz);
    return c;
  }

  function readTheme() {
    dark = document.documentElement.dataset.theme === 'dark';
    ink = dark ? '236,232,225' : '30,28,25';
    RAMP = dark
      ? ['#ffd2a8', '#ffb07a', '#ff8a52', '#f0622e', '#d43d1c', '#b02012', '#8e0e0c', '#6e0606']
      : ['#f2b38a', '#e88a5a', '#d9603a', '#c03a22', '#9c1d14', '#7a0f0c', '#5c0808', '#3f0404'];
    RAMP_A = dark ? [0.42, 0.5, 0.58, 0.66, 0.74, 0.82, 0.9, 0.96] : [0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88, 0.96];
    buildStars();
    // Hot inner accretion flow around the black hole
    accSprite = glow([
      [0, 'rgba(255,255,255,1)'],
      [0.12, dark ? 'rgba(205,220,255,0.95)' : 'rgba(120,150,230,0.95)'],
      [0.35, dark ? 'rgba(160,185,255,0.25)' : 'rgba(80,110,210,0.2)'],
      [1, 'rgba(160,185,255,0)']
    ]);
  }

  function buildStars() {
    if (!W) return;
    starsCanvas = document.createElement('canvas');
    starsCanvas.width = Math.round(W * dpr); starsCanvas.height = Math.round(H * dpr);
    var g = starsCanvas.getContext('2d');
    g.scale(dpr, dpr);
    var n = Math.round(W * H / 9000);
    for (var i = 0; i < n; i++) {
      var r = Math.random() < 0.9 ? 0.6 : 1.1;
      g.globalAlpha = 0.12 + Math.random() * (dark ? 0.4 : 0.25);
      g.fillStyle = 'rgb(' + ink + ')';
      g.beginPath(); g.arc(Math.random() * W, Math.random() * H, r, 0, 6.2832); g.fill();
    }
  }

  // Keep 10–15% of the donor off the left edge at its leftmost point
  // (worst case over the slow rotation, so the system doesn't drift)
  function placeSystem() {
    if (mobile) return;
    CX = S * (XD + 0.75 * R_EGG_D);
  }

  function layout() {
    W = canvas.clientWidth || window.innerWidth;
    H = canvas.clientHeight || window.innerHeight;
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    mobile = W < 820;
    if (mobile) {
      S = Math.min(W * 0.37, H * 0.32); CX = W * 0.42; CY = H * 0.28;
      fadeStart = W * 2; fadeEnd = W * 3;
    } else {
      S = Math.min(W * 0.184, H * 0.37); CY = H * 0.67;
      fadeStart = W * 0.4; fadeEnd = W * 0.6;
    }
    placeSystem();
    readTheme();
  }

  // ---- Drawing ----
  // corotating (x, y, z) → screen, using the current orbital phase
  function projX(x, y, z) {
    var xr = x * cph - y * sph, yr = x * sph + y * cph;
    var sy = yr * cosI + z * sinI;
    return CX + (xr * cosP - sy * sinP) * S;
  }
  function projY(x, y, z) {
    var xr = x * cph - y * sph, yr = x * sph + y * cph;
    var sy = yr * cosI + z * sinI;
    return CY - (xr * sinP + sy * cosP) * S;
  }

  // Lobes are close to bodies of revolution about the binary axis, so the
  // silhouette = outline with the axis foreshortened by the viewing geometry.
  function lobePath(pts, xc, dir) {
    var ox = projX(xc, 0, 0), oy = projY(xc, 0, 0);
    var ax = projX(xc + dir, 0, 0) - ox, ay = projY(xc + dir, 0, 0) - oy;
    var L = Math.sqrt(ax * ax + ay * ay), ux = ax / L, uy = ay / L, nx = -uy, ny = ux;
    ctx.beginPath();
    for (var k = 0; k < pts.length; k += 2) {
      var px = ox + ux * pts[k] * L + nx * pts[k + 1] * S;
      var py = oy + uy * pts[k] * L + ny * pts[k + 1] * S;
      if (k) ctx.lineTo(px, py); else ctx.moveTo(px, py);
    }
    ctx.closePath();
    return { ox: ox, oy: oy, ux: ux, uy: uy, L: L };
  }

  function drawDonor() {
    var g = lobePath(LOBE_D, XD, -1);
    // Hotter progenitor → cool K giant as it is stripped
    var k = 1 - (Q - M_DONOR_NOW / M_BH) / ((M_DONOR_0 - M_DONOR_NOW) / M_BH);
    var hx = g.ox + g.ux * R_DONOR * 0.45 * g.L, hy = g.oy + g.uy * R_DONOR * 0.45 * g.L;
    var gr = ctx.createRadialGradient(hx, hy, 0, g.ox, g.oy, R_DONOR * S * 1.05);
    gr.addColorStop(0, '#fff0d8');
    gr.addColorStop(0.35, k < 0.5 ? '#ffc98f' : '#f59a5c');
    gr.addColorStop(0.75, k < 0.5 ? '#e0844a' : '#c2461f');
    gr.addColorStop(1, '#7a1e10');
    ctx.globalAlpha = dark ? 0.5 : 0.45;
    ctx.fillStyle = gr;
    ctx.fill();
  }

  function drawLobes() {
    ctx.globalAlpha = dark ? 0.22 : 0.2;
    ctx.strokeStyle = 'rgb(' + ink + ')';
    ctx.lineWidth = 0.8;
    ctx.setLineDash([3, 4]);
    lobePath(LOBE_A, XA, 1); ctx.stroke();
    lobePath(LOBE_D, XD, -1); ctx.stroke();
    ctx.setLineDash([]);
    ctx.globalAlpha = 0.5;
    ctx.fillStyle = 'rgb(' + ink + ')';
    ctx.beginPath(); ctx.arc(projX(L1, 0, 0), projY(L1, 0, 0), 1.6, 0, 6.2832); ctx.fill();
  }

  function drawGas() {
    var span = 1.4 * S;
    sgx0 = CX - span; sgy0 = CY - span;
    var need = Math.ceil(2 * span / CELLPX) + 1;
    if (need !== sgw) { sgw = sgh = need; SGRID = new Uint16Array(sgw * sgh); } else SGRID.fill(0);
    var i;
    for (i = 0; i < NP; i++) {
      if (!ACT[i]) continue;
      var x = projX(PX[i], PY[i], PZ[i]), y = projY(PX[i], PY[i], PZ[i]);
      SXs[i] = x; SYs[i] = y;
      var gx = ((x - sgx0) / CELLPX) | 0, gy = ((y - sgy0) / CELLPX) | 0;
      if (gx >= 0 && gy >= 0 && gx < sgw && gy < sgh) SGRID[gy * sgw + gx]++;
    }
    // Column density Σ: 3×3 cell sum, log-scaled into 8 colour buckets
    var smax = 1;
    for (i = 0; i < NP; i++) {
      if (!ACT[i]) continue;
      var cx = ((SXs[i] - sgx0) / CELLPX) | 0, cy = ((SYs[i] - sgy0) / CELLPX) | 0, sum = 0;
      for (var oy = -1; oy <= 1; oy++) {
        var yy = cy + oy; if (yy < 0 || yy >= sgh) continue;
        for (var ox = -1; ox <= 1; ox++) {
          var xx = cx + ox; if (xx < 0 || xx >= sgw) continue;
          sum += SGRID[yy * sgw + xx];
        }
      }
      BK[i] = sum > 255 ? 255 : sum;
      if (sum > smax) smax = sum;
    }
    var norm = 1 / Math.log(1 + Math.max(10, smax * 0.85));
    for (i = 0; i < NP; i++) if (ACT[i]) BK[i] = Math.min(7, (Math.log(1 + BK[i]) * norm * 8) | 0);
    var sz = mobile ? 1.6 : 2.2, half = sz / 2;
    for (var b = 0; b < 8; b++) {
      ctx.globalAlpha = RAMP_A[b];
      ctx.fillStyle = RAMP[b];
      ctx.beginPath();
      for (i = 0; i < NP; i++) {
        if (ACT[i] && BK[i] === b) ctx.rect(SXs[i] - half, SYs[i] - half, sz, sz);
      }
      ctx.fill();
    }
    // Black hole + hot inner flow
    var ax = projX(XA, 0, 0), ay = projY(XA, 0, 0), s = mobile ? 17 : 25;
    ctx.globalAlpha = 1;
    ctx.drawImage(accSprite, ax - s / 2, ay - s / 2, s, s);
    ctx.fillStyle = '#000';
    ctx.beginPath(); ctx.arc(ax, ay, 1.7, 0, 6.2832); ctx.fill();
  }

  function draw() {
    cph = Math.cos(phase); sph = Math.sin(phase);
    placeSystem();
    ctx.globalAlpha = 1;
    ctx.globalCompositeOperation = 'source-over';
    ctx.clearRect(0, 0, W, H);

    drawLobes();
    // Depth order: toward-viewer coordinate is −y_inertial·sin i, so the donor
    // (x > 0) is in front of the accretor when sin(phase) < 0
    if (sph < 0) { drawGas(); drawDonor(); }
    else { drawDonor(); drawGas(); }

    // Fade out toward the middle of the page
    if (!mobile) {
      ctx.globalAlpha = 1;
      ctx.globalCompositeOperation = 'destination-out';
      var fg = ctx.createLinearGradient(fadeStart, 0, fadeEnd, 0);
      fg.addColorStop(0, 'rgba(0,0,0,0)');
      fg.addColorStop(1, 'rgba(0,0,0,1)');
      ctx.fillStyle = fg;
      ctx.fillRect(fadeStart, 0, W - fadeStart, H);
    }
    // Faint star backdrop behind everything
    ctx.globalCompositeOperation = 'destination-over';
    if (starsCanvas) ctx.drawImage(starsCanvas, 0, 0, W, H);
    ctx.globalCompositeOperation = 'source-over';
    ctx.globalAlpha = 1;
    updateLabels(false);
  }

  // ---- Loop ----
  function frame(now) {
    var dt = last ? (now - last) / 1000 : 0.016;
    last = now;
    if (dt > 0.05) dt = 0.05;
    t += dt;
    step(dt);
    draw();
    if (running) rafId = requestAnimationFrame(frame);
  }
  function play() {
    wantRunning = true;
    if (running || document.hidden) return;
    running = true; last = 0;
    rafId = requestAnimationFrame(frame);
  }
  function stop() { running = false; cancelAnimationFrame(rafId); }
  function pause() { wantRunning = false; stop(); }
  function staticFrame() { draw(); updateLabels(true); }

  // ---- Events ----
  var resizeTimer;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () { layout(); if (!running) staticFrame(); }, 120);
  });
  document.addEventListener('visibilitychange', function () {
    if (document.hidden) stop(); else if (wantRunning) play();
  });
  new MutationObserver(function () { readTheme(); if (!running) staticFrame(); })
    .observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  var btn = document.getElementById('bgToggle');
  var paused = false;
  try { paused = localStorage.getItem('bgPaused') === '1'; } catch (e) {}
  function syncBtn() { if (btn) btn.textContent = wantRunning ? 'Pause background' : 'Play background'; }
  if (btn) btn.addEventListener('click', function () {
    if (wantRunning) { pause(); staticFrame(); } else { play(); }
    try { localStorage.setItem('bgPaused', wantRunning ? '0' : '1'); } catch (e) {}
    syncBtn();
  });

  // Test hook: fast-forward the simulation
  window.__roche = {
    advance: function (sec) {
      for (var i = 0; i < sec * 60; i++) { t += 1 / 60; step(1 / 60); }
      draw(); updateLabels(true);
      var a = 0, d = 0; for (var k = 0; k < NP; k++) { a += ACT[k]; d += ACT[k] && DISK[k]; }
      return { tMyr: Math.round(tMyr), q: +Q.toFixed(4), active: a, disk: d, stream: a - d, CX: Math.round(CX), S: Math.round(S) };
    }
  };

  // ---- Start ----
  init();
  layout();
  if (reduceMotion || paused) {
    // Still frame: show the present-day system with a settled flow
    tMyr = T_END; setGeometry(tMyr); geomAt = tMyr;
    for (var i = 0; i < 150; i++) step(1 / 60);
    staticFrame();
  } else play();
  updateLabels(true);
  syncBtn();
})();
