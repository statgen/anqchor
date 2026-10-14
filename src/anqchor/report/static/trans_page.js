/* anqchor trans report: hotspot controls, bins along the genome, the hotspot grid and
   table, and the cross-chromosome QQ plot, drawn from ANQ_DATA.panels (one panel per
   dataset, normally one), or, for a comparison of runs, from ANQ_DATA.comparison: shared
   hotspot controls with the hotspot count of every run, the bins, grid and table of one
   selected run, λGC per run and the QQ plot of every run. Figures go through
   ANQ.renderFigure, which re-themes them for dark mode; the grid and tables use CSS colors
   that follow the theme. */
(function () {
  "use strict";
  var D = window.ANQ_DATA || {};
  var PANELS = D.panels || [];
  var CMP = D.comparison || null;

  /* ------------------------------------------------------------- helpers */
  function isNum(v) { return typeof v === "number" && isFinite(v); }
  function n0(v) { return Math.round(v).toLocaleString("en-US"); }
  function clone(o) { return JSON.parse(JSON.stringify(o)); }
  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined && text !== null) e.textContent = String(text);
    return e;
  }
  function phone() { return window.innerWidth < 640; }
  function fmtLp(v) { return v >= 100 ? v.toFixed(0) : v.toFixed(1); }
  function mb(bp) { var v = bp / 1e6; return Math.abs(v - Math.round(v)) < 1e-9 ? String(Math.round(v)) : String(Number(v.toPrecision(4))); }
  function chromName(c) { c = String(c); return /^chr/i.test(c) ? c : "chr" + c; }
  function binName(P, i) {
    var H = P.hot, s = H.bins.start[i];
    return chromName(H.bins.chrom[i]) + ": " + mb(s) + "–" + mb(s + H.binBp) + " Mb";
  }
  function layoutFor(P, extra) {
    var lay = clone(P.layout);
    Object.keys(extra).forEach(function (k) {
      var v = extra[k];
      if (v && typeof v === "object" && !Array.isArray(v) && lay[k] && typeof lay[k] === "object") {
        var m = {};
        Object.keys(lay[k]).forEach(function (q) { m[q] = lay[k][q]; });
        Object.keys(v).forEach(function (q) { m[q] = v[q]; });
        lay[k] = m;
      } else { lay[k] = v; }
    });
    return lay;
  }

  /* one tooltip element (the cis report's style) */
  var tip = null;
  function showTip(evt, title, rows) {
    if (!tip) { tip = el("div", "cis-tip"); tip.setAttribute("role", "status"); document.body.appendChild(tip); }
    tip.innerHTML = "";
    if (title) tip.appendChild(el("div", "tt-title", title));
    rows.forEach(function (r) {
      var row = el("div", "tt-row");
      if (r.color) {
        var k = el("span", "tt-key");
        k.style.setProperty("--c", ANQ.color(r.color));
        row.appendChild(k);
      }
      row.appendChild(el("b", "tt-v", r.value));
      if (r.sub) row.appendChild(el("span", "tt-sub", r.sub));
      if (r.label) row.appendChild(el("span", "tt-l", r.label));
      tip.appendChild(row);
    });
    tip.style.display = "block";
    var e = evt || {}, x = e.clientX || 0, y = e.clientY || 0;
    var w = tip.offsetWidth, h = tip.offsetHeight;
    var left = x + 16 + w > window.innerWidth - 8 ? x - 16 - w : x + 16;
    tip.style.left = Math.max(8, Math.min(window.innerWidth - w - 8, left)) + "px";
    tip.style.top = Math.max(8, Math.min(window.innerHeight - h - 8, y - h / 2)) + "px";
  }
  function hideTip() { if (tip) tip.style.display = "none"; }
  function onHover(id, fn) {
    var gd = document.getElementById(id);
    if (!gd || !gd.on || gd._trHover) return;
    gd._trHover = true;
    gd.on("plotly_hover", function (ev) { try { fn(ev); } catch (e) { console.error(e); } });
    gd.on("plotly_unhover", hideTip);
    gd.addEventListener("mouseleave", hideTip);
  }

  /* ------------------------------------------------------------ state */
  /* index of the selected threshold in P.hot (-1: not in this run's grid) */
  function tIndex(P) {
    var ths = P.hot.thresholds;
    for (var i = 0; i < ths.length; i++) if (ths[i].key === P.tKey) return i;
    return P.cmp ? -1 : 0;
  }
  function curTh(P) {
    var i = tIndex(P);
    return i >= 0 ? P.hot.thresholds[i] : { key: P.tKey, label: P.tKey, lp: Infinity };
  }
  function counts(P) { return P.hot.bins.n[P.tKey] || []; }
  function nHot(P) {
    var i = tIndex(P), j = P.hot.k.indexOf(P.k);
    var row = P.hot.counts[i] || [];
    return j >= 0 && isNum(row[j]) ? row[j] : 0;
  }
  function select(P, tKey, k) {
    if (tKey !== null) P.tKey = tKey;
    if (k !== null) P.k = k;
    document.querySelectorAll("[data-hot-t='" + P.sfx + "'] button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.getAttribute("data-v") === P.tKey));
    });
    document.querySelectorAll("[data-hot-k='" + P.sfx + "'] button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(Number(b.getAttribute("data-v")) === P.k));
    });
    P.page = 0;
    drawHot(P);
  }
  function drawHot(P) {
    if (P.after) P.after(P);
    drawHero(P);
    drawMan(P);
    drawGrid(P);
    drawTable(P);
  }

  /* ------------------------------------------------------------ hero */
  function drawHero(P) {
    var n = nHot(P), th = curTh(P);
    var bin = mb(P.hot.binBp) + "-Mb";
    var num = document.querySelector("[data-hot-n='" + P.sfx + "']");
    var say = document.querySelector("[data-hot-say='" + P.sfx + "']");
    if (num) num.textContent = n0(n);
    if (say && P.cmp) {
      say.textContent = "Hotspots: " + (P.sameBin === false ? "bins (width per run)" : bin + " bins") + " with " +
        (P.k === 1 ? "a trans signal" : "trans signals for ≥ " + P.k + " genes") + " at " + th.label;
      return;
    }
    if (say) {
      say.textContent = (n === 1 ? "hotspot" : "hotspots") + ": " + bin + " bins with " +
        (P.k === 1 ? "a trans signal" : "trans signals for ≥ " + P.k + " genes") + " at " + th.label;
    }
  }

  /* ------------------------------------------- bins along the genome */
  function drawMan(P) {
    var id = "fig-man" + P.sfx, H = P.hot, B = H.bins, N = counts(P), C = P.colors;
    if (!document.getElementById(id)) return;
    var hot = { x: [], y: [], i: [] }, rest = { x: [], y: [], i: [] };
    var maxN = 0;
    for (var i = 0; i < B.x.length; i++) {
      var v = N[i] || 0;
      maxN = Math.max(maxN, v);
      var g = v >= P.k ? hot : rest;
      g.x.push(B.x[i]); g.y.push(v); g.i.push(i);
    }
    function bars(g, color, name) {
      return { type: "bar", x: g.x, y: g.y, customdata: g.i, width: H.binBp * 0.9, name: name,
               marker: { color: color, line: { width: 0 } }, hoverinfo: "none", showlegend: false };
    }
    var top = Math.max(maxN, P.k, 1) * 1.12;
    var shapes = [], ticks = [], labels = [];
    var narrow = phone();
    H.chroms.forEach(function (c, j) {
      if (j % 2 === 1) {
        shapes.push({ type: "rect", xref: "x", yref: "paper", x0: c.start, x1: c.end, y0: 0, y1: 1,
                      fillcolor: C.grid, opacity: 0.45, line: { width: 0 }, layer: "below" });
      }
      ticks.push((c.start + c.end) / 2);
      labels.push(narrow && H.chroms.length > 12 ? String(c.name).replace(/^chr/i, "") : chromName(c.name));
    });
    shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: P.k, y1: P.k,
                  line: { color: C.muted, width: 1, dash: "dot" } });
    var last = H.chroms.length ? H.chroms[H.chroms.length - 1].end : 1;
    var lay = layoutFor(P, {
      height: 360, showlegend: false, hovermode: "closest", bargap: 0, shapes: shapes,
      margin: { l: 64, r: 16, t: 12, b: 48 },
      xaxis: { range: [0, last], tickvals: ticks, ticktext: labels, showgrid: false, zeroline: false,
               title: { text: "" } },
      yaxis: { range: [0, top], rangemode: "tozero", title: { text: "Genes with a trans signal" } }
    });
    if (maxN === 0) {
      lay.annotations = [{ xref: "paper", yref: "paper", x: 0.5, y: 0.5, showarrow: false,
                           text: "No trans signal at " + curTh(P).label,
                           font: { color: C.ink2, size: 16 } }];
    }
    /* other bins: muted on a run page; in a comparison a lighter neutral that no run uses
       (runs past the eighth are drawn in the greys) */
    var restColor = P.cmp ? C.axis : C.muted;
    ANQ.renderFigure(id, { data: [bars(rest, restColor, "other bins"), bars(hot, C.trans, "hotspots")], layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points[0];
      if (!pt) return;
      var i = pt.customdata, th = curTh(P), n = counts(P)[i] || 0;
      var rows = [{ color: n >= P.k ? C.trans : restColor, value: n0(n), sub: (n === 1 ? "gene" : "genes") + " · " + th.label, label: "" }];
      var genes = (B.genes[i] || []).filter(function (g) { return g[1] > th.lp; });
      genes.slice(0, 8).forEach(function (g) { rows.push({ value: fmtLp(g[1]), sub: "", label: g[0] }); });
      if (genes.length > 8) rows.push({ value: "+" + n0(n - 8), sub: "more", label: "" });
      showTip(ev.event, binName(P, i) + " · " + n0(B.nVar[i]) + (B.nVar[i] === 1 ? " variant" : " variants"), rows);
    });
  }

  /* ------------------------------------------------------ hotspot grid */
  function drawGrid(P) {
    var tbl = document.querySelector("[data-hot-grid='" + P.sfx + "']");
    if (!tbl) return;
    var H = P.hot, maxV = 1;
    H.counts.forEach(function (r) { r.forEach(function (v) { maxV = Math.max(maxV, v || 0); }); });
    tbl.innerHTML = "";
    var thead = el("thead"), tr = el("tr");
    tr.appendChild(el("th", "corner", "p <"));
    H.k.forEach(function (k) { var th = el("th", "num", "≥ " + k); th.scope = "col"; tr.appendChild(th); });
    thead.appendChild(tr);
    tbl.appendChild(thead);
    var tb = el("tbody");
    H.thresholds.forEach(function (t, i) {
      var row = el("tr");
      var th = el("th", "rowh", t.short);
      th.scope = "row"; th.title = t.label;
      row.appendChild(th);
      H.k.forEach(function (k, j) {
        var v = (H.counts[i] || [])[j] || 0;
        var td = el("td", "num", n0(v));
        /* magnitude: one hue mixed into the surface (follows the theme) */
        var f = v > 0 ? 8 + 52 * Math.log(1 + v) / Math.log(1 + maxV) : 0;
        td.style.background = f ? "color-mix(in oklab, " + (P.gridVar || "var(--series-1)") + " " + f.toFixed(0) + "%, var(--surface))" : "var(--surface)";
        if (!v) td.classList.add("zero");
        if (t.key === P.tKey && k === P.k) { td.classList.add("sel"); td.setAttribute("aria-current", "true"); }
        td.tabIndex = 0;
        td.title = n0(v) + (v === 1 ? " bin" : " bins") + " with ≥ " + k + (k === 1 ? " gene" : " genes") + " at " + t.label;
        var act = function () { select(P, t.key, k); };
        td.addEventListener("click", act);
        td.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); } });
        row.appendChild(td);
      });
      tb.appendChild(row);
    });
    tbl.appendChild(tb);
  }

  /* ----------------------------------------------------- hotspot table */
  var CMP_TOP = 10; /* comparison: bins shown before "Show all" */
  function drawTable(P) {
    var tbl = document.querySelector("[data-hot-tbl='" + P.sfx + "']");
    if (!tbl) return;
    var H = P.hot, B = H.bins, N = counts(P), th = curTh(P), p = phone();
    var rows = [];
    for (var i = 0; i < N.length; i++) if ((N[i] || 0) >= P.k) rows.push(i);
    rows.sort(function (a, b) { return (N[b] - N[a]) || (a - b); });
    var cnt = document.querySelector("[data-hot-cnt='" + P.sfx + "']");
    if (cnt) cnt.textContent = "(" + n0(rows.length) + ")";
    var cols = [{ k: "bin", label: "Bin", w: 200, pw: 128 },
                { k: "nv", label: "Variants", num: true, w: 96, pw: 72 },
                { k: "ng", label: "Genes", num: true, w: 80, pw: 64 },
                { k: "lp", label: "Top −log10 p", plabel: "Top", num: true, w: 120, pw: 60 },
                { k: "genes", label: "Genes with a trans signal", w: 520, pw: 260 }];
    var cg = tbl.querySelector("colgroup"), total = 0;
    cg.innerHTML = "";
    cols.forEach(function (c) { var col = el("col"), w = p ? c.pw : c.w; col.style.width = w + "px"; total += w; cg.appendChild(col); });
    tbl.style.minWidth = total + "px";
    var thead = tbl.querySelector("thead");
    thead.innerHTML = "";
    var tr = el("tr");
    cols.forEach(function (c) { var h = el("th", c.num ? "num" : "", p && c.plabel ? c.plabel : c.label); h.scope = "col"; tr.appendChild(h); });
    thead.appendChild(tr);
    var collapsed = P.cmp && !P.showAll && rows.length > CMP_TOP;
    var size = collapsed ? CMP_TOP : (P.size || 100), pages = Math.max(1, Math.ceil(rows.length / (size || rows.length || 1)));
    if (P.page >= pages) P.page = pages - 1;
    var start = size ? P.page * size : 0, end = size ? Math.min(rows.length, start + size) : rows.length;
    var tbody = tbl.querySelector("tbody");
    tbody.innerHTML = "";
    if (!rows.length) {
      var r0 = el("tr"), td0 = el("td", "muted", "No bin reaches ≥ " + P.k + (P.k === 1 ? " gene" : " genes") + " at " + th.label + ".");
      td0.colSpan = cols.length; r0.appendChild(td0); tbody.appendChild(r0);
    }
    for (var q = start; q < end; q++) {
      var i2 = rows[q], genes = (B.genes[i2] || []).filter(function (g) { return g[1] > th.lp; });
      var row = el("tr");
      row.appendChild(el("td", "", binName(P, i2)));
      row.appendChild(el("td", "num", n0(B.nVar[i2])));
      row.appendChild(el("td", "num", n0(N[i2])));
      row.appendChild(el("td", "num", genes.length ? fmtLp(genes[0][1]) : "–"));
      var names = genes.map(function (g) { return g[0]; });
      var more = N[i2] - genes.length;
      var td = el("td", "glist", names.join(", ") + (more > 0 ? " … +" + n0(more) : ""));
      td.title = genes.map(function (g) { return g[3] + "  " + fmtLp(g[1]) + "  (" + g[2] + ")"; }).join("\n");
      row.appendChild(td);
      tbody.appendChild(row);
    }
    var foot = document.querySelector("[data-hot-foot='" + P.sfx + "']");
    if (foot) {
      foot.innerHTML = "";
      if (collapsed) {
        var more = el("button", "show-all", "Show all " + n0(rows.length) + " bins");
        more.type = "button";
        more.addEventListener("click", function () { P.showAll = true; P.page = 0; drawTable(P); });
        foot.appendChild(more);
      } else if (rows.length > 100) {
        var pg = el("span", "tpager");
        var prev = el("button", "", "‹"), info = el("span", "num", n0(start + 1) + "–" + n0(end) + " of " + n0(rows.length)), next = el("button", "", "›");
        prev.type = next.type = "button";
        prev.setAttribute("aria-label", "Previous page"); next.setAttribute("aria-label", "Next page");
        prev.disabled = P.page === 0; next.disabled = P.page >= pages - 1;
        prev.addEventListener("click", function () { P.page = Math.max(0, P.page - 1); drawTable(P); });
        next.addEventListener("click", function () { P.page += 1; drawTable(P); });
        pg.appendChild(prev); pg.appendChild(info); pg.appendChild(next);
        var sz = el("span", "tsize");
        sz.appendChild(document.createTextNode("Rows"));
        [[100, "100"], [500, "500"], [0, "All"]].forEach(function (o) {
          var b = el("button", "", o[1]);
          b.type = "button";
          b.setAttribute("aria-pressed", String((P.size === undefined ? 100 : P.size) === o[0]));
          b.addEventListener("click", function () { P.size = o[0]; P.page = 0; drawTable(P); });
          sz.appendChild(b);
        });
        foot.appendChild(sz); foot.appendChild(pg);
      }
    }
  }

  /* --------------------------------------------------------------- QQ */
  function drawQQ(P) {
    var id = "fig-qq" + P.sfx, Q = P.qq, C = P.colors;
    if (!document.getElementById(id) || !Q || !Q.exp) return;
    var mx = 0, my = 0;
    Q.exp.forEach(function (v) { if (isNum(v)) mx = Math.max(mx, v); });
    Q.obs.forEach(function (v) { if (isNum(v)) my = Math.max(my, v); });
    Q.hi.forEach(function (v) { if (isNum(v)) my = Math.max(my, v); });
    var xr = mx * 1.04, yr = Math.max(my, mx) * 1.04;
    var data = [
      { type: "scatter", mode: "lines", x: Q.exp, y: Q.lo, line: { width: 0, color: C.grid }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", x: Q.exp, y: Q.hi, line: { width: 0, color: C.grid }, fill: "tonexty",
        fillcolor: C.grid, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "markers", x: Q.exp, y: Q.obs, hoverinfo: "none", showlegend: false,
        marker: { color: C.trans, size: 7, line: { width: 0 } } }
    ];
    var lay = layoutFor(P, {
      height: 460, showlegend: false, hovermode: "closest",
      margin: { l: 64, r: 16, t: 12, b: 56 },
      shapes: [{ type: "line", x0: 0, y0: 0, x1: Math.min(xr, yr), y1: Math.min(xr, yr),
                 line: { color: C.muted, width: 1, dash: "dot" } }],
      xaxis: { range: [0, xr], title: { text: "Expected −log10 p" }, zeroline: false },
      yaxis: { range: [0, yr], title: { text: "Observed −log10 p" }, zeroline: false }
    });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points.filter(function (x) { return x.curveNumber === 2; })[0];
      if (!pt) return;
      var i = pt.pointIndex;
      showTip(ev.event, "Cross-chromosome cells", [
        { color: C.trans, value: fmtLp(Q.obs[i]), sub: "observed −log10 p", label: "" },
        { value: fmtLp(Q.exp[i]), sub: "expected", label: "" },
        { value: fmtLp(Q.lo[i]) + "–" + fmtLp(Q.hi[i]), sub: "95% range", label: "" }
      ]);
    });
  }


  /* ========================================================== comparison */
  function runCount(C, r) {
    var H = r.hot, i = -1, j = H.k.indexOf(C.P.k);
    for (var q = 0; q < H.thresholds.length; q++) if (H.thresholds[q].key === C.P.tKey) i = q;
    if (i < 0 || j < 0) return null;
    var v = (H.counts[i] || [])[j];
    return isNum(v) ? v : null;
  }
  /* hotspot count per run at the shared threshold and minimum (rows = runs) */
  function cmpRunBars(C) {
    var box = document.querySelector("[data-hot-runs='" + C.sfx + "']");
    if (!box) return;
    box.innerHTML = "";
    box.style.setProperty("--ncol", "1");
    C.runs.forEach(function (r, i) {
      var n = runCount(C, r), row = el("div", "hits-row hr-row" + (i === C.sel ? " is-sel" : ""));
      row.style.setProperty("--c", r.css);
      row.tabIndex = 0;
      row.setAttribute("role", "button");
      row.setAttribute("aria-pressed", String(i === C.sel));
      row.title = r.label + ": show its bins below";
      var hn = el("span", "hn"), key = el("span", "key");
      key.setAttribute("aria-hidden", "true");
      hn.appendChild(key);
      hn.appendChild(el("span", "hn-l", r.label));
      row.appendChild(hn);
      var hc = el("span", "hc"), hv = el("span", "hv");
      hv.appendChild(el("b", "", n === null ? "–" : n0(n)));
      hv.appendChild(el("span", "hp", n === null ? "not computed for this run" : "of " + n0(r.nBins) + " bins"));
      hc.appendChild(hv);
      var tr = el("span", "track"), fill = el("span", "fill");
      fill.style.width = (n && r.nBins ? Math.min(100, 100 * n / r.nBins) : 0).toFixed(1) + "%";
      tr.appendChild(fill);
      hc.appendChild(tr);
      row.appendChild(hc);
      var act = function () { cmpSelect(C, i); };
      row.addEventListener("click", act);
      row.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); } });
      box.appendChild(row);
    });
  }
  function cmpSelect(C, i) {
    C.sel = i;
    var r = C.runs[i], P = C.P;
    P.hot = r.hot;
    P.colors.trans = r.color;
    P.gridVar = r.css;
    P.page = 0;
    P.showAll = false;
    var seg = document.querySelector("[data-seg='hot" + C.sfx + "']");
    if (seg) seg.querySelectorAll("button").forEach(function (b) { b.setAttribute("aria-pressed", String(Number(b.getAttribute("data-v")) === i)); });
    var pick = document.querySelector("[data-pick='hot" + C.sfx + "']");
    if (pick) pick.value = String(i);
    document.querySelectorAll("[data-hot-run='" + C.sfx + "']").forEach(function (e) { e.textContent = r.label; });
    if (!C.sameBin) document.querySelectorAll("[data-hot-bin='" + C.sfx + "']").forEach(function (e) { e.textContent = mb(r.hot.binBp) + "-Mb bins"; });
    drawHot(P);
  }

  /* λGC per run: one dot per run (rows), dashed line at 1 */
  function drawLambda(C) {
    var id = "fig-lambda" + C.sfx;
    if (!document.getElementById(id)) return;
    var runs = C.runs, keys = runs.map(function (r, i) { return "r" + i; }), vals = [1];
    runs.forEach(function (r) { if (isNum(r.lambda)) vals.push(r.lambda); });
    var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
    var pad = Math.max(0.05, (hi - lo) * 0.08);
    var data = runs.map(function (r, i) {
      return { type: "scatter", mode: "markers", y: [keys[i]], x: [isNum(r.lambda) ? r.lambda : null],
               marker: { color: r.color, size: 15, line: { color: C.colors.surface, width: 1.5 } },
               hoverinfo: "none", showlegend: false };
    });
    var n = phone() ? 12 : 28;
    var lay = layoutFor(C, {
      height: 72 + 48 * runs.length, hovermode: "closest", showlegend: false,
      margin: { l: 16, t: 8, r: 64, b: 40 },
      xaxis: { range: [Math.min(0.95, lo - pad), hi + pad], title: { text: "" }, zeroline: false },
      yaxis: { type: "category", categoryorder: "array", categoryarray: keys, autorange: "reversed",
               tickvals: keys, ticktext: runs.map(function (r) { return r.label.length > n ? r.label.slice(0, n - 1) + "…" : r.label; }),
               tickfont: { size: 15, color: C.colors.ink }, showgrid: false, title: { text: "" } },
      shapes: [{ type: "line", xref: "x", yref: "paper", x0: 1, x1: 1, y0: 0, y1: 1,
                 line: { color: C.colors.muted, width: 1, dash: "dot" } }],
      annotations: runs.map(function (r, i) {
        return { x: 1, y: keys[i], xref: "paper", yref: "y", text: isNum(r.lambda) ? r.lambda.toFixed(2) : "–",
                 showarrow: false, xanchor: "left", xshift: 10, font: { color: C.colors.ink, size: 15 } };
      })
    });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points[0];
      if (!pt) return;
      var r = runs[pt.curveNumber];
      showTip(ev.event, r.label, [{ color: r.color, value: isNum(r.lambda) ? r.lambda.toFixed(3) : "–", sub: "λGC", label: "" }]);
    });
  }

  /* QQ plot of every run: one line per run over the band of the run with the most cells */
  function drawQQcmp(C) {
    var id = "fig-qq" + C.sfx;
    if (!document.getElementById(id)) return;
    var runs = C.runs.filter(function (r) { return r.qq && r.qq.exp; });
    if (!runs.length) return;
    var band = runs.reduce(function (a, b) { return (b.qq.n > a.qq.n) ? b : a; }, runs[0]).qq;
    var mx = 0, my = 0;
    runs.forEach(function (r) {
      r.qq.exp.forEach(function (v) { if (isNum(v)) mx = Math.max(mx, v); });
      r.qq.obs.forEach(function (v) { if (isNum(v)) my = Math.max(my, v); });
    });
    var bandTop = 0;
    band.hi.forEach(function (v) { if (isNum(v)) bandTop = Math.max(bandTop, v); });
    my = Math.max(my, bandTop);
    /* zoom: the y axis up to 1.5 × the band (runs above it are clipped and named at the top) */
    var cap = 1.5 * Math.max(bandTop, mx);
    var clipped = runs.filter(function (r) { return r.qq.obs.some(function (v) { return isNum(v) && v > cap; }); });
    var seg = document.querySelector("[data-seg='qq" + C.sfx + "']");
    if (seg) seg.style.display = clipped.length ? "" : "none";
    var zoom = clipped.length > 0 && C.qqScale !== "full";
    var xr = mx * 1.04, yr = zoom ? cap : Math.max(my, mx) * 1.04;
    var notes = zoom ? clipped.map(function (r, j) {
      var top = 0;
      r.qq.obs.forEach(function (v) { if (isNum(v)) top = Math.max(top, v); });
      return { x: 0.02, y: 1, xref: "paper", yref: "paper", showarrow: false, xanchor: "left", yanchor: "top",
               yshift: -4 - 20 * j, text: r.label + " ↑ " + fmtLp(top) + " (off scale)", font: { color: r.color, size: 14 } };
    }) : [];
    var data = [
      { type: "scatter", mode: "lines", x: band.exp, y: band.lo, line: { width: 0, color: C.colors.grid }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", x: band.exp, y: band.hi, line: { width: 0, color: C.colors.grid }, fill: "tonexty",
        fillcolor: C.colors.grid, hoverinfo: "skip", showlegend: false }
    ];
    runs.forEach(function (r) {
      data.push({ type: "scatter", mode: "lines", x: r.qq.exp, y: r.qq.obs, name: r.label, hoverinfo: "none",
                  line: { color: r.color, width: 2.5, dash: r.dash } });
    });
    var p = phone();
    var lay = layoutFor(C, {
      height: 480, showlegend: true, hovermode: "closest",
      legend: { orientation: "h", traceorder: "normal", x: 0, y: p ? -0.22 : 1.02, xanchor: "left", yanchor: p ? "top" : "bottom",
                font: { size: 14, color: C.colors.ink2 } },
      margin: { l: 64, r: 16, t: p ? 12 : 40, b: p ? 110 : 56 },
      annotations: notes,
      shapes: [{ type: "line", x0: 0, y0: 0, x1: Math.min(xr, yr), y1: Math.min(xr, yr),
                 line: { color: C.colors.muted, width: 1, dash: "dot" } }],
      xaxis: { range: [0, xr], title: { text: "Expected −log10 p" }, zeroline: false },
      yaxis: { range: [0, yr], title: { text: "Observed −log10 p" }, zeroline: false }
    });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points.filter(function (x) { return x.curveNumber >= 2; })[0];
      if (!pt) return;
      var r = runs[pt.curveNumber - 2], i = pt.pointIndex;
      showTip(ev.event, r.label, [
        { color: r.color, value: fmtLp(r.qq.obs[i]), sub: "observed −log10 p", label: "" },
        { value: fmtLp(r.qq.exp[i]), sub: "expected", label: "" }
      ]);
    });
  }

  function initComparison(C) {
    var d = C.default || {};
    C.sel = 0;
    C.P = { cmp: true, sfx: C.sfx, hot: C.runs[0].hot, layout: C.layout, page: 0,
            colors: JSON.parse(JSON.stringify(C.colors)),
            tKey: d.threshold || (C.thresholds[0] || {}).key,
            k: isNum(d.k) ? d.k : (C.k[0] || 1),
            after: function () { cmpRunBars(C); } };
    C.P.colors.trans = C.runs[0].color;
    C.P.gridVar = C.runs[0].css;
    C.P.sameBin = C.sameBin !== false;
    if (!C.P.sameBin) document.querySelectorAll("[data-hot-bin='" + C.sfx + "']").forEach(function (e) { e.textContent = mb(C.runs[0].hot.binBp) + "-Mb bins"; });
    wire(C.P);
    var seg = document.querySelector("[data-seg='hot" + C.sfx + "']");
    if (seg) seg.querySelectorAll("button").forEach(function (b) {
      b.addEventListener("click", function () { cmpSelect(C, Number(b.getAttribute("data-v"))); });
    });
    var pick = document.querySelector("[data-pick='hot" + C.sfx + "']");
    if (pick) pick.addEventListener("change", function () { cmpSelect(C, Number(pick.value)); });
    var qseg = document.querySelector("[data-seg='qq" + C.sfx + "']");
    if (qseg) qseg.querySelectorAll("button").forEach(function (b) {
      b.addEventListener("click", function () {
        C.qqScale = b.getAttribute("data-v");
        qseg.querySelectorAll("button").forEach(function (x) { x.setAttribute("aria-pressed", String(x === b)); });
        drawQQcmp(C);
      });
    });
    drawHot(C.P);
    drawLambda(C);
    drawQQcmp(C);
  }

  /* ------------------------------------------------------------- wiring */
  function wire(P) {
    document.querySelectorAll("[data-hot-t='" + P.sfx + "'] button").forEach(function (b) {
      b.addEventListener("click", function () { select(P, b.getAttribute("data-v"), null); });
    });
    document.querySelectorAll("[data-hot-k='" + P.sfx + "'] button").forEach(function (b) {
      b.addEventListener("click", function () { select(P, null, Number(b.getAttribute("data-v"))); });
    });
  }
  ANQ.onReady(function () {
    if (CMP) initComparison(CMP);
    PANELS.forEach(function (P) {
      var d = (P.hot && P.hot.default) || {};
      P.tKey = d.threshold || (P.hot.thresholds[0] || {}).key;
      P.k = isNum(d.k) ? d.k : (P.hot.k[0] || 1);
      P.page = 0;
      wire(P);
      drawHot(P);
      drawQQ(P);
    });
    window.addEventListener("scroll", hideTip, { passive: true });
    var lastW = window.innerWidth, timer = null;
    window.addEventListener("resize", function () {
      if (timer) clearTimeout(timer);
      timer = setTimeout(function () {
        if (Math.abs(window.innerWidth - lastW) < 40) return;
        lastW = window.innerWidth;
        PANELS.forEach(function (P) { drawMan(P); drawTable(P); drawQQ(P); });
        if (CMP) { drawMan(CMP.P); drawTable(CMP.P); drawLambda(CMP); drawQQcmp(CMP); }
      }, 200);
    });
  });
})();
