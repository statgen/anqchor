/* anqchor cis report: power curve, signal strength, effective sample size, agreement
   charts and the "All pairs" table, drawn from ANQ_DATA.panels (one panel per dataset,
   normally one). Figures are registered with ANQ.renderFigure, which re-themes them for
   dark mode. */
(function () {
  "use strict";
  var D = window.ANQ_DATA || {};
  var PANELS = D.panels || [];
  var TICKS = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
  var SPARSE = [0, 1, 10, 100, 1000, 10000, 100000];
  var RATIO_TICKS = [0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 4, 5, 10, 20];

  /* ------------------------------------------------------------- helpers */
  function isNum(v) { return typeof v === "number" && isFinite(v); }
  function l1p(v) { return Math.log(1 + Math.max(0, v)); }
  function n0(v) { return Math.round(v).toLocaleString("en-US"); }
  function pct(v) { return (100 * v).toFixed(v >= 0.995 || v < 0.1 ? 1 : 0) + "%"; }
  function g3(v) { return isNum(v) ? String(Number(v.toPrecision(3))) : "–"; }
  function clone(o) { return JSON.parse(JSON.stringify(o)); }
  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined && text !== null) e.textContent = String(text);
    return e;
  }
  function phone() { return window.innerWidth < 640; }
  function short(s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + "…" : s; }
  /* p-value text of a -log10 p threshold; the levels and the curve count p < threshold */
  function pText(t) {
    if (t <= 0) return "any p";
    var p = Math.pow(10, -t);
    if (p >= 0.001) return "p < " + String(Number(p.toPrecision(2)));
    var e = p.toExponential(1).split("e");
    return "p < " + e[0] + "e" + Number(e[1]);
  }
  /* a log1p axis with tick labels in natural units (sparser ticks for narrow plots) */
  function log1pAxis(maxv, title, sparse) {
    var top = Math.max(1, maxv);
    var vals = (sparse ? SPARSE : TICKS).filter(function (t) { return t <= top * 1.05; });
    return {
      range: [0, l1p(top) * 1.04],
      tickvals: vals.map(l1p),
      ticktext: vals.map(String),
      title: { text: title }
    };
  }
  function layoutFor(P, extra) {
    var lay = clone(P.layout);
    Object.keys(extra).forEach(function (k) {
      if (extra[k] && typeof extra[k] === "object" && !Array.isArray(extra[k]) && lay[k] &&
          typeof lay[k] === "object") {
        var m = {};
        Object.keys(lay[k]).forEach(function (q) { m[q] = lay[k][q]; });
        Object.keys(extra[k]).forEach(function (q) {
          if (extra[k][q] && typeof extra[k][q] === "object" && m[q] && typeof m[q] === "object" &&
              !Array.isArray(m[q])) {
            var mm = {};
            Object.keys(m[q]).forEach(function (z) { mm[z] = m[q][z]; });
            Object.keys(extra[k][q]).forEach(function (z) { mm[z] = extra[k][q][z]; });
            m[q] = mm;
          } else { m[q] = extra[k][q]; }
        });
        lay[k] = m;
      } else { lay[k] = extra[k]; }
    });
    return lay;
  }
  /* threshold and reference lines: one style in every chart, clear of the gridlines */
  function vline(x, P) {
    return { type: "line", xref: "x", yref: "paper", x0: x, x1: x, y0: 0, y1: 1,
             line: { color: P.colors.ink2, width: 1.5, dash: "dot" } };
  }
  function seg(x0, y0, x1, y1, color, dash) {
    return { type: "line", x0: x0, y0: y0, x1: x1, y1: y1, line: { color: color, width: 1.5, dash: dash || "dash" } };
  }
  function figWidth(id) {
    var gd = document.getElementById(id);
    return gd && gd.clientWidth ? gd.clientWidth : 600;
  }
  /* category ticks for horizontal charts: short labels on narrow screens */
  function catAxis(keys, labels) {
    var n = phone() ? 12 : 28;
    return { type: "category", categoryorder: "array", categoryarray: keys, autorange: "reversed",
             tickvals: keys, ticktext: labels.map(function (l) { return short(l, n); }),
             showgrid: false, title: { text: "" } };
  }

  /* one tooltip element for the custom-hover charts (values first, then the run) */
  var tip = null;
  function showTip(evt, title, rows) {
    if (!tip) { tip = el("div", "cis-tip"); tip.setAttribute("role", "status"); document.body.appendChild(tip); }
    tip.innerHTML = "";
    if (title) tip.appendChild(el("div", "tt-title", title));
    rows.forEach(function (r) {
      var row = el("div", "tt-row");
      var k = el("span", "tt-key" + (r.dash && r.dash !== "solid" ? " dashed" : ""));
      k.style.setProperty("--c", ANQ.color(r.color));
      row.appendChild(k);
      row.appendChild(el("b", "tt-v", r.value));
      if (r.sub) row.appendChild(el("span", "tt-sub", r.sub));
      var l = el("span", "tt-l", r.label);
      if (r.label) l.title = r.label;
      row.appendChild(l);
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
    if (!gd || !gd.on || gd._cisHover) return;
    gd._cisHover = true;
    gd.on("plotly_hover", function (ev) { try { fn(ev); } catch (e) { console.error(e); } });
    gd.on("plotly_unhover", hideTip);
    gd.addEventListener("mouseleave", hideTip);
  }

  /* estimated rows of a horizontal Plotly legend in the element's width */
  function legendRows(id, runs) {
    var w = figWidth(id), x = 0, rows = 1, avail = Math.max(120, w - 80);
    runs.forEach(function (r) {
      var iw = 48 + 8 * String(r.label).length;
      if (x > 0 && x + iw > avail) { rows++; x = 0; }
      x += iw;
    });
    return rows;
  }

  /* ---------------------------------------------------------- power curve */
  function drawPower(P, mode) {
    var id = "fig-power" + P.sfx, C = P.curve, narrow = phone();
    var xs = C.grid.map(l1p), top = C.grid[C.grid.length - 1] || 1;
    var many = P.runs.length > 3;
    /* later runs drawn first, so that the earlier ones stay visible where lines overlap */
    var data = P.runs.slice().reverse().map(function (r) {
      var n = C.counts[r.run] || [];
      return {
        type: "scatter", mode: "lines", name: r.label, x: xs, opacity: many ? 0.85 : 1,
        y: n.map(function (v) { return mode === "f" ? (r.nTested ? v / r.nTested : null) : v; }),
        line: { color: r.color, width: many ? 1.5 : 2, shape: "hv", dash: r.dash || "solid" }, hoverinfo: "none"
      };
    });
    var tested = P.runs.map(function (r) { return r.nTested; });
    var maxN = Math.max.apply(null, tested.concat([1]));
    var same = tested.every(function (v) { return v === tested[0]; });
    var shapes = [], notes = [];
    P.cutoffs.forEach(function (c, j) {
      shapes.push(vline(l1p(c.t), P));
      notes.push({ x: l1p(c.t), y: 0, xref: "x", yref: "paper",
                   text: narrow ? c.label.replace(/^p < /, "") : c.label, showarrow: false,
                   xanchor: "left", yanchor: "bottom", xshift: 4, yshift: 4 + (narrow ? 22 * j : 0),
                   font: { color: P.colors.ink2, size: 13 } });
    });
    if (same) {
      var yT = mode === "f" ? 1 : maxN;
      shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: yT, y1: yT,
                    line: { color: P.colors.muted, width: 1, dash: "dot" } });
      notes.push({ x: 1, y: yT, xref: "paper", yref: "y", text: n0(maxN) + " tested", showarrow: false,
                   xanchor: "right", yanchor: "bottom", font: { color: P.colors.ink2, size: 13 } });
    }
    /* the legend sits above the plot: grow the figure by its rows so the plot keeps its height */
    var multi = P.runs.length > 1, rowsL = multi ? legendRows(id, P.runs) : 0;
    P.powerMode = mode;
    var lay = layoutFor(P, {
      height: 400 + 24 * rowsL, showlegend: multi, hovermode: "x", shapes: shapes, annotations: notes,
      legend: { traceorder: "reversed" },
      margin: { t: multi ? 12 + 24 * rowsL : 16 },
      xaxis: Object.assign(log1pAxis(top, "−log10 p threshold", narrow), {
        showspikes: true, spikemode: "across", spikesnap: "data", spikecolor: P.colors.axis,
        spikethickness: 1, spikedash: "solid" }),
      yaxis: mode === "f"
        ? { title: { text: "Reference pairs replicated" }, tickformat: ".0%", range: [0, 1.06] }
        : { title: { text: "Reference pairs replicated" }, range: [0, maxN * 1.06], rangemode: "tozero" }
    });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points[0];
      if (!pt) return;
      var i = pt.pointIndex, t = C.grid[i];
      var rows = P.runs.map(function (r) {
        var n = (C.counts[r.run] || [])[i];
        return { color: r.color, dash: r.dash, n: isNum(n) ? n : -1, value: isNum(n) ? n0(n) : "–",
                 sub: r.nTested && isNum(n) ? pct(n / r.nTested) : "", label: r.label };
      });
      rows.sort(function (a, b) { return b.n - a.n; });
      showTip(ev.event, pText(t) + "  (−log10 p > " + g3(t) + ")", rows);
    });
  }

  /* ------------------------------------------------------- signal strength */
  function drawDist(P, mode) {
    var id = "fig-dist" + P.sfx;
    var rows = P.runs.map(function (r, i) {
      return { key: "r" + i, label: r.label, color: r.color, s: P.dist[mode][i] };
    });
    if (P.ref[mode]) rows.push({ key: "ref", label: "Reference study", color: P.colors.ink2, s: P.ref[mode], ref: true });
    rows = rows.filter(function (r) { return r.s; });
    var tf = mode === "lp" ? l1p : function (v) { return v; };
    var data = rows.map(function (r) {
      var s = r.s;
      return {
        type: "box", orientation: "h", name: r.label, y: [r.key],
        q1: [tf(s.q25)], median: [tf(s.med)], q3: [tf(s.q75)],
        lowerfence: [tf(s.q05)], upperfence: [tf(s.q95)],
        fillcolor: r.ref ? "rgba(0,0,0,0)" : r.color, line: { color: r.color, width: 2 },
        opacity: r.ref ? 1 : 0.9, whiskerwidth: 0.4, width: 0.56, hoverinfo: "none", showlegend: false
      };
    });
    /* the runs set the range; a larger reference study is clipped at the right edge */
    var runRows = rows.filter(function (r) { return !r.ref; });
    var maxv = Math.max.apply(null, runRows.map(function (r) { return r.s.q95; }).concat([1]));
    var xaxis = mode === "lp" ? log1pAxis(maxv, "−log10 p", phone()) :
      { title: { text: "|z|" }, range: [0, maxv * 1.05], rangemode: "tozero" };
    var keys = rows.map(function (r) { return r.key; });
    var lay = layoutFor(P, {
      height: 120 + 52 * rows.length, showlegend: false, hovermode: "closest",
      margin: { l: 16, t: 8 }, xaxis: xaxis,
      yaxis: catAxis(keys, rows.map(function (r) { return r.label; }))
    });
    /* medians: a surface-colored tick inside each filled box (ink in the hollow one) */
    lay.shapes = rows.map(function (r) {
      return { type: "line", xref: "x", yref: "y", x0: tf(r.s.med), x1: tf(r.s.med),
               y0: keys.indexOf(r.key) - 0.28, y1: keys.indexOf(r.key) + 0.28,
               line: { color: r.ref ? P.colors.ink2 : P.colors.surface, width: 3 } };
    });
    if (mode === "lp") P.cutoffs.forEach(function (c) { lay.shapes.push(vline(l1p(c.t), P)); });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points[0];
      if (!pt) return;
      var r = rows[pt.curveNumber];
      if (!r) return;
      var s = r.s, u = mode === "lp" ? "−log10 p" : "|z|";
      var lines = [
        { color: r.color, value: g3(s.med), sub: "median " + u, label: "" },
        { color: r.color, value: g3(s.q25) + "–" + g3(s.q75), sub: "middle half", label: "" },
        { color: r.color, value: g3(s.q05) + "–" + g3(s.q95), sub: "5th–95th", label: "" }
      ];
      if (mode === "z") lines.push({ color: r.color, value: g3(s.msq), sub: "mean z²", label: "" });
      showTip(ev.event, r.label + " · " + n0(s.n) + " pairs", lines);
    });
  }

  /* effective sample size: dots and 95% intervals on a log axis centred on 1 */
  function drawNeff(P) {
    var id = "fig-neff" + P.sfx;
    if (!document.getElementById(id)) return;
    var runs = P.runs, keys = runs.map(function (r, i) { return "r" + i; });
    var data = runs.map(function (r, i) {
      var v = P.neff[i] || {}, ok = isNum(v.v) && v.v > 0;
      var err = ok && isNum(v.lo) && isNum(v.hi)
        ? { type: "data", symmetric: false, array: [v.hi - v.v], arrayminus: [v.v - v.lo],
            color: r.color, thickness: 2.5, width: 0 }
        : { visible: false };
      return {
        type: "scatter", mode: "markers", name: r.label, y: [keys[i]], x: [ok ? v.v : null],
        marker: i === 0 ? { color: P.colors.surface, size: 13, symbol: "circle", line: { color: r.color, width: 2.5 } }
                        : { color: r.color, size: 14, line: { color: P.colors.surface, width: 1.5 } },
        error_x: err, hoverinfo: "none", showlegend: false
      };
    });
    var lg = [Math.log10(2)];
    P.neff.forEach(function (v) {
      [v.v, v.lo, v.hi].forEach(function (x) { if (isNum(x) && x > 0) lg.push(Math.abs(Math.log10(x))); });
    });
    var m = Math.max.apply(null, lg) * 1.12;
    var ticks = RATIO_TICKS.filter(function (t) { return Math.abs(Math.log10(t)) <= m; });
    if (ticks.indexOf(4) >= 0 && ticks.indexOf(5) >= 0) ticks.splice(ticks.indexOf(5), 1);
    if (ticks.indexOf(0.25) >= 0 && ticks.indexOf(0.2) >= 0) ticks.splice(ticks.indexOf(0.2), 1);
    var lay = layoutFor(P, {
      height: 72 + 48 * runs.length, hovermode: "closest", showlegend: false,
      margin: { l: 16, t: 8, r: 76, b: 40 },
      xaxis: { type: "log", range: [-m, m], tickvals: ticks,
               ticktext: ticks.map(function (t) { return "×" + t; }), title: { text: "" }, zeroline: false },
      yaxis: catAxis(keys, runs.map(function (r) { return r.label; })),
      shapes: [vline(1, P)],
      /* values in a column right of the plot, clear of the dots and intervals */
      annotations: runs.map(function (r, i) {
        var v = P.neff[i] || {};
        return { x: 1, y: keys[i], xref: "paper", yref: "y", text: isNum(v.v) ? "×" + v.v.toFixed(2) : "NA",
                 showarrow: false, xanchor: "left", xshift: 10, font: { color: P.colors.ink, size: 15 } };
      })
    });
    ANQ.renderFigure(id, { data: data, layout: lay });
    onHover(id, function (ev) {
      var pt = ev.points && ev.points[0];
      if (!pt) return;
      var r = runs[pt.curveNumber], v = P.neff[pt.curveNumber] || {};
      var rows = [{ color: r.color, value: isNum(v.v) ? "×" + v.v.toFixed(2) : "NA", sub: "", label: r.label }];
      if (isNum(v.lo)) rows.push({ color: r.color, value: v.lo.toFixed(2) + "–" + v.hi.toFixed(2), sub: "95% interval", label: "" });
      showTip(ev.event, pt.curveNumber === 0 ? "Baseline" : "vs " + runs[0].label, rows);
    });
  }

  /* --------------------------------------------------- agreement scatters */
  function pearson(x, y) {
    var n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    for (var i = 0; i < x.length; i++) {
      var a = x[i], b = y[i];
      if (!isNum(a) || !isNum(b)) continue;
      n++; sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
    }
    if (n < 3) return null;
    var cov = sxy - sx * sy / n, vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
    return vx > 0 && vy > 0 ? cov / Math.sqrt(vx * vy) : null;
  }
  function slope0(x, y) {
    var sxy = 0, sxx = 0;
    for (var i = 0; i < x.length; i++) {
      if (!isNum(x[i]) || !isNum(y[i])) continue;
      sxy += x[i] * y[i]; sxx += x[i] * x[i];
    }
    return sxx > 0 ? sxy / sxx : null;
  }
  function sgn(v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }
  function scatterFig(P, key, run, narrow) {
    var pr = P.pairs, v = P.vals[run.run], base = P.vals[P.baseline];
    var xs, ys, xt, yt, logScale = false;
    if (key === "beta") { xs = pr.beta; ys = v.beta; xt = "Reference effect"; yt = "Effect, " + run.label; }
    else if (key === "z") { xs = pr.z; ys = v.z; xt = "Reference z"; yt = "z, " + run.label; }
    else if (key === "lp") { xs = pr.lp; ys = v.lp; xt = "Reference −log10 p"; yt = "−log10 p, " + run.label; logScale = true; }
    else {
      xs = base.z; ys = v.z;
      xt = "z, " + P.runs[0].label; yt = "z, " + run.label;
    }
    var tf = logScale ? l1p : function (a) { return a; };
    var groups = [[], []];  /* 0: same sign, 1: opposite */
    for (var i = 0; i < xs.length; i++) {
      if (!isNum(xs[i]) || !isNum(ys[i])) continue;
      /* against the baseline: the two runs disagree in sign; otherwise: opposite to the reference */
      var opp = key === "zb" ? sgn(xs[i]) * sgn(ys[i]) < 0 : v.sign[i] === 0;
      groups[opp ? 1 : 0].push(i);
    }
    var many = xs.length > 4000;
    var data = groups.map(function (idx, g) {
      return {
        type: many ? "scattergl" : "scatter", mode: "markers",
        name: g ? "opposite sign" : "same sign",
        x: idx.map(function (i) { return tf(xs[i]); }), y: idx.map(function (i) { return tf(ys[i]); }),
        customdata: idx.map(function (i) { return [pr.trait[i], pr.variant[i], g3(xs[i]), g3(ys[i])]; }),
        marker: g ? { color: P.colors.ink, symbol: "diamond-open", size: 11, line: { color: P.colors.ink, width: 2 } }
                  : { color: run.color, size: many ? 5 : 8, opacity: 0.8, line: { color: P.colors.surface, width: many ? 0 : 1 } },
        hovertemplate: "<b>%{customdata[0]}</b><br>%{customdata[1]}<br>" + xt + " %{customdata[2]}<br>" +
                       yt + " %{customdata[3]}<extra></extra>",
        showlegend: false
      };
    });
    var idxAll = groups[0].concat(groups[1]);
    function absMax(arr) {
      var m = 0;
      idxAll.forEach(function (i) { m = Math.max(m, Math.abs(arr[i])); });
      return m;
    }
    var xl = Math.max(absMax(xs), logScale ? 1 : 1e-9) * 1.05, yl = Math.max(absMax(ys), logScale ? 1 : 1e-9) * 1.05;
    if (key === "zb") { xl = yl = Math.max(xl, yl); }  /* same scale: y = x is the reference line */
    var axX, axY, shapes = [];
    var zl = { zeroline: true, zerolinecolor: P.colors.axis, zerolinewidth: 1 };
    if (logScale) {
      axX = log1pAxis(xl, "", narrow); axY = log1pAxis(yl, "", narrow);
      /* a little room below 0, so that markers at p ~ 1 are not cut */
      [axX, axY].forEach(function (a) { a.range[0] = -0.03 * a.range[1]; });
    } else {
      axX = Object.assign({ range: [-xl, xl] }, zl); axY = Object.assign({ range: [-yl, yl] }, zl);
    }
    var b = null;
    if (key === "zb") shapes.push(seg(-xl, -xl, xl, xl, P.colors.ink2));
    else if (!logScale) {
      b = slope0(xs, ys);
      if (isNum(b)) {
        /* the fit through 0, clipped to the plot */
        var x1 = Math.min(xl, Math.abs(b) > 0 ? yl / Math.abs(b) : xl);
        shapes.push(seg(-x1, -b * x1, x1, b * x1, P.colors.ink2));
      }
    }
    var lay = layoutFor(P, {
      height: 360, showlegend: false, hovermode: "closest", shapes: shapes,
      margin: { l: 60, r: 12, t: 8, b: 52 },
      xaxis: Object.assign(clone(axX), { title: { text: xt } }),
      yaxis: Object.assign(clone(axY), { title: { text: yt } })
    });
    var r = pearson(xs, ys), nOpp = groups[1].length, n = groups[0].length + nOpp;
    var stat = "r " + (r === null ? "–" : r.toFixed(3));
    if (isNum(b)) stat += " · slope " + g3(b);
    if (key === "beta" && n) stat += " · opposite sign " + n0(nOpp);
    if (key === "zb") {
      var k = P.runs.map(function (x) { return x.run; }).indexOf(run.run), ne = P.neff[k] || {};
      if (isNum(ne.v)) stat += " · effective size ×" + ne.v.toFixed(2);
      if (nOpp) stat += " · sign differs " + n0(nOpp);
    }
    return { spec: { data: data, layout: lay }, stat: stat };
  }
  function drawConc(P, runKey) {
    var run = P.runs.filter(function (r) { return r.run === runKey; })[0] || P.runs[0];
    if (!run) return;
    P.selected = run.run;
    var root = document.getElementById("agreement" + P.sfx);
    if (!root) return;
    var grid = root.querySelector(".conc-grid");
    var hideZb = run.run === P.baseline;
    /* three plots in one row when "z vs baseline" is hidden (the baseline is selected) */
    var relayout = false;
    if (P.comparison && grid && grid.classList.contains("n3") !== hideZb) {
      grid.classList.toggle("n3", hideZb);
      relayout = true;
    }
    root.querySelectorAll(".conc-card").forEach(function (card) {
      var k = card.getAttribute("data-k");
      var hide = k === "zb" && hideZb;
      card.hidden = hide;
      if (hide) return;
      var gd = document.getElementById("fig-conc-" + k + P.sfx);
      var f = scatterFig(P, k, run, figWidth("fig-conc-" + k + P.sfx) < 380 || phone());
      card.querySelector(".cc-stat").textContent = f.stat;
      if (gd) ANQ.renderFigure(gd.id, f.spec);
      /* the cards changed width with the grid: fit the plots to them */
      if (gd && relayout && window.Plotly && window.Plotly.Plots) window.Plotly.Plots.resize(gd);
    });
  }

  /* ------------------------------------------------------ all pairs table */
  function tableCols(P) {
    var cols = [{ k: "pair", label: "Pair", w: 300, pw: 150, sticky: true },
                { k: "ref", label: "Reference", plabel: "Ref.", num: true, w: 112, pw: 64,
                  title: "Reference study −log10 p" }];
    if (P.comparison) cols.push({ k: "opp", label: "Opp.", num: true, w: 72, pw: 56,
                                  title: "Runs with the opposite sign to the reference" });
    P.runs.forEach(function (r, j) {
      cols.push({ k: "r" + j, run: r, label: P.comparison ? r.label : "−log10 p", num: true,
                  w: 116, pw: 84, title: (P.comparison ? r.label + ": " : "") + "−log10 p" });
    });
    return cols;
  }
  function tableRows(P) {
    var T = P.table, n = T.trait.length, out = [];
    for (var i = 0; i < n; i++) {
      var opp = 0;
      P.runs.forEach(function (r) { if (T.sign[r.run][i] === 0) opp++; });
      out.push({ i: i, opp: opp });
    }
    return out;
  }
  function cellVal(P, row, c) {
    var T = P.table, i = row.i;
    if (c.k === "pair") return T.trait[i];
    if (c.k === "ref") return T.ref[i];
    if (c.k === "opp") return row.opp;
    return T.lp[c.run.run][i];
  }
  function tView(P) {
    var S = P.tstate, q = S.q, T = P.table;
    var rows = P.trows.filter(function (r) {
      return !q || (T.trait[r.i] + "\u0000" + T.variant[r.i]).toLowerCase().indexOf(q) >= 0;
    });
    var col = P.tcols.filter(function (c) { return c.k === S.col; })[0] || P.tcols[1], dir = S.dir;
    rows.sort(function (a, b) {
      var x = cellVal(P, a, col), y = cellVal(P, b, col);
      var xm = !isNum(x) && typeof x !== "string", ym = !isNum(y) && typeof y !== "string";
      if (xm || ym) return xm === ym ? a.i - b.i : (xm ? 1 : -1);
      if (typeof x === "number" && typeof y === "number") return (x - y) * dir || a.i - b.i;
      return String(x).localeCompare(String(y), undefined, { numeric: true }) * dir;
    });
    return rows;
  }
  function tPager(P) {
    var box = el("span", "tpager");
    var prev = el("button", "", "‹"), info = el("span", "num"), next = el("button", "", "›");
    prev.type = next.type = "button";
    prev.setAttribute("aria-label", "Previous page"); next.setAttribute("aria-label", "Next page");
    prev.addEventListener("click", function () { P.tstate.page = Math.max(0, P.tstate.page - 1); drawTable(P); });
    next.addEventListener("click", function () { P.tstate.page += 1; drawTable(P); });
    box.appendChild(prev); box.appendChild(info); box.appendChild(next);
    return { box: box, prev: prev, next: next, info: info };
  }
  function tToolbar(P) {
    var bar = document.querySelector("[data-ptbl='" + P.sfx + "']");
    var foot = document.querySelector("[data-ptfoot='" + P.sfx + "']");
    if (!bar) return false;
    var row = el("div", "tbar-row");
    var inp = el("input");
    inp.type = "search"; inp.placeholder = "Search trait or variant…";
    inp.setAttribute("aria-label", "Search trait or variant");
    inp.addEventListener("input", function () {
      P.tstate.q = inp.value.trim().toLowerCase(); P.tstate.page = 0; drawTable(P);
    });
    row.appendChild(inp);
    var top = tPager(P);
    row.appendChild(top.box);
    var size = el("span", "tsize");
    size.setAttribute("role", "group"); size.setAttribute("aria-label", "Rows per page");
    size.appendChild(document.createTextNode("Rows"));
    [[100, "100"], [500, "500"], [0, "All"]].forEach(function (o) {
      var b = el("button", "", o[1]);
      b.type = "button"; b._size = o[0];
      b.addEventListener("click", function () { P.tstate.size = o[0]; P.tstate.page = 0; drawTable(P); });
      size.appendChild(b);
    });
    row.appendChild(size);
    var dl = el("button", "btn", "Download TSV");
    dl.type = "button";
    dl.addEventListener("click", function () { tDownload(P); });
    row.appendChild(dl);
    bar.appendChild(row);
    var bottom = tPager(P);
    if (foot) foot.appendChild(bottom.box);
    P.tui = { pagers: [top, bottom], size: size };
    return true;
  }
  function fmtLp(v) { return v >= 100 ? v.toFixed(0) : v.toFixed(1); }
  function drawTable(P) {
    if (!P.tui && !tToolbar(P)) return;
    var wrap = document.getElementById("ptbl" + P.sfx);
    var table = wrap.querySelector("table"), S = P.tstate, T = P.table, p = phone();
    var cols = P.tcols, total = 0;
    var cg = table.querySelector("colgroup");
    cg.innerHTML = "";
    cols.forEach(function (c) {
      var col = el("col"), w = p ? c.pw : c.w;
      col.style.width = w + "px"; total += w;
      cg.appendChild(col);
    });
    table.style.width = total + "px";
    var thead = table.querySelector("thead");
    thead.innerHTML = "";
    var tr = el("tr");
    cols.forEach(function (c) {
      var th = el("th", (c.num ? "num" : "") + (c.sticky ? " pin" : ""));
      th.scope = "col"; th.tabIndex = 0;
      th.title = c.title || c.label;
      if (c.run) {
        th.style.setProperty("--c", ANQ.color(c.run.color));
        if (P.comparison) th.classList.add("runh");
      }
      var lab = el("span", "th-l", p && c.plabel ? c.plabel : c.label);
      lab.appendChild(el("span", "sorti"));
      th.appendChild(lab);
      th.setAttribute("aria-sort", S.col === c.k ? (S.dir > 0 ? "ascending" : "descending") : "none");
      var act = function () {
        if (S.col === c.k) S.dir = -S.dir; else { S.col = c.k; S.dir = c.num ? -1 : 1; }
        S.page = 0; drawTable(P);
      };
      th.addEventListener("click", act);
      th.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); } });
      tr.appendChild(th);
    });
    thead.appendChild(tr);
    var rows = tView(P), n = rows.length;
    var ps = S.size > 0 ? S.size : Math.max(1, n), pages = Math.max(1, Math.ceil(n / ps));
    if (S.page >= pages) S.page = pages - 1;
    var start = S.page * ps, end = Math.min(n, start + ps);
    P.tui.pagers.forEach(function (pg, j) {
      pg.info.textContent = n ? n0(start + 1) + "–" + n0(end) + " of " + n0(n) : "No rows";
      pg.prev.disabled = S.page === 0;
      pg.next.disabled = S.page >= pages - 1;
      if (j === 1) pg.box.hidden = pages <= 1;
    });
    P.tui.size.querySelectorAll("button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b._size === S.size));
    });
    var tbody = table.querySelector("tbody");
    tbody.innerHTML = "";
    var frag = document.createDocumentFragment();
    for (var k = start; k < end; k++) {
      var row = rows[k], i = row.i, trr = el("tr");
      cols.forEach(function (c) {
        var td;
        if (c.k === "pair") {
          td = el("td", "pin");
          var two = el("span", "two");
          two.appendChild(el("span", "", T.trait[i]));
          two.appendChild(el("span", "sub idc", T.variant[i]));
          td.appendChild(two);
          td.title = T.trait[i] + "\n" + T.variant[i];
        } else if (c.k === "ref") {
          var v = T.ref[i];
          td = el("td", "num" + (isNum(v) ? "" : " muted"), isNum(v) ? fmtLp(v) : "–");
        } else if (c.k === "opp") {
          td = el("td", "num" + (row.opp ? " bad" : " muted"), row.opp ? row.opp : "0");
        } else {
          var lp = T.lp[c.run.run][i], sg = T.sign[c.run.run][i];
          if (!isNum(lp)) {
            td = el("td", "num muted", "–");
            td.title = T.why[c.run.run][i] || "not tested";
          } else if (sg === 0) {
            td = el("td", "num bad", "−" + fmtLp(lp));
            td.title = "opposite sign to the reference";
          } else td = el("td", "num", fmtLp(lp));
        }
        trr.appendChild(td);
      });
      frag.appendChild(trr);
    }
    tbody.appendChild(frag);
  }
  function tDownload(P) {
    var T = P.table;
    var head = ["trait", "variant", "ref_log10p", "n_opposite_sign"];
    P.runs.forEach(function (r) { head.push("log10p:" + r.label, "sign:" + r.label); });
    var lines = [head.join("\t")];
    function v(x) { return isNum(x) ? String(x) : "NA"; }
    tView(P).forEach(function (row) {
      var i = row.i, f = [T.trait[i], T.variant[i], v(T.ref[i]), String(row.opp)];
      P.runs.forEach(function (r) {
        var s = T.sign[r.run][i];
        f.push(v(T.lp[r.run][i]), s === 1 ? "same" : (s === 0 ? "opposite" : "NA"));
      });
      lines.push(f.join("\t"));
    });
    var blob = new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = P.file || "cis_pairs.tsv";
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }

  /* ------------------------------------------------------------- wiring */
  function wireSeg(name, fn) {
    var seg = document.querySelector("[data-seg='" + name + "']");
    if (!seg) return;
    seg.querySelectorAll("button").forEach(function (b) {
      b.addEventListener("click", function () {
        seg.querySelectorAll("button").forEach(function (o) { o.setAttribute("aria-pressed", String(o === b)); });
        fn(b.getAttribute("data-v"));
      });
    });
  }
  function wirePick(name, fn) {
    var s = document.querySelector("[data-pick='" + name + "']");
    if (s) s.addEventListener("change", function () { fn(s.value); });
  }
  ANQ.onReady(function () {
    PANELS.forEach(function (P) {
      P.colors = P.colors || {};
      drawPower(P, "n");
      drawDist(P, "z");
      drawNeff(P);
      drawConc(P, P.selected);
      wireSeg("power" + P.sfx, function (v) { drawPower(P, v); });
      wireSeg("dist" + P.sfx, function (v) { P.distMode = v; drawDist(P, v); });
      wireSeg("conc" + P.sfx, function (v) { drawConc(P, v); });
      wirePick("conc" + P.sfx, function (v) { drawConc(P, v); });
      if (P.table) {
        P.tcols = tableCols(P);
        P.trows = tableRows(P);
        P.tstate = { q: "", col: "ref", dir: -1, page: 0, size: 100 };
        drawTable(P);
      }
    });
    window.addEventListener("scroll", hideTip, { passive: true });
    /* legends, short labels and tick density follow the width: redraw when it changes */
    var lastW = window.innerWidth, timer = null;
    window.addEventListener("resize", function () {
      if (timer) clearTimeout(timer);
      timer = setTimeout(function () {
        if (Math.abs(window.innerWidth - lastW) < 40) return;
        lastW = window.innerWidth;
        PANELS.forEach(function (P) {
          drawPower(P, P.powerMode || "n");
          drawDist(P, P.distMode || "z");
          drawNeff(P);
          drawConc(P, P.selected);
          if (P.table) drawTable(P);
        });
      }, 200);
    });
  });
  ANQ.CIS = { pText: pText, log1pAxis: log1pAxis, phone: phone };
})();
