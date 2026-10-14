/* anqchor match report: dashboard, threshold sliders, samples to check, match chart,
   duplicates and the sample table. Statuses are recomputed in the browser with ANQ_MATCH
   (a port of qpgen's rules); the files in the run directory keep the run's.
   "Match score" on the page is qpgen's Z-score (z_self, z_1..z_5). */
(function () {
  "use strict";
  var D = window.ANQ_DATA || {};
  var M = window.ANQ_MATCH || ANQ_MATCH;
  var C = D.colors;
  var RUN = D.config;
  var st = { zt: RUN.z_threshold, zd: RUN.z_diff };
  var CATS = D.categories, CAT = {};
  CATS.forEach(function (c, i) { c.order = i; CAT[c.key] = c; });
  var MAX_CARDS = 12, MAX_LABELS = 8, MAX_UNASSIGNED_CARDS = 6, PHONE = 700;
  /* chart drawing order: common results first, so the rare ones sit on top */
  var CHART_ORDER = ["strict", "noid", "newid", "lenient", "nomatch", "several", "taken", "twoids", "different"];
  var SYMBOL = { strict: "circle", lenient: "diamond", different: "x", twoids: "star", several: "square",
                 nomatch: "triangle-up", newid: "circle", taken: "x", noid: "circle" };

  /* ------------------------------------------------------------- helpers */
  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined && text !== null) e.textContent = String(text);
    return e;
  }
  function isNum(v) { return typeof v === "number" && isFinite(v); }
  function f2(v) { return isNum(v) ? v.toFixed(2) : "–"; }
  function n0(v) { return v.toLocaleString("en-US"); }
  function rankText(n) { return "rank " + n0(n); }
  function code(id) { return el("code", "", id); }
  function icon(c) {
    var i = el("span", "ic", c.icon);
    i.setAttribute("aria-hidden", "true");
    return i;
  }
  function chip(key, count, title) {
    var c = CAT[key], s = el("span", "chip-status r-" + c.role);
    s.appendChild(icon(c));
    s.appendChild(el("span", "cl", c.label + (count !== undefined ? " " + n0(count) : "")));
    if (title) s.title = title;
    return s;
  }
  function plainChip(role, ic, label) {
    var s = el("span", "chip-status r-" + role);
    s.appendChild(icon({ icon: ic }));
    s.appendChild(document.createTextNode(label));
    return s;
  }
  function phone() { return window.innerWidth < PHONE; }
  function changed() { return Math.abs(st.zt - RUN.z_threshold) > 1e-9 || Math.abs(st.zd - RUN.z_diff) > 1e-9; }
  function otherOf(r) {
    /* best genotype other than the assigned one (the best genotype when none is assigned) */
    if (r.geno_id_assigned && r.rank_self === 1) return { id: r.id_2, z: r.z_2 };
    return { id: r.id_1, z: r.z_1 };
  }
  function ownAhead(r) {
    var o = otherOf(r);
    return r.geno_id_assigned && isNum(r.z_self) && isNum(o.z) ? r.z_self - o.z : null;
  }
  function top5(r) {
    var out = [];
    for (var k = 1; k <= 5; k++) if (r["id_" + k]) out.push(r["id_" + k] + " " + f2(r["z_" + k]));
    return out.join(", ");
  }

  /* ------------------------------------------------------------- datasets */
  D.datasets.forEach(function (ds) {
    ds.objs = ds.rows.map(function (a) {
      var o = {};
      D.fields.forEach(function (f, i) { o[f] = a[i]; });
      return o;
    });
    ds.mapped = ds.objs.map(function (r) { return !!r.geno_id_assigned; });
    ds.nUnassigned = ds.mapped.filter(function (m) { return !m; }).length;
    ds.claimed = {};
    ds.objs.forEach(function (r) { if (r.geno_id_assigned) ds.claimed[r.geno_id_assigned] = true; });
    ds.runEvals = M.evaluateAll(ds.objs, RUN.z_threshold, RUN.z_diff);
    ds.runCats = ds.runEvals.map(function (e) { return e.cat; });
    ds.runTiles = M.tileCounts(ds.runCats, ds.mapped, ds.runEvals.map(function (e) { return e.dup; }));
    ds.panel = document.getElementById("panel-" + ds.slug);
    ds.hidden = {};
    ds.filter = null;
    ds.query = "";
    ds.sort = { col: "result", dir: 1 };
    ds.page = 0;
    ds.pageSize = window.innerWidth < 600 ? 100 : 500;
    ds.sig = "";
  });
  var active = D.datasets[0];

  /* ------------------------------------------------------------ dashboard */
  var TILE_ROLE = { strict: "good", lenient: "warning", different: "critical", twoids: "dup", several: "serious",
                    none: "serious", nodata: "muted", unassigned: "info" };
  var TILE_OPTIONAL = { twoids: true, several: true, nodata: true };
  var TILE_VALUE = { twoids: "dups" };
  var TILE_FILTER = {
    examined: null,
    assigned: null,
    strict: { key: "strict", test: function (ds, i) { return ds.cats[i] === "strict"; } },
    lenient: { key: "lenient", test: function (ds, i) { return ds.mapped[i] && (ds.cats[i] === "strict" || ds.cats[i] === "lenient"); } },
    different: { key: "different", test: function (ds, i) { return ds.cats[i] === "different"; } },
    twoids: { key: "twoids", test: function (ds, i) { return ds.evals[i].dup; } },
    several: { key: "several", test: function (ds, i) { return ds.mapped[i] && ds.cats[i] === "several"; } },
    none: { key: "none", test: function (ds, i) { return ds.mapped[i] && ds.cats[i] === "nomatch"; } },
    nodata: { key: "nodata", test: function (ds, i) { return ds.mapped[i] && ds.cats[i] === "nodata"; } },
    unassigned: { key: "unassigned", test: function (ds, i) { return !ds.mapped[i]; } }
  };
  function pct(n, d) { return d ? (100 * n / d).toFixed(1) + "%" : ""; }
  /* Mirrors unassigned_sub in anqchor.report.match. */
  function unassignedSub(t) {
    var p = [];
    [["newid", "new ID found"], ["taken", "ID already used"], ["un_twoids", "two genotype IDs"],
     ["unclear", "no clear match"], ["un_nodata", "not scored"]].forEach(function (k) {
      if (t[k[0]]) p.push(n0(t[k[0]]) + " " + k[1]);
    });
    return p.join(" · ");
  }
  function shortList(ids) {
    return ids.length <= 3 ? ids.join(", ") : ids[0] + ", " + ids[1] + " and " + (ids.length - 2) + " more";
  }
  /* Mirrors tile_subs in anqchor.report.match. */
  function drawTiles(ds) {
    var t = M.tileCounts(ds.cats, ds.mapped, ds.evals.map(function (e) { return e.dup; })), r = ds.runTiles;
    var off = changed();
    var diffIds = [], dupIds = [], weakIds = [];
    ds.objs.forEach(function (o, i) {
      if (ds.cats[i] === "different") diffIds.push(o.pheno_id + " → " + o.id_1);
      if (ds.evals[i].dup) dupIds.push(o.pheno_id);
      if (ds.cats[i] === "nomatch" && ds.mapped[i]) weakIds.push(o.pheno_id);
    });
    var subs = {
      examined: "",
      assigned: t.examined > t.assigned ? n0(t.examined - t.assigned) + " without" : "",
      strict: pct(t.strict, t.assigned),
      lenient: t.lenient_only ? "strict + " + n0(t.lenient_only) : "",
      different: t.different ? shortList(diffIds) : "",
      twoids: t.dups ? shortList(dupIds) : "",
      several: "",
      none: t.none ? shortList(weakIds) : "",
      nodata: "",
      unassigned: unassignedSub(t)
    };
    var dash = ds.panel.querySelector(".dash");
    dash.querySelectorAll(".dash-tile").forEach(function (b) {
      var k = b.getAttribute("data-tile"), vk = TILE_VALUE[k] || k, v = t[vk];
      if (v === undefined) return;
      if (TILE_OPTIONAL[k]) b.hidden = !(v || r[vk]);
      b.querySelector(".t-value").textContent = n0(v);
      var sub = b.querySelector(".t-sub");
      sub.textContent = subs[k];
      sub.title = subs[k];
      var role = TILE_ROLE[k] && (v || k === "unassigned") ? TILE_ROLE[k] : "";
      b.className = "dash-tile" + (role ? " has-role r-" + role : "");
      var delta = off && r[vk] !== v, runEl = b.querySelector(".t-run");
      runEl.hidden = !delta;
      runEl.textContent = delta ? "run: " + n0(r[vk]) : "";
      if (delta) b.classList.add("changed");
      b.setAttribute("aria-pressed", String(!!(ds.filter && ds.filter.key === k)));
    });
  }

  /* ----------------------------------------------------- samples to check */
  function sevSort(ds) {
    return function (a, b) {
      var d = CAT[ds.cats[a]].order - CAT[ds.cats[b]].order;
      if (d) return d;
      var x = ds.objs[a].z_self, y = ds.objs[b].z_self;
      if (!isNum(x) && !isNum(y)) { x = ds.objs[a].z_1; y = ds.objs[b].z_1; }
      if (!isNum(x)) return isNum(y) ? 1 : 0;
      if (!isNum(y)) return -1;
      return x - y;
    };
  }
  /* every sample that needs a look: assigned and not a strict match, matching two
     genotype IDs, or without an assigned ID */
  function toCheck(ds, i) { return (ds.mapped[i] && ds.cats[i] !== "strict") || ds.evals[i].dup || !ds.mapped[i]; }
  function barRow(grid, label, z, max, own, extra, cut) {
    var l = el("span", "bl");
    if (typeof label === "string") l.textContent = label; else l.appendChild(label);
    var tr = el("span", "track"), fill = el("span", "fill" + (own ? " own" : ""));
    fill.style.width = (isNum(z) ? Math.max(0, Math.min(1, z / max)) * 100 : 0).toFixed(1) + "%";
    tr.appendChild(fill);
    if (cut !== undefined) {
      var c = el("span", "cut");
      c.style.left = (Math.max(0, Math.min(1, cut / max)) * 100).toFixed(1) + "%";
      c.title = "Lenient cutoff " + cut.toFixed(2);
      tr.appendChild(c);
    }
    grid.appendChild(l); grid.appendChild(tr);
    grid.appendChild(el("span", "bv", f2(z)));
    grid.appendChild(el("span", "bx", extra || ""));
  }
  function strictOwner(ds, geno, cats) {
    cats = cats || ds.cats;
    for (var j = 0; j < ds.objs.length; j++) {
      if (ds.objs[j].geno_id_assigned === geno && cats[j] === "strict") return j;
    }
    return -1;
  }
  function idLabel(r, id) {
    var c = code(id || "–");
    if (id && id === r.geno_id_assigned) c.className = "own-id";
    return c;
  }
  /* "Best: X, ahead of next by 1.97 (needs more than 2.00)." */
  function leadNote(r) {
    if (!r.id_1 || !isNum(r.z_1) || !isNum(r.z_2)) return "";
    return "Best: " + r.id_1 + ", ahead of next by " + (r.z_1 - r.z_2).toFixed(2) +
      " (needs more than " + st.zd.toFixed(2) + ").";
  }
  function noData(r) { return r.qpgen_status === D.notScored ? ": not scored (no covariates)" : ": no data"; }
  function card(ds, i, max) {
    var r = ds.objs[i], c = ds.cats[i], e = ds.evals[i];
    var a = el("article", "check-card r-" + CAT[c].role);
    a.appendChild(el("div")).appendChild(chip(c, undefined, e.status));
    var title = el("div", "cc-title"), bars = el("div", "bars"), note = "";
    function arrowTitle(to) {
      title.appendChild(code(r.pheno_id));
      title.appendChild(document.createTextNode(" → "));
      title.appendChild(code(to || "–"));
    }
    if (e.dup) {
      title.appendChild(code(r.pheno_id));
      title.appendChild(document.createTextNode(": two genotype IDs"));
      barRow(bars, idLabel(r, r.id_1), r.z_1, max, r.id_1 === r.geno_id_assigned, "");
      barRow(bars, idLabel(r, r.id_2), r.z_2, max, r.id_2 === r.geno_id_assigned, "");
      barRow(bars, "next", r.z_3, max, false, "");
      note = r.id_1 + " and " + r.id_2 + " are probably one person.";
      if (ds.mapped[i] && !e.twoIds) note += " Own ID: " + rankText(r.rank_self) + ".";
    } else if (c === "different") {
      arrowTitle(r.id_1);
      barRow(bars, "own ID", r.z_self, max, true, isNum(r.rank_self) ? rankText(r.rank_self) : "", st.zt);
      barRow(bars, code(r.id_1), r.z_1, max, false, "");
      var j = r.id_1 ? strictOwner(ds, r.id_1) : -1;
      if (j >= 0) note = r.id_1 + " also matches its own sample" + (ds.objs[j].pheno_id !== r.id_1 ? " " + ds.objs[j].pheno_id : "") + ".";
    } else if (c === "newid" || c === "taken") {
      arrowTitle(r.id_1);
      barRow(bars, code(r.id_1), r.z_1, max, false, "best");
      if (r.id_2) barRow(bars, code(r.id_2), r.z_2, max, false, "next");
      if (c === "taken") {
        var k = strictOwner(ds, r.id_1);
        note = r.id_1 + " already matches sample " + (k >= 0 ? ds.objs[k].pheno_id : r.id_1) + ".";
      }
    } else if (!ds.mapped[i]) {
      title.appendChild(code(r.pheno_id));
      if (c === "nodata") {
        title.appendChild(document.createTextNode(noData(r)));
      } else {
        title.appendChild(document.createTextNode(": no assigned ID"));
        barRow(bars, code(r.id_1), r.z_1, max, false, "best");
        if (r.id_2) barRow(bars, code(r.id_2), r.z_2, max, false, "next");
        if (c === "noid") note = leadNote(r);
        if (c === "noid" && r.id_1 && !ds.claimed[r.id_1]) note += " No other sample has " + r.id_1 + ".";
      }
    } else {
      title.appendChild(code(r.pheno_id));
      if (c === "nodata") {
        title.appendChild(document.createTextNode(noData(r)));
      } else {
        title.appendChild(document.createTextNode(isNum(r.rank_self) ? ": own ID at " + rankText(r.rank_self) : ""));
        var o = otherOf(r);
        barRow(bars, "own ID", r.z_self, max, true, "", st.zt);
        if (o.id) barRow(bars, code(o.id), o.z, max, false, "best");
        if (c === "nomatch" && o.id && o.id === r.id_1) {
          if (r.id_2) barRow(bars, code(r.id_2), r.z_2, max, false, "next");
          note = leadNote(r);
          var h = strictOwner(ds, r.id_1);
          if (h >= 0) note += " " + r.id_1 + " also matches its own sample" +
            (ds.objs[h].pheno_id !== r.id_1 ? " " + ds.objs[h].pheno_id : "") + ".";
        }
      }
    }
    a.appendChild(title);
    if (bars.childNodes.length) a.appendChild(bars);
    if (note) a.appendChild(el("div", "cc-note", note));
    a.title = "Top 5: " + top5(r);
    return a;
  }
  function unassignedCard(ds, idx) {
    var a = el("article", "check-card r-info");
    a.appendChild(el("div")).appendChild(plainChip("info", "+", "No assigned ID"));
    a.appendChild(el("div", "cc-title", n0(idx.length) + " samples without an assigned ID"));
    var counts = {}, unclaimed = 0;
    idx.forEach(function (i) {
      var c = ds.cats[i], r = ds.objs[i];
      counts[c] = (counts[c] || 0) + 1;
      if (c === "noid" && r.id_1 && !ds.claimed[r.id_1]) unclaimed++;
    });
    var list = el("div", "cc-counts");
    CATS.forEach(function (c) {
      if (!counts[c.key]) return;
      var b = el("button");
      b.type = "button";
      b.appendChild(chip(c.key, counts[c.key]));
      b.addEventListener("click", function () {
        setFilter(ds, { key: "un:" + c.key, test: function (d, i) { return !d.mapped[i] && d.cats[i] === c.key; } });
        openTable(ds);
      });
      list.appendChild(b);
    });
    a.appendChild(list);
    if (unclaimed) {
      a.appendChild(el("div", "cc-note", n0(unclaimed) + " of " + n0(counts.noid) + " with no clear match: " +
        "their best genotype isn't used by any other sample, so it is probably theirs (lead below " +
        st.zd.toFixed(2) + ", so no ID is suggested)."));
    }
    var b2 = el("button", "linkbtn", "Show in table →");
    b2.type = "button";
    b2.addEventListener("click", function () { setFilter(ds, TILE_FILTER.unassigned); openTable(ds); });
    a.appendChild(b2);
    return a;
  }
  function notScoredCard(ds, idx) {
    var a = el("article", "check-card r-muted");
    a.appendChild(el("div")).appendChild(chip("nodata", undefined));
    var ids = idx.map(function (i) { return ds.objs[i].pheno_id; });
    var why = idx.every(function (i) { return ds.objs[i].qpgen_status === D.notScored; }) ? " (no covariates)" : "";
    a.appendChild(el("div", "cc-title", n0(idx.length) + " samples not scored" + why));
    a.appendChild(el("div", "cc-note", ids.length <= 8 ? ids.join(", ") : ids.slice(0, 8).join(", ") + " …"));
    var b = el("button", "linkbtn", "Show in table →");
    b.type = "button";
    b.addEventListener("click", function () {
      setFilter(ds, { key: "cat:nodata", test: function (d, i) { return d.cats[i] === "nodata"; } });
      openTable(ds);
    });
    a.appendChild(b);
    return a;
  }
  function drawCards(ds) {
    var sec = document.getElementById("check-" + ds.slug);
    var grid = sec.querySelector(".check-grid"), more = sec.querySelector(".more-link");
    var assigned = [], un = [], unAll = [];
    ds.cats.forEach(function (c, i) {
      if (!ds.mapped[i]) unAll.push(i);
      if (!toCheck(ds, i)) return;
      (ds.mapped[i] || ds.evals[i].dup ? assigned : un).push(i);
    });
    var summarize = unAll.length > MAX_UNASSIGNED_CARDS;
    /* two or more not scored: one card (samples without an ID are in the summary card) */
    var nsIdx = assigned.concat(summarize ? [] : un).filter(function (i) { return ds.cats[i] === "nodata"; });
    var foldNoData = nsIdx.length >= 2;
    function keep(i) { return !(foldNoData && ds.cats[i] === "nodata"); }
    assigned = assigned.filter(keep); un = un.filter(keep);
    var idx = summarize ? assigned : assigned.concat(un);
    idx.sort(sevSort(ds));
    var total = assigned.length + un.length + (foldNoData ? nsIdx.length : 0);
    sec.querySelector(".cnt").textContent = "(" + n0(total) + ")";
    grid.innerHTML = ""; more.innerHTML = "";
    if (!total) {
      var e = el("div", "check-card empty r-good");
      e.appendChild(icon({ icon: "✓" }));
      e.appendChild(document.createTextNode("All " + n0(ds.mapped.filter(Boolean).length) + " match their own ID."));
      grid.appendChild(e);
      return;
    }
    if (summarize) grid.appendChild(unassignedCard(ds, unAll));
    if (foldNoData) grid.appendChild(notScoredCard(ds, nsIdx));
    var shown = idx.slice(0, MAX_CARDS - (summarize ? 1 : 0) - (foldNoData ? 1 : 0)), max = Math.max(1, st.zt);
    shown.forEach(function (i) {
      var r = ds.objs[i];
      [r.z_self, r.z_1, r.z_2].forEach(function (z) { if (isNum(z) && z > max) max = z; });
    });
    max *= 1.05;
    shown.forEach(function (i) { grid.appendChild(card(ds, i, max)); });
    if (shown.length < idx.length) {
      var b = el("button", "linkbtn", "Show all " + n0(total) + " in table →");
      b.type = "button";
      b.addEventListener("click", function () {
        setFilter(ds, { key: "check", test: toCheck });
        openTable(ds);
      });
      more.appendChild(b);
    }
  }

  /* ---------------------------------------------------------------- chart */
  /* x: best match score (z_1), y: runner-up score (z_2). Below the dashed line the best
     genotype leads the runner-up by more than the clear lead; on the diagonal they tie. */
  function drawLegend(ds, counts, anyOpen) {
    var row = document.querySelector("#chart-" + ds.slug + " .legend-row");
    row.innerHTML = "";
    CATS.forEach(function (c) {
      var k = c.key;
      if (!counts[k] && k !== "strict") return;
      var b = el("button");
      b.type = "button";
      b.setAttribute("aria-pressed", String(!ds.hidden[k]));
      b.title = (ds.hidden[k] ? "Show " : "Hide ") + c.label.toLowerCase();
      b.appendChild(chip(k, counts[k] || 0));
      b.addEventListener("click", function () { ds.hidden[k] = !ds.hidden[k]; drawChart(ds); });
      row.appendChild(b);
    });
    if (anyOpen) {
      var key = el("span", "open-key");
      key.appendChild(el("span", "ring"));
      key.appendChild(document.createTextNode("hollow: no assigned ID"));
      row.appendChild(key);
    }
  }
  function labelled(ds, i) {
    var c = ds.cats[i];
    return ds.evals[i].dup || c === "different" || c === "several" || c === "nomatch" || c === "taken" ||
      (c === "noid" && ds.nUnassigned <= MAX_UNASSIGNED_CARDS) || (c === "newid" && ds.nUnassigned <= 3);
  }
  function hoverText(r, label) {
    function who(id) { return id && id === r.geno_id_assigned ? " (own ID)" : ""; }
    return "<b>" + r.pheno_id + "</b> · " + label + (r.geno_id_assigned ? "" : " · no assigned ID") +
      "<br>best " + (r.id_1 || "–") + " " + f2(r.z_1) + who(r.id_1) +
      "<br>runner-up " + (r.id_2 || "–") + " " + f2(r.z_2) + who(r.id_2) +
      (r.geno_id_assigned && r.rank_self > 2 ? "<br>own ID " + f2(r.z_self) + " (" + rankText(r.rank_self) + ")" : "");
  }
  function drawChart(ds) {
    var fig = document.getElementById(ds.fig);
    if (!fig) return;
    var counts = {}, pts = {}, xs = [], ys = [], anyOpen = false;
    ds.objs.forEach(function (r, i) {
      var c = ds.cats[i];
      if (!isNum(r.z_1) || !isNum(r.z_2) || c === "nodata") return;
      counts[c] = (counts[c] || 0) + 1;
      var key = c + (ds.mapped[i] ? "" : ":open");
      (pts[key] = pts[key] || []).push(i);
      if (!ds.mapped[i]) anyOpen = true;
      xs.push(r.z_1); ys.push(r.z_2);
    });
    drawLegend(ds, counts, anyOpen);
    /* axes start just below the lowest score (at 0 or below when scores reach it) */
    function lowOf(v) {
      if (!v.length) return 0;
      var m = Math.min.apply(null, v);
      return m > 1 ? Math.floor((m - 0.3) * 2) / 2 : Math.min(0, Math.floor(m - 0.3));
    }
    var xlo = lowOf(xs), ylo = lowOf(ys), lo = Math.min(xlo, ylo);
    var xhi = xs.length ? Math.max.apply(null, xs) + 0.5 : 1;
    var yhi = ys.length ? Math.max.apply(null, ys) + 0.5 : 1;
    var traces = [];
    var labels = [];
    CHART_ORDER.forEach(function (k) {
      [false, true].forEach(function (open) {
        var sel = pts[k + (open ? ":open" : "")];
        if (!sel || ds.hidden[k]) return;
        var meta = CAT[k], color = C.role[meta.role], sym = SYMBOL[k] + (open ? "-open" : "");
        var marker;
        /* routine results small and light, so problem results stand out */
        if (open && k === "newid") marker = { color: color, size: 7, opacity: 0.7, symbol: sym, line: { color: color, width: 1.5 } };
        else if (open) marker = { color: color, size: k === "noid" ? 9 : 12, symbol: sym, line: { color: color, width: 2 } };
        else if (k === "strict") marker = { color: color, size: 6, opacity: 0.5, symbol: sym };
        else if (k === "lenient") marker = { color: color, size: 11, symbol: sym, line: { color: C.surface, width: 1.5 } };
        else marker = { color: color, size: k === "twoids" ? 16 : 13, symbol: sym, line: { color: C.surface, width: 1.5 } };
        traces.push({
          type: "scatter", mode: "markers", name: meta.label,
          x: sel.map(function (i) { return ds.objs[i].z_1; }),
          y: sel.map(function (i) { return ds.objs[i].z_2; }),
          customdata: sel.map(function (i) { return hoverText(ds.objs[i], meta.label); }),
          meta: sel.map(function (i) { return ds.objs[i].pheno_id; }),
          hovertemplate: "%{customdata}<extra></extra>",
          marker: marker
        });
        sel.forEach(function (i) { if (labelled(ds, i)) labels.push(i); });
      });
    });
    /* a sample that matches two genotype IDs but keeps its own ID first (a strict match) gets
       a ring, so that a duplicated genotype stands out */
    var ring = [];
    ds.objs.forEach(function (r, i) {
      if (ds.evals[i].dup && ds.cats[i] !== "twoids" && !ds.hidden[ds.cats[i]] && isNum(r.z_1) && isNum(r.z_2)) ring.push(i);
    });
    if (ring.length) traces.push({ type: "scatter", mode: "markers", name: "Matches two genotype IDs", hoverinfo: "skip",
      x: ring.map(function (i) { return ds.objs[i].z_1; }), y: ring.map(function (i) { return ds.objs[i].z_2; }),
      marker: { symbol: "star-open", size: 18, color: C.role.dup, line: { color: C.role.dup, width: 2 } } });
    if (ring.length) {
      var key = el("span", "open-key", "☆ own ID ties with a second genotype ID");
      document.querySelector("#chart-" + ds.slug + " .legend-row").appendChild(key);
    }
    labels.sort(sevSort(ds));
    var ph = window.innerWidth < 600;
    labels = labels.slice(0, ph ? 3 : MAX_LABELS);
    var height = ph ? 380 : 460, mg = D.chartLayout.margin;
    var plotW = Math.max(200, (fig.clientWidth || (ph ? 340 : 1100)) - mg.l - mg.r), plotH = height - mg.t - mg.b;
    var zd = st.zd;
    var ann = [
      { x: Math.min(xhi, yhi) - 0.2, y: Math.min(xhi, yhi) - 0.2, xref: "x", yref: "y", text: "tie", showarrow: false,
        xanchor: "right", yanchor: "bottom", xshift: -4, font: { color: C.ink2, size: 14 } },
      { x: xhi, y: ylo, xref: "x", yref: "y", text: "one genotype clearly best", showarrow: false, xanchor: "right",
        yanchor: "bottom", xshift: -8, yshift: 6, font: { color: C.ink2, size: 14 } }
    ];
    var pxX = plotW / Math.max(1, xhi - xlo), pxY = plotH / Math.max(1, yhi - ylo);
    /* the dashed line's label, left of the line and inside the plot (below the top edge) */
    var leadY = Math.min(yhi - 26 / pxY, xhi - zd - 0.15);
    /* the "tie" label goes when it would collide with this one (narrow plots) */
    var tieX = (Math.min(xhi, yhi) - 0.2) * pxX, leadX = (leadY + zd) * pxX;
    if (leadY > ylo + 0.5 && Math.abs(leadX - tieX) < 150 && (Math.min(xhi, yhi) - 0.2 - leadY) * pxY < 30) ann.shift();
    if (leadY > ylo + 0.5) ann.push({ x: leadY + zd, y: leadY, xref: "x", yref: "y",
      text: "clear lead " + zd.toFixed(2), showarrow: false, xanchor: "right", yanchor: "bottom", xshift: -6,
      font: { color: C.ink2, size: 14 } });
    /* label offsets (pixels): pick the spot whose box overlaps no other label or the plot
       edge, and the fewest points */
    var X0 = xlo * pxX, X1 = xhi * pxX, Y0 = -yhi * pxY, Y1 = -ylo * pxY;
    var hard = labels.map(function (i) { return [ds.objs[i].z_1 * pxX - 8, -ds.objs[i].z_2 * pxY - 8, 16, 16]; });
    hard.push([X0 - 400, Y0 - 400, X1 - X0 + 800, 400]);       /* above the plot */
    hard.push([X1, Y0 - 400, 400, Y1 - Y0 + 800]);             /* right of the plot */
    hard.push([X0 - 400, Y0 - 400, 400, Y1 - Y0 + 800]);       /* left of the plot */
    hard.push([X0 - 400, Y1, X1 - X0 + 800, 400]);             /* below the plot */
    hard.push([X1 - 210, Y1 - 28, 210, 28]);                   /* the shaded area's label */
    var tieAt = Math.min(xhi, yhi) - 0.2;
    hard.push([tieAt * pxX - 40, -tieAt * pxY - 22, 36, 22]);   /* "tie" */
    if (leadY > ylo + 0.5) hard.push([(leadY + zd) * pxX - 130, -leadY * pxY - 22, 124, 22]);
    var dots = [];
    ds.objs.forEach(function (r, i) {
      if (isNum(r.z_1) && isNum(r.z_2) && ds.cats[i] !== "nodata" && !ds.hidden[ds.cats[i]]) {
        dots.push([r.z_1 * pxX - 4, -r.z_2 * pxY - 4, 8, 8]);
      }
    });
    function overlaps(b, x0, y0, w, h) { return x0 < b[0] + b[2] && b[0] < x0 + w && y0 < b[1] + b[3] && b[1] < y0 + h; }
    var SPOTS = [[34, -22], [-34, -22], [34, 22], [-34, 22], [0, -34], [0, 34], [60, 0], [-60, 0],
                 [56, -44], [-56, -44], [56, 44], [-56, 44]];
    labels.forEach(function (i) {
      var r = ds.objs[i], x = r.z_1, y = r.z_2, w = 8.5 * String(r.pheno_id).length + 10, h = 20;
      var px = x * pxX, py = -y * pxY, pick = SPOTS[0], best = Infinity;
      SPOTS.forEach(function (sp) {
        var cx = px + sp[0], cy = py + sp[1];
        var x0 = sp[0] > 0 ? cx : (sp[0] < 0 ? cx - w : cx - w / 2), y0 = cy - h / 2, score = 0;
        hard.forEach(function (b) { if (overlaps(b, x0, y0, w, h)) score += 1000; });
        dots.forEach(function (b) { if (overlaps(b, x0, y0, w, h)) score += 1; });
        if (score < best) { best = score; pick = sp; }
      });
      var bx = px + pick[0];
      hard.push([pick[0] > 0 ? bx : (pick[0] < 0 ? bx - w : bx - w / 2), py + pick[1] - h / 2, w, h]);
      ann.push({ x: x, y: y, xref: "x", yref: "y", text: r.pheno_id, showarrow: true, arrowhead: 0, arrowwidth: 1,
        arrowcolor: C.ink2, standoff: 7, ax: pick[0], ay: pick[1],
        xanchor: pick[0] > 0 ? "left" : (pick[0] < 0 ? "right" : "center"), yanchor: "middle",
        font: { color: C.ink, size: 14 }, bgcolor: C.surface, borderpad: 2 });
    });
    var base = JSON.parse(JSON.stringify(D.chartLayout));
    base.height = height;
    base.xaxis.range = [xlo, xhi];
    base.yaxis.range = [ylo, yhi];
    var far = Math.max(xhi, yhi) + 50;
    base.shapes = [
      { type: "path", xref: "x", yref: "y", layer: "below", line: { width: 0 }, fillcolor: C.muted, opacity: 0.1,
        path: "M " + (lo + zd) + " " + lo + " L " + far + " " + lo + " L " + far + " " + (far - zd) + " Z" },
      { type: "line", xref: "x", yref: "y", x0: lo, y0: lo, x1: far, y1: far, layer: "below",
        line: { color: C.axis, width: 1, dash: "dot" } },
      { type: "line", xref: "x", yref: "y", x0: lo + zd, y0: lo, x1: far, y1: far - zd, layer: "below",
        line: { color: C.ink2, width: 1.5, dash: "dash" } }
    ];
    base.annotations = ann;
    fig.style.minHeight = base.height + "px";
    ANQ.renderFigure(ds.fig, { data: traces, layout: base });
    if (!fig._anqClick && fig.on) {
      fig._anqClick = true;
      fig.on("plotly_click", function (ev) {
        var p = ev && ev.points && ev.points[0];
        if (!p || !p.data || !p.data.meta) return;
        var id = p.data.meta[p.pointIndex];
        if (!id) return;
        ds.query = String(id).toLowerCase();
        ds.queryRaw = id;
        ds.filter = null;
        var inp = ds.panel.querySelector(".tbar input[type=search]");
        if (inp) inp.value = id;
        ds.page = 0;
        drawTable(ds, true);
        drawTiles(ds);
        openTable(ds);
      });
    }
  }

  /* ------------------------------------------------------------ duplicates */
  /* Samples that match two genotype IDs about equally (ANQ_MATCH.isDuplicate, at the
     sliders' values): the two genotypes are probably one person. */
  function drawDups(ds) {
    var sec = document.getElementById("dups-" + ds.slug);
    var list = sec.querySelector(".dup-list");
    var idx = [];
    ds.evals.forEach(function (e, i) { if (e.dup) idx.push(i); });
    sec.querySelector(".cnt").textContent = idx.length ? "(" + idx.length + ")" : "";
    list.innerHTML = "";
    if (!idx.length) {
      list.className = "dup-list is-empty";
      var e = el("div", "dup-empty");
      e.appendChild(icon({ icon: "✓" }));
      e.appendChild(document.createTextNode("No sample matches two genotype IDs."));
      e.className = "dup-empty r-good";
      list.appendChild(e);
      return;
    }
    list.className = "dup-list";
    var max = 1;
    idx.forEach(function (i) { var z = ds.objs[i].z_1; if (isNum(z) && z > max) max = z; });
    max *= 1.05;
    idx.forEach(function (i) {
      var r = ds.objs[i], row = el("div", "dup-row");
      var head = el("div", "dup-ids");
      head.appendChild(code(r.pheno_id));
      head.appendChild(el("span", "dup-sub", r.geno_id_assigned ? "" : "no assigned ID"));
      var bars = el("div", "bars");
      barRow(bars, idLabel(r, r.id_1), r.z_1, max, r.id_1 === r.geno_id_assigned, "");
      barRow(bars, idLabel(r, r.id_2), r.z_2, max, r.id_2 === r.geno_id_assigned, "");
      barRow(bars, "next", r.z_3, max, false, "");
      var flag = el("div", "dup-flag");
      flag.appendChild(chip(ds.cats[i]));
      row.appendChild(head); row.appendChild(bars); row.appendChild(flag);
      list.appendChild(row);
    });
  }

  /* ---------------------------------------------------------- all samples */
  /* lead of the row: own ID over the best other genotype, or (no assigned ID) the best
     genotype over the next */
  function leadOf(r) {
    if (r.geno_id_assigned) return ownAhead(r);
    return isNum(r.z_1) && isNum(r.z_2) ? r.z_1 - r.z_2 : null;
  }
  function signed(m) { return (m > 0 ? "+" : m < 0 ? "−" : "") + Math.abs(m).toFixed(2); }
  /* w: fixed width in px (none: share the rest); ph: hidden at phone width */
  var COLS = [
    { k: "sample", label: "Sample",
      val: function (ds, i) { return ds.objs[i].pheno_id; },
      cell: function (ds, i) {
        var r = ds.objs[i], c = ds.cats[i];
        if (!phone() || (c === "strict" && ds.mapped[i]) || !r.id_1) return [r.pheno_id, "idc"];
        /* phone: the genotype the result is about, under the sample */
        var o = otherOf(r), box = el("span", "two");
        box.appendChild(el("span", "", r.pheno_id));
        var sub = el("span", "sub", "→ " + o.id);
        sub.title = (ds.mapped[i] ? "best other " : "best ") + o.id + " " + f2(o.z);
        box.appendChild(sub);
        return [box, "idc"];
      } },
    { k: "assigned", label: "Assigned ID", ph: true, show: function (ds) { return !ds.sameIds; },
      val: function (ds, i) { return ds.objs[i].geno_id_assigned; },
      cell: function (ds, i) {
        var r = ds.objs[i], same = r.geno_id_assigned === r.pheno_id;
        return [same || !r.geno_id_assigned ? "" : r.geno_id_assigned, "idc"];
      } },
    { k: "result", label: "Result", w: 250, pw: 118,
      val: function (ds, i) { return CAT[ds.cats[i]].order; },
      cell: function (ds, i) {
        var c = ds.cats[i], e = ds.evals[i];
        if (c === "strict") {
          var s = el("span", "res-strict r-good");
          s.appendChild(icon(CAT.strict));
          s.appendChild(el("span", "cl", phone() ? "Strict" : CAT.strict.label));
          s.title = CAT.strict.label + " (" + e.status + ")";
          return [s, "res"];
        }
        return [chip(c, undefined, CAT[c].label + " (" + e.status + ")"), "res"];
      } },
    { k: "own", label: "Own score", plabel: "Own", num: true, w: 124, pw: 64,
      val: function (ds, i) { return ds.objs[i].z_self; },
      cell: function (ds, i) {
        return ds.mapped[i] ? [f2(ds.objs[i].z_self), "num"] : ["no ID", "num muted"];
      } },
    { k: "rank", label: "Own rank", num: true, w: 116, ph: true,
      val: function (ds, i) { return ds.objs[i].rank_self; },
      cell: function (ds, i) {
        var v = ds.objs[i].rank_self;
        return [isNum(v) ? n0(v) : "–", "num" + (isNum(v) && v > 1 ? " warn" : " muted")];
      } },
    { k: "other", label: "Best other", ph: true, gap: true,
      val: function (ds, i) { return otherOf(ds.objs[i]).id; },
      cell: function (ds, i) {
        var r = ds.objs[i];
        return [otherOf(r).id || "–", "idc gapl" + (ds.cats[i] === "different" ? " bad" : "")];
      } },
    { k: "otherz", label: "Its score", num: true, w: 112, ph: true,
      val: function (ds, i) { return otherOf(ds.objs[i]).z; },
      cell: function (ds, i) { return [f2(otherOf(ds.objs[i]).z), "num"]; } },
    { k: "ahead", label: "Lead", num: true, w: 92, pw: 62,
      title: "Own score minus the best other's (negative: another genotype fits better). " +
        "Without an assigned ID: the best genotype's lead over the next, in italics.",
      val: function (ds, i) { return leadOf(ds.objs[i]); },
      cell: function (ds, i) {
        var r = ds.objs[i], m = leadOf(r);
        if (m === null) return ["–", "num muted"];
        if (!ds.mapped[i]) return [signed(m), "num best-lead"];
        /* red only where another genotype wins (not for accepted lenient matches) */
        var c = ds.cats[i], lost = c === "different" || c === "nomatch" || c === "several";
        return [signed(m), "num" + (m < 0 && lost ? " bad" : "")];
      } },
    { k: "suggest", label: "Suggested ID", ph: true, gap: true,
      show: function (ds) { return ds.evals.some(function (e) { return !!e.proposed; }); },
      val: function (ds, i) { return ds.evals[i].proposed; },
      cell: function (ds, i) { return [ds.evals[i].proposed || "", "idc gapl"]; } }
  ];
  function visibleCols(ds) {
    var p = phone();
    return COLS.filter(function (c) { return (!c.show || c.show(ds)) && !(p && c.ph); });
  }
  function setFilter(ds, f) {
    ds.filter = f;
    ds.page = 0;
    drawTable(ds, true);
    drawTiles(ds);
  }
  function openTable(ds) {
    var d = document.getElementById("all-" + ds.slug);
    d.open = true;
    d.scrollIntoView({ behavior: "smooth", block: "start" });
  }
  function view(ds) {
    var q = ds.query, idx = [];
    for (var i = 0; i < ds.objs.length; i++) {
      if (ds.filter && !ds.filter.test(ds, i)) continue;
      if (q) {
        var r = ds.objs[i];
        if ([r.pheno_id, r.geno_id_assigned, r.id_1, r.id_2].join("\u0000").toLowerCase().indexOf(q) < 0) continue;
      }
      idx.push(i);
    }
    var col = COLS.filter(function (c) { return c.k === ds.sort.col; })[0], dir = ds.sort.dir;
    var sev = sevSort(ds);
    idx.sort(function (a, b) {
      if (col.k === "result") return sev(a, b) * dir;
      var x = col.val(ds, a), y = col.val(ds, b);
      var xm = x === null || x === undefined || x === "", ym = y === null || y === undefined || y === "";
      if (xm || ym) return xm === ym ? 0 : (xm ? 1 : -1);
      if (typeof x === "number" && typeof y === "number") return (x - y) * dir;
      return String(x).localeCompare(String(y), undefined, { numeric: true }) * dir;
    });
    return idx;
  }
  function pager(ds) {
    var box = el("span", "tpager");
    var prev = el("button", "", "‹"), info = el("span", "num"), next = el("button", "", "›");
    prev.type = next.type = "button";
    prev.setAttribute("aria-label", "Previous page"); next.setAttribute("aria-label", "Next page");
    prev.addEventListener("click", function () { ds.page = Math.max(0, ds.page - 1); drawTable(ds, true); });
    next.addEventListener("click", function () { ds.page += 1; drawTable(ds, true); });
    box.appendChild(prev); box.appendChild(info); box.appendChild(next);
    return { box: box, prev: prev, next: next, info: info };
  }
  function buildToolbar(ds) {
    var bar = ds.panel.querySelector(".tbar");
    bar.innerHTML = "";
    var row1 = el("div", "tbar-row");
    var inp = el("input");
    inp.type = "search"; inp.placeholder = "Search ID…"; inp.setAttribute("aria-label", "Search ID");
    inp.value = ds.queryRaw || "";
    inp.addEventListener("input", function () {
      ds.queryRaw = inp.value; ds.query = inp.value.trim().toLowerCase(); ds.page = 0; drawTable(ds, true);
    });
    row1.appendChild(inp);
    var top = pager(ds);
    row1.appendChild(top.box);
    var size = el("span", "tsize");
    size.setAttribute("role", "group"); size.setAttribute("aria-label", "Rows per page");
    size.appendChild(document.createTextNode("Rows"));
    [[100, "100"], [500, "500"], [0, "All"]].forEach(function (o) {
      var b = el("button", "", o[1]);
      b.type = "button"; b._size = o[0];
      b.addEventListener("click", function () { ds.pageSize = o[0]; ds.page = 0; drawTable(ds, true); });
      size.appendChild(b);
    });
    row1.appendChild(size);
    var dl = el("button", "btn", "Download TSV");
    dl.type = "button";
    dl.addEventListener("click", function () { download(ds); });
    row1.appendChild(dl);
    bar.appendChild(row1);
    var chips = el("div", "fchips");
    chips.setAttribute("role", "group"); chips.setAttribute("aria-label", "Show");
    bar.appendChild(chips);
    var foot = ds.panel.querySelector(".tfoot");
    foot.innerHTML = "";
    var bottom = pager(ds);
    foot.appendChild(bottom.box);
    ds.ui = { pagers: [top, bottom], size: size, chips: chips };
  }
  function filterButton(ds, key, content, test) {
    var b = el("button");
    b.type = "button";
    if (typeof content === "string") { b.className = "all"; b.textContent = content; } else b.appendChild(content);
    b.setAttribute("aria-pressed", String(key === null ? !ds.filter : !!(ds.filter && ds.filter.key === key)));
    b.addEventListener("click", function () { setFilter(ds, key === null ? null : { key: key, test: test }); });
    return b;
  }
  function drawChips(ds) {
    var box = ds.ui.chips, counts = {};
    ds.cats.forEach(function (c) { counts[c] = (counts[c] || 0) + 1; });
    box.innerHTML = "";
    box.appendChild(filterButton(ds, null, "All " + n0(ds.objs.length)));
    CATS.forEach(function (c) {
      if (!counts[c.key]) return;
      box.appendChild(filterButton(ds, "cat:" + c.key, chip(c.key, counts[c.key] || 0),
        function (d, i) { return d.cats[i] === c.key; }));
    });
    if (ds.nUnassigned) {
      box.appendChild(filterButton(ds, "unassigned", plainChip("info", "+", "No assigned ID " + n0(ds.nUnassigned)),
        TILE_FILTER.unassigned.test));
    }
  }
  function drawHead(ds, cols) {
    var table = ds.panel.querySelector("table.mtable");
    var cg = table.querySelector("colgroup"), thead = table.querySelector("thead");
    var p = phone();
    cg.innerHTML = "";
    cols.forEach(function (c) {
      var col = el("col"), w = p ? (c.pw || c.w) : c.w;
      if (w) col.style.width = w + "px";
      cg.appendChild(col);
    });
    thead.innerHTML = "";
    var tr = el("tr");
    cols.forEach(function (c) {
      var th = el("th", (c.num ? "num" : "") + (c.gap ? " gapl" : ""));
      th.scope = "col";
      if (c.title) th.title = c.title;
      th.tabIndex = 0;
      var lab = el("span", "th-l", p && c.plabel ? c.plabel : c.label);
      lab.appendChild(el("span", "sorti"));
      th.appendChild(lab);
      th.setAttribute("aria-sort", ds.sort.col === c.k ? (ds.sort.dir > 0 ? "ascending" : "descending") : "none");
      var act = function () {
        if (ds.sort.col === c.k) ds.sort.dir = -ds.sort.dir; else ds.sort = { col: c.k, dir: c.num ? -1 : 1 };
        ds.page = 0;
        drawTable(ds, true);
      };
      th.addEventListener("click", act);
      th.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); } });
      tr.appendChild(th);
    });
    thead.appendChild(tr);
  }
  function drawTable(ds, force) {
    var sig = ds.cats.join(",") + "|" + ds.evals.map(function (e) { return e.proposed; }).join(",") + "|" +
      ds.query + "|" + (ds.filter ? ds.filter.key : "") + "|" + ds.sort.col + ds.sort.dir + "|" + ds.page + "|" +
      ds.pageSize + "|" + phone();
    if (!force && sig === ds.sig) return;
    ds.sig = sig;
    if (!ds.ui) buildToolbar(ds);
    var cols = visibleCols(ds);
    drawChips(ds);
    drawHead(ds, cols);
    var idx = view(ds), n = idx.length;
    var ps = ds.pageSize > 0 ? ds.pageSize : Math.max(1, n);
    var pages = Math.max(1, Math.ceil(n / ps));
    if (ds.page >= pages) ds.page = pages - 1;
    var start = ds.page * ps, end = Math.min(n, start + ps);
    ds.ui.pagers.forEach(function (pg, j) {
      pg.info.textContent = n ? n0(start + 1) + "–" + n0(end) + " of " + n0(n) : "No rows";
      pg.prev.disabled = ds.page === 0;
      pg.next.disabled = ds.page >= pages - 1;
      if (j === 1) pg.box.hidden = pages <= 1;
    });
    ds.ui.size.querySelectorAll("button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b._size === ds.pageSize));
    });
    var tbody = ds.panel.querySelector("table.mtable tbody");
    tbody.innerHTML = "";
    var frag = document.createDocumentFragment();
    for (var k = start; k < end; k++) {
      var i = idx[k], r = ds.objs[i], c = ds.cats[i], e = ds.evals[i];
      var tr = el("tr", c === "strict" && !e.dup ? "" : "flag r-" + (e.dup ? "dup" : CAT[c].role));
      var runCode = r.qpgen_status;
      tr.title = "qpgen: " + e.status + (e.status !== runCode ? " (saved: " + runCode + ")" : "") +
        " · traits " + (isNum(r.n_traits) ? r.n_traits : "–") + " · cor " + f2(r.cor_self) + "\nTop 5: " + top5(r);
      cols.forEach(function (col) {
        var cell = col.cell(ds, i), td = el("td", cell[1]);
        if (typeof cell[0] === "string") {
          td.textContent = cell[0];
          if (col.k === "sample" || col.k === "assigned" || col.k === "other" || col.k === "suggest") td.title = cell[0];
        } else td.appendChild(cell[0]);
        tr.appendChild(td);
      });
      frag.appendChild(tr);
    }
    tbody.appendChild(frag);
  }
  function download(ds) {
    var head = ["sample", "assigned_id", "result", "status_now", "status_run", "own_score", "own_rank", "best_other_id",
                "best_other_score", "own_ahead_by", "best_id", "best_score", "best_lead_over_next",
                "two_genotype_ids", "suggested_id"];
    var lines = [head.join("\t")];
    function v(x) { return x === null || x === undefined || x === "" ? "NA" : String(x); }
    view(ds).forEach(function (i) {
      var r = ds.objs[i], e = ds.evals[i], m = ownAhead(r), o = otherOf(r);
      lines.push([r.pheno_id, v(r.geno_id_assigned), CAT[ds.cats[i]].label, e.status, r.qpgen_status, v(r.z_self),
                  v(r.rank_self), v(o.id), v(o.z), m === null ? "NA" : m.toFixed(6), v(r.id_1), v(r.z_1),
                  e.gap === null ? "NA" : e.gap.toFixed(6), e.dup ? "yes" : "no", v(e.proposed)].join("\t"));
    });
    var blob = new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "match_samples" + (D.datasets.length > 1 ? "." + ds.slug : "") + ".tsv";
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }

  /* ----------------------------------------------------------- thresholds */
  var zt = document.getElementById("anq-zt"), ztR = document.getElementById("anq-zt-r");
  var zd = document.getElementById("anq-zd"), zdR = document.getElementById("anq-zd-r");
  var pending = false;
  function schedule() {
    if (pending) return;
    pending = true;
    (window.requestAnimationFrame || function (f) { setTimeout(f, 16); })(function () { pending = false; recompute(); });
  }
  function bind(num, range, key) {
    range.addEventListener("input", function () {
      st[key] = parseFloat(range.value);
      num.value = st[key].toFixed(2);
      schedule();
    });
    num.addEventListener("input", function () {
      var v = parseFloat(num.value), bad = !isFinite(v) || v < 0;
      num.classList.toggle("invalid", bad);
      num.setCustomValidity(bad ? "Enter a number of 0 or more" : "");
      if (bad) return;
      st[key] = v;
      if (v > parseFloat(range.max)) {
        /* widen the slider so the thumb shows the typed value; keep the run mark in place */
        range.max = String(Math.ceil(v));
        var mark = range.parentNode.querySelector(".th-mark");
        if (mark) mark.style.left = "calc(8px + (100% - 16px) * " + (RUN[key === "zt" ? "z_threshold" : "z_diff"] /
          parseFloat(range.max)).toFixed(4) + ")";
      }
      range.value = String(v);
      schedule();
    });
    num.addEventListener("change", function () {
      num.value = st[key].toFixed(2);
      num.classList.remove("invalid");
      num.setCustomValidity("");
    });
  }
  function drawThresholdState() {
    var off = changed();
    document.getElementById("anq-reset").disabled = !off;
    document.getElementById("anq-th-note").hidden = !off;
    ztR.setAttribute("aria-valuetext", st.zt.toFixed(2));
    zdR.setAttribute("aria-valuetext", st.zd.toFixed(2));
  }
  function yaml() {
    return "pipelines:\n  match:\n    z_threshold: " + Number(st.zt.toFixed(2)) + "\n    z_diff: " +
      Number(st.zd.toFixed(2)) + "\n";
  }
  function copyText(text, done) {
    function fallback() {
      var ta = document.createElement("textarea");
      ta.value = text; ta.setAttribute("readonly", ""); ta.style.position = "absolute"; ta.style.left = "-9999px";
      document.body.appendChild(ta); ta.select();
      var ok = false;
      try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
      ta.remove(); done(ok);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(function () { done(true); }, fallback);
    } else { fallback(); }
  }

  /* ------------------------------------------------------------ recompute */
  function recompute() {
    D.datasets.forEach(function (ds) {
      ds.evals = M.evaluateAll(ds.objs, st.zt, st.zd);
      ds.cats = ds.evals.map(function (e) { return e.cat; });
      drawTiles(ds);
      drawCards(ds);
      drawChart(ds);
      drawDups(ds);
      drawTable(ds, false);
    });
    drawThresholdState();
  }

  function showDataset(slug) {
    D.datasets.forEach(function (ds) {
      var on = ds.slug === slug;
      ds.panel.hidden = !on;
      var ctx = document.querySelector("[data-ctx='" + ds.slug + "']");
      if (ctx) ctx.hidden = !on;
      if (on) {
        active = ds;
        ds.panel.querySelector(".th-slot").appendChild(document.getElementById("anq-th"));
      }
    });
    document.querySelectorAll("[data-ds-set]").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.getAttribute("data-ds-set") === slug));
    });
    if (window.Plotly) active.panel.querySelectorAll(".plotly-fig").forEach(function (f) {
      if (f.offsetParent) window.Plotly.Plots.resize(f);
    });
  }

  ANQ.onReady(function () {
    bind(zt, ztR, "zt");
    bind(zd, zdR, "zd");
    document.getElementById("anq-reset").addEventListener("click", function () {
      st.zt = RUN.z_threshold; st.zd = RUN.z_diff;
      zt.value = st.zt.toFixed(2); zd.value = st.zd.toFixed(2);
      ztR.value = String(st.zt); zdR.value = String(st.zd);
      recompute();
    });
    document.getElementById("anq-copy").addEventListener("click", function () {
      var msg = document.getElementById("anq-copy-msg");
      copyText(yaml(), function (ok) {
        msg.textContent = ok ? "Copied." : "Copy failed.";
        setTimeout(function () { msg.textContent = ""; }, 2500);
      });
    });
    D.datasets.forEach(function (ds) {
      ds.panel.querySelectorAll(".dash-tile").forEach(function (b) {
        b.addEventListener("click", function () {
          var f = TILE_FILTER[b.getAttribute("data-tile")];
          setFilter(ds, f && ds.filter && ds.filter.key === f.key ? null : f);
          openTable(ds);
        });
      });
      ds.panel.querySelectorAll("details").forEach(function (d) {
        d.addEventListener("toggle", function () {
          if (d.open && window.Plotly) d.querySelectorAll(".plotly-fig").forEach(function (f) { window.Plotly.Plots.resize(f); });
        });
      });
    });
    document.querySelectorAll("[data-ds-set]").forEach(function (b) {
      b.addEventListener("click", function () { showDataset(b.getAttribute("data-ds-set")); });
    });
    ANQ.onTheme(function () { D.datasets.forEach(function (ds) { if (ds.evals) drawChart(ds); }); });
    var resizeTimer = null, lastPhone = window.innerWidth < 600, lastNarrow = phone();
    window.addEventListener("resize", function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () {
        var ph = window.innerWidth < 600;
        if (ph !== lastPhone) { lastPhone = ph; D.datasets.forEach(drawChart); }
        if (phone() !== lastNarrow) { lastNarrow = phone(); D.datasets.forEach(function (ds) { drawTable(ds, true); }); }
      }, 150);
    });
    showDataset(D.datasets[0].slug);
    D.datasets.forEach(function (ds) {
      /* the page's rules at the run's thresholds against the saved statuses */
      var n = 0;
      ds.objs.forEach(function (r, i) { if (r.qpgen_status && ds.runEvals[i].status !== r.qpgen_status) n++; });
      var al = ds.panel.querySelector(".js-disagree");
      if (al && n) {
        al.textContent = "⚠ " + n + " sample" + (n === 1 ? "" : "s") +
          ": the page's result differs from the saved result (hover the table row).";
        al.hidden = false;
      }
    });
    recompute();
  });
})();
