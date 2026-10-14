/* anqchor report runtime: theme, Plotly figures, tables. No external dependencies
   besides plotly.js. Exposes window.ANQ for page-specific scripts. */
(function () {
  "use strict";
  var T = window.ANQ_THEME || {};
  var ANQ = (window.ANQ = { figures: {}, tables: {}, _ready: [] });
  /* Page scripts register callbacks with ANQ.onReady(fn); they run after tables and
     figures are rendered (immediately when that has already happened). */
  ANQ._booted = false;
  ANQ.onReady = function (fn) {
    if (ANQ._booted) {
      try { fn(); } catch (e) { console.error("anqchor: page script", e); }
    } else {
      ANQ._ready.push(fn);
    }
  };

  /* ---------------------------------------------------------------- theme */
  function storedTheme() {
    try { return localStorage.getItem("anq-theme") || "auto"; } catch (e) { return "auto"; }
  }
  function isDark() {
    var t = document.documentElement.getAttribute("data-theme");
    if (t === "dark") return true;
    if (t === "light") return false;
    return !!(window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  }
  function setTheme(mode) {
    if (mode === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", mode);
    try { localStorage.setItem("anq-theme", mode); } catch (e) { /* storage unavailable */ }
    document.querySelectorAll("[data-theme-set]").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.getAttribute("data-theme-set") === mode));
    });
    rethemeAll();
  }
  ANQ.isDark = isDark;
  ANQ.setTheme = setTheme;

  var colorMap = null;
  function buildColorMap() {
    var m = {};
    var L = T.light || {}, D = T.dark || {};
    Object.keys(L).forEach(function (k) { if (D[k]) m[L[k].toLowerCase()] = D[k]; });
    var s = T.seriesMap || {};
    Object.keys(s).forEach(function (k) { m[k.toLowerCase()] = s[k]; });
    return m;
  }
  function mapColors(obj) {
    if (Array.isArray(obj)) return obj.map(mapColors);
    if (obj && typeof obj === "object") {
      var o = {};
      Object.keys(obj).forEach(function (k) { o[k] = mapColors(obj[k]); });
      return o;
    }
    if (typeof obj === "string" && obj.charAt(0) === "#") {
      var r = colorMap[obj.toLowerCase()];
      return r || obj;
    }
    return obj;
  }
  function themed(spec) {
    if (!isDark()) return spec;
    if (!colorMap) colorMap = buildColorMap();
    return mapColors(spec);
  }

  /* -------------------------------------------------------------- figures */
  var PLOT_CONFIG = { responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toggleSpikelines"] };

  var warnedNoPlotly = false;
  function renderFigure(id, spec) {
    var el = document.getElementById(id);
    if (!el) return;
    if (!window.Plotly) {
      /* plotlyjs: cdn and no network (HPC node, offline viewer, e-mailed page) */
      ANQ.figures[id] = spec;
      if (!el.querySelector(".plotly-missing")) {
        var msg = document.createElement("p");
        msg.className = "plotly-missing";
        msg.textContent = "This chart needs plotly.js from the CDN, which could not be loaded (offline?). " +
          "Re-render the report with --plotlyjs embed to view it offline.";
        el.appendChild(msg);
      }
      if (!warnedNoPlotly && window.console) {
        warnedNoPlotly = true;
        console.warn("anqchor: plotly.js is not loaded; charts are not drawn");
      }
      return;
    }
    ANQ.figures[id] = spec;
    var s = themed(spec);
    window.Plotly.react(el, s.data, s.layout, PLOT_CONFIG);
  }
  var themeListeners = [];
  function rethemeAll() {
    themeListeners.forEach(function (fn) {
      try { fn(isDark()); } catch (e) { console.error("anqchor: theme listener", e); }
    });
    Object.keys(ANQ.figures).forEach(function (id) { renderFigure(id, ANQ.figures[id]); });
    Object.keys(ANQ.tables).forEach(function (id) { ANQ.tables[id].draw(); });
  }
  /* Page scripts that build figure specs themselves register fn(isDark) here; it runs
     before the stored figures are re-themed. */
  ANQ.onTheme = function (fn) { themeListeners.push(fn); };
  ANQ.renderFigure = renderFigure;
  ANQ.rethemeAll = rethemeAll;
  ANQ.color = function (hex) {
    if (!isDark()) return hex;
    if (!colorMap) colorMap = buildColorMap();
    return colorMap[String(hex).toLowerCase()] || hex;
  };

  /* --------------------------------------------------------------- tables */
  function fmt(v, digits) {
    if (v === null || v === undefined || v === "") return "";
    if (typeof v === "number") {
      if (Number.isInteger(v)) return v.toLocaleString("en-US");
      var a = Math.abs(v);
      if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(Math.max(1, digits - 1));
      return String(Number(v.toPrecision(digits)));
    }
    return String(v);
  }
  var STATUS_GROUP = { PASS: "good", WARN: "warning", FAIL: "critical", NO_TEST: "info" };
  function cellNode(v, col, spec) {
    var td = document.createElement("td");
    if (typeof v === "number") td.className = "num";
    if (spec.status.indexOf(col) >= 0 && v) {
      var chip = document.createElement("span");
      chip.className = "chip";
      var dot = document.createElement("span");
      dot.className = "status-dot status-" + (STATUS_GROUP[v] || "info");
      dot.setAttribute("aria-hidden", "true");
      chip.appendChild(dot);
      chip.appendChild(document.createTextNode(String(v)));
      td.appendChild(chip);
    } else {
      td.textContent = fmt(v, spec.digits);
    }
    return td;
  }
  function toTsv(spec, rows) {
    var lines = [spec.columns.join("\t")];
    rows.forEach(function (r) {
      lines.push(r.map(function (v) { return v === null || v === undefined ? "NA" : String(v); }).join("\t"));
    });
    return lines.join("\n") + "\n";
  }

  function Table(container, spec) {
    this.el = container;
    this.spec = spec;
    this.rows = spec.rows.slice();
    this.view = this.rows;
    this.sortCol = -1;
    this.sortDir = 1;
    this.page = 0;
    this.query = "";
    this.pageSize = spec.pageSize;
    this.build();
  }
  Table.prototype.build = function () {
    var self = this, spec = this.spec;
    this.el.innerHTML = "";
    var controls = document.createElement("div");
    controls.className = "tbl-controls";
    if (spec.searchable) {
      var input = document.createElement("input");
      input.type = "search";
      input.placeholder = "Filter rows";
      input.setAttribute("aria-label", "Filter rows");
      input.addEventListener("input", function () { self.query = input.value.toLowerCase(); self.page = 0; self.refresh(); });
      controls.appendChild(input);
    }
    this.info = document.createElement("span");
    controls.appendChild(this.info);
    this.prev = document.createElement("button");
    this.prev.type = "button"; this.prev.textContent = "Previous";
    this.prev.addEventListener("click", function () { self.page = Math.max(0, self.page - 1); self.draw(); });
    this.next = document.createElement("button");
    this.next.type = "button"; this.next.textContent = "Next";
    this.next.addEventListener("click", function () { self.page += 1; self.draw(); });
    var dl = document.createElement("button");
    dl.type = "button"; dl.textContent = "Download TSV";
    dl.addEventListener("click", function () {
      var blob = new Blob([toTsv(spec, self.view)], { type: "text/tab-separated-values" });
      var a = document.createElement("a");
      a.href = URL.createObjectURL(blob); a.download = spec.filename;
      document.body.appendChild(a); a.click();
      setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
    });
    controls.appendChild(this.prev); controls.appendChild(this.next);
    if (spec.rows.length > 100) {
      /* rows per page: 100 | 500 | All */
      var seg = document.createElement("span");
      seg.className = "tbl-pagesize";
      seg.setAttribute("role", "group");
      seg.setAttribute("aria-label", "Rows per page");
      seg.appendChild(document.createTextNode("Rows "));
      this.sizeButtons = [];
      [[100, "100"], [500, "500"], [0, "All"]].forEach(function (o) {
        var b = document.createElement("button");
        b.type = "button"; b.textContent = o[1];
        b.addEventListener("click", function () { self.pageSize = o[0]; self.page = 0; self.draw(); });
        b._size = o[0];
        self.sizeButtons.push(b);
        seg.appendChild(b);
      });
      controls.appendChild(seg);
    }
    controls.appendChild(dl);
    this.el.appendChild(controls);
    var wrap = document.createElement("div");
    wrap.className = "tbl-wrap";
    this.table = document.createElement("table");
    var thead = document.createElement("thead"), tr = document.createElement("tr");
    spec.columns.forEach(function (c, i) {
      var th = document.createElement("th");
      th.textContent = c;
      th.scope = "col";
      th.tabIndex = 0;
      if (spec.help && spec.help[c]) th.title = spec.help[c];
      var act = function () { self.sortBy(i); };
      th.addEventListener("click", act);
      th.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); act(); } });
      tr.appendChild(th);
    });
    thead.appendChild(tr);
    this.thead = thead;
    this.tbody = document.createElement("tbody");
    this.table.appendChild(thead); this.table.appendChild(this.tbody);
    wrap.appendChild(this.table);
    this.el.appendChild(wrap);
    this.refresh();
  };
  Table.prototype.setRows = function (rows) { this.rows = rows.slice(); this.refresh(); };
  Table.prototype.sortBy = function (i) {
    if (this.sortCol === i) this.sortDir = -this.sortDir; else { this.sortCol = i; this.sortDir = 1; }
    this.refresh();
  };
  Table.prototype.refresh = function () {
    var q = this.query, i = this.sortCol, d = this.sortDir;
    var rows = q ? this.rows.filter(function (r) {
      return r.some(function (v) { return v !== null && String(v).toLowerCase().indexOf(q) >= 0; });
    }) : this.rows.slice();
    if (i >= 0) {
      rows.sort(function (a, b) {
        var x = a[i], y = b[i];
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        if (typeof x === "number" && typeof y === "number") return (x - y) * d;
        return String(x).localeCompare(String(y), undefined, { numeric: true }) * d;
      });
    }
    var self = this;
    Array.prototype.forEach.call(this.thead.querySelectorAll("th"), function (th, k) {
      th.setAttribute("aria-sort", k === self.sortCol ? (self.sortDir > 0 ? "ascending" : "descending") : "none");
    });
    this.view = rows;
    this.draw();
  };
  Table.prototype.draw = function () {
    var spec = this.spec, n = this.view.length;
    var ps = this.pageSize > 0 ? this.pageSize : Math.max(1, n);
    var cur = this.pageSize;
    (this.sizeButtons || []).forEach(function (b) { b.setAttribute("aria-pressed", String(b._size === cur)); });
    var pages = Math.max(1, Math.ceil(n / ps));
    if (this.page >= pages) this.page = pages - 1;
    var start = this.page * ps, end = Math.min(n, start + ps);
    this.tbody.innerHTML = "";
    var frag = document.createDocumentFragment();
    for (var r = start; r < end; r++) {
      var tr = document.createElement("tr");
      var row = this.view[r];
      for (var c = 0; c < spec.columns.length; c++) tr.appendChild(cellNode(row[c], spec.columns[c], spec));
      frag.appendChild(tr);
    }
    this.tbody.appendChild(frag);
    this.info.textContent = n === 0 ? "No rows" : (start + 1) + "–" + end + " of " + n.toLocaleString("en-US") + " rows";
    this.prev.disabled = this.page === 0;
    this.next.disabled = this.page >= pages - 1;
  };
  ANQ.renderTable = function (id, spec) {
    var el = document.getElementById(id);
    if (!el) return null;
    var t = new Table(el, spec);
    ANQ.tables[id] = t;
    return t;
  };
  ANQ.fmt = fmt;

  /* ----------------------------------------------------------------- boot */
  function boot() {
    var mode = storedTheme();
    if (mode !== "auto") document.documentElement.setAttribute("data-theme", mode);
    document.querySelectorAll("[data-theme-set]").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b.getAttribute("data-theme-set") === mode));
      b.addEventListener("click", function () { setTheme(b.getAttribute("data-theme-set")); });
    });
    document.querySelectorAll("script.table-data").forEach(function (s) {
      try { ANQ.renderTable(s.getAttribute("data-target"), JSON.parse(s.textContent)); }
      catch (e) { console.error("anqchor: table", e); }
    });
    document.querySelectorAll("script.fig-data").forEach(function (s) {
      try { renderFigure(s.getAttribute("data-target"), JSON.parse(s.textContent)); }
      catch (e) { console.error("anqchor: figure", e); }
    });
    if (window.matchMedia) {
      var mq = window.matchMedia("(prefers-color-scheme: dark)");
      var h = function () { if (!document.documentElement.getAttribute("data-theme")) rethemeAll(); };
      if (mq.addEventListener) mq.addEventListener("change", h); else if (mq.addListener) mq.addListener(h);
    }
    ANQ._booted = true;
    ANQ._ready.forEach(function (fn) {
      try { fn(); } catch (e) { console.error("anqchor: page script", e); }
    });
    ANQ._ready = [];
  }
  ANQ.boot = boot;
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
