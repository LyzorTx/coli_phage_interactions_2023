// main.js — vanilla JS + Alpine.js glue for the CHISEL explainability UI.
//
// Always fetches from the sibling ./data/ directory (same-origin). GitHub Release asset
// URLs do not serve Access-Control-Allow-Origin headers, so a cross-origin fetch from
// Pages is blocked by the browser; the publish-explainability-ui workflow bakes the
// release's JSONs + Parquets into the Pages deploy artifact so `./data/...` just works.
//
// Owner/repo are still derived from window.location when available — used only for the
// optional footer release-metadata fetch against api.github.com (which DOES send CORS
// headers). Falls back gracefully to "local" mode with no release footer.
function resolveDataBase() {
  const loc = window.location;
  const host = loc.hostname;
  if (host.endsWith('.github.io')) {
    const owner = host.split('.')[0];
    const repo = (loc.pathname.split('/').filter(Boolean)[0] ?? '').trim();
    if (owner && repo) {
      return {
        mode: 'pages',
        owner,
        repo,
        base: './data',
        apiBase: `https://api.github.com/repos/${owner}/${repo}`,
      };
    }
  }
  return { mode: 'local', base: './data' };
}

const DATA_SOURCE = resolveDataBase();

const SLOT_COLORS = {
  host_defense: '#8264c9',
  host_surface: '#2aa198',
  host_typing: '#268bd2',
  host_stats: '#b58900',
  phage_projection: '#d33682',
  phage_stats: '#cb4b16',
  phage_kmer: '#6c71c4',
  pair_depo_capsule: '#859900',
  pair_receptor_omp: '#dc322f',
  pair_concentration: '#073642',
  other: '#839496',
};

const SLOT_PREFIXES = Object.keys(SLOT_COLORS)
  .filter((slot) => slot !== 'other')
  .map((slot) => ({ slot, prefix: slot + '__' }));

function featureToSlot(feature) {
  for (const { slot, prefix } of SLOT_PREFIXES) {
    if (feature.startsWith(prefix)) return slot;
  }
  return 'other';
}

const RELIABILITY_VARIANT_ORDER = [
  { key: 'guelin_raw', label: 'Guelin raw', color: '#268bd2', dash: 'solid' },
  { key: 'guelin_loof_calibrated', label: 'Guelin LOOF-calibrated', color: '#268bd2', dash: 'dot' },
  { key: 'basel_raw', label: 'BASEL raw', color: '#dc322f', dash: 'solid' },
  { key: 'basel_guelin_calibrator_applied', label: 'BASEL (Guelin calibrator applied)', color: '#dc322f', dash: 'dot' },
];

function ui() {
  return {
    tabs: {
      overview: { label: 'Overview' },
      crossSource: { label: 'Cross-source' },
      calibration: { label: 'Calibration' },
      featureImportance: { label: 'Feature importance' },
      slots: { label: 'Per-slot breakdown' },
      predictions: { label: 'Predictions table' },
      pairExplorer: { label: 'Pair explorer' },
    },
    active: 'overview',
    loading: true,
    error: null,

    summary: null,
    crossSource: null,
    featureImportance: null,
    reliability: null,
    slotManifest: null,
    predictionsBact: null,
    predictionsPhage: null,

    reliabilityAxis: 'bacteria',
    fiAxis: 'bacteria_axis',
    slotAxis: 'bacteria_axis',
    predAxis: 'bacteria',
    predSource: 'all',
    predSort: 'absErrDesc',
    predRowsVisible: [],
    predPageSize: 500,

    dataSource: DATA_SOURCE,
    releaseMeta: null,

    shapAxis: 'bacteria_axis',
    shapBact: '',
    shapPhage: '',
    shapPair: null,
    shapLoading: false,
    shapError: null,
    duckdbConn: null,
    duckdbInitPromise: null,

    async init() {
      try {
        this.fetchReleaseMeta(); // fire-and-forget; footer updates when it resolves
        const [summary, crossSource, featureImportance, predBact, predPhage, reliability, slots] =
          await Promise.all([
            this.fetchJSON('ch05_summary.json'),
            this.fetchJSON('ch05_cross_source.json'),
            this.fetchJSON('ch05_feature_importance.json'),
            this.fetchJSON('ch05_predictions_bact.json'),
            this.fetchJSON('ch05_predictions_phage.json'),
            this.fetchJSON('ch05_reliability.json'),
            this.fetchJSON('slot_manifest.json'),
          ]);
        this.summary = summary;
        this.crossSource = crossSource;
        this.featureImportance = featureImportance;
        this.predictionsBact = predBact;
        this.predictionsPhage = predPhage;
        this.reliability = reliability;
        this.slotManifest = slots;
        this.loading = false;
        await this.$nextTick();
        this.renderAll();
      } catch (err) {
        console.error(err);
        this.error = `Failed to load snapshot data: ${err.message}. Run build_snapshot.py and serve this directory.`;
        this.loading = false;
      }
    },

    async fetchJSON(name) {
      const res = await fetch(`${this.dataSource.base}/${name}`);
      if (!res.ok) throw new Error(`${name}: HTTP ${res.status}`);
      return res.json();
    },

    async fetchReleaseMeta() {
      if (!this.dataSource.apiBase) return;
      const cacheKey = `release-meta:${this.dataSource.apiBase}`;
      const cached = sessionStorage.getItem(cacheKey);
      if (cached) {
        try { this.releaseMeta = JSON.parse(cached); return; } catch (_) { /* fall through */ }
      }
      try {
        const res = await fetch(`${this.dataSource.apiBase}/releases/latest`, { cache: 'no-cache' });
        if (!res.ok) return;
        const body = await res.json();
        this.releaseMeta = { tag: body.tag_name, published_at: body.published_at, html_url: body.html_url };
        sessionStorage.setItem(cacheKey, JSON.stringify(this.releaseMeta));
      } catch (err) {
        console.warn('Release metadata fetch failed (cosmetic footer only):', err);
      }
    },

    renderAll() {
      this.renderOverview();
      this.renderCrossSource();
      this.renderReliability();
      this.renderFeatureImportance();
      this.renderSlots();
      this.rebuildPredRows();
    },

    onTabClick(key) {
      this.active = key;
      if (key === 'pairExplorer' && !this.duckdbConn && !this.duckdbInitPromise) {
        // Kick off DuckDB-WASM init on first visit. Does not block the tab switch;
        // selectCustomPair() will await it.
        this.initDuckDB().catch((err) => {
          console.warn('DuckDB-WASM init failed', err);
          this.shapError = `DuckDB-WASM unavailable: ${err.message}`;
        });
      }
    },

    drillInto(row) {
      this.shapAxis = this.predAxis === 'bacteria' ? 'bacteria_axis' : 'phage_axis';
      this.shapBact = row.bacteria;
      this.shapPhage = row.phage;
      this.active = 'pairExplorer';
      this.selectCustomPair();
    },

    selectCustomPair() {
      if (!this.shapBact || !this.shapPhage) {
        this.shapError = 'Enter both a bacterium and a phage.';
        return;
      }
      this.shapError = null;
      this.loadShapForPair();
    },

    // ---- formatting helpers ----
    axisLabel(axis) { return axis === 'bacteria' ? 'Bacteria-axis' : 'Phage-axis'; },
    fmtCI(ciBlock) {
      if (!ciBlock) return '—';
      const { point, lo, hi } = ciBlock;
      if (point == null) return '—';
      const body = point.toFixed(4);
      if (lo == null || hi == null) return body;
      return `${body} <span class="ci">[${lo.toFixed(4)}, ${hi.toFixed(4)}]</span>`;
    },
    fmtInt(n) { return (n ?? 0).toLocaleString(); },

    // ---- Overview ----
    renderOverview() {
      if (!this.summary) return;
      const axes = ['bacteria', 'phage'];
      const aucTraces = [];
      const brierTraces = [];
      axes.forEach((axis) => {
        const a = this.summary.axes[axis].auc;
        aucTraces.push({
          name: this.axisLabel(axis),
          y: [this.axisLabel(axis)],
          x: [a.point],
          error_x: { type: 'data', symmetric: false, array: [a.hi - a.point], arrayminus: [a.point - a.lo] },
          type: 'scatter', mode: 'markers',
          marker: { size: 12, color: '#0b6cba' },
          xaxis: 'x1', yaxis: 'y1',
          showlegend: false,
          hovertemplate: `${this.axisLabel(axis)} AUC: %{x:.4f} [${a.lo.toFixed(4)}, ${a.hi.toFixed(4)}]<extra></extra>`,
        });
        const b = this.summary.axes[axis].brier;
        brierTraces.push({
          name: this.axisLabel(axis),
          y: [this.axisLabel(axis)],
          x: [b.point],
          error_x: { type: 'data', symmetric: false, array: [b.hi - b.point], arrayminus: [b.point - b.lo] },
          type: 'scatter', mode: 'markers',
          marker: { size: 12, color: '#cb4b16' },
          xaxis: 'x2', yaxis: 'y2',
          showlegend: false,
          hovertemplate: `${this.axisLabel(axis)} Brier: %{x:.4f} [${b.lo.toFixed(4)}, ${b.hi.toFixed(4)}]<extra></extra>`,
        });
      });
      Plotly.react(
        'overview-lollipop',
        [...aucTraces, ...brierTraces],
        {
          grid: { rows: 1, columns: 2, pattern: 'independent' },
          xaxis: { title: 'AUC', range: [0.7, 0.95], domain: [0, 0.48] },
          xaxis2: { title: 'Brier', range: [0.1, 0.22], domain: [0.52, 1.0] },
          yaxis: { automargin: true, anchor: 'x1' },
          yaxis2: { automargin: true, anchor: 'x2' },
          margin: { l: 120, r: 20, t: 10, b: 40 },
          height: 220,
        },
        { displayModeBar: false, responsive: true },
      );
    },

    // ---- Cross-source ----
    crossSourceRows() {
      if (!this.crossSource) return [];
      const rows = [];
      for (const src of (this.crossSource.bacteria ?? [])) rows.push({ axis: 'bacteria', ...src });
      for (const src of (this.crossSource.phage ?? [])) rows.push({ axis: 'phage', ...src });
      return rows;
    },
    renderCrossSource() {
      if (!this.crossSource) return;
      const axes = ['bacteria', 'phage'];
      const sources = ['guelin', 'basel'];
      const traces = [];
      sources.forEach((source) => {
        const xs = [];
        const ys = [];
        const errs = [];
        const hover = [];
        axes.forEach((axis) => {
          const entry = (this.crossSource[axis] ?? []).find((r) => r.source === source);
          if (!entry) return;
          xs.push(entry.auc.point);
          ys.push(this.axisLabel(axis));
          const hi = entry.auc.hi == null ? 0 : entry.auc.hi - entry.auc.point;
          const lo = entry.auc.lo == null ? 0 : entry.auc.point - entry.auc.lo;
          errs.push({ hi, lo });
          hover.push(`${source} ${axis}: ${entry.auc.point.toFixed(4)} (n=${entry.n_pairs})`);
        });
        traces.push({
          name: source,
          type: 'scatter',
          mode: 'markers',
          x: xs,
          y: ys,
          error_x: {
            type: 'data',
            symmetric: false,
            array: errs.map((e) => e.hi),
            arrayminus: errs.map((e) => e.lo),
          },
          marker: { size: 12, color: source === 'guelin' ? '#268bd2' : '#dc322f' },
          text: hover,
          hovertemplate: '%{text}<extra></extra>',
        });
      });
      Plotly.react(
        'cross-source-plot',
        traces,
        {
          xaxis: { title: 'AUC', range: [0.65, 0.98] },
          yaxis: { automargin: true },
          margin: { l: 140, r: 20, t: 10, b: 40 },
          height: 260,
          legend: { orientation: 'h', x: 0, y: -0.2 },
        },
        { displayModeBar: false, responsive: true },
      );
    },

    // ---- Calibration ----
    reliabilityVariants() {
      if (!this.reliability) return [];
      return RELIABILITY_VARIANT_ORDER.map((v) => {
        const key = `${this.reliabilityAxis}_axis_${v.key}`;
        const body = this.reliability.variants[key];
        return { key, ece: body?.ece ?? null };
      });
    },
    renderReliability() {
      if (!this.reliability) return;
      const traces = [];
      RELIABILITY_VARIANT_ORDER.forEach((v) => {
        const key = `${this.reliabilityAxis}_axis_${v.key}`;
        const body = this.reliability.variants[key];
        if (!body) return;
        traces.push({
          name: v.label,
          type: 'scatter',
          mode: 'lines+markers',
          x: body.bins.map((b) => b.mean_predicted),
          y: body.bins.map((b) => b.observed_rate),
          line: { color: v.color, dash: v.dash },
          marker: { color: v.color, size: 7 },
          hovertemplate: `${v.label}<br>bin [%{customdata[0]:.2f}, %{customdata[1]:.2f}]` +
            `<br>predicted %{x:.3f} → observed %{y:.3f}<extra></extra>`,
          customdata: body.bins.map((b) => [b.bin_lo, b.bin_hi]),
        });
      });
      traces.push({
        name: 'perfect',
        x: [0, 1], y: [0, 1],
        mode: 'lines',
        line: { color: '#aaa', dash: 'dash' },
        hoverinfo: 'skip',
        showlegend: true,
      });
      Plotly.react(
        'reliability-plot',
        traces,
        {
          xaxis: { title: 'mean predicted probability', range: [0, 1] },
          yaxis: { title: 'observed fraction positive', range: [0, 1] },
          margin: { l: 60, r: 20, t: 10, b: 50 },
          height: 420,
          legend: { orientation: 'h', x: 0, y: -0.15 },
        },
        { displayModeBar: false, responsive: true },
      );
    },

    // ---- Feature importance ----
    renderFeatureImportance() {
      if (!this.featureImportance) return;
      const rows = this.featureImportance[this.fiAxis].slice(0, 30);
      rows.reverse(); // Plotly bar charts grow bottom-up.
      const trace = {
        type: 'bar',
        orientation: 'h',
        x: rows.map((r) => r.mean_importance),
        y: rows.map((r) => r.feature),
        marker: {
          color: rows.map((r) => SLOT_COLORS[r.slot] ?? SLOT_COLORS.other),
        },
        customdata: rows.map((r) => [r.slot, r.n_folds_selected, r.is_concentration_feature]),
        hovertemplate: '%{y}<br>slot: %{customdata[0]}<br>importance: %{x:.1f}' +
          '<br>selected in %{customdata[1]}/30 fits<extra></extra>',
      };
      Plotly.react(
        'fi-plot',
        [trace],
        {
          margin: { l: 320, r: 40, t: 10, b: 40 },
          xaxis: { title: 'mean LightGBM importance' },
          yaxis: { automargin: true, tickfont: { size: 11 } },
          height: Math.max(640, rows.length * 22),
          showlegend: false,
        },
        { displayModeBar: false, responsive: true },
      );
    },

    // ---- Per-slot breakdown ----
    renderSlots() {
      if (!this.slotManifest) return;
      const container = document.getElementById('slot-cards');
      if (!container) return;
      const rows = [...this.slotManifest.slots];
      rows.sort((a, b) => (b.axes?.[this.slotAxis]?.cumulative_importance ?? 0) -
                         (a.axes?.[this.slotAxis]?.cumulative_importance ?? 0));
      container.innerHTML = '';
      for (const slot of rows) {
        const axisInfo = slot.axes?.[this.slotAxis] ?? { feature_count: 0, cumulative_importance: 0 };
        if (axisInfo.feature_count === 0) continue;
        const card = document.createElement('div');
        card.className = 'slot-card';
        const color = SLOT_COLORS[slot.slot_name] ?? SLOT_COLORS.other;
        card.innerHTML = `
          <h3 style="color: ${color}">${slot.slot_name}
            <span class="role-pill">${slot.block_role}</span>
          </h3>
          <div class="description">${slot.description}</div>
          <div class="stat-row">
            <span>feature count</span><span>${axisInfo.feature_count}</span>
            <span>cumulative importance</span><span>${axisInfo.cumulative_importance.toFixed(0)}</span>
            <span>entity key</span><span><code>${slot.entity_key}</code></span>
          </div>
        `;
        container.appendChild(card);
      }
    },

    // ---- Pair explorer (DuckDB-WASM + SHAP waterfall) ----

    async initDuckDB() {
      if (this.duckdbInitPromise) return this.duckdbInitPromise;
      this.shapLoading = true;
      this.duckdbInitPromise = (async () => {
        const duckdb = await import('https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm');
        const bundles = duckdb.getJsDelivrBundles();
        const bundle = await duckdb.selectBundle(bundles);
        const workerUrl = URL.createObjectURL(
          new Blob([`importScripts("${bundle.mainWorker}");`], { type: 'text/javascript' }),
        );
        const worker = new Worker(workerUrl);
        const logger = new duckdb.ConsoleLogger(duckdb.LogLevel.WARNING);
        const db = new duckdb.AsyncDuckDB(logger, worker);
        await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
        URL.revokeObjectURL(workerUrl);
        // Register the six Parquet assets as virtual files. The DATA_SOURCE.base URL is
        // either the release-download root (Pages) or ./data (local preview).
        const parquetNames = [
          'bacteria_axis_shap_values.parquet',
          'bacteria_axis_shap_base_values.parquet',
          'bacteria_axis_feature_values.parquet',
          'phage_axis_shap_values.parquet',
          'phage_axis_shap_base_values.parquet',
          'phage_axis_feature_values.parquet',
        ];
        for (const name of parquetNames) {
          await db.registerFileURL(
            name,
            `${this.dataSource.base}/${name}`,
            4, // duckdb.DuckDBDataProtocol.HTTP
            false,
          );
        }
        this.duckdbConn = await db.connect();
        this.shapLoading = false;
      })();
      return this.duckdbInitPromise;
    },

    async loadShapForPair() {
      try {
        this.shapError = null;
        this.shapLoading = true;
        await this.initDuckDB();
        const conn = this.duckdbConn;
        if (!conn) throw new Error('DuckDB connection not initialised');
        const axis = this.shapAxis; // 'bacteria_axis' | 'phage_axis'
        const bact = this.shapBact.replace(/'/g, "''");
        const phage = this.shapPhage.replace(/'/g, "''");
        const shapQ = `SELECT * FROM '${axis}_shap_values.parquet' ` +
          `WHERE bacteria = '${bact}' AND phage = '${phage}' LIMIT 1`;
        const valuesQ = `SELECT * FROM '${axis}_feature_values.parquet' ` +
          `WHERE bacteria = '${bact}' AND phage = '${phage}' LIMIT 1`;
        const baseQ = `SELECT base_value FROM '${axis}_shap_base_values.parquet' ` +
          `WHERE pair_id = '${bact}__${phage}' LIMIT 1`;
        const [shapRes, valuesRes, baseRes] = await Promise.all([
          conn.query(shapQ),
          conn.query(valuesQ),
          conn.query(baseQ),
        ]);
        const shapRow = shapRes.numRows > 0 ? shapRes.get(0) : null;
        const valuesRow = valuesRes.numRows > 0 ? valuesRes.get(0) : null;
        const baseRow = baseRes.numRows > 0 ? baseRes.get(0) : null;
        if (!shapRow || !valuesRow || !baseRow) {
          throw new Error(`Pair not found in ${axis}: ${this.shapBact} × ${this.shapPhage}`);
        }
        // Build the top-20-contribution list from the SHAP row. Arrow rows are dict-like
        // but we convert to plain object for Alpine reactivity.
        const shapObj = Object.fromEntries(
          shapRow.toArray().map((v, i) => [shapRes.schema.fields[i].name, v]),
        );
        const valuesObj = Object.fromEntries(
          valuesRow.toArray().map((v, i) => [valuesRes.schema.fields[i].name, v]),
        );
        const baseValue = Number(baseRow.toArray()[0]);

        const shapFeatures = [];
        for (const [col, value] of Object.entries(shapObj)) {
          if (!col.startsWith('shap__')) continue;
          const featureName = col.slice('shap__'.length);
          const numericShap = Number(value);
          if (!Number.isFinite(numericShap) || numericShap === 0) continue;
          shapFeatures.push({
            feature: featureName,
            slot: featureToSlot(featureName),
            shap: numericShap,
            raw: valuesObj[featureName] ?? null,
          });
        }
        shapFeatures.sort((a, b) => Math.abs(b.shap) - Math.abs(a.shap));
        const top = shapFeatures.slice(0, 20);
        const shapSum = shapFeatures.reduce((acc, f) => acc + f.shap, 0);
        const logit = baseValue + shapSum;
        const labelValue = valuesObj.label ?? valuesObj.label_row_binary ?? null;

        this.shapPair = {
          bacteria: this.shapBact,
          phage: this.shapPhage,
          label: labelValue,
          base: baseValue,
          shapSum,
          logit,
          topFeatures: top,
        };
        this.shapLoading = false;
        await this.$nextTick();
        this.renderWaterfall(top);
      } catch (err) {
        console.error(err);
        this.shapError = err.message;
        this.shapLoading = false;
        this.shapPair = null;
      }
    },

    renderWaterfall(top) {
      if (!top || !top.length) return;
      const reversed = [...top].reverse(); // Plotly stacks bottom-up
      const trace = {
        type: 'bar',
        orientation: 'h',
        x: reversed.map((r) => r.shap),
        y: reversed.map((r) => r.feature),
        marker: {
          color: reversed.map((r) => (r.shap > 0 ? '#0b6cba' : '#dc322f')),
        },
        customdata: reversed.map((r) => [r.slot, r.raw]),
        hovertemplate: '%{y}<br>slot: %{customdata[0]}<br>raw value: %{customdata[1]}' +
          '<br>SHAP: %{x:.4f}<extra></extra>',
      };
      Plotly.react(
        'waterfall-plot',
        [trace],
        {
          xaxis: { title: 'SHAP contribution (+/− log-odds)' },
          yaxis: { automargin: true, tickfont: { size: 11 } },
          margin: { l: 320, r: 40, t: 10, b: 40 },
          height: Math.max(480, reversed.length * 24),
          showlegend: false,
        },
        { displayModeBar: false, responsive: true },
      );
    },

    // ---- Predictions table ----
    rebuildPredRows() {
      const src = this.predAxis === 'bacteria' ? this.predictionsBact : this.predictionsPhage;
      if (!src) { this.predRowsVisible = []; return; }
      let rows = src;
      if (this.predSource !== 'all') rows = rows.filter((r) => r.source === this.predSource);
      rows = rows.map((r) => ({ ...r, absErr: Math.abs(r.predicted - r.label) }));
      const sorts = {
        absErrDesc: (a, b) => b.absErr - a.absErr,
        predictedDesc: (a, b) => b.predicted - a.predicted,
        predictedAsc: (a, b) => a.predicted - b.predicted,
        bacteria: (a, b) => a.bacteria.localeCompare(b.bacteria),
        phage: (a, b) => a.phage.localeCompare(b.phage),
      };
      rows.sort(sorts[this.predSort] ?? sorts.absErrDesc);
      this.predRowsVisible = rows;
    },
  };
}
