const state = {
  bulkPreviewRows: [],
  jobPoller: null,
  currentJobId: null,
  snapshotNames: [],
  boltzConfigs: [],
  epitopePlots: [],
  diversityPlots: [],
  lastJobStatus: null,
  clusterStatus: null,
  boltzLocalRuns: {},
  boltzLocalRunsBySpec: {},
  expandedTargets: new Set(),
  diversityCsv: null,
  diversityHtml: null,
  diversityMessage: null,
  diversityFiles: [],
  diversityOutputDir: null,
  epitopeCsvName: null,
  epitopeHotspotCsvName: null,
  binderRows: [],
  binderTotal: 0,
  binderPage: 1,
  binderPageSize: 100,
  binderCsvName: null,
  binderMessage: null,
  binderPdbIds: [],
  boltzDesignCounts: {},
  bulkPreviewTimer: null,
  bulkPreviewSig: null,
  binderFilterTimer: null,
  rfaConfigs: {},
  configContext: null,
  runContext: null,
  runMode: null,
};

const PIPELINE_RERUN_DELAY_MS = 3000; // 3 seconds

// Keep in sync with scripts/pymol_utils.py base_palette.
const EPITOPE_COLORS = [
  [0.90, 0.40, 0.00],
  [0.00, 0.55, 0.60],
  [0.60, 0.20, 0.75],
  [0.10, 0.60, 0.10],
  [0.20, 0.40, 0.85],
  [0.80, 0.00, 0.20],
  [0.75, 0.55, 0.10],
  [0.00, 0.70, 0.95],
  [0.95, 0.55, 0.80],
  [0.30, 0.30, 0.30],
];

function epitopeColor(label) {
  const match = String(label || '').match(/epitope[_-]?(\d+)/i);
  if (!match) return null;
  const idx = Number(match[1]);
  if (!Number.isFinite(idx) || idx <= 0) return null;
  const base = EPITOPE_COLORS[(idx - 1) % EPITOPE_COLORS.length];
  const rgb = base.map((val) => Math.round(val * 255));
  return `rgb(${rgb.join(',')})`;
}

function formatEngineLabel(engine) {
  const raw = String(engine || '').trim().toLowerCase();
  if (!raw) return '—';
  if (raw.includes('boltz')) return 'BoltzGen';
  if (raw.includes('rfa')) return 'RFAntibody';
  return engine;
}

function bulkFileSrc(value) {
  if (!value) return null;
  const text = String(value);
  const baseName = text.split('/').pop();
  const name = baseName || text;
  return `/api/bulk/file?name=${encodeURIComponent(name)}&t=${Date.now()}`;
}

function isBoltzEngine(engine) {
  const raw = String(engine || '').trim().toLowerCase();
  return !raw || raw.includes('boltz');
}

const el = {
  bulkStatus: document.querySelector('#bulk-status'),
  bulkAlgorithmBtn: document.querySelector('#bulk-algorithm-btn'),
  bulkCsvInput: document.querySelector('#bulk-csv-input'),
  bulkEpitopes: document.querySelector('#bulk-epitopes'),
  bulkEpitopePrompt: document.querySelector('#bulk-epitope-prompt'),
  bulkLaunchPymol: document.querySelector('#bulk-launch-pymol'),
  bulkExportInsights: document.querySelector('#bulk-export-insights'),
  bulkExportDesigns: document.querySelector('#bulk-export-designs'),
  bulkForceInit: document.querySelector('#bulk-force-init'),
  designEngine: document.querySelector('#design-engine'),
  designTotal: document.querySelector('#design-total'),
  bulkRunPrefix: document.querySelector('#bulk-run-prefix'),
  bulkThrottle: document.querySelector('#bulk-throttle'),
  bulkPreviewBtn: document.querySelector('#bulk-preview'),
  bulkRunBtn: document.querySelector('#bulk-run'),
  bulkPreviewPanel: document.querySelector('#bulk-preview-panel'),
  bulkPreviewTable: document.querySelector('#bulk-preview-table tbody'),
  bulkPreviewSummary: document.querySelector('#bulk-preview-summary'),
  bulkPreviewRefresh: document.querySelector('#bulk-preview-refresh'),
  bulkVisualizeEpitopes: document.querySelector('#bulk-visualize-epitopes'),
  jobAlert: document.querySelector('#job-alert'),
  jobLog: document.querySelector('#job-log'),
  jobMeta: document.querySelector('#job-meta'),
  snapshotSection: document.querySelector('#snapshot-section'),
  snapshotGrid: document.querySelector('#snapshot-grid'),
  snapshotDownloadPdf: document.querySelector('#snapshot-download-pdf'),
  epitopePlotSection: document.querySelector('#epitope-plot-section'),
  epitopePlotGrid: document.querySelector('#epitope-plot-grid'),
  epitopeReportDownload: document.querySelector('#epitope-report-download'),
  epitopeCsvDownload: document.querySelector('#epitope-csv-download'),
  epitopeHotspotCsvDownload: document.querySelector('#epitope-hotspot-csv-download'),
  epitopeMetricsSection: document.querySelector('#epitope-metrics-section'),
  diversitySection: document.querySelector('#diversity-section'),
  diversityGrid: document.querySelector('#diversity-grid'),
  diversityRefresh: document.querySelector('#diversity-refresh'),
  diversityRebuild: document.querySelector('#diversity-rebuild'),
  diversityDownloadCsv: document.querySelector('#diversity-download-csv'),
  diversityDownloadHtml: document.querySelector('#diversity-download-html'),
  boltzPanel: document.querySelector('#boltz-config-panel'),
  boltzTable: document.querySelector('#boltz-config-table tbody'),
  boltzSummary: document.querySelector('#boltz-config-summary'),
  boltzDesignCount: document.querySelector('#boltz-design-count'),
  boltzTimeHours: document.querySelector('#boltz-time-hours'),
  boltzCropRadius: document.querySelector('#boltz-crop-radius'),
  boltzRegenerate: document.querySelector('#boltz-config-regenerate'),
  boltzRefresh: document.querySelector('#boltz-config-refresh'),
  boltzPlotDiversity: document.querySelector('#boltz-plot-diversity'),
  epitopeDiversitySelection: document.querySelector('#epitope-diversity-selection'),
  epitopeDiversityPlot: document.querySelector('#epitope-diversity-plot'),
  boltzShowRunSelection: document.querySelector('#boltz-show-run-selection'),
  boltzShowRunAll: document.querySelector('#boltz-show-run-all'),
  boltzRerunRange: document.querySelector('#boltz-rerun-range'),
  boltzShowRunRange: document.querySelector('#boltz-show-run-range'),
  boltzRerunRangeModal: document.querySelector('#boltz-rerun-range-modal'),
  boltzRerunRangeStart: document.querySelector('#boltz-rerun-range-start'),
  boltzRerunRangeEnd: document.querySelector('#boltz-rerun-range-end'),
  boltzRerunRangeExpected: document.querySelector('#boltz-rerun-range-expected'),
  boltzRerunRangeAttempts: document.querySelector('#boltz-rerun-range-attempts'),
  boltzRerunRangeForce: document.querySelector('#boltz-rerun-range-force'),
  boltzRerunRangePrompt: document.querySelector('#boltz-rerun-range-prompt'),
  boltzRerunRangeConfirm: document.querySelector('#boltz-rerun-range-confirm'),
  boltzRerunRangeCancel: document.querySelector('#boltz-rerun-range-cancel'),
  boltzRerunRangeClose: document.querySelector('#boltz-rerun-range-close'),
  bulkAlgorithmModal: document.querySelector('#bulk-algorithm-modal'),
  bulkAlgorithmClose: document.querySelector('#bulk-algorithm-close'),
  boltzRunRangeModal: document.querySelector('#boltz-run-range-modal'),
  boltzRunRangeStart: document.querySelector('#boltz-run-range-start'),
  boltzRunRangeEnd: document.querySelector('#boltz-run-range-end'),
  boltzRunRangeMinAllowed: document.querySelector('#boltz-run-range-min-allowed'),
  boltzRunRangeMinEpitopes: document.querySelector('#boltz-run-range-min-epitopes'),
  boltzRunRangeFilterNote: document.querySelector('#boltz-run-range-filter-note'),
  boltzRunRangeAllowedNote: document.querySelector('#boltz-run-range-allowed-note'),
  boltzRunRangeCountNote: document.querySelector('#boltz-run-range-count-note'),
  boltzRunRangeConfirm: document.querySelector('#boltz-run-range-confirm'),
  boltzRunRangeCancel: document.querySelector('#boltz-run-range-cancel'),
  boltzRunRangeClose: document.querySelector('#boltz-run-range-close'),
  boltzRegenerateRangeModal: document.querySelector('#boltz-regenerate-range-modal'),
  boltzRegenerateRangeStart: document.querySelector('#boltz-regenerate-range-start'),
  boltzRegenerateRangeEnd: document.querySelector('#boltz-regenerate-range-end'),
  boltzRegenerateRangeConfirm: document.querySelector('#boltz-regenerate-range-confirm'),
  boltzRegenerateRangeCancel: document.querySelector('#boltz-regenerate-range-cancel'),
  boltzRegenerateRangeClose: document.querySelector('#boltz-regenerate-range-close'),
  binderPanel: document.querySelector('#boltz-binders-panel'),
  binderTable: document.querySelector('#boltz-binders-table tbody'),
  binderSummary: document.querySelector('#boltz-binders-summary'),
  binderCsvNote: document.querySelector('#boltz-binders-csv'),
  binderRefresh: document.querySelector('#boltz-binders-refresh'),
  binderExportBtn: document.querySelector('#boltz-binders-export'),
  binderDownload: document.querySelector('#boltz-binders-download'),
  binderPagination: document.querySelector('#boltz-binders-pagination'),
  binderPageLabel: document.querySelector('#boltz-binders-page-label'),
  binderFilterPdb: document.querySelector('#binder-filter-pdb'),
  binderFilterEpitope: document.querySelector('#binder-filter-epitope'),
  binderFilterEngine: document.querySelector('#binder-filter-engine'),
  binderOrderBy: document.querySelector('#binder-order-by'),
  binderExportModal: document.querySelector('#binder-export-modal'),
  binderExportClose: document.querySelector('#binder-export-close'),
  binderExportCancel: document.querySelector('#binder-export-cancel'),
  binderExportConfirm: document.querySelector('#binder-export-confirm'),
  binderExportSelections: document.querySelector('#binder-export-selections'),
  binderExportCount: document.querySelector('#binder-export-count'),
  binderExportSummary: document.querySelector('#binder-export-summary'),
  boltzConfigModal: document.querySelector('#boltz-config-modal'),
  boltzConfigTitle: document.querySelector('#boltz-config-title'),
  boltzConfigBody: document.querySelector('#boltz-config-body'),
  boltzConfigEngine: document.querySelector('#bulk-config-engine'),
  boltzConfigClose: document.querySelector('#boltz-config-close'),
  boltzLogModal: document.querySelector('#boltz-log-modal'),
  boltzLogTitle: document.querySelector('#boltz-log-title'),
  boltzLogBody: document.querySelector('#boltz-log-body'),
  boltzLogClose: document.querySelector('#boltz-log-close'),
  boltzRunModal: document.querySelector('#boltz-run-modal'),
  boltzRunTitle: document.querySelector('#boltz-run-title'),
  boltzRunBody: document.querySelector('#boltz-run-body'),
  boltzRunEngine: document.querySelector('#bulk-run-engine'),
  boltzRunClose: document.querySelector('#boltz-run-close'),
  pipelineRerunModal: document.querySelector('#pipeline-rerun-modal'),
  pipelineRerunTitle: document.querySelector('#pipeline-rerun-title'),
  pipelineRerunExpected: document.querySelector('#pipeline-rerun-expected'),
  pipelineRerunAttempts: document.querySelector('#pipeline-rerun-attempts'),
  pipelineRerunForce: document.querySelector('#pipeline-rerun-force'),
  pipelineRerunPrompt: document.querySelector('#pipeline-rerun-prompt'),
  pipelineRerunConfirm: document.querySelector('#pipeline-rerun-confirm'),
  pipelineRerunCancel: document.querySelector('#pipeline-rerun-cancel'),
  pipelineRerunClose: document.querySelector('#pipeline-rerun-close'),
};

function formatRangeLabel(value) {
  if (Array.isArray(value) && value.length) {
    const start = value[0];
    const end = value.length > 1 ? value[1] : value[0];
    if (start !== undefined && end !== undefined) return `${start}-${end}`;
  }
  if (value && typeof value === 'object' && 'start' in value && 'end' in value) {
    return `${value.start}-${value.end}`;
  }
  if (value === 0) return '0';
  if (value === null || value === undefined) return null;
  const text = String(value).trim();
  return text || null;
}

function formatPercent(value, decimals = 1) {
  const num = Number(value);
  if (!Number.isFinite(num)) return null;
  return `${(num * 100).toFixed(decimals)}%`;
}

function formatNumberValue(value, decimals = 1, fallback = '—') {
  const num = Number(value);
  if (!Number.isFinite(num)) return fallback;
  return num.toFixed(decimals);
}

function asNumber(value) {
  const num = Number(value);
  return Number.isFinite(num) ? num : null;
}

function describeSeries(values = []) {
  const arr = values.filter((v) => Number.isFinite(v));
  if (!arr.length) return null;
  const sorted = [...arr].sort((a, b) => a - b);
  const count = sorted.length;
  const sum = sorted.reduce((acc, val) => acc + val, 0);
  const mid = Math.floor(count / 2);
  const median = count % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
  return {
    count,
    min: sorted[0],
    max: sorted[count - 1],
    mean: sum / count,
    median,
  };
}

async function copyTextToClipboard(text, btn = null) {
  if (!text) return;
  try {
    if (navigator?.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
    } else {
      const ta = document.createElement('textarea');
      ta.value = text;
      ta.style.position = 'fixed';
      ta.style.top = '-1000px';
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      ta.remove();
    }
    if (btn) {
      const prev = btn.textContent;
      btn.textContent = 'Copied';
      setTimeout(() => { btn.textContent = prev; }, 1200);
    }
  } catch (err) {
    showAlert('Copy failed. Try selecting the text manually.');
  }
}

function renderRunCommandBlocks(text, options = {}) {
  if (!el.boltzRunBody) return;
  el.boltzRunBody.innerHTML = '';
  const singleBlock = Boolean(options.singleBlock);
  if (!text) {
    el.boltzRunBody.textContent = 'No commands available.';
    return;
  }
  if (singleBlock) {
    const wrapper = document.createElement('div');
    wrapper.style.marginBottom = '12px';

    const header = document.createElement('div');
    header.style.display = 'flex';
    header.style.justifyContent = 'space-between';
    header.style.alignItems = 'center';
    header.style.gap = '8px';
    header.style.marginBottom = '6px';

    const heading = document.createElement('div');
    heading.style.fontWeight = '600';
    heading.textContent = 'Commands';
    header.appendChild(heading);

    const copyBtn = document.createElement('button');
    copyBtn.type = 'button';
    copyBtn.className = 'ghost';
    copyBtn.textContent = 'Copy';
    copyBtn.addEventListener('click', () => copyTextToClipboard(String(text), copyBtn));
    header.appendChild(copyBtn);

    const block = document.createElement('pre');
    block.textContent = String(text).trim() || '—';
    block.style.background = '#f1f5f9';
    block.style.border = '1px solid #e2e8f0';
    block.style.borderRadius = '8px';
    block.style.padding = '10px';
    block.style.margin = '0';
    block.style.whiteSpace = 'pre-wrap';

    wrapper.appendChild(header);
    wrapper.appendChild(block);
    el.boltzRunBody.appendChild(wrapper);
    return;
  }
  const lines = String(text).split('\n');
  const sections = [];
  let current = { title: null, lines: [] };
  lines.forEach((line) => {
    const trimmed = line.trim();
    if (trimmed.startsWith('#')) {
      if (current.title || current.lines.length) sections.push(current);
      current = { title: trimmed.replace(/^#+\s*/, ''), lines: [] };
    } else {
      current.lines.push(line);
    }
  });
  if (current.title || current.lines.length) sections.push(current);
  if (!sections.length) {
    sections.push({ title: 'Commands', lines });
  }

  sections.forEach((section, idx) => {
    const title = section.title || `Step ${idx + 1}`;
    const commandText = section.lines.join('\n').trim();
    const wrapper = document.createElement('div');
    wrapper.style.marginBottom = '12px';

    const header = document.createElement('div');
    header.style.display = 'flex';
    header.style.justifyContent = 'space-between';
    header.style.alignItems = 'center';
    header.style.gap = '8px';
    header.style.marginBottom = '6px';

    const heading = document.createElement('div');
    heading.style.fontWeight = '600';
    heading.textContent = title;
    header.appendChild(heading);

    const copyBtn = document.createElement('button');
    copyBtn.type = 'button';
    copyBtn.className = 'ghost';
    copyBtn.textContent = 'Copy';
    copyBtn.addEventListener('click', () => copyTextToClipboard(commandText, copyBtn));
    header.appendChild(copyBtn);

    const block = document.createElement('pre');
    block.textContent = commandText || '—';
    block.style.background = '#f1f5f9';
    block.style.border = '1px solid #e2e8f0';
    block.style.borderRadius = '8px';
    block.style.padding = '10px';
    block.style.margin = '0';
    block.style.whiteSpace = 'pre-wrap';

    wrapper.appendChild(header);
    wrapper.appendChild(block);
    el.boltzRunBody.appendChild(wrapper);
  });
}

function bulkPreviewSignature(text) {
  const trimmed = String(text || '').trim();
  if (!trimmed) return null;
  const head = trimmed.slice(0, 96);
  const tail = trimmed.slice(-96);
  return `${trimmed.length}:${head}:${tail}`;
}

function scheduleBulkPreview() {
  if (!el.bulkCsvInput) return;
  if (state.bulkPreviewTimer) clearTimeout(state.bulkPreviewTimer);
  state.bulkPreviewTimer = setTimeout(() => {
    const csvText = el.bulkCsvInput?.value || '';
    const sig = bulkPreviewSignature(csvText);
    if (!sig) {
      state.bulkPreviewSig = null;
      state.bulkPreviewRows = [];
      renderBulkPreview([]);
      state.boltzConfigs = [];
      renderBoltzConfigs();
      return;
    }
    if (sig === state.bulkPreviewSig) return;
    state.bulkPreviewSig = sig;
    previewBulkCsv({ silent: true });
  }, 650);
}

function getBinderFilters() {
  return {
    pdb: (el.binderFilterPdb?.value || '').trim(),
    epitope: (el.binderFilterEpitope?.value || '').trim(),
    engine: (el.binderFilterEngine?.value || '').trim(),
    orderBy: (el.binderOrderBy?.value || '').trim(),
  };
}

function scheduleBinderRefresh({ silent = true } = {}) {
  if (state.binderFilterTimer) clearTimeout(state.binderFilterTimer);
  state.binderFilterTimer = setTimeout(() => {
    state.binderPage = 1;
    refreshDiversity({ silent, page: 1 });
  }, 300);
}

function renderHistogram(container, values = [], { bins = 8, label = 'values' } = {}) {
  if (!container) return;
  container.innerHTML = '';
  container.classList.remove('empty');
  const list = (values || []).filter((v) => Number.isFinite(v));
  if (!list.length) {
    container.classList.add('empty');
    container.textContent = 'No data yet.';
    return;
  }
  const min = Math.min(...list);
  const max = Math.max(...list);
  const span = max - min || 1;
  const bucketCount = Math.max(3, Math.min(bins, list.length));
  const binWidth = span / bucketCount || 1;
  const counts = Array.from({ length: bucketCount }, () => 0);
  list.forEach((val) => {
    const idx = Math.min(bucketCount - 1, Math.floor((val - min) / binWidth));
    counts[idx] += 1;
  });
  const peak = Math.max(...counts) || 1;
  counts.forEach((count, idx) => {
    const bar = document.createElement('div');
    bar.className = 'mini-bar';
    const pct = Math.max(6, Math.round((count / peak) * 100));
    bar.style.height = `${pct}%`;
    bar.dataset.count = count;
    const start = min + idx * binWidth;
    const end = idx === bucketCount - 1 ? max : start + binWidth;
    const lo = Math.round(start);
    const hi = Math.round(end);
    bar.title = `${count} ${label}: ${lo}-${hi}`;
    container.appendChild(bar);
  });
}

function renderStatPills(container, rows = []) {
  if (!container) return;
  container.innerHTML = '';
  rows.filter((row) => row && row.label && row.value).forEach((row) => {
    const pill = document.createElement('div');
    pill.className = 'metric-pill';

    const head = document.createElement('div');
    head.className = 'pill-head';
    const label = document.createElement('div');
    label.className = 'pill-label';
    label.textContent = row.label;
    const value = document.createElement('div');
    value.className = 'pill-value';
    value.textContent = row.value;
    head.appendChild(label);
    head.appendChild(value);
    pill.appendChild(head);

    if (row.hint) {
      const hint = document.createElement('div');
      hint.className = 'pill-hint';
      hint.textContent = row.hint;
      pill.appendChild(hint);
    }
    if (row.extra instanceof HTMLElement) {
      pill.appendChild(row.extra);
    }

    container.appendChild(pill);
  });
}

function buildSplitBar(fraction) {
  const frac = Math.min(1, Math.max(0, Number(fraction) || 0));
  const bar = document.createElement('div');
  bar.className = 'metric-split';
  const a = document.createElement('div');
  a.className = 'hydro-a';
  a.style.width = `${(frac * 100).toFixed(1)}%`;
  const b = document.createElement('div');
  b.className = 'hydro-b';
  b.style.width = `${(100 - frac * 100).toFixed(1)}%`;
  bar.appendChild(a);
  bar.appendChild(b);
  return bar;
}

function setBadge(badge, text, color = null) {
  if (!badge) return;
  if (!text) {
    badge.hidden = true;
    return;
  }
  badge.hidden = false;
  badge.textContent = text;
  badge.style.background = color || 'rgba(37, 99, 235, 0.12)';
  badge.style.color = color ? '#0f172a' : 'var(--color-accent)';
}

function stylePymolButton(btn, hasPrep) {
  if (!btn) return;
  if (hasPrep) {
    btn.disabled = false;
    btn.style.background = '';
    btn.style.borderColor = '';
    btn.style.color = '';
    btn.title = '';
    return;
  }
  btn.disabled = true;
  btn.style.background = 'rgba(248, 113, 113, 0.15)';
  btn.style.borderColor = 'rgba(248, 113, 113, 0.35)';
  btn.style.color = '#b91c1c';
  btn.title = 'Hotspot bundle missing; re-run pipeline to enable PyMOL.';
}

function normalizeStatus(status) {
  const value = (status || '').toLowerCase();
  if (value === 'running' || value === 'pending') return { text: value === 'pending' ? 'Queued' : 'Running', tone: 'info' };
  if (value === 'success') return { text: 'Complete', tone: 'success' };
  if (value === 'failed' || value === 'canceled') return { text: 'Failed', tone: 'danger' };
  return { text: 'No result', tone: 'muted' };
}

function buildStatusBadge(status) {
  const meta = normalizeStatus(status);
  const pill = document.createElement('span');
  pill.className = 'badge';
  pill.textContent = meta.text;
  if (meta.tone === 'success') {
    pill.style.background = 'rgba(134, 239, 172, 0.25)';
    pill.style.color = '#14532d';
  } else if (meta.tone === 'danger') {
    pill.style.background = 'rgba(248, 113, 113, 0.25)';
    pill.style.color = '#7f1d1d';
  } else if (meta.tone === 'info') {
    pill.style.background = 'rgba(191, 219, 254, 0.4)';
    pill.style.color = '#1d4ed8';
  } else {
    pill.style.background = 'rgba(148, 163, 184, 0.25)';
    pill.style.color = '#475569';
  }
  return pill;
}

function buildLocalResultBadge(count = null) {
  const pill = document.createElement('span');
  pill.className = 'badge';
  const num = Number(count);
  const hasCount = Number.isFinite(num);
  if (!hasCount) {
    pill.textContent = '—';
    pill.style.background = 'rgba(148, 163, 184, 0.25)';
    pill.style.color = '#475569';
    return pill;
  }
  pill.textContent = `${num} design${num === 1 ? '' : 's'}`;
  if (num > 0) {
    pill.style.background = 'rgba(134, 239, 172, 0.25)';
    pill.style.color = '#14532d';
  } else {
    pill.style.background = 'rgba(148, 163, 184, 0.25)';
    pill.style.color = '#475569';
  }
  return pill;
}

function normalizeBoltzCounts(counts = null) {
  const output = {};
  if (!counts || typeof counts !== 'object') return output;
  Object.entries(counts).forEach(([key, value]) => {
    const upper = (key || '').trim().toUpperCase();
    if (!upper) return;
    const num = Number(value);
    output[upper] = Number.isFinite(num) ? num : 0;
  });
  return output;
}

function getBoltzDesignCount() {
  const raw = Number(el.boltzDesignCount?.value || 0);
  if (!Number.isFinite(raw) || raw <= 0) return null;
  return Math.max(1, Math.round(raw));
}

function getBoltzTimeHours() {
  const raw = Number(el.boltzTimeHours?.value || 0);
  if (!Number.isFinite(raw) || raw <= 0) return null;
  return Math.min(240, Math.max(1, Math.round(raw)));
}

function getBoltzCropRadius() {
  const raw = Number(el.boltzCropRadius?.value || 0);
  if (!Number.isFinite(raw) || raw <= 0) return null;
  return raw;
}

function toggleModal(modalEl, isOpen) {
  if (!modalEl) return;
  modalEl.hidden = !isOpen;
  if (isOpen) {
    document.body.classList.add('modal-open');
  } else {
    const openModals = document.querySelectorAll('.modal:not([hidden])');
    if (!openModals.length) {
      document.body.classList.remove('modal-open');
    }
  }
}

function epitopeLabel(cfg = {}, fallbackIndex = null) {
  if (cfg.epitope_name) return cfg.epitope_name;
  if (cfg.epitope_id) return cfg.epitope_id;
  if (fallbackIndex !== null) return `Epitope ${fallbackIndex}`;
  return 'Epitope';
}

function normalizeEpitopeToken(value) {
  const raw = String(value || '').trim();
  if (!raw) return '';
  const match = raw.match(/epitope[_\-\s]*(\d+)/i);
  if (match) return `epitope_${Number(match[1])}`;
  if (/^\d+$/.test(raw)) return `epitope_${Number(raw)}`;
  return raw.toLowerCase().replace(/\s+/g, '');
}

function isSurfaceAllowed(cfg) {
  return cfg?.hotspot_surface_ok !== false;
}

function surfaceEpitopeLabel(cfg, fallbackIndex = null) {
  const base = epitopeLabel(cfg, fallbackIndex);
  if (cfg?.hotspot_surface_ok === false) {
    return `${base} (SASA filtered)`;
  }
  return base;
}

function availableEpitopeLabels(configs = []) {
  const labels = [];
  const seen = new Set();
  configs.forEach((cfg, idx) => {
    if (!isSurfaceAllowed(cfg)) return;
    const label = epitopeLabel(cfg, idx + 1);
    const key = String(label || '').trim().toLowerCase();
    if (!key || seen.has(key)) return;
    seen.add(key);
    labels.push(label);
  });
  return labels;
}

function epitopeKeysForConfig(cfg = {}, fallbackIndex = null) {
  const keys = new Set();
  const candidates = [
    cfg.epitope_name,
    cfg.epitope_id,
    epitopeLabel(cfg, fallbackIndex),
  ];
  candidates.forEach((candidate) => {
    const normalized = normalizeEpitopeToken(candidate);
    if (normalized) keys.add(normalized);
  });
  return keys;
}

function parseEpitopeSelectionInput(raw) {
  const text = String(raw || '').trim();
  if (!text) return [];
  return text
    .split(/[\n,]+/)
    .map((part) => part.trim())
    .filter(Boolean)
    .map((part) => {
      const match = part.match(/^([0-9A-Za-z]{4})\s*[:/\-]\s*(.+)$/);
      if (!match) return null;
      return {
        pdbId: (match[1] || '').trim().toUpperCase(),
        epitope: (match[2] || '').trim(),
        raw: part,
      };
    })
    .filter(Boolean);
}

function parseBinderExportSelections(raw) {
  const tokens = String(raw || '')
    .split(/[\n,]+/)
    .map((part) => part.trim())
    .filter(Boolean);
  const entries = [];
  const invalid = [];
  tokens.forEach((token) => {
    const match = token.match(/^([0-9A-Za-z]{4})\s*[:/\-]\s*(.+)$/);
    if (!match) {
      invalid.push(token);
      return;
    }
    const pdbId = (match[1] || '').trim().toUpperCase();
    const epRaw = (match[2] || '').trim();
    const epitope = normalizeEpitopeToken(epRaw);
    if (!pdbId || !epitope) {
      invalid.push(token);
      return;
    }
    entries.push({ pdbId, epitope, raw: token });
  });
  return { entries, invalid };
}

async function downloadNamedFile(name) {
  if (!name) return;
  const url = `/api/bulk/file?name=${encodeURIComponent(name)}`;
  const res = await fetch(url);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || `Failed to download ${name} (${res.status})`);
  }
  const blob = await res.blob();
  const cd = res.headers.get('content-disposition') || '';
  const match = cd.match(/filename=\"?([^\";]+)\"?/i);
  const filename = match?.[1] || name;
  const objUrl = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = objUrl;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(objUrl), 5000);
}

function summarizeTargetStatus(target) {
  if (target.target_job_status) return target.target_job_status;
  const configs = Array.isArray(target.configs) ? target.configs : [];
  if (!configs.length) return null;
  const statuses = configs.map((cfg) => (cfg.job_status || '').toLowerCase());
  if (statuses.some((s) => s === 'running' || s === 'pending')) return 'running';
  if (statuses.some((s) => s === 'failed')) return 'failed';
  if (statuses.some((s) => s === 'success')) return 'success';
  return null;
}

function hasLocalBoltzResults(pdbId) {
  const key = (pdbId || '').toUpperCase();
  const count = state.boltzLocalRuns?.[key];
  const num = Number(count);
  if (Number.isFinite(num)) return num > 0;
  return Boolean(count);
}

function getBoltzDesignTotal(pdbId) {
  const key = (pdbId || '').toUpperCase();
  return state.boltzLocalRuns?.[key] ?? null;
}

function getBoltzSpecKey(cfg) {
  const direct = cfg?.epitope_id;
  if (direct) return String(direct).trim().toUpperCase();
  const path = cfg?.config_path ? String(cfg.config_path) : '';
  if (!path) return null;
  const parts = path.split('/').filter(Boolean);
  const hit = parts.find((part) => /epitope/i.test(part));
  return hit ? hit.toUpperCase() : null;
}

function getBoltzDesignTotalForSpec(pdbId, cfg) {
  const key = (pdbId || '').toUpperCase();
  const specKey = getBoltzSpecKey(cfg);
  if (!specKey) return null;
  const specMap = state.boltzLocalRunsBySpec?.[key] || null;
  if (!specMap) return null;
  const num = Number(specMap[specKey]);
  return Number.isFinite(num) ? num : null;
}

function setActionButtonsDisabled(disabled) {
  [el.bulkRunBtn, el.bulkVisualizeEpitopes].forEach((btn) => {
    if (btn) btn.disabled = disabled;
  });
}

function showAlert(message, isError = true) {
  if (!el.jobAlert) return;
  el.jobAlert.textContent = message;
  el.jobAlert.hidden = false;
  el.jobAlert.classList.toggle('success', !isError);
}

function resetJobLog(message = 'Waiting for commands...') {
  if (el.jobLog) {
    el.jobLog.textContent = message;
  }
}

function appendLog(line) {
  if (!el.jobLog || !line) return;
  const atBottom = Math.abs(el.jobLog.scrollHeight - el.jobLog.scrollTop - el.jobLog.clientHeight) < 8;
  el.jobLog.textContent += `\n${line}`;
  if (atBottom) {
    el.jobLog.scrollTop = el.jobLog.scrollHeight;
  }
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function snapshotsChanged(nextNames = []) {
  const prev = state.snapshotNames || [];
  if (prev.length !== nextNames.length) return true;
  for (let i = 0; i < prev.length; i += 1) {
    if (prev[i] !== nextNames[i]) return true;
  }
  return false;
}

function renderEpitopeMetrics(meta = []) {
  if (!el.epitopeMetricsSection) return;
  el.epitopeMetricsSection.hidden = true;
}

function renderSnapshots(meta = []) {
  if (!el.snapshotGrid || !el.snapshotSection) return;
  el.snapshotGrid.innerHTML = '';
  const rawList = Array.isArray(meta) ? meta : [];
  const byPdb = new Map();
  rawList.forEach((item) => {
    const key = (item?.pdb_id || '').toUpperCase();
    if (!key) return;
    const prev = byPdb.get(key);
    if (!prev || (item.created_at || 0) > (prev.created_at || 0)) {
      byPdb.set(key, item);
    }
  });
  const list = byPdb.size ? Array.from(byPdb.values()) : rawList;
  renderEpitopeMetrics(list);
  if (!list.length) {
    state.snapshotNames = [];
    el.snapshotSection.hidden = true;
    if (el.snapshotDownloadPdf) {
      el.snapshotDownloadPdf.hidden = true;
      el.snapshotDownloadPdf.disabled = true;
    }
    return;
  }
  state.snapshotNames = list.map((item) => item.filename).filter(Boolean);
  if (el.snapshotDownloadPdf) {
    el.snapshotDownloadPdf.hidden = false;
    el.snapshotDownloadPdf.disabled = false;
  }
  const palette = ['#2563eb', '#f97316', '#22c55e', '#a855f7', '#ef4444', '#0ea5e9', '#14b8a6', '#facc15'];
  list.forEach((item) => {
    const card = document.createElement('div');
    card.className = 'snapshot-card';
    const header = document.createElement('div');
    header.className = 'snapshot-header';
    const title = document.createElement('div');
    title.textContent = `${item.pdb_id || 'Target'} · ${item.epitopes?.length || 0} epitopes`;
    header.appendChild(title);
    card.appendChild(header);

    const warnings = Array.isArray(item.warnings) ? item.warnings.filter(Boolean) : [];
    if (warnings.length) {
      const warnBox = document.createElement('div');
      warnBox.className = 'snapshot-warnings';
      warnings.forEach((msg) => {
        const warnLine = document.createElement('div');
        warnLine.textContent = msg;
        warnBox.appendChild(warnLine);
      });
      card.appendChild(warnBox);
    }

    const img = document.createElement('img');
    img.src = item.url;
    img.alt = `${item.pdb_id || 'Target'} hotspot snapshot`;
    img.style.maxWidth = '100%';
    img.style.border = '1px solid #e2e8f0';
    img.style.borderRadius = '8px';
    img.loading = 'lazy';
    card.appendChild(img);

    const epList = document.createElement('ul');
    epList.className = 'snapshot-ep-list';
    const epitopes = Array.isArray(item.epitopes) ? item.epitopes : [];
    epitopes.forEach((ep, idx) => {
      const li = document.createElement('li');
      const color = palette[idx % palette.length];
      const name = ep.name || `Epitope ${idx + 1}`;
      const hotspots = (ep.hotspots || []).join(', ');
      const masks = (ep.mask_residues || []).join(', ');
      const residueText = hotspots || masks || 'No residues recorded';
      li.innerHTML = `<span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${color};margin-right:6px;"></span><strong>${name}</strong>: ${residueText}`;
      epList.appendChild(li);
    });
    card.appendChild(epList);

    const alignment = item.alignment || {};
    const vendorRange = formatRangeLabel(alignment.vendor_range);
    const chainRanges = Array.isArray(alignment.chain_ranges) ? alignment.chain_ranges : [];
    const epCoverage = Array.isArray(alignment.epitope_coverage) ? alignment.epitope_coverage : [];
    const alignNote = alignment.note;
    const hasAlignment = vendorRange || chainRanges.length || epCoverage.length || alignNote;
    if (hasAlignment) {
      const alignBox = document.createElement('div');
      alignBox.className = 'snapshot-align';

      const vendorLine = document.createElement('div');
      vendorLine.innerHTML = `<strong>Product range:</strong> ${vendorRange || 'Not recorded'}`;
      alignBox.appendChild(vendorLine);

      if (chainRanges.length) {
        const chainLine = document.createElement('div');
        const chainText = chainRanges.map((cr) => {
          const chainLabel = cr.chain || '?';
          const pdbRange = formatRangeLabel(cr.range || cr.chain_range) || '—';
          const vendorOverlap = formatRangeLabel(cr.vendor_overlap);
          const stats = [];
          const idPct = formatPercent(cr.identity);
          const covPct = formatPercent(cr.coverage);
          if (idPct) stats.push(`${idPct} id`);
          if (covPct) stats.push(`${covPct} cov`);
          const statsText = stats.length ? ` [${stats.join(' · ')}]` : '';
          const overlapText = vendorOverlap ? ` vs ${vendorOverlap}` : '';
          return `${chainLabel}:${pdbRange}${overlapText}${statsText}`;
        }).join(' · ');
        chainLine.innerHTML = `<strong>PDB overlap:</strong> ${chainText || 'Not mapped'}`;
        alignBox.appendChild(chainLine);
      }

      if (epCoverage.length) {
        const covList = document.createElement('ul');
        covList.className = 'snapshot-coverage';
        epCoverage.forEach((cov) => {
          const li = document.createElement('li');
          const totalNum = Number(cov.total);
          const coveredNum = Number(cov.covered);
          const total = Number.isFinite(totalNum) ? totalNum : null;
          const covered = Number.isFinite(coveredNum) ? coveredNum : null;
          const status = cov.status === 'ok'
            ? 'within product range'
            : cov.status === 'outside'
              ? 'outside product range'
              : 'range unknown';
          const countText = total ? ` (${covered ?? '?'} / ${total})` : '';
          const outsideText = Array.isArray(cov.outside) && cov.outside.length
            ? ` — outside: ${cov.outside.join(', ')}`
            : '';
          li.textContent = `${cov.name || 'Epitope'}: ${status}${countText}${outsideText}`;
          if (cov.status === 'outside') li.classList.add('warn');
          covList.appendChild(li);
        });
        alignBox.appendChild(covList);
      }

      if (alignNote) {
        const noteEl = document.createElement('div');
        noteEl.className = 'snapshot-note';
        noteEl.textContent = alignNote;
        alignBox.appendChild(noteEl);
      }

      card.appendChild(alignBox);
    }
    el.snapshotGrid.appendChild(card);
  });
  el.snapshotSection.hidden = false;
}

async function hydrateSnapshotsFromCache(options = {}) {
  const {
    pdbIds = [],
    limit = 12,
    force = false,
    silent = false,
  } = options;
  if (!el.snapshotGrid || !el.snapshotSection) return;
  if (!force) {
    if (state.jobPoller) return;
    if (el.snapshotGrid.childElementCount > 0) return;
  }
  const params = new URLSearchParams();
  const ids = Array.isArray(pdbIds)
    ? Array.from(new Set(pdbIds.map((id) => (id || '').trim().toUpperCase()).filter(Boolean)))
    : [];
  if (ids.length) params.append('pdb_ids', ids.join(','));
  const capped = Number(limit);
  params.append('limit', String(Number.isFinite(capped) && capped > 0 ? Math.min(200, Math.round(capped)) : 12));
  try {
    const res = await fetch(`/api/pymol/snapshots/latest?${params.toString()}`, { cache: 'no-store' });
    if (!res.ok) {
      if (!silent) appendLog(`Cached snapshot fetch failed (${res.status})`);
      return;
    }
    const data = await res.json();
    renderSnapshots(Array.isArray(data) ? data : []);
    if (!silent && Array.isArray(data) && data.length) {
      appendLog(`Loaded ${data.length} cached hotspot snapshot${data.length === 1 ? '' : 's'}.`);
    }
  } catch (err) {
    if (!silent) appendLog(`Cached snapshot fetch error: ${err.message || err}`);
  }
}

function renderEpitopePlots(items = []) {
  if (!el.epitopePlotGrid || !el.epitopePlotSection) return;
  el.epitopePlotGrid.innerHTML = '';
  const list = Array.isArray(items) ? items.filter(Boolean) : [];
  if (!list.length) {
    state.epitopePlots = [];
    el.epitopePlotSection.hidden = true;
    state.epitopeCsvName = null;
    state.epitopeHotspotCsvName = null;
    if (el.epitopeCsvDownload) el.epitopeCsvDownload.hidden = true;
    if (el.epitopeHotspotCsvDownload) el.epitopeHotspotCsvDownload.hidden = true;
    return;
  }
  state.epitopePlots = list;
  list.forEach((item) => {
    const label = item.label || (item.src || '').split('/').pop() || 'Epitope plot';
    const src = item.src;
    if (!src) return;
    const card = document.createElement('div');
    card.className = 'snapshot-card';
    const header = document.createElement('div');
    header.className = 'snapshot-header';
    const title = document.createElement('div');
    title.textContent = label;
    header.appendChild(title);
    card.appendChild(header);

    const imgBox = document.createElement('div');
    imgBox.className = 'snapshot-image';
    const img = document.createElement('img');
    img.src = src;
    img.alt = label;
    img.loading = 'lazy';
    img.style.maxWidth = '100%';
    img.style.border = '1px solid #e2e8f0';
    img.style.borderRadius = '8px';
    imgBox.appendChild(img);
    card.appendChild(imgBox);

    el.epitopePlotGrid.appendChild(card);
  });
  if (el.epitopeCsvDownload) {
    el.epitopeCsvDownload.hidden = !state.epitopeCsvName;
  }
  if (el.epitopeHotspotCsvDownload) {
    el.epitopeHotspotCsvDownload.hidden = !state.epitopeHotspotCsvName;
  }
  el.epitopePlotSection.hidden = false;
}

function renderDiversityPlots(items = []) {
  if (!el.diversityGrid || !el.diversitySection) return;
  el.diversityGrid.innerHTML = '';
  const outputDir = (state.diversityOutputDir || '').trim();
  if (outputDir) {
    const dirNote = document.createElement('div');
    dirNote.className = 'help-text';
    dirNote.style.margin = '0 0 8px 0';
    dirNote.textContent = `Output directory: ${outputDir}`;
    el.diversityGrid.appendChild(dirNote);
  }
  const files = Array.isArray(state.diversityFiles) ? state.diversityFiles : [];
  if (files.length) {
    const details = document.createElement('details');
    details.style.margin = '0 0 10px 0';
    const summary = document.createElement('summary');
    summary.textContent = `Detected ${files.length} all_designs_metrics.csv file${files.length === 1 ? '' : 's'}`;
    summary.style.cursor = 'pointer';
    summary.style.userSelect = 'none';
    details.appendChild(summary);

    const pre = document.createElement('pre');
    pre.style.marginTop = '8px';
    pre.style.maxHeight = '220px';
    pre.style.overflow = 'auto';
    pre.style.fontSize = '12px';
    pre.style.background = '#0f172a';
    pre.style.color = '#e2e8f0';
    pre.style.padding = '10px 12px';
    pre.style.borderRadius = '8px';
    pre.textContent = files
      .map((f) => (f && (f.path || f.metrics_path)) ? `${f.pdb_id || ''}\t${f.epitope_name || ''}\t${f.datetime || ''}\t${f.path || f.metrics_path}` : String(f))
      .join('\n');
    details.appendChild(pre);
    el.diversityGrid.appendChild(details);
  }
  const list = Array.isArray(items) ? items.filter((p) => p && (p.png_name || p.svg_name)) : [];
  if (!list.length) {
    state.diversityPlots = [];
    const empty = document.createElement('div');
    empty.className = 'help-text';
    empty.style.margin = '0';
    const message = (state.diversityMessage || '').trim();
    empty.textContent =
      message ||
      'No diversity plots yet. For BoltzGen, download metrics to targets/<PDB>/designs/boltzgen/epitope_*/final_ranked_designs/all_designs_metrics.csv (or designs/boltzgen/<run_label>/epitope_*/final_ranked_designs). For RFAntibody, sync af3_rankings.tsv under targets/<PDB>/designs/rfantibody/_assessments/<run_label>/ (legacy: designs/_assessments), then click Refresh to rebuild all_design_metrics.';
    el.diversityGrid.appendChild(empty);
  }
  state.diversityPlots = list;
  list.forEach((plot) => {
    const card = document.createElement('div');
    card.className = 'snapshot-card';
    const header = document.createElement('div');
    header.className = 'snapshot-header';
    const title = document.createElement('div');
    title.textContent = plot.pdb_id || 'BoltzGen';
    header.appendChild(title);
    const linkBar = document.createElement('div');
    linkBar.className = 'snapshot-links';
    if (plot.png_name) {
      const pngLink = document.createElement('a');
      pngLink.href = `/api/bulk/file?name=${encodeURIComponent(plot.png_name)}`;
      pngLink.textContent = 'PNG';
      pngLink.target = '_blank';
      linkBar.appendChild(pngLink);
    }
    if (plot.svg_name) {
      const svgLink = document.createElement('a');
      svgLink.href = `/api/bulk/file?name=${encodeURIComponent(plot.svg_name)}`;
      svgLink.textContent = 'SVG';
      svgLink.target = '_blank';
      if (linkBar.childElementCount) linkBar.appendChild(document.createTextNode(' · '));
      linkBar.appendChild(svgLink);
    }
    header.appendChild(linkBar);
    card.appendChild(header);

    const pathRow = document.createElement('div');
    pathRow.className = 'snapshot-note';
    const pngPath = plot.png_path || (plot.png_name ? `${outputDir}/${plot.png_name}` : '');
    const svgPath = plot.svg_path || (plot.svg_name ? `${outputDir}/${plot.svg_name}` : '');
    const parts = [];
    if (pngPath) parts.push(`PNG: ${pngPath}`);
    if (svgPath) parts.push(`SVG: ${svgPath}`);
    pathRow.textContent = parts.length ? parts.join(' · ') : 'Plot paths unavailable';
    card.appendChild(pathRow);

    if (plot.png_name) {
      const imgBox = document.createElement('div');
      imgBox.className = 'snapshot-image';
      const img = document.createElement('img');
      img.src = `/api/bulk/file?name=${encodeURIComponent(plot.png_name)}`;
      img.alt = plot.pdb_id || 'BoltzGen diversity';
      img.loading = 'lazy';
      img.style.maxWidth = '100%';
      img.style.border = '1px solid #e2e8f0';
      img.style.borderRadius = '8px';
      imgBox.appendChild(img);
      card.appendChild(imgBox);
    }

    const legend = document.createElement('div');
    legend.className = 'snapshot-note';
    const colors = plot.epitope_colors || {};
    if (Object.keys(colors).length) {
      legend.innerHTML = Object.entries(colors)
        .map(([name, color]) => `<span style="display:inline-flex;align-items:center;gap:6px;margin-right:8px;"><span style="width:12px;height:12px;border-radius:999px;background:${color};display:inline-block;"></span>${name}</span>`)
        .join('');
    } else {
      legend.textContent = 'Epitope colors not available';
    }
    card.appendChild(legend);

    el.diversityGrid.appendChild(card);
  });

  if (el.diversityDownloadCsv) {
    el.diversityDownloadCsv.hidden = !state.diversityCsv;
  }
  if (el.diversityDownloadHtml) {
    el.diversityDownloadHtml.hidden = !state.diversityHtml;
  }
  el.diversitySection.hidden = false;
}

async function downloadEpitopeReportHtml() {
  if (el.epitopeReportDownload) el.epitopeReportDownload.disabled = true;
  try {
    if (!state.currentJobId) {
      showAlert('Run a bulk job first to generate a report.');
      return;
    }
    const res = await fetch(`/api/jobs/${encodeURIComponent(state.currentJobId)}/epitope-report`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to generate report (${res.status})`);
    }
    const blob = await res.blob();
    const cd = res.headers.get('content-disposition') || '';
    const match = cd.match(/filename=\"?([^\";]+)\"?/i);
    const filename = match?.[1] || 'epitope_report.html';
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  } catch (err) {
    showAlert(err.message || 'Unable to generate report.');
  } finally {
    if (el.epitopeReportDownload) el.epitopeReportDownload.disabled = false;
  }
}

async function downloadEpitopeCsv() {
  if (el.epitopeCsvDownload) el.epitopeCsvDownload.disabled = true;
  try {
    if (!state.epitopeCsvName) {
      showAlert('No epitope diversity CSV available.');
      return;
    }
    const src = bulkFileSrc(state.epitopeCsvName);
    if (!src) {
      showAlert('CSV not available.');
      return;
    }
    const res = await fetch(src);
    if (!res.ok) {
      const detail = await res.text().catch(() => '');
      throw new Error(detail || `Failed to download CSV (${res.status})`);
    }
    const blob = await res.blob();
    const filename = state.epitopeCsvName || 'epitope_diversity.csv';
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  } catch (err) {
    showAlert(err.message || 'Unable to download CSV.');
  } finally {
    if (el.epitopeCsvDownload) el.epitopeCsvDownload.disabled = false;
  }
}

async function downloadEpitopeHotspotCsv() {
  if (el.epitopeHotspotCsvDownload) el.epitopeHotspotCsvDownload.disabled = true;
  try {
    if (!state.epitopeHotspotCsvName) {
      showAlert('No hotspot diversity CSV available.');
      return;
    }
    const src = bulkFileSrc(state.epitopeHotspotCsvName);
    if (!src) {
      showAlert('CSV not available.');
      return;
    }
    const res = await fetch(src);
    if (!res.ok) {
      const detail = await res.text().catch(() => '');
      throw new Error(detail || `Failed to download CSV (${res.status})`);
    }
    const blob = await res.blob();
    const filename = state.epitopeHotspotCsvName || 'hotspot_diversity.csv';
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  } catch (err) {
    showAlert(err.message || 'Unable to download hotspot CSV.');
  } finally {
    if (el.epitopeHotspotCsvDownload) el.epitopeHotspotCsvDownload.disabled = false;
  }
}

async function refreshDiversity({ silent = false, page = null, force = false } = {}) {
  try {
    const params = new URLSearchParams();
    const nextPage = page ?? state.binderPage ?? 1;
    params.append('page', String(nextPage));
    params.append('page_size', String(state.binderPageSize || 100));
    const filters = getBinderFilters();
    if (filters.pdb) params.append('filter_pdb', filters.pdb);
    if (filters.epitope) params.append('filter_epitope', filters.epitope);
    if (filters.engine) params.append('filter_engine', filters.engine);
    if (filters.orderBy) params.append('order_by', filters.orderBy);
    const baseUrl = force ? '/api/bulk/boltzgen/diversity/refresh' : '/api/bulk/boltzgen/diversity';
    const res = await fetch(`${baseUrl}?${params.toString()}`, { method: force ? 'POST' : 'GET' });
    if (!res.ok) throw new Error('Unable to load diversity report');
    const data = await res.json();
    state.diversityCsv = data.csv_name || null;
    state.diversityHtml = data.html_name || null;
    state.diversityMessage = data.message || null;
    state.diversityFiles = Array.isArray(data.metrics_files) ? data.metrics_files : [];
    state.diversityOutputDir = data.output_dir || null;
    state.diversityPlots = Array.isArray(data.plots) ? data.plots : [];
    state.boltzDesignCounts = normalizeBoltzCounts(data.binder_counts);
    const pdbIds = new Set();
    const binderRows = Array.isArray(data.binder_rows) ? data.binder_rows : [];
    (state.diversityPlots || []).forEach((plot) => {
      const val = (plot.pdb_id || '').trim().toUpperCase();
      if (val) pdbIds.add(val);
    });
    (state.diversityFiles || []).forEach((file) => {
      const val = (file.pdb_id || '').trim().toUpperCase();
      if (val) pdbIds.add(val);
    });
    binderRows.forEach((row) => {
      const val = (row.pdb_id || '').trim().toUpperCase();
      if (val) pdbIds.add(val);
    });
    state.binderPdbIds = Array.from(pdbIds);
    if (state.diversityCsv) {
      state.binderCsvName = state.diversityCsv;
    }
    renderDiversityPlots(state.diversityPlots);
    if (!silent && data.message) {
      showAlert(data.message, false);
    }
    if (binderRows.length || typeof data.binder_total === 'number') {
      state.binderRows = binderRows;
      state.binderTotal = Number(data.binder_total || binderRows.length || 0);
      state.binderPage = Number(data.binder_page || nextPage || 1);
      state.binderPageSize = Number(data.binder_page_size || state.binderPageSize || 100);
      state.binderMessage = data.binder_message || data.message || null;
      renderBinderRows(state.binderRows);
    } else {
      state.binderRows = [];
      state.binderTotal = 0;
      state.binderMessage = data.message || null;
      await loadBinderTable({ page: nextPage, silent: true });
    }
    if (Object.keys(state.boltzDesignCounts || {}).length) {
      state.boltzLocalRuns = { ...state.boltzDesignCounts };
      state.boltzLocalRunsBySpec = {};
      if (Array.isArray(state.boltzConfigs) && state.boltzConfigs.length) {
        renderBoltzConfigs();
      }
    }
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
  }
}

async function forceRefreshDiversity() {
  const ok = window.confirm('Clear cached diversity plots and rebuild now? This can take a while.');
  if (!ok) return;
  if (el.diversityRebuild) el.diversityRebuild.disabled = true;
  try {
    await refreshDiversity({ silent: false, force: true });
  } finally {
    if (el.diversityRebuild) el.diversityRebuild.disabled = false;
  }
}

async function downloadDiversityFile(kind) {
  const name = kind === 'csv' ? state.diversityCsv : state.diversityHtml;
  if (!name) {
    showAlert(`No ${kind.toUpperCase()} available yet.`);
    return;
  }

  const button = kind === 'csv' ? el.diversityDownloadCsv : el.diversityDownloadHtml;
  if (button) button.disabled = true;

  try {
    const url = `/api/bulk/file?name=${encodeURIComponent(name)}`;
    const res = await fetch(url);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to download (${res.status})`);
    }
    const blob = await res.blob();
    const cd = res.headers.get('content-disposition') || '';
    const match = cd.match(/filename=\"?([^\";]+)\"?/i);
    const filename = match?.[1] || name;

    const objectUrl = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = objectUrl;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(objectUrl), 5000);
  } catch (err) {
    showAlert(err.message || `Unable to download ${kind.toUpperCase()}.`);
  } finally {
    if (button) button.disabled = false;
  }
}

function renderBinderRows(rows = []) {
  if (!el.binderTable || !el.binderPanel) return;
  const tbody = el.binderTable;
  tbody.innerHTML = '';
  const list = Array.isArray(rows) ? rows.filter(Boolean) : [];

  const noteRow = () => {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = 9;
    td.style.color = '#64748b';
    td.textContent = state.binderTotal ? 'No binders on this page.' : (state.binderMessage || 'No binders found.');
    tr.appendChild(td);
    return tr;
  };

  if (!list.length) {
    tbody.appendChild(noteRow());
  } else {
    list.forEach((row, idx) => {
      const epLabel = row.epitope_id || row.epitope || '—';
      const tr = document.createElement('tr');
      const engineLabel = formatEngineLabel(row.engine);
      const values = [
        row.pdb_id || '—',
        epLabel,
        engineLabel,
        row.iptm !== null && row.iptm !== undefined ? Number(row.iptm).toFixed(3) : '—',
        row.rmsd !== null && row.rmsd !== undefined ? Number(row.rmsd).toFixed(3) : '—',
        row.hotspot_dist !== null && row.hotspot_dist !== undefined ? Number(row.hotspot_dist).toFixed(2) : '—',
        row.ipsae_min !== null && row.ipsae_min !== undefined ? Number(row.ipsae_min).toFixed(3) : '—',
      ];
      values.forEach((val, colIdx) => {
        const td = document.createElement('td');
        td.textContent = val;
        if (colIdx === 1) {
          const color = epitopeColor(val);
          if (color) td.style.color = color;
        }
        tr.appendChild(td);
      });

      const pathTd = document.createElement('td');
      if (row.design_path) {
        const link = document.createElement('a');
        link.href = `file://${row.design_path}`;
        link.textContent = row.design_path.split('/').pop() || row.design_path;
        link.title = row.design_path;
        link.target = '_blank';
        link.rel = 'noopener noreferrer';
        pathTd.appendChild(link);
      } else {
        pathTd.textContent = '—';
      }
      tr.appendChild(pathTd);

      const pymolTd = document.createElement('td');
      const pymolBtn = document.createElement('button');
      pymolBtn.type = 'button';
      pymolBtn.textContent = 'PyMOL';
      pymolBtn.dataset.action = 'pymol-binder';
      pymolBtn.dataset.index = String(idx);
      if (!row.design_path || !isBoltzEngine(row.engine)) pymolBtn.disabled = true;
      pymolTd.appendChild(pymolBtn);
      tr.appendChild(pymolTd);

      tbody.appendChild(tr);
    });
  }

  if (el.binderSummary) {
    const msg = state.binderMessage || '';
    const summaryText = msg || `Showing ${list.length} of ${state.binderTotal || rows.length} binders`;
    el.binderSummary.textContent = summaryText;
    el.binderSummary.hidden = false;
  }

  if (el.binderCsvNote) {
    const csvName = state.binderCsvName || state.diversityCsv;
    if (csvName) {
      const outputDir = (state.diversityOutputDir || '').trim();
      const fullPath = outputDir ? `${outputDir}/${csvName}` : csvName;
      el.binderCsvNote.textContent = `Detected CSV: ${fullPath}`;
      el.binderCsvNote.hidden = false;
    } else {
      el.binderCsvNote.hidden = true;
      el.binderCsvNote.textContent = '';
    }
  }

  const totalPages = Math.max(1, Math.ceil((state.binderTotal || 0) / (state.binderPageSize || 1)));
  if (el.binderPagination && el.binderPageLabel) {
    el.binderPagination.hidden = totalPages <= 1;
    el.binderPageLabel.textContent = `Page ${state.binderPage || 1} of ${totalPages}`;
    const buttons = el.binderPagination.querySelectorAll('button[data-page-nav]');
    buttons.forEach((btn) => {
      const dir = btn.dataset.pageNav;
      if (dir === 'prev') btn.disabled = (state.binderPage || 1) <= 1;
      if (dir === 'next') btn.disabled = (state.binderPage || 1) >= totalPages;
    });
  }

  if (el.binderDownload) {
    const csvName = state.binderCsvName || state.diversityCsv;
    el.binderDownload.hidden = !csvName;
  }

  el.binderPanel.hidden = false;
}

async function loadBinderTable({ page = 1, silent = false, force = false } = {}) {
  if (!el.binderPanel) return;
  const ids = Array.from(
    new Set([
      ...(Array.isArray(state.bulkPreviewRows)
        ? state.bulkPreviewRows.map((row) => (row.resolved_pdb_id || row.pdb_id || '').trim().toUpperCase())
        : []),
      ...(Array.isArray(state.binderPdbIds) ? state.binderPdbIds : []),
    ].filter(Boolean)),
  );
  if (!ids.length) {
    el.binderPanel.hidden = true;
    return;
  }
  const params = new URLSearchParams();
  params.append('pdb_ids', ids.join(','));
  params.append('page', String(page));
  params.append('page_size', String(state.binderPageSize || 100));
  const filters = getBinderFilters();
  if (filters.pdb) params.append('filter_pdb', filters.pdb);
  if (filters.epitope) params.append('filter_epitope', filters.epitope);
  if (filters.engine) params.append('filter_engine', filters.engine);
  if (filters.orderBy) params.append('order_by', filters.orderBy);
  if (force) params.append('_ts', String(Date.now()));
  try {
    if (el.binderRefresh) el.binderRefresh.disabled = true;
    const res = await fetch(`/api/bulk/boltzgen/binders?${params.toString()}`, { cache: 'no-store' });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to load binders (${res.status})`);
    }
    const body = await res.json();
    state.binderRows = Array.isArray(body.rows) ? body.rows : [];
    state.binderTotal = Number(body.total_rows || state.binderRows.length || 0);
    state.binderPage = Number(body.page || page || 1);
    state.binderPageSize = Number(body.page_size || state.binderPageSize || 100);
    state.binderCsvName = body.csv_name || state.diversityCsv || null;
    state.binderMessage = body.message || null;
    renderBinderRows(state.binderRows);
    if (!silent && body.message) {
      showAlert(body.message, false);
    }
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
    if (el.binderPanel) el.binderPanel.hidden = true;
  } finally {
    if (el.binderRefresh) el.binderRefresh.disabled = false;
  }
}

async function downloadBinderCsv() {
  const csvName = state.binderCsvName || state.diversityCsv;
  if (!csvName) {
    showAlert('No binder CSV available yet.');
    return;
  }
  if (el.binderDownload) el.binderDownload.disabled = true;
  try {
    const url = `/api/bulk/file?name=${encodeURIComponent(csvName)}`;
    const res = await fetch(url);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to download binder CSV (${res.status})`);
    }
    const blob = await res.blob();
    const cd = res.headers.get('content-disposition') || '';
    const match = cd.match(/filename=\"?([^\";]+)\"?/i);
    const filename = match?.[1] || csvName;
    const objUrl = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = objUrl;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(objUrl), 5000);
  } catch (err) {
    showAlert(err.message || 'Unable to download binder CSV.');
  } finally {
    if (el.binderDownload) el.binderDownload.disabled = false;
  }
}

async function exportSelectedBinders() {
  const selectionRaw = el.binderExportSelections?.value || '';
  const { entries, invalid } = parseBinderExportSelections(selectionRaw);
  if (!entries.length) {
    showAlert('Enter at least one PDB:epitope selection to export.');
    return;
  }
  if (invalid.length) {
    showAlert(`Invalid selections skipped: ${invalid.join(', ')}`);
  }
  const seen = new Set();
  const uniqueEntries = entries.filter((entry) => {
    const key = `${entry.pdbId}:${entry.epitope}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });
  const countRaw = Number(el.binderExportCount?.value || 0);
  const perGroup = Number.isFinite(countRaw) && countRaw > 0 ? Math.floor(countRaw) : 0;
  if (!perGroup) {
    showAlert('Enter a valid number of binders per antigen:epitope.');
    return;
  }

  if (el.binderExportConfirm) el.binderExportConfirm.disabled = true;
  try {
    const selections = uniqueEntries.map((entry) => `${entry.pdbId}:${entry.epitope}`);
    const payload = {
      selections,
      per_group: perGroup,
      include_summary: Boolean(el.binderExportSummary?.checked),
    };
    const res = await fetch('/api/bulk/boltzgen/binders/export', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Export failed (${res.status})`);
    }
    const body = await res.json();
    if (!body.csv_name) {
      showAlert(body.message || 'No binder export produced.');
      return;
    }
    await downloadNamedFile(body.csv_name);
    if (body.summary_csv_name) {
      await downloadNamedFile(body.summary_csv_name);
    }
    showAlert(body.message || `Export ready: ${body.csv_name}`, false);
  } catch (err) {
    showAlert(err.message || 'Failed to export selected binders.');
  } finally {
    if (el.binderExportConfirm) el.binderExportConfirm.disabled = false;
  }
}

async function launchBinderPymol(index) {
  const row = Array.isArray(state.binderRows) ? state.binderRows[index] : null;
  if (!row) {
    showAlert('No binder selected.');
    return;
  }
  if (!isBoltzEngine(row.engine)) {
    showAlert('PyMOL launch is only available for BoltzGen binders.');
    return;
  }
  if (!row.design_path) {
    showAlert('Design path unavailable for this binder.');
    return;
  }
  try {
    const epLabel = row.epitope_id || row.epitope || null;
    const payload = {
      pdb_id: row.pdb_id,
      design_path: row.design_path,
      epitope_label: epLabel,
      binding_label: row.binding_label || null,
      include_label: row.include_label || null,
      target_path: row.target_path || null,
    };
    const res = await fetch('/api/bulk/boltzgen/binder/pymol', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `PyMOL launch failed (${res.status})`);
    }
    const body = await res.json();
    showAlert(`PyMOL script ready: ${body.script_path}`, false);
  } catch (err) {
    showAlert(err.message || 'Failed to launch PyMOL.');
  }
}

async function launchBoltzPymol(pdbId, epitopeName = null, triggerBtn = null) {
  if (!pdbId) {
    showAlert('Missing PDB ID for PyMOL.');
    return;
  }
  if (triggerBtn) triggerBtn.disabled = true;
  try {
    const payload = {
      launch: true,
      bundle_only: false,
      epitope_name: epitopeName || null,
    };
    const res = await fetch(`/api/targets/${encodeURIComponent(pdbId)}/pymol/hotspots`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `PyMOL launch failed (${res.status})`);
    }
    const body = await res.json();
    const scopeText = epitopeName ? `epitope ${epitopeName}` : 'all epitopes';
    showAlert(`PyMOL hotspot bundle ready (${pdbId} · ${scopeText})`, false);
    appendLog(body.message || `PyMOL bundle ready for ${pdbId}`);
  } catch (err) {
    showAlert(err.message || 'Failed to launch PyMOL.');
  } finally {
    if (triggerBtn) triggerBtn.disabled = false;
  }
}

function openPipelineRerunModal(pdbId, antigenUrl = null) {
  if (!el.pipelineRerunModal) return;
  if (!pdbId) {
    showAlert('Missing PDB ID for pipeline rerun.');
    return;
  }
  el.pipelineRerunModal.dataset.pdbId = pdbId;
  if (antigenUrl) {
    el.pipelineRerunModal.dataset.antigenUrl = antigenUrl;
  } else {
    delete el.pipelineRerunModal.dataset.antigenUrl;
  }
  if (el.pipelineRerunTitle) el.pipelineRerunTitle.textContent = `Select epitopes (LLM) · ${pdbId}`;
  if (el.pipelineRerunExpected) el.pipelineRerunExpected.value = '10';
  if (el.pipelineRerunAttempts) el.pipelineRerunAttempts.value = '3';
  if (el.pipelineRerunForce) el.pipelineRerunForce.checked = true;
  if (el.pipelineRerunPrompt) el.pipelineRerunPrompt.value = '';
  toggleModal(el.pipelineRerunModal, true);
}

function openPipelineRerunModalBulk(targets = []) {
  if (!el.pipelineRerunModal) return;
  const ids = targets.map((t) => (t?.pdb_id || '').trim().toUpperCase()).filter(Boolean);
  if (!ids.length) {
    showAlert('No targets selected for pipeline rerun.');
    return;
  }
  el.pipelineRerunModal.dataset.bulkIds = ids.join(',');
  el.pipelineRerunModal.dataset.pdbId = '';
  el.pipelineRerunModal.dataset.antigenUrl = '';
  if (el.pipelineRerunTitle) {
    el.pipelineRerunTitle.textContent = `Select epitopes (LLM) · ${ids.length} targets`;
  }
  if (el.pipelineRerunExpected) el.pipelineRerunExpected.value = '10';
  if (el.pipelineRerunAttempts) el.pipelineRerunAttempts.value = '3';
  if (el.pipelineRerunForce) el.pipelineRerunForce.checked = true;
  if (el.pipelineRerunPrompt) el.pipelineRerunPrompt.value = '';
  toggleModal(el.pipelineRerunModal, true);
}

function findBulkRowForPdb(pdbId) {
  const key = (pdbId || '').trim().toUpperCase();
  if (!key) return null;
  const rows = Array.isArray(state.bulkPreviewRows) ? state.bulkPreviewRows : [];
  return rows.find((row) => {
    const resolved = (row.resolved_pdb_id || '').trim().toUpperCase();
    const raw = (row.pdb_id || '').trim().toUpperCase();
    return resolved === key || raw === key;
  }) || null;
}

function findAccessionsFromInput(pdbId) {
  const key = (pdbId || '').trim().toUpperCase();
  const rawText = el.bulkCsvInput?.value || '';
  if (!key || !rawText.trim()) return null;
  const lines = rawText.split(/\r?\n/).filter((line) => line.trim());
  if (!lines.length) return null;
  const headerLine = lines[0];
  const delimiter = headerLine.includes('\t') && !headerLine.includes(',') ? '\t' : ',';
  const header = headerLine.split(delimiter).map((cell) => cell.trim().toLowerCase());
  const pdbIdx = header.findIndex((cell) => (
    ['chosen_pdb', 'pdb', 'pdb_id', 'pdbid'].includes(cell)
  ));
  const accIdx = header.findIndex((cell) => (
    ['vendor_accession'].includes(cell)
  ));
  const rangeIdx = header.findIndex((cell) => (
    ['vendor_range', 'pdb_vendor_intersection', 'vendor_overlap_range'].includes(cell)
  ));
  if (pdbIdx < 0 || accIdx < 0) return null;
  for (let i = 1; i < lines.length; i += 1) {
    const cols = lines[i].split(delimiter);
    if (cols.length <= Math.max(pdbIdx, accIdx)) continue;
    const pdbVal = (cols[pdbIdx] || '').trim().toUpperCase();
    if (pdbVal === key) {
      const accVal = (cols[accIdx] || '').trim();
      const rangeVal = rangeIdx >= 0 && cols.length > rangeIdx ? (cols[rangeIdx] || '').trim() : '';
      return {
        accession: accVal || null,
        vendor_range: rangeVal || null,
      };
    }
  }
  return null;
}

async function queuePipelineRerunTargets(list, triggerBtn = null, overrides = {}) {
  if (!Array.isArray(list) || !list.length) {
    showAlert('No targets selected for pipeline rerun.');
    return;
  }
  if (triggerBtn) triggerBtn.disabled = true;
  try {
    const expectedRaw = overrides.expected ?? el.pipelineRerunExpected?.value;
    const attemptsRaw = overrides.attempts ?? el.pipelineRerunAttempts?.value;
    const force = typeof overrides.force === 'boolean' ? overrides.force : Boolean(el.pipelineRerunForce?.checked);
    const promptRaw = overrides.prompt ?? el.pipelineRerunPrompt?.value ?? '';
    const expected = expectedRaw ? Number(expectedRaw) || null : null;
    const attempts = attemptsRaw ? Number(attemptsRaw) || 3 : 3;
    const decideScopePrompt = String(promptRaw || '').trim() || null;
    const designCountRaw = overrides.designCount ?? getBoltzDesignCount();
    const designCount = designCountRaw ? Number(designCountRaw) || null : null;
    for (let i = 0; i < list.length; i += 1) {
      const item = list[i];
      const id = (item?.pdb_id || item || '').trim();
      if (!id) continue;
      const row = findBulkRowForPdb(id);
      const fallback = row ? null : findAccessionsFromInput(id);
      const antigenUrl = (item?.antigen_url || row?.antigen_url || '').trim() || null;
      const payload = {
        pdb_id: id,
        force,
        expected_epitopes: expected,
        decide_scope_attempts: attempts,
        decide_scope_prompt: decideScopePrompt,
        design_count: designCount,
        antigen_url: antigenUrl,
        target_accession: row?.accession || fallback?.accession || null,
        target_vendor_range: row?.vendor_range || fallback?.vendor_range || null,
      };
      const res = await fetch('/api/targets/pipeline/refresh', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });
      if (!res.ok) {
        const detail = await res.json().catch(() => ({}));
        throw new Error(detail.detail || `Pipeline refresh failed (${res.status})`);
      }
      const body = await res.json();
      if (body.job_id) {
        startJobPolling(body.job_id);
      }
      if (i < list.length - 1 && PIPELINE_RERUN_DELAY_MS > 0) {
        appendLog(`Waiting ${Math.round(PIPELINE_RERUN_DELAY_MS / 1000)}s before next pipeline refresh...`);
        await sleep(PIPELINE_RERUN_DELAY_MS);
      }
    }
    showAlert(`Pipeline refresh queued for ${list.length} targets.`, false);
  } catch (err) {
    showAlert(err.message || 'Failed to rerun pipeline.');
  } finally {
    if (triggerBtn) triggerBtn.disabled = false;
  }
}

function openPipelineRerunRangeModal() {
  if (!el.boltzRerunRangeModal) return;
  if (el.boltzRerunRangeExpected) el.boltzRerunRangeExpected.value = '10';
  if (el.boltzRerunRangeAttempts) el.boltzRerunRangeAttempts.value = '3';
  if (el.boltzRerunRangeForce) el.boltzRerunRangeForce.checked = true;
  if (el.boltzRerunRangePrompt) el.boltzRerunRangePrompt.value = '';
  toggleModal(el.boltzRerunRangeModal, true);
}

function getBoltzRowRange(startEl, endEl, targets, label) {
  const startRaw = Number(startEl?.value || 0);
  const endRaw = Number(endEl?.value || 0);
  if (!Number.isFinite(startRaw) || !Number.isFinite(endRaw) || startRaw <= 0 || endRaw <= 0) {
    showAlert(`Enter a valid ${label} row range (start/end).`);
    return null;
  }
  const start = Math.min(startRaw, endRaw);
  const end = Math.max(startRaw, endRaw);
  if (start < 1 || end > targets.length) {
    showAlert(`Row range must be between 1 and ${targets.length}.`);
    return null;
  }
  return { start, end, slice: targets.slice(start - 1, end) };
}

function getRunCommandFilters() {
  const allowedRaw = Number(el.boltzRunRangeMinAllowed?.value || 0);
  const epitopesRaw = Number(el.boltzRunRangeMinEpitopes?.value || 0);
  const minAllowedLength = Number.isFinite(allowedRaw) && allowedRaw > 0 ? Math.round(allowedRaw) : null;
  const minEpitopes = Number.isFinite(epitopesRaw) && epitopesRaw > 0 ? Math.round(epitopesRaw) : null;
  return { minAllowedLength, minEpitopes };
}

function getTargetEpitopeCount(target) {
  const direct = asNumber(target?.epitope_count);
  if (direct !== null) return direct;
  const configs = Array.isArray(target?.configs) ? target.configs : [];
  return configs.length;
}

function targetPassesRunFilters(target, filters) {
  if (!filters) return true;
  if (filters.minAllowedLength) {
    const allowed = asNumber(target?.allowed_epitope_length);
    if (allowed === null || allowed < filters.minAllowedLength) return false;
  }
  if (filters.minEpitopes) {
    const count = getTargetEpitopeCount(target);
    if (!Number.isFinite(count) || count < filters.minEpitopes) return false;
  }
  return true;
}

function applyRunCommandFilters(targets, filters) {
  const pass = [];
  const fail = [];
  (targets || []).forEach((target) => {
    if (targetPassesRunFilters(target, filters)) {
      pass.push(target);
    } else {
      fail.push(target);
    }
  });
  return { pass, fail };
}

function runCommandFilterSection(filters, totalCount, passCount) {
  const lines = [];
  if (filters?.minAllowedLength) {
    lines.push(`allowed_epitope_range length >= ${filters.minAllowedLength}`);
  }
  if (filters?.minEpitopes) {
    lines.push(`epitopes selected >= ${filters.minEpitopes}`);
  }
  if (!lines.length) return '';
  const summary = [`# Filter criteria (target.yaml)`, ...lines];
  if (Number.isFinite(totalCount) && Number.isFinite(passCount)) {
    summary.push(`targets included: ${passCount} of ${totalCount}`);
  }
  return summary.join('\n');
}

function updateRunCommandFilterNote() {
  if (!el.boltzRunRangeFilterNote) return;
  const filters = getRunCommandFilters();
  const parts = [];
  if (filters.minAllowedLength) {
    parts.push(`allowed_epitope_range length >= ${filters.minAllowedLength}`);
  }
  if (filters.minEpitopes) {
    parts.push(`epitopes selected >= ${filters.minEpitopes}`);
  }
  el.boltzRunRangeFilterNote.textContent = parts.length
    ? `Filters applied to commands (target.yaml): ${parts.join('; ')}.`
    : 'No filters applied to manual commands.';
  if (el.boltzRunRangeAllowedNote) {
    el.boltzRunRangeAllowedNote.textContent = (
      'allowed_epitope_range is read from target.yaml (allowed_epitope_range or allowed_range). '
      + 'Length is computed by summing inclusive residue spans per token '
      + '(e.g. A:10-20,B:1-5 -> 16).'
    );
  }
  updateRunCommandCountNote(filters);
}

function previewRunCommandRange(targets) {
  const startRaw = Number(el.boltzRunRangeStart?.value || 0);
  const endRaw = Number(el.boltzRunRangeEnd?.value || 0);
  if (!Number.isFinite(startRaw) || !Number.isFinite(endRaw) || startRaw <= 0 || endRaw <= 0) {
    return null;
  }
  const start = Math.min(startRaw, endRaw);
  const end = Math.max(startRaw, endRaw);
  if (!targets.length || start < 1 || end > targets.length) return null;
  const slice = targets.slice(start - 1, end);
  return { start, end, slice };
}

function updateRunCommandCountNote(filters) {
  if (!el.boltzRunRangeCountNote) return;
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    el.boltzRunRangeCountNote.textContent = 'Load BoltzGen configs to preview retained targets.';
    return;
  }
  const range = previewRunCommandRange(targets);
  if (!range) {
    el.boltzRunRangeCountNote.textContent = 'Enter a valid row range to preview retained targets.';
    return;
  }
  const { pass } = applyRunCommandFilters(range.slice, filters);
  el.boltzRunRangeCountNote.textContent = `Retained targets: ${pass.length} of ${range.slice.length} (rows ${range.start}-${range.end}).`;
}

function prependRunCommandFilters(text, filters, totalCount, passCount) {
  const section = runCommandFilterSection(filters, totalCount, passCount);
  if (!section) return text;
  return `${section}\n\n${text || ''}`.trim();
}

async function submitPipelineRerunRangeModal(triggerBtn = null) {
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    showAlert('No BoltzGen configs loaded.');
    return;
  }
  const range = getBoltzRowRange(el.boltzRerunRangeStart, el.boltzRerunRangeEnd, targets, 'pipeline');
  if (!range) return;
  const slice = range.slice;
  if (!slice.length) {
    showAlert('No targets found for that row range.');
    return;
  }
  const expectedOverride = el.boltzRerunRangeExpected?.value || null;
  const attemptsOverride = el.boltzRerunRangeAttempts?.value || null;
  const forceOverride = Boolean(el.boltzRerunRangeForce?.checked);
  const promptOverride = String(el.boltzRerunRangePrompt?.value || '').trim() || null;
  await queuePipelineRerunTargets(slice, triggerBtn, {
    expected: expectedOverride,
    attempts: attemptsOverride,
    force: forceOverride,
    prompt: promptOverride,
  });
  toggleModal(el.boltzRerunRangeModal, false);
}

function openRunCommandRangeModal() {
  if (!el.boltzRunRangeModal) return;
  if (el.boltzRunRangeMinAllowed && !el.boltzRunRangeMinAllowed.value) {
    el.boltzRunRangeMinAllowed.value = '50';
  }
  if (el.boltzRunRangeMinEpitopes && !el.boltzRunRangeMinEpitopes.value) {
    el.boltzRunRangeMinEpitopes.value = '10';
  }
  updateRunCommandFilterNote();
  toggleModal(el.boltzRunRangeModal, true);
}

function openRegenerateRangeModal() {
  if (!el.boltzRegenerateRangeModal) return;
  toggleModal(el.boltzRegenerateRangeModal, true);
}

async function submitRegenerateRangeModal(triggerBtn = null) {
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    showAlert('No BoltzGen configs loaded.');
    return;
  }
  const range = getBoltzRowRange(
    el.boltzRegenerateRangeStart,
    el.boltzRegenerateRangeEnd,
    targets,
    'rebuild',
  );
  if (!range) return;
  const slice = range.slice;
  if (!slice.length) {
    showAlert('No targets found for that row range.');
    return;
  }
  const ids = slice
    .map((row) => (row?.pdb_id || '').trim().toUpperCase())
    .filter(Boolean);
  await regenerateBoltzConfigs({ pdbIds: ids, triggerBtn });
  toggleModal(el.boltzRegenerateRangeModal, false);
}

async function showRunCommandRange() {
  state.runContext = null;
  state.runMode = 'range';
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    showAlert('No targets loaded.');
    return;
  }
  const range = getBoltzRowRange(el.boltzRunRangeStart, el.boltzRunRangeEnd, targets, 'command');
  if (!range) return;
  const { start, end, slice } = range;
  if (!slice.length) {
    showAlert('No targets found for that row range.');
    return;
  }
  const filters = getRunCommandFilters();
  const { pass } = applyRunCommandFilters(slice, filters);
  const filteredSlice = pass;
  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing cluster commands...';
  if (el.boltzRunTitle) {
    const countLabel = filteredSlice.length !== slice.length
      ? ` (${filteredSlice.length}/${slice.length})`
      : '';
    el.boltzRunTitle.textContent = `Manual run · rows ${start}-${end}${countLabel}`;
  }
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }
  const engine = getRunEngine();
  if (engine === 'rfantibody') {
    const entries = await fetchRfaLauncherEntries(filteredSlice);
    const text = filteredSlice.length
      ? buildRfaRunCommandsTextForTargets(entries)
      : '# Manual RFA pipeline submission\nNo targets passed filters.';
    renderRunCommandBlocks(prependRunCommandFilters(text, filters, slice.length, filteredSlice.length));
  } else {
    const text = filteredSlice.length
      ? buildAllRunCommandsText(filteredSlice)
      : '# Manual BoltzGen submission\nNo targets passed filters.';
    renderRunCommandBlocks(prependRunCommandFilters(text, filters, slice.length, filteredSlice.length));
  }
  toggleModal(el.boltzRunRangeModal, false);
}

async function submitPipelineRerun(triggerBtn = null) {
  const bulkIds = el.pipelineRerunModal?.dataset?.bulkIds;
  const pdbId = el.pipelineRerunModal?.dataset?.pdbId;
  const expectedRaw = el.pipelineRerunExpected?.value;
  const attemptsRaw = el.pipelineRerunAttempts?.value;
  const expected = expectedRaw ? Number(expectedRaw) || null : null;
  const attempts = attemptsRaw ? Number(attemptsRaw) || 3 : 3;
  const decideScopePrompt = String(el.pipelineRerunPrompt?.value || '').trim() || null;
  const designCount = getBoltzDesignCount();
  const force = Boolean(el.pipelineRerunForce?.checked);
  if (triggerBtn) triggerBtn.disabled = true;
  try {
    if (bulkIds) {
      const list = bulkIds.split(',').map((id) => id.trim()).filter(Boolean);
      await queuePipelineRerunTargets(list, null, { prompt: decideScopePrompt });
    } else if (pdbId) {
      const row = findBulkRowForPdb(pdbId);
      const fallback = row ? null : findAccessionsFromInput(pdbId);
      const antigenUrl = (el.pipelineRerunModal?.dataset?.antigenUrl || row?.antigen_url || '').trim() || null;
      const payload = {
        pdb_id: pdbId,
        force,
        expected_epitopes: expected,
        decide_scope_attempts: attempts,
        decide_scope_prompt: decideScopePrompt,
        design_count: designCount,
        antigen_url: antigenUrl,
        target_accession: row?.accession || fallback?.accession || null,
        target_vendor_range: row?.vendor_range || fallback?.vendor_range || null,
      };
      const res = await fetch('/api/targets/pipeline/refresh', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });
      if (!res.ok) {
        const detail = await res.json().catch(() => ({}));
        throw new Error(detail.detail || `Pipeline refresh failed (${res.status})`);
      }
      const body = await res.json();
      showAlert(body.message || 'Pipeline refresh queued.', false);
      if (body.job_id) {
        startJobPolling(body.job_id);
      }
    } else {
      showAlert('Missing PDB ID.');
    }
    toggleModal(el.pipelineRerunModal, false);
  } catch (err) {
    showAlert(err.message || 'Failed to rerun pipeline.');
  } finally {
    if (triggerBtn) triggerBtn.disabled = false;
  }
}

async function loadSnapshotMetadata(names) {
  if (!names || !names.length) {
    renderSnapshots([]);
    return;
  }
  const params = new URLSearchParams();
  params.append('names', names.join(','));
  try {
    const res = await fetch(`/api/pymol/snapshots/metadata?${params.toString()}`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to load snapshot metadata (${res.status})`);
    }
    const body = await res.json();
    renderSnapshots(body);
  } catch (err) {
    renderSnapshots([]);
    appendLog(`Snapshot metadata error: ${err.message || err}`);
  }
}

function renderBoltzConfigs() {
  if (!el.boltzPanel || !el.boltzTable) return;
  const tbody = el.boltzTable;
  tbody.innerHTML = '';
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    const emptyRow = document.createElement('tr');
    const emptyCell = document.createElement('td');
    emptyCell.colSpan = 8;
    emptyCell.className = 'empty-note';
    emptyCell.textContent = 'No BoltzGen configs yet. Paste CSV/TSV above to auto-detect targets.';
    emptyRow.appendChild(emptyCell);
    tbody.appendChild(emptyRow);
    if (el.boltzSummary) el.boltzSummary.hidden = true;
    el.boltzPanel.hidden = false;
    return;
  }
  let configCount = 0;
  targets.forEach((target, idx) => {
    const configs = Array.isArray(target.configs) ? target.configs : [];
    configCount += configs.length;
    const key = (target.pdb_id || '').toUpperCase();
    const isExpanded = state.expandedTargets instanceof Set ? state.expandedTargets.has(key) : false;
    const hasPrep = Boolean(target.has_prep);
    const availableLabels = availableEpitopeLabels(configs);
    let epitopeCell = '—';
    let epitopeTitle = '';
    if (configs.length) {
      epitopeCell = `${availableLabels.length} epitope${availableLabels.length === 1 ? '' : 's'}`;
      if (availableLabels.length) {
        epitopeTitle = availableLabels.join(', ');
      }
    }
    const tr = document.createElement('tr');
    tr.className = `target-row ${isExpanded ? 'expanded' : 'collapsed'}`;
    tr.dataset.pdbId = key;
    const allowedLength = asNumber(target.allowed_epitope_length);
    const cells = [
      idx + 1,
      target.preset_name || '—',
      target.pdb_id || '',
      target.allowed_epitope_range || '—',
      allowedLength === null ? '—' : String(Math.round(allowedLength)),
      epitopeCell,
    ];
    cells.forEach((val, cellIdx) => {
      const td = document.createElement('td');
      if (cellIdx === 0) {
        td.innerHTML = `<span class="toggle-icon">▸</span>${val}`;
      } else {
        td.textContent = val;
      }
      if (cellIdx === 5 && epitopeTitle) {
        td.title = epitopeTitle;
      }
      tr.appendChild(td);
    });
    const statusCell = document.createElement('td');
    statusCell.className = 'boltz-status-cell';
    statusCell.appendChild(buildLocalResultBadge(getBoltzDesignTotal(target.pdb_id)));
    tr.appendChild(statusCell);

    const cmdCell = document.createElement('td');
    cmdCell.className = 'boltz-cmd-cell';
    const cfgBtn = document.createElement('button');
    cfgBtn.type = 'button';
    cfgBtn.textContent = 'Config';
    cfgBtn.className = 'ghost';
    cfgBtn.dataset.action = 'show-config-target';
    cfgBtn.dataset.pdbId = target.pdb_id || '';
    cmdCell.appendChild(cfgBtn);

    const debugBtn = document.createElement('button');
    debugBtn.type = 'button';
    debugBtn.textContent = 'Debug';
    debugBtn.className = 'ghost';
    debugBtn.dataset.action = 'show-target-yaml';
    debugBtn.dataset.pdbId = target.pdb_id || '';
    debugBtn.title = 'Show target.yaml';
    cmdCell.appendChild(debugBtn);

    const manualBtn = document.createElement('button');
    manualBtn.type = 'button';
    manualBtn.textContent = 'Command';
    manualBtn.className = 'ghost';
    manualBtn.dataset.action = 'show-run';
    manualBtn.dataset.pdbId = target.pdb_id || '';
    manualBtn.dataset.scope = 'target';
    cmdCell.appendChild(manualBtn);

    const rerunBtn = document.createElement('button');
    rerunBtn.type = 'button';
    rerunBtn.textContent = 'Select epitopes (LLM)';
    rerunBtn.className = 'ghost';
    rerunBtn.dataset.action = 'rerun-pipeline';
    rerunBtn.dataset.pdbId = target.pdb_id || '';
    if (target.antigen_url) rerunBtn.dataset.antigenUrl = target.antigen_url;
    cmdCell.appendChild(rerunBtn);

    const pymolBtn = document.createElement('button');
    pymolBtn.type = 'button';
    pymolBtn.textContent = 'PyMOL';
    pymolBtn.className = 'ghost';
    pymolBtn.dataset.action = 'pymol-target';
    pymolBtn.dataset.pdbId = target.pdb_id || '';
    stylePymolButton(pymolBtn, hasPrep);
    cmdCell.appendChild(pymolBtn);

    tr.appendChild(cmdCell);
    tbody.appendChild(tr);

    configs.forEach((cfg, cfgIdx) => {
      const epRow = document.createElement('tr');
      epRow.className = 'epitope-row';
      epRow.dataset.parent = key;
      epRow.hidden = !isExpanded;
      const epCells = [
        '',
        '',
        '',
        '',
        '',
        surfaceEpitopeLabel(cfg, cfgIdx + 1),
      ];
      epCells.forEach((val, cellIdx) => {
        const td = document.createElement('td');
        td.textContent = val || '';
        if (cellIdx === 5) td.style.paddingLeft = '14px';
        epRow.appendChild(td);
      });
      if (cfg?.hotspot_surface_ok === false) {
        epRow.style.opacity = '0.65';
        epRow.title = 'Filtered: hotspot below SASA cutoff.';
      }
      const epStatus = document.createElement('td');
      epStatus.className = 'boltz-status-cell';
      epStatus.appendChild(buildLocalResultBadge(getBoltzDesignTotalForSpec(target.pdb_id, cfg)));
      epRow.appendChild(epStatus);

      const epCmd = document.createElement('td');
      epCmd.className = 'boltz-cmd-cell';
      const epCfgBtn = document.createElement('button');
      epCfgBtn.type = 'button';
      epCfgBtn.textContent = 'Config';
      epCfgBtn.className = 'ghost';
      epCfgBtn.dataset.action = 'show-config-epitope';
      epCfgBtn.dataset.pdbId = target.pdb_id || '';
      epCfgBtn.dataset.configPath = cfg.config_path || '';
      epCfgBtn.dataset.epitopeName = epitopeLabel(cfg, cfgIdx + 1);
      epCmd.appendChild(epCfgBtn);

      const epDebug = document.createElement('button');
      epDebug.type = 'button';
      epDebug.textContent = 'Debug';
      epDebug.className = 'ghost';
      epDebug.dataset.action = 'show-target-yaml';
      epDebug.dataset.pdbId = target.pdb_id || '';
      epDebug.title = 'Show target.yaml';
      epCmd.appendChild(epDebug);

      const epManual = document.createElement('button');
      epManual.type = 'button';
      epManual.textContent = 'Command';
      epManual.className = 'ghost';
      epManual.dataset.action = 'show-run';
      epManual.dataset.pdbId = target.pdb_id || '';
      epManual.dataset.configPath = cfg.config_path || '';
      epManual.dataset.epitopeName = epitopeLabel(cfg, cfgIdx + 1);
      if (cfg?.hotspot_surface_ok === false) {
        epManual.disabled = true;
        epManual.title = 'Filtered by SASA: hotspot not surface exposed.';
      }
      epCmd.appendChild(epManual);

      const epPymol = document.createElement('button');
      epPymol.type = 'button';
      epPymol.textContent = 'PyMOL';
      epPymol.className = 'ghost';
      epPymol.dataset.action = 'pymol-epitope';
      epPymol.dataset.pdbId = target.pdb_id || '';
      epPymol.dataset.epitopeName = epitopeLabel(cfg, cfgIdx + 1);
      stylePymolButton(epPymol, hasPrep);
      epCmd.appendChild(epPymol);

      epRow.appendChild(epCmd);
      tbody.appendChild(epRow);
    });
  });

  if (el.boltzSummary) {
    el.boltzSummary.textContent = `${configCount} config${configCount === 1 ? '' : 's'} across ${targets.length} target${targets.length === 1 ? '' : 's'}`;
    el.boltzSummary.hidden = false;
  }
  el.boltzPanel.hidden = false;
}

async function loadClusterStatus() {
  try {
    const res = await fetch('/api/cluster/status');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    state.clusterStatus = await res.json();
  } catch (err) {
    state.clusterStatus = null;
  }
}

async function loadBoltzRunStatuses(pdbIds = []) {
  const ids = Array.from(new Set(
    (pdbIds || []).map((id) => (id || '').toUpperCase()).filter(Boolean),
  ));
  const results = {};
  const resultsBySpec = {};
  if (!ids.length) {
    state.boltzLocalRuns = results;
    state.boltzLocalRunsBySpec = resultsBySpec;
    return results;
  }
  const counts = state.boltzDesignCounts || {};
  if (Object.keys(counts).length) {
    ids.forEach((id) => {
      const num = Number(counts[id]);
      results[id] = Number.isFinite(num) ? num : 0;
    });
    state.boltzLocalRuns = results;
    state.boltzLocalRunsBySpec = resultsBySpec;
    return results;
  }
  await Promise.all(ids.map(async (id) => {
    try {
      const res = await fetch(`/api/targets/${id}/boltzgen/runs`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const body = await res.json();
      const runs = Array.isArray(body.runs) ? body.runs : [];
      const specCounts = {};
      let totalDesigns = 0;
      let sawCount = false;
      let sawRun = false;
      runs.forEach((run) => {
        sawRun = true;
        const specs = Array.isArray(run.specs) ? run.specs : [];
        specs.forEach((spec) => {
          const specName = String(spec?.name || '').trim();
          const count = Number(spec?.design_count);
          if (!Number.isFinite(count)) return;
          if (specName) {
            const specKey = specName.toUpperCase();
            specCounts[specKey] = Number(specCounts[specKey] || 0) + count;
          }
          totalDesigns += count;
          sawCount = true;
        });
      });
      if (!sawRun) {
        results[id] = 0;
        resultsBySpec[id] = {};
      } else if (sawCount) {
        results[id] = totalDesigns;
        resultsBySpec[id] = specCounts;
      } else {
        results[id] = null;
        resultsBySpec[id] = {};
      }
    } catch (err) {
      results[id] = null;
      resultsBySpec[id] = {};
    }
  }));
  state.boltzLocalRuns = results;
  state.boltzLocalRunsBySpec = resultsBySpec;
  return results;
}

async function loadBoltzConfigs(options = {}) {
  const { silent = false } = options;
  if (!el.boltzPanel || !el.boltzTable) return;
  const ids = Array.from(new Set(
    (state.bulkPreviewRows || [])
      .map((row) => (row.resolved_pdb_id || row.pdb_id || '').trim().toUpperCase())
      .filter(Boolean),
  ));
  if (!ids.length) {
    state.boltzConfigs = [];
    renderBoltzConfigs();
    return;
  }
  try {
    const params = new URLSearchParams();
    params.append('pdb_ids', ids.join(','));
    const res = await fetch(`/api/bulk/boltzgen/configs?${params.toString()}`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to load BoltzGen configs (${res.status})`);
    }
    const body = await res.json();
    const targets = Array.isArray(body.targets) ? body.targets : [];
    const presetMap = new Map();
    (state.bulkPreviewRows || []).forEach((row) => {
      const key = (row.resolved_pdb_id || '').toUpperCase();
      if (key) presetMap.set(key, row.preset_name || row.gene || '');
    });
    targets.forEach((target) => {
      const key = (target.pdb_id || '').toUpperCase();
      if (presetMap.has(key)) {
        target.preset_name = presetMap.get(key);
      }
    });
    state.expandedTargets = new Set();
    state.boltzConfigs = targets;
    await loadBoltzRunStatuses(targets.map((t) => t.pdb_id));
    renderBoltzConfigs();
    await loadBinderTable({ page: state.binderPage || 1, silent: true });
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
    renderBoltzConfigs();
  }
}

async function regenerateBoltzConfigs(options = {}) {
  const { pdbIds = null, triggerBtn = null } = options;
  const candidates = Array.isArray(state.boltzConfigs) && state.boltzConfigs.length
    ? state.boltzConfigs
    : (state.bulkPreviewRows || []);
  const ids = Array.from(new Set(
    (pdbIds && pdbIds.length ? pdbIds : candidates)
      .map((row) => (row?.pdb_id || row?.resolved_pdb_id || row || '').trim().toUpperCase())
      .filter(Boolean),
  ));
  if (!ids.length) {
    showAlert('No BoltzGen configs loaded.');
    return;
  }
  const designCount = getBoltzDesignCount() || 100;
  const cropRadius = getBoltzCropRadius();
  const activeBtn = triggerBtn || el.boltzRegenerate;
  if (activeBtn) activeBtn.disabled = true;
  try {
    const res = await fetch('/api/bulk/boltzgen/configs/regenerate', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        pdb_ids: ids,
        design_count: designCount,
        boltzgen_crop_radius: cropRadius,
      }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to rebuild configs (${res.status})`);
    }
    const body = await res.json();
    const results = Array.isArray(body.results) ? body.results : [];
    const ok = results.filter((row) => row && row.status === 'ok');
    const skipped = results.filter((row) => row && row.status === 'skipped');
    const failed = results.filter((row) => row && row.status === 'error');
    let message = `Rebuilt configs for ${ok.length} target${ok.length === 1 ? '' : 's'}.`;
    if (skipped.length) message += ` Skipped ${skipped.length}.`;
    if (failed.length) message += ` Failed ${failed.length}.`;
    showAlert(message, failed.length > 0);
    if (ok.length) {
      await loadBoltzConfigs({ silent: true });
    }
  } catch (err) {
    showAlert(err.message || 'Failed to rebuild BoltzGen configs.');
  } finally {
    if (activeBtn) activeBtn.disabled = false;
  }
}

async function plotBoltzAntigenDiversity() {
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  const ids = Array.from(new Set(targets.map((t) => (t?.pdb_id || '').trim().toUpperCase()).filter(Boolean)));
  if (!ids.length) {
    showAlert('No BoltzGen configs loaded.');
    return;
  }
  const btn = el.boltzPlotDiversity;
  if (btn) btn.disabled = true;
  try {
    const res = await fetch('/api/bulk/boltzgen/antigen-diversity', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ pdb_ids: ids }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to plot diversity (${res.status})`);
    }
    const body = await res.json();
    const count = Array.isArray(body.plots) ? body.plots.length : 0;
    const outDir = body.output_dir ? ` Output dir: ${body.output_dir}` : '';
    const message = body.message || `Generated ${count} antigen diversity plot${count === 1 ? '' : 's'}.${outDir}`;
    showAlert(message, false);
  } catch (err) {
    showAlert(err.message || 'Failed to plot diversity.');
  } finally {
    if (btn) btn.disabled = false;
  }
}

async function plotEpitopeDiversitySelection() {
  const raw = (el.epitopeDiversitySelection?.value || '').trim();
  if (!raw) {
    showAlert('Enter selections like 5WT9:epitope_1, 3J8F:epitope_4.');
    return;
  }
  const selections = raw
    .split(/[\n,]+/)
    .map((part) => part.trim())
    .filter(Boolean);
  if (!selections.length) {
    showAlert('No valid selections found.');
    return;
  }
  const btn = el.epitopeDiversityPlot;
  if (btn) btn.disabled = true;
  try {
    const res = await fetch('/api/bulk/boltzgen/epitope-diversity', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ selections }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to plot epitope diversity (${res.status})`);
    }
    const body = await res.json();
    const plots = Array.isArray(body.plots) ? body.plots : [];
    state.epitopeCsvName = body.csv_name || null;
    state.epitopeHotspotCsvName = body.hotspot_csv_name || null;
    const items = plots
      .map((plot) => {
        const src = bulkFileSrc(plot.png_name || plot.png_path || plot.svg_name || plot.svg_path);
        const label = plot.title || plot.png_name || plot.svg_name || 'Epitope diversity';
        if (!src) return null;
        return { src, label };
      })
      .filter(Boolean);
    renderEpitopePlots(items);
    const message = body.message || `Generated ${items.length} epitope diversity plot(s).`;
    showAlert(message, false);
  } catch (err) {
    showAlert(err.message || 'Failed to plot epitope diversity.');
  } finally {
    if (btn) btn.disabled = false;
  }
}

async function fetchBoltzConfig(pdbId, configPath) {
  const params = new URLSearchParams();
  params.append('pdb_id', pdbId);
  params.append('config_path', configPath);
  const res = await fetch(`/api/bulk/boltzgen/config?${params.toString()}`);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || `Unable to load config (${res.status})`);
  }
  return res.json();
}

async function fetchTargetYaml(pdbId) {
  const res = await fetch(`/api/targets/${encodeURIComponent(pdbId)}/target-yaml`);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || `Unable to load target.yaml (${res.status})`);
  }
  return res.json();
}

function getConfigEngine() {
  const value = (el.boltzConfigEngine?.value || 'boltzgen').trim().toLowerCase();
  return value || 'boltzgen';
}

function getRunEngine() {
  const control = el.boltzRunEngine;
  if (!control) return 'boltzgen';
  if (control.tagName === 'SELECT') {
    const value = (control.value || 'boltzgen').trim().toLowerCase();
    return value || 'boltzgen';
  }
  const active = control.querySelector('button.active');
  const value = active?.dataset?.engine || control.querySelector('button[data-engine]')?.dataset?.engine || 'boltzgen';
  return (value || 'boltzgen').trim().toLowerCase();
}

function setRunEngine(engine) {
  const control = el.boltzRunEngine;
  const value = (engine || 'boltzgen').trim().toLowerCase() || 'boltzgen';
  if (!control) return value;
  if (control.tagName === 'SELECT') {
    control.value = value;
    return value;
  }
  const buttons = Array.from(control.querySelectorAll('button[data-engine]'));
  let matched = false;
  buttons.forEach((btn) => {
    const isActive = (btn.dataset.engine || '').trim().toLowerCase() === value;
    btn.classList.toggle('active', isActive);
    btn.setAttribute('aria-selected', isActive ? 'true' : 'false');
    if (isActive) matched = true;
  });
  if (!matched && buttons.length) {
    buttons.forEach((btn, idx) => {
      const isActive = idx === 0;
      btn.classList.toggle('active', isActive);
      btn.setAttribute('aria-selected', isActive ? 'true' : 'false');
    });
    return (buttons[0].dataset.engine || value).trim().toLowerCase();
  }
  return value;
}

function sanitizeRfaEpitopeName(name) {
  return String(name || '').trim().replace(/[ /]+/g, '_');
}

async function fetchRfaConfigs(pdbId, options = {}) {
  const { force = false } = options;
  const key = (pdbId || '').toUpperCase();
  if (!key) return null;
  if (!force && state.rfaConfigs && state.rfaConfigs[key]) {
    return state.rfaConfigs[key];
  }
  const params = new URLSearchParams();
  params.append('pdb_ids', key);
  const res = await fetch(`/api/bulk/rfa/configs?${params.toString()}`);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || `Unable to load RFA configs (${res.status})`);
  }
  const payload = await res.json();
  const targets = Array.isArray(payload.targets) ? payload.targets : [];
  const target = targets.find((item) => (item?.pdb_id || '').toUpperCase() === key) || {
    pdb_id: key,
    scripts: [],
  };
  if (!state.rfaConfigs) state.rfaConfigs = {};
  state.rfaConfigs[key] = target;
  return target;
}

async function fetchRfaScript(pdbId, scriptPath) {
  const params = new URLSearchParams();
  params.append('pdb_id', pdbId);
  params.append('script_path', scriptPath);
  const res = await fetch(`/api/bulk/rfa/config?${params.toString()}`);
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail || `Unable to load script (${res.status})`);
  }
  return res.json();
}

async function showRfaConfig(pdbId, epitopeName = null) {
  if (!el.boltzConfigBody || !el.boltzConfigTitle) return;
  try {
    const target = await fetchRfaConfigs(pdbId, { force: true });
    const scripts = Array.isArray(target?.scripts) ? target.scripts : [];
    const epKey = epitopeName ? sanitizeRfaEpitopeName(epitopeName) : '';
    const filtered = epKey
      ? scripts.filter((script) => script.stage !== 'launcher' && script.epitope === epKey)
      : scripts;
    if (!filtered.length) {
      el.boltzConfigTitle.textContent = epKey
        ? `RFA pipeline scripts · ${pdbId} · ${epitopeName}`
        : `RFA pipeline scripts · ${pdbId}`;
      el.boltzConfigBody.textContent = 'No RFA pipeline scripts found.';
      toggleModal(el.boltzConfigModal, true);
      return;
    }
    const texts = [];
    for (const script of filtered) {
      try {
        const content = await fetchRfaScript(pdbId, script.path);
        const label = script.stage ? `${script.stage}` : 'script';
        texts.push(`# ${label} · ${content.script_name}\n${content.script_text || ''}`);
      } catch (err) {
        const label = script.name || script.path || 'script';
        texts.push(`# ${label}\nError loading script: ${err.message || err}`);
      }
    }
    el.boltzConfigTitle.textContent = epKey
      ? `RFA pipeline scripts · ${pdbId} · ${epitopeName}`
      : `RFA pipeline scripts · ${pdbId}`;
    el.boltzConfigBody.textContent = texts.join('\n\n');
    toggleModal(el.boltzConfigModal, true);
  } catch (err) {
    showAlert(err.message || String(err));
  }
}

async function renderConfigModal() {
  const ctx = state.configContext;
  if (!ctx || !ctx.pdbId) return;
  const engine = getConfigEngine();
  if (engine === 'rfantibody') {
    await showRfaConfig(ctx.pdbId, ctx.epitopeName || null);
  } else {
    await showBoltzConfig(ctx.pdbId, ctx.configPath || null);
  }
}

function setConfigContext(pdbId, configPath = null, epitopeName = null) {
  state.configContext = {
    pdbId: pdbId || '',
    configPath: configPath || null,
    epitopeName: epitopeName || null,
  };
}

async function showBoltzConfig(pdbId, configPath = null) {
  if (!el.boltzConfigBody || !el.boltzConfigTitle) return;
  try {
    if (configPath) {
      const cfg = await fetchBoltzConfig(pdbId, configPath);
      el.boltzConfigTitle.textContent = `BoltzGen config · ${cfg.epitope_name || configPath}`;
      el.boltzConfigBody.textContent = cfg.yaml_text || 'No config content found.';
    } else {
      const target = (state.boltzConfigs || []).find((t) => (t.pdb_id || '').toUpperCase() === (pdbId || '').toUpperCase());
      const configs = Array.isArray(target?.configs) ? target.configs : [];
      if (!configs.length) throw new Error('No configs detected for this target.');
      const texts = [];
      for (const cfg of configs) {
        try {
          const content = await fetchBoltzConfig(pdbId, cfg.config_path);
          texts.push(`# ${content.epitope_name || cfg.config_path}\n${content.yaml_text || ''}`);
        } catch (err) {
          texts.push(`# ${epitopeLabel(cfg)}\nError loading config: ${err.message || err}`);
        }
      }
      el.boltzConfigTitle.textContent = `BoltzGen configs · ${pdbId}`;
      el.boltzConfigBody.textContent = texts.join('\n\n');
    }
    toggleModal(el.boltzConfigModal, true);
  } catch (err) {
    showAlert(err.message || String(err));
  }
}

async function showTargetYaml(pdbId) {
  if (!el.boltzConfigBody || !el.boltzConfigTitle) return;
  try {
    state.configContext = null;
    const payload = await fetchTargetYaml(pdbId);
    el.boltzConfigTitle.textContent = `target.yaml · ${payload.pdb_id || pdbId}`;
    el.boltzConfigBody.textContent = payload.yaml_text || 'No target.yaml content found.';
    toggleModal(el.boltzConfigModal, true);
  } catch (err) {
    showAlert(err.message || String(err));
  }
}

async function showBoltzLog(jobId, title = 'Job log') {
  if (!jobId) {
    showAlert('No job ID available yet.');
    return;
  }
  try {
    const job = await fetchJob(jobId);
    if (el.boltzLogTitle) el.boltzLogTitle.textContent = title;
    if (el.boltzLogBody) {
      el.boltzLogBody.textContent = (job.logs && job.logs.length)
        ? job.logs.join('\n')
        : 'No logs available yet.';
    }
    toggleModal(el.boltzLogModal, true);
  } catch (err) {
    showAlert(err.message || String(err));
  }
}

function clusterDefaults() {
  const info = state.clusterStatus || {};
  const boltz = info.boltzgen || {};
  const baseRootRaw = (info.target_root || info.remote_root || '<remote_root>').toString().replace(/\/$/, '');
  const workspaceRoot = baseRootRaw.toLowerCase().endsWith('/targets')
    ? baseRootRaw.replace(/\/+targets$/i, '')
    : baseRootRaw;
  const targetBaseRaw = (info.target_root || `${workspaceRoot}/targets`).toString().replace(/\/$/, '');
  const lastSegment = targetBaseRaw.split('/').filter(Boolean).pop() || '';
  const targetRoot = lastSegment.toLowerCase() === 'targets' ? targetBaseRaw : `${targetBaseRaw}/targets`;
  const toolsRoot = workspaceRoot ? `${workspaceRoot}/lib` : '<remote_root>/lib';
  const localRootRaw = (info.local_root || '').toString().replace(/\/$/, '');
  const localWorkspaceRoot = localRootRaw || '';
  const localTargetsRoot = localWorkspaceRoot
    ? (localWorkspaceRoot.toLowerCase().endsWith('/targets') ? localWorkspaceRoot : `${localWorkspaceRoot}/targets`)
    : '';
  return {
    remoteRoot: workspaceRoot,
    targetRoot,
    toolsRoot,
    localWorkspaceRoot,
    localTargetsRoot,
    sshTarget: info.ssh_target || '<cluster>',
    condaActivate: boltz.conda_activate || info.conda_activate || '',
    boltz,
  };
}

function runLabelFor(pdbId) {
  const prefix = (el.bulkRunPrefix?.value || 'boltz').trim() || 'boltz';
  const cleanedPrefix = prefix.replace(/\s+/g, '_');
  const now = new Date();
  const stamp = `${now.getFullYear()}${String(now.getMonth() + 1).padStart(2, '0')}${String(now.getDate()).padStart(2, '0')}_${String(now.getHours()).padStart(2, '0')}${String(now.getMinutes()).padStart(2, '0')}`;
  const upper = (pdbId || '').toUpperCase();
  return `${cleanedPrefix}_${upper || 'target'}_${stamp}`;
}

function outputRootFor(baseRoot, runLabel) {
  const cleaned = `${baseRoot || ''}`.replace(/\/$/, '');
  if (!cleaned || !runLabel) return cleaned;
  const parts = cleaned.split('/').filter(Boolean);
  const last = parts[parts.length - 1] || '';
  if (last === runLabel) return cleaned;
  return `${cleaned}/${runLabel}`;
}

function buildRunCommandText(pdbId, specPaths = [], epitopeName = null) {
  const upper = (pdbId || '').toUpperCase();
  const {
    remoteRoot,
    targetRoot,
    toolsRoot,
    sshTarget,
    condaActivate,
    boltz,
    localTargetsRoot,
    localWorkspaceRoot,
  } = clusterDefaults();
  const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${upper}`.replace(/\/+/g, '/');
  const localTarget = localTargetsRoot ? `${localTargetsRoot}/${upper}` : `targets/${upper}`;
  const localTools = localWorkspaceRoot ? `${localWorkspaceRoot}/lib/tools` : 'lib/tools';
  const runLabel = runLabelFor(upper);
  const designCount = getBoltzDesignCount() || 100;
  const timeHours = getBoltzTimeHours() || boltz.time_hours || 12;
  const partition = boltz.partition || '<partition>';
  const gpus = boltz.gpus || '<gpus>';
  const cpus = boltz.cpus || '<cpus>';
  const memVal = boltz.mem_gb ?? '<mem>';
  const memLabel = Number.isFinite(memVal) ? `${memVal}G` : `${memVal}`;
  const outputRootBase = boltz.output_root
    ? `${boltz.output_root}`.replace(/\/$/, '')
    : `${remoteTarget}/designs/boltzgen`;
  const outputRoot = outputRootFor(outputRootBase, runLabel);
  const localOutputDir = `${localTarget}/designs/boltzgen/${runLabel}`;
  const cacheDir = boltz.cache_dir ? `${boltz.cache_dir}` : null;
  const extraArgs = Array.isArray(boltz.extra_args) ? boltz.extra_args : (boltz.extra_args ? [boltz.extra_args] : []);
  const specs = (specPaths && specPaths.length) ? specPaths : ['configs/*/boltzgen_config.yaml'];
  const remoteSpecs = specs.map((spec) => `${remoteTarget}/${spec}`.replace(/\/+/g, '/'));
  const specLines = remoteSpecs.map((spec) => `  --spec ${spec} \\`).join('\n');
  const cacheLine = cacheDir ? `  --cache_dir ${cacheDir} \\` : null;
  const extraLine = extraArgs.length ? `  --extra_run_args ${extraArgs.join(' ')} \\` : null;
  const condaLine = condaActivate ? condaActivate : '# conda activate bg';
  const verifyLine = `ssh ${sshTarget} "ls ${remoteSpecs.join(' ')}"`;
  const envRoot = remoteRoot || '<remote_root>';
  const envTargetRoot = targetRoot || `${envRoot}/targets`;
  const pipelinePath = `${toolsRoot || '<remote_root>/lib'}/tools/boltzgen/pipeline.py`;
  return [
    `# Manual BoltzGen submission${epitopeName ? ` (${epitopeName})` : ''}`,
    '',
    '# 1) Sync target configs + tools to cluster',
    `rsync -az ${localTarget} ${sshTarget}:${targetRoot}`,
    `rsync -az ${localTools} ${sshTarget}:${toolsRoot || '<remote_root>/lib'}`,
    '',
    '# 2) Verify boltzgen configs exist on cluster',
    verifyLine,
    '',
    '# 3) Launch BoltzGen pipeline',
    condaLine,
    `INITBINDER_ROOT=${envRoot} INITBINDER_TARGET_ROOT=${envTargetRoot} \\`,
    `python ${pipelinePath} pipeline ${upper} \\`,
    `  --run_label ${runLabel} \\`,
    `  --num_designs ${designCount} \\`,
    `  --output_root ${outputRoot} \\`,
    `  --partition ${partition} --gpus ${gpus} --cpus ${cpus} --mem ${memLabel} --time_h ${timeHours} \\`,
    specLines,
    cacheLine,
    extraLine,
    '  --submit',
    '',
    '# 4) Monitor',
    `ssh ${sshTarget} "squeue -u $USER | grep ${runLabel}"`,
    `ssh ${sshTarget} "tail -f ${outputRoot}/launcher.log"`,
    '',
    '# 5) Pull results back',
    `mkdir -p ${localOutputDir}`,
    `rsync -az ${sshTarget}:${outputRoot}/ ${localOutputDir}/`,
  ].filter(Boolean).join('\n');
}

function buildRfaRunCommandText(pdbId, launcherPath) {
  const upper = (pdbId || '').toUpperCase();
  const {
    remoteRoot, targetRoot, sshTarget, localTargetsRoot, localWorkspaceRoot,
  } = clusterDefaults();
  const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${upper}`.replace(/\/+/g, '/');
  const localTarget = localTargetsRoot ? `${localTargetsRoot}/${upper}` : `targets/${upper}`;
  const remoteRootClean = remoteRoot ? remoteRoot.replace(/\/$/, '') : '<remote_root>';
  const remoteLauncher = remoteRoot
    ? `${remoteRootClean}/${launcherPath}`.replace(/\/+/g, '/')
    : launcherPath;
  const localTools = localWorkspaceRoot ? `${localWorkspaceRoot}/tools` : 'tools';

  return [
    '# Manual RFA pipeline submission',
    '',
    '# 1) Sync target + launcher scripts to cluster',
    `rsync -az ${localTarget} ${sshTarget}:${targetRoot}`,
    `rsync -az ${localTools} ${sshTarget}:${remoteRootClean}`,
    '',
    '# 2) Verify launcher exists',
    `ssh ${sshTarget} "ls ${remoteLauncher}"`,
    '',
    '# 3) Launch pipeline',
    `ssh ${sshTarget} "cd ${remoteRootClean} && bash ${remoteLauncher}"`,
    '',
    '# 4) Monitor',
    `ssh ${sshTarget} "squeue -u $USER | grep rfa-"`,
    `ssh ${sshTarget} "ls ${remoteRootClean}/slurm_logs | tail -n 5"`,
    '',
    '# 5) Pull results back (optional)',
    `mkdir -p ${localTarget}/designs/rfantibody/`,
    `rsync -az ${sshTarget}:${remoteTarget}/designs/rfantibody/ ${localTarget}/designs`,
  ].filter(Boolean).join('\n');
}

function buildRfaRunCommandsTextForTargets(entries = []) {
  const {
    remoteRoot, targetRoot, sshTarget, localTargetsRoot, localWorkspaceRoot,
  } = clusterDefaults();
  const remoteRootClean = remoteRoot ? remoteRoot.replace(/\/$/, '') : '<remote_root>';
  const localTools = localWorkspaceRoot ? `${localWorkspaceRoot}/tools` : 'tools';

  const lines = [
    '# Manual RFA pipeline submission (all targets)',
    '',
    '# 1) Sync targets + tools to cluster',
  ];

  entries.forEach((entry) => {
    const pdb = (entry?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    lines.push(`rsync -az ${localTarget} ${sshTarget}:${targetRoot}`);
  });
  lines.push(`rsync -az ${localTools} ${sshTarget}:${remoteRootClean}`);

  lines.push('', '# 2) Verify launchers exist');
  entries.forEach((entry) => {
    const pdb = (entry?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const launcherPath = entry?.launcherPath || '';
    if (!launcherPath) {
      lines.push(`# ${pdb}: launcher missing (rebuild configs)`);
      return;
    }
    const remoteLauncher = remoteRoot
      ? `${remoteRootClean}/${launcherPath}`.replace(/\/+/g, '/')
      : launcherPath;
    lines.push(`ssh ${sshTarget} "ls ${remoteLauncher}"`);
  });

  lines.push('', '# 3) Launch pipelines');
  entries.forEach((entry) => {
    const pdb = (entry?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const launcherPath = entry?.launcherPath || '';
    if (!launcherPath) {
      lines.push(`# ${pdb}: launcher missing (rebuild configs)`);
      return;
    }
    const remoteLauncher = remoteRoot
      ? `${remoteRootClean}/${launcherPath}`.replace(/\/+/g, '/')
      : launcherPath;
    lines.push(`ssh ${sshTarget} "cd ${remoteRootClean} && bash ${remoteLauncher}"`);
  });

  lines.push('', '# 4) Monitor');
  lines.push(`ssh ${sshTarget} "squeue -u $USER | grep rfa-"`);
  lines.push(`ssh ${sshTarget} "ls ${remoteRootClean}/slurm_logs | tail -n 5"`);

  lines.push('', '# 5) Pull results back (optional)');
  entries.forEach((entry) => {
    const pdb = (entry?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    lines.push(`mkdir -p ${localTarget}/designs/rfantibody/`);
    lines.push(`rsync -az ${sshTarget}:${remoteTarget}/designs/rfantibody/ ${localTarget}/designs`);
  });

  return lines.filter(Boolean).join('\n');
}

async function fetchRfaLauncherEntries(targets = []) {
  const entries = [];
  for (const target of targets) {
    const pdbId = (target?.pdb_id || '').toUpperCase();
    if (!pdbId) continue;
    try {
      const rfaTarget = await fetchRfaConfigs(pdbId, { force: true });
      entries.push({ pdbId, launcherPath: rfaTarget?.launcher_path || '' });
    } catch (err) {
      entries.push({ pdbId, launcherPath: '' });
    }
  }
  return entries;
}

function buildAllRunCommandsText(targets = []) {
  const {
    remoteRoot,
    targetRoot,
    toolsRoot,
    sshTarget,
    condaActivate,
    boltz,
    localTargetsRoot,
    localWorkspaceRoot,
  } = clusterDefaults();
  const envRoot = remoteRoot || '<remote_root>';
  const envTargetRoot = targetRoot || `${envRoot}/targets`;
  const pipelinePath = `${toolsRoot || '<remote_root>/lib'}/tools/boltzgen/pipeline.py`;
  const condaLine = condaActivate ? condaActivate : '# conda activate bg';
  const localTools = localWorkspaceRoot ? `${localWorkspaceRoot}/lib/tools` : 'lib/tools';

  const designCount = getBoltzDesignCount() || 100;
  const timeHours = getBoltzTimeHours() || boltz.time_hours || 12;
  const partition = boltz.partition || '<partition>';
  const gpus = boltz.gpus || '<gpus>';
  const cpus = boltz.cpus || '<cpus>';
  const memVal = boltz.mem_gb ?? '<mem>';
  const memLabel = Number.isFinite(memVal) ? `${memVal}G` : `${memVal}`;
  const cacheDir = boltz.cache_dir ? `${boltz.cache_dir}` : null;
  const extraArgs = Array.isArray(boltz.extra_args) ? boltz.extra_args : (boltz.extra_args ? [boltz.extra_args] : []);

  const lines = [
    '# Manual BoltzGen submission (all targets)',
    '',
    '# 1) Sync target configs + tools to cluster',
  ];
  targets.forEach((t) => {
    const pdb = (t?.pdb_id || '').toUpperCase();
    if (!pdb) return;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    lines.push(`rsync -az ${localTarget} ${sshTarget}:${targetRoot}`);
  });
  lines.push(`rsync -az ${localTools} ${sshTarget}:${toolsRoot || '<remote_root>/lib'}`);

  lines.push('', '# 2) Verify boltzgen configs exist on cluster');
  targets.forEach((t) => {
    const pdb = (t?.pdb_id || '').toUpperCase();
    if (!pdb) return;
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const hasConfigs = Array.isArray(t?.configs) && t.configs.length > 0;
    const allowedConfigs = hasConfigs ? t.configs.filter(isSurfaceAllowed) : [];
    const specs = allowedConfigs.map((cfg) => cfg.config_path).filter(Boolean);
    const remoteSpecs = specs.length
      ? specs
      : (hasConfigs ? [] : ['configs/*/boltzgen_config.yaml']);
    if (!remoteSpecs.length) {
      lines.push(`# ${pdb}: no surface-exposed epitopes (SASA filtered)`);
      return;
    }
    lines.push(`ssh ${sshTarget} "ls ${remoteSpecs.map((spec) => `${remoteTarget}/${spec}`.replace(/\/+/g, '/')).join(' ')}"`);
  });

  lines.push('', '# 3) Launch BoltzGen pipeline', condaLine);
  lines.push(`export INITBINDER_ROOT=${envRoot}`);
  lines.push(`export INITBINDER_TARGET_ROOT=${envTargetRoot}`);

  const monitorLines = [];
  const pullLines = [];
  targets.forEach((t) => {
    const pdb = (t?.pdb_id || '').toUpperCase();
    if (!pdb) return;
    const runLabel = runLabelFor(pdb);
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const outputRootBase = boltz.output_root
      ? `${boltz.output_root}`.replace(/\/$/, '')
      : `${remoteTarget}/designs/boltzgen`;
    const outputRoot = outputRootFor(outputRootBase, runLabel);
    const hasConfigs = Array.isArray(t?.configs) && t.configs.length > 0;
    const allowedConfigs = hasConfigs ? t.configs.filter(isSurfaceAllowed) : [];
    const specs = allowedConfigs.map((cfg) => cfg.config_path).filter(Boolean);
    const remoteSpecs = specs.length
      ? specs
      : (hasConfigs ? [] : ['configs/*/boltzgen_config.yaml']);
    if (!remoteSpecs.length) {
      lines.push('', `# ${pdb}: no surface-exposed epitopes (SASA filtered)`);
      return;
    }
    const specLines = remoteSpecs.map((spec) => `  --spec ${spec} \\`);
    const cacheLine = cacheDir ? `  --cache_dir ${cacheDir} \\` : null;
    const extraLine = extraArgs.length ? `  --extra_run_args ${extraArgs.join(' ')} \\` : null;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    const localOutputDir = `${localTarget}/designs/boltzgen/${runLabel}`;

    lines.push(
      '',
      `python ${pipelinePath} pipeline ${pdb} \\`,
      `  --run_label ${runLabel} \\`,
      `  --num_designs ${designCount} \\`,
      `  --output_root ${outputRoot} \\`,
      `  --partition ${partition} --gpus ${gpus} --cpus ${cpus} --mem ${memLabel} --time_h ${timeHours} \\`,
      ...specLines,
      cacheLine,
      extraLine,
      '  --submit',
    );

    monitorLines.push(`ssh ${sshTarget} "squeue -u $USER | grep ${runLabel}"`);
    monitorLines.push(`ssh ${sshTarget} "tail -f ${outputRoot}/launcher.log"`);
    pullLines.push(`mkdir -p ${localOutputDir}`);
    pullLines.push(`rsync -az ${sshTarget}:${outputRoot}/ ${localOutputDir}/`);
  });

  if (monitorLines.length) {
    lines.push('', '# Monitor', ...monitorLines);
  }
  if (pullLines.length) {
    lines.push('', '# Pull results back', ...pullLines);
  }

  return lines.filter(Boolean).join('\n');
}

function buildSelectedRunCommandsText(selections = []) {
  const {
    remoteRoot,
    targetRoot,
    toolsRoot,
    sshTarget,
    condaActivate,
    boltz,
    localTargetsRoot,
    localWorkspaceRoot,
  } = clusterDefaults();
  const envRoot = remoteRoot || '<remote_root>';
  const envTargetRoot = targetRoot || `${envRoot}/targets`;
  const pipelinePath = `${toolsRoot || '<remote_root>/lib'}/tools/boltzgen/pipeline.py`;
  const condaLine = condaActivate ? condaActivate : '# conda activate bg';
  const localTools = localWorkspaceRoot ? `${localWorkspaceRoot}/lib/tools` : 'lib/tools';

  const grouped = new Map();
  selections.forEach((sel) => {
    const pdb = (sel?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const entry = grouped.get(pdb) || { pdbId: pdb, specs: new Set(), epitopeNames: new Set() };
    if (sel?.configPath) entry.specs.add(sel.configPath);
    if (sel?.epitopeName) entry.epitopeNames.add(sel.epitopeName);
    grouped.set(pdb, entry);
  });
  const targets = Array.from(grouped.values());
  if (!targets.length) {
    return '# Manual BoltzGen submission (selection)\nNo matching targets.';
  }

  const designCount = getBoltzDesignCount() || 100;
  const timeHours = getBoltzTimeHours() || boltz.time_hours || 12;
  const partition = boltz.partition || '<partition>';
  const gpus = boltz.gpus || '<gpus>';
  const cpus = boltz.cpus || '<cpus>';
  const memVal = boltz.mem_gb ?? '<mem>';
  const memLabel = Number.isFinite(memVal) ? `${memVal}G` : `${memVal}`;
  const cacheDir = boltz.cache_dir ? `${boltz.cache_dir}` : null;
  const extraArgs = Array.isArray(boltz.extra_args) ? boltz.extra_args : (boltz.extra_args ? [boltz.extra_args] : []);

  const lines = [
    '# Manual BoltzGen submission (selection)',
    '',
    '# 1) Sync target configs + tools to cluster',
  ];
  targets.forEach((t) => {
    const pdb = (t?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    lines.push(`rsync -az ${localTarget} ${sshTarget}:${targetRoot}`);
  });
  lines.push(`rsync -az ${localTools} ${sshTarget}:${toolsRoot || '<remote_root>/lib'}`);

  lines.push('', '# 2) Verify boltzgen configs exist on cluster');
  targets.forEach((t) => {
    const pdb = (t?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const specs = Array.from(t?.specs || []);
    const remoteSpecs = (specs.length ? specs : ['configs/*/boltzgen_config.yaml'])
      .map((spec) => `${remoteTarget}/${spec}`.replace(/\/+/g, '/'));
    lines.push(`ssh ${sshTarget} "ls ${remoteSpecs.join(' ')}"`);
  });

  lines.push('', '# 3) Launch BoltzGen pipeline', condaLine);
  lines.push(`export INITBINDER_ROOT=${envRoot}`);
  lines.push(`export INITBINDER_TARGET_ROOT=${envTargetRoot}`);

  const monitorLines = [];
  const pullLines = [];
  targets.forEach((t) => {
    const pdb = (t?.pdbId || '').toUpperCase();
    if (!pdb) return;
    const epNames = Array.from(t?.epitopeNames || []).filter(Boolean);
    if (epNames.length) {
      lines.push('', `# ${pdb} (${epNames.join(', ')})`);
    } else {
      lines.push('', `# ${pdb}`);
    }
    const runLabel = runLabelFor(pdb);
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const outputRootBase = boltz.output_root
      ? `${boltz.output_root}`.replace(/\/$/, '')
      : `${remoteTarget}/designs/boltzgen`;
    const outputRoot = outputRootFor(outputRootBase, runLabel);
    const specs = Array.from(t?.specs || []);
    const remoteSpecs = (specs.length ? specs : ['configs/*/boltzgen_config.yaml'])
      .map((spec) => `${remoteTarget}/${spec}`.replace(/\/+/g, '/'));
    const specLines = remoteSpecs.map((spec) => `  --spec ${spec} \\`);
    const cacheLine = cacheDir ? `  --cache_dir ${cacheDir} \\` : null;
    const extraLine = extraArgs.length ? `  --extra_run_args ${extraArgs.join(' ')} \\` : null;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;
    const localOutputDir = `${localTarget}/designs/boltzgen/${runLabel}`;

    lines.push(
      `python ${pipelinePath} pipeline ${pdb} \\`,
      `  --run_label ${runLabel} \\`,
      `  --num_designs ${designCount} \\`,
      `  --output_root ${outputRoot} \\`,
      `  --partition ${partition} --gpus ${gpus} --cpus ${cpus} --mem ${memLabel} --time_h ${timeHours} \\`,
      ...specLines,
      cacheLine,
      extraLine,
      '  --submit',
    );

    monitorLines.push(`ssh ${sshTarget} "squeue -u $USER | grep ${runLabel}"`);
    monitorLines.push(`ssh ${sshTarget} "tail -f ${outputRoot}/launcher.log"`);
    pullLines.push(`mkdir -p ${localOutputDir}`);
    pullLines.push(`rsync -az ${sshTarget}:${outputRoot}/ ${localOutputDir}/`);
  });

  if (monitorLines.length) {
    lines.push('', '# Monitor', ...monitorLines);
  }
  if (pullLines.length) {
    lines.push('', '# Pull results back', ...pullLines);
  }

  return lines.filter(Boolean).join('\n');
}

function matchEpitopeConfig(target, epitopeToken) {
  const configs = Array.isArray(target?.configs) ? target.configs : [];
  if (!configs.length) return null;
  const raw = String(epitopeToken || '').trim();
  const normalized = normalizeEpitopeToken(raw);
  const normalizedUpper = normalized ? normalized.toUpperCase() : '';
  for (let idx = 0; idx < configs.length; idx += 1) {
    const cfg = configs[idx];
    const keys = epitopeKeysForConfig(cfg, idx + 1);
    if (normalized && keys.has(normalized)) return cfg;
    if (raw) {
      const candidates = [
        cfg?.epitope_name,
        cfg?.epitope_id,
        epitopeLabel(cfg, idx + 1),
      ].map((val) => String(val || '').trim().toLowerCase());
      if (candidates.includes(raw.toLowerCase())) return cfg;
    }
    const cfgPath = String(cfg?.config_path || '').toUpperCase();
    if (normalizedUpper && cfgPath.includes(normalizedUpper)) return cfg;
  }
  return null;
}

async function showBoltzRunCommand(pdbId, configPath = null, epitopeName = null) {
  if (!pdbId) {
    showAlert('Missing PDB ID for manual command.');
    return;
  }
  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing cluster commands...';
  if (el.boltzRunTitle) {
    el.boltzRunTitle.textContent = configPath
      ? `Manual run · ${epitopeName || configPath}`
      : `Manual run · ${pdbId}`;
  }
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }
  const target = (state.boltzConfigs || []).find((t) => (t.pdb_id || '').toUpperCase() === (pdbId || '').toUpperCase());
  let specs = [];
  if (configPath) {
    if (target && Array.isArray(target.configs)) {
      const cfg = target.configs.find((entry) => entry?.config_path === configPath);
      if (cfg && cfg.hotspot_surface_ok === false) {
        showAlert('Epitope filtered by SASA; no run command generated.');
        return;
      }
    }
    specs = [configPath];
  } else if (target && Array.isArray(target.configs)) {
    specs = target.configs.filter(isSurfaceAllowed).map((cfg) => cfg.config_path).filter(Boolean);
    if (!specs.length) {
      renderRunCommandBlocks('# No surface-exposed epitopes available for this target.');
      return;
    }
  }
  const text = buildRunCommandText(pdbId, specs, epitopeName);
  renderRunCommandBlocks(text, { singleBlock: true });
}

async function showRfaRunCommand(pdbId) {
  if (!pdbId) {
    showAlert('Missing PDB ID for manual command.');
    return;
  }
  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing RFA pipeline commands...';
  if (el.boltzRunTitle) el.boltzRunTitle.textContent = `RFA pipeline · ${pdbId}`;
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }
  try {
    const target = await fetchRfaConfigs(pdbId, { force: true });
    const launcherPath = target?.launcher_path || '';
    if (!launcherPath) {
      renderRunCommandBlocks('# No RFA pipeline launcher found. Use Rebuild configs to generate scripts.');
      return;
    }
    const text = buildRfaRunCommandText(pdbId, launcherPath);
    renderRunCommandBlocks(text);
  } catch (err) {
    showAlert(err.message || String(err));
  }
}

async function renderRunCommandModal() {
  const ctx = state.runContext;
  if (!ctx || !ctx.pdbId) {
    if (state.runMode === 'range') {
      await showRunCommandRange();
    } else {
      await showRunCommandAll();
    }
    return;
  }
  const engine = getRunEngine();
  if (engine === 'rfantibody') {
    await showRfaRunCommand(ctx.pdbId);
  } else {
    await showBoltzRunCommand(ctx.pdbId, ctx.configPath || null, ctx.epitopeName || null);
  }
}

function setRunContext(pdbId, configPath = null, epitopeName = null) {
  state.runContext = {
    pdbId: pdbId || '',
    configPath: configPath || null,
    epitopeName: epitopeName || null,
  };
  state.runMode = 'single';
}

async function showRunCommand(pdbId, configPath = null, epitopeName = null) {
  setRunContext(pdbId, configPath, epitopeName);
  await renderRunCommandModal();
}

async function showRfaRunCommandAll() {
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    renderRunCommandBlocks('# No targets detected.');
    return;
  }
  const filters = getRunCommandFilters();
  const { pass } = applyRunCommandFilters(targets, filters);
  const filteredTargets = pass;
  if (el.boltzRunTitle && filteredTargets.length !== targets.length) {
    el.boltzRunTitle.textContent = `Manual run · all targets (${filteredTargets.length}/${targets.length})`;
  }
  if (!filteredTargets.length) {
    const text = '# Manual RFA pipeline submission\nNo targets passed filters.';
    renderRunCommandBlocks(prependRunCommandFilters(text, filters, targets.length, filteredTargets.length));
    return;
  }
  const entries = await fetchRfaLauncherEntries(filteredTargets);
  const text = entries.length
    ? buildRfaRunCommandsTextForTargets(entries)
    : '# Manual RFA pipeline submission\nNo targets passed filters.';
  renderRunCommandBlocks(prependRunCommandFilters(text, filters, targets.length, filteredTargets.length));
}

async function showRunCommandAll() {
  state.runContext = null;
  state.runMode = 'all';
  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing cluster commands...';
  if (el.boltzRunTitle) el.boltzRunTitle.textContent = 'Manual run · all targets';
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }
  const engine = getRunEngine();
  if (engine === 'rfantibody') {
    await showRfaRunCommandAll();
    return;
  }
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    renderRunCommandBlocks('# No BoltzGen targets detected.');
    return;
  }
  const filters = getRunCommandFilters();
  const { pass } = applyRunCommandFilters(targets, filters);
  const filteredTargets = pass;
  if (el.boltzRunTitle && filteredTargets.length !== targets.length) {
    el.boltzRunTitle.textContent = `Manual run · all targets (${filteredTargets.length}/${targets.length})`;
  }
  const text = filteredTargets.length
    ? buildAllRunCommandsText(filteredTargets)
    : '# Manual BoltzGen submission\nNo targets passed filters.';
  renderRunCommandBlocks(prependRunCommandFilters(text, filters, targets.length, filteredTargets.length));
}

async function showRunCommandSelection() {
  state.runContext = null;
  state.runMode = 'selection';
  const raw = (el.epitopeDiversitySelection?.value || '').trim();
  if (!raw) {
    showAlert('Enter selections like 5WT9:epitope_1, 3J8F:epitope_4.');
    return;
  }
  const selections = parseEpitopeSelectionInput(raw);
  if (!selections.length) {
    showAlert('No valid selections found. Use PDB:epitope_#.');
    return;
  }
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    showAlert('No BoltzGen configs loaded.');
    return;
  }
  const engine = getRunEngine();
  if (engine === 'rfantibody') {
    showAlert('Selection-based commands are only available for BoltzGen.');
    return;
  }

  const targetMap = new Map(
    targets.map((t) => [(t?.pdb_id || '').toUpperCase(), t]),
  );
  const matched = [];
  const missingTargets = [];
  const missingEpitopes = [];
  const filteredEpitopes = [];
  selections.forEach((sel) => {
    const target = targetMap.get(sel.pdbId);
    if (!target) {
      missingTargets.push(sel.raw);
      return;
    }
    const cfg = matchEpitopeConfig(target, sel.epitope);
    if (!cfg) {
      missingEpitopes.push(sel.raw);
      return;
    }
    if (cfg.hotspot_surface_ok === false) {
      filteredEpitopes.push(sel.raw);
      return;
    }
    matched.push({
      pdbId: sel.pdbId,
      configPath: cfg.config_path || '',
      epitopeName: sel.epitope || epitopeLabel(cfg),
    });
  });

  if (!matched.length) {
    showAlert('No matching configs found for that selection.');
    return;
  }

  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing cluster commands...';
  if (el.boltzRunTitle) {
    el.boltzRunTitle.textContent = `Manual run · selection (${matched.length})`;
  }
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }

  let text = buildSelectedRunCommandsText(matched);
  const notes = [];
  if (missingTargets.length) {
    notes.push(`# Missing targets: ${missingTargets.join(', ')}`);
  }
  if (missingEpitopes.length) {
    notes.push(`# Missing epitopes: ${missingEpitopes.join(', ')}`);
  }
  if (filteredEpitopes.length) {
    notes.push(`# SASA filtered epitopes: ${filteredEpitopes.join(', ')}`);
  }
  if (notes.length) {
    text = `${notes.join('\n')}\n\n${text}`;
  }
  renderRunCommandBlocks(text, { singleBlock: true });
}

async function runBoltzTarget(pdbId, triggerEl = null) {
  const designs = getBoltzDesignCount();
  if (!designs) {
    showAlert('Enter a positive epitope design count.');
    return;
  }
  const payload = {
    pdb_id: pdbId,
    design_count: designs,
    run_label_prefix: (el.bulkRunPrefix?.value || 'boltz').trim() || 'boltz',
    time_hours: getBoltzTimeHours(),
  };
  if (triggerEl) triggerEl.disabled = true;
  try {
    const res = await fetch('/api/bulk/boltzgen/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to submit BoltzGen run (${res.status})`);
    }
    const body = await res.json();
    showAlert(body.message || 'BoltzGen config run queued.', false);
    await loadBoltzConfigs({ silent: true });
  } catch (err) {
    showAlert(err.message || String(err));
  } finally {
    if (triggerEl) triggerEl.disabled = false;
  }
}

async function runBoltzEpitope(pdbId, configPath, triggerEl = null) {
  const designs = getBoltzDesignCount();
  if (!designs) {
    showAlert('Enter a positive epitope design count.');
    return;
  }
  const payload = {
    pdb_id: pdbId,
    config_path: configPath,
    design_count: designs,
    run_label_prefix: (el.bulkRunPrefix?.value || 'boltz').trim() || 'boltz',
    time_hours: getBoltzTimeHours(),
  };
  if (triggerEl) triggerEl.disabled = true;
  try {
    const res = await fetch('/api/bulk/boltzgen/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed to submit epitope run (${res.status})`);
    }
    const body = await res.json();
    showAlert(body.message || 'Epitope run queued.', false);
    await loadBoltzConfigs({ silent: true });
  } catch (err) {
    showAlert(err.message || String(err));
  } finally {
    if (triggerEl) triggerEl.disabled = false;
  }
}

function handleBinderTableClick(event) {
  const btn = event.target.closest('button[data-action="pymol-binder"]');
  if (!btn) return;
  const idx = Number(btn.dataset.index || -1);
  if (Number.isFinite(idx) && idx >= 0) {
    launchBinderPymol(idx);
  }
}

function handleBinderPagination(event) {
  const btn = event.target.closest('button[data-page-nav]');
  if (!btn) return;
  const dir = btn.dataset.pageNav;
  const current = state.binderPage || 1;
  const totalPages = Math.max(1, Math.ceil((state.binderTotal || 0) / (state.binderPageSize || 1)));
  if (dir === 'prev' && current > 1) {
    refreshDiversity({ silent: true, page: current - 1 });
  } else if (dir === 'next' && current < totalPages) {
    refreshDiversity({ silent: true, page: current + 1 });
  }
}

function handleBoltzTableClick(event) {
  const btn = event.target.closest('button[data-action]');
  if (!btn) {
    const row = event.target.closest('tr.target-row');
    if (row && row.dataset.pdbId) {
      toggleEpitopeRows(row.dataset.pdbId);
    }
    return;
  }
  const action = btn.dataset.action;
  const pdbId = btn.dataset.pdbId;
  const configPath = btn.dataset.configPath;
  if (action === 'run-target') {
    runBoltzTarget(pdbId, btn);
  } else if (action === 'run-epitope') {
    runBoltzEpitope(pdbId, configPath, btn);
  } else if (action === 'show-config-target') {
    setConfigContext(pdbId, null, null);
    renderConfigModal();
  } else if (action === 'show-config-epitope') {
    setConfigContext(pdbId, configPath, btn.dataset.epitopeName || null);
    renderConfigModal();
  } else if (action === 'show-target-yaml') {
    showTargetYaml(pdbId);
  } else if (action === 'show-log') {
    showBoltzLog(btn.dataset.jobId, btn.dataset.logTitle || 'Job log');
  } else if (action === 'show-run') {
    showRunCommand(pdbId, configPath, btn.dataset.epitopeName || null);
  } else if (action === 'pymol-target') {
    launchBoltzPymol(pdbId, null, btn);
  } else if (action === 'pymol-epitope') {
    launchBoltzPymol(pdbId, btn.dataset.epitopeName || null, btn);
  } else if (action === 'rerun-pipeline') {
    openPipelineRerunModal(pdbId, btn.dataset.antigenUrl || null);
  }
}

function toggleEpitopeRows(pdbId) {
  const key = (pdbId || '').toUpperCase();
  if (!(state.expandedTargets instanceof Set)) {
    state.expandedTargets = new Set();
  }
  if (state.expandedTargets.has(key)) {
    state.expandedTargets.delete(key);
  } else {
    state.expandedTargets.add(key);
  }
  const isExpanded = state.expandedTargets.has(key);
  const rows = el.boltzTable?.querySelectorAll(`tr.epitope-row[data-parent="${key}"]`);
  if (rows) {
    rows.forEach((row) => {
      row.hidden = !isExpanded;
    });
  }
  const parentRow = el.boltzTable?.querySelector(`tr.target-row[data-pdb-id="${key}"]`);
  if (parentRow) {
    parentRow.classList.toggle('expanded', isExpanded);
    parentRow.classList.toggle('collapsed', !isExpanded);
  }
}

function renderBulkPreview(rows = [], summary = '') {
  if (!el.bulkPreviewPanel || !el.bulkPreviewTable) return;
  const tbody = el.bulkPreviewTable;
  tbody.innerHTML = '';
  const list = Array.isArray(rows) ? rows : [];
  if (!list.length) {
    el.bulkPreviewPanel.hidden = true;
    return;
  }
  list.forEach((row) => {
    const tr = document.createElement('tr');
    const values = [
      row.raw_index ?? '',
      row.preset_name || '',
      row.protein_name || '',
      row.antigen_url || '',
      row.accession || '',
      row.resolved_pdb_id || row.pdb_id || '—',
      Array.isArray(row.warnings) && row.warnings.length ? row.warnings.join('; ') : '',
    ];
    values.forEach((value) => {
      const td = document.createElement('td');
      td.textContent = value;
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
  if (el.bulkPreviewSummary) {
    el.bulkPreviewSummary.textContent = summary || `${list.length} rows parsed`;
    el.bulkPreviewSummary.hidden = false;
  }
  el.bulkPreviewPanel.hidden = false;
}

async function previewBulkCsv(options = {}) {
  const { silent = false } = options;
  const csvText = el.bulkCsvInput?.value || '';
  if (!csvText.trim()) {
    if (!silent) showAlert('Paste a CSV/TSV payload first.');
    return;
  }
  const numEpitopes = (() => {
    const raw = Number(el.bulkEpitopes?.value || 0);
    if (!Number.isFinite(raw) || raw <= 0) return null;
    return Math.min(32, Math.max(1, Math.round(raw)));
  })();
  const payload = {
    csv_text: csvText,
    num_epitopes: numEpitopes,
    decide_scope_prompt: (el.bulkEpitopePrompt?.value || '').trim() || null,
  };
  try {
    const res = await fetch('/api/bulk/preview', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Preview failed with ${res.status}`);
    }
    const body = await res.json();
    state.bulkPreviewRows = Array.isArray(body.rows) ? body.rows : [];
    state.binderPage = 1;
    renderBulkPreview(state.bulkPreviewRows, body.message || '');
    loadBoltzConfigs({ silent: true });
    loadBinderTable({ page: 1, silent: true });
    const previewPdbIds = Array.from(
      new Set(
        state.bulkPreviewRows
          .map((row) => (row.resolved_pdb_id || row.pdb_id || '').trim().toUpperCase())
          .filter(Boolean),
      ),
    );
    hydrateSnapshotsFromCache({
      pdbIds: previewPdbIds,
      limit: Math.max(12, Math.max(previewPdbIds.length, 1) * 3),
      force: true,
      silent: true,
    });
    setBadge(el.bulkStatus, body.message || 'Preview ready', 'rgba(134, 239, 172, 0.25)');
    if (!silent) showAlert(body.message || 'Preview ready.', false);
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
    setBadge(el.bulkStatus, 'Preview failed', 'rgba(248, 113, 113, 0.25)');
  }
}

function collectDesignSettings() {
  const settings = {
    model_engine: (el.designEngine?.value || 'boltzgen').trim() || 'boltzgen',
    total_designs: Number(el.designTotal?.value || 0) || 90,
    run_assess: true,
    run_label_prefix: (el.bulkRunPrefix?.value || '').trim() || null,
    boltz_time_hours: getBoltzTimeHours(),
    boltzgen_crop_radius: getBoltzCropRadius(),
  };
  return settings;
}

async function startBulkRun(options = {}) {
  const {
    submitDesignsOverride = null,
    triggerEl = null,
    llmDelaySeconds = null,
    decideScopeAttempts = null,
    launchPymolOverride = null,
    renderSnapshotsOverride = null,
  } = options;
  const csvText = el.bulkCsvInput?.value || '';
  if (!csvText.trim()) {
    showAlert('Paste a CSV/TSV payload first.');
    return;
  }
  const numEpitopes = (() => {
    const raw = Number(el.bulkEpitopes?.value || 0);
    if (!Number.isFinite(raw) || raw <= 0) return null;
    return Math.min(32, Math.max(1, Math.round(raw)));
  })();
  const throttle = (() => {
    const rawVal = el.bulkThrottle?.value;
    const raw = Number(rawVal);
    if (!Number.isFinite(raw)) return 3;
    if (raw < 0) return 0;
    return raw;
  })();
  const launchPymol = typeof launchPymolOverride === 'boolean'
    ? launchPymolOverride
    : (el.bulkLaunchPymol ? el.bulkLaunchPymol.checked : true);
  const renderSnapshots = typeof renderSnapshotsOverride === 'boolean'
    ? renderSnapshotsOverride
    : launchPymol;
  const payload = {
    csv_text: csvText,
    num_epitopes: numEpitopes,
    decide_scope_prompt: (el.bulkEpitopePrompt?.value || '').trim() || null,
    launch_pymol: launchPymol,
    render_pymol_snapshots: renderSnapshots,
    export_insights: el.bulkExportInsights ? el.bulkExportInsights.checked : true,
    export_designs: el.bulkExportDesigns ? el.bulkExportDesigns.checked : true,
    submit_designs: submitDesignsOverride ?? true,
    force_init: el.bulkForceInit ? el.bulkForceInit.checked : false,
    prepare_targets: true,
    design_settings: collectDesignSettings(),
    throttle_seconds: throttle,
    llm_delay_seconds: Number.isFinite(llmDelaySeconds) && llmDelaySeconds > 0 ? llmDelaySeconds : 0,
    decide_scope_attempts: Number.isFinite(decideScopeAttempts) && decideScopeAttempts > 0
      ? Math.min(5, Math.max(1, Math.round(decideScopeAttempts)))
      : 1,
  };
  setActionButtonsDisabled(true);
  setBadge(el.bulkStatus, 'Queuing…');
  resetJobLog('Waiting for commands...');
  try {
    const res = await fetch('/api/bulk/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Bulk submission failed with ${res.status}`);
    }
    const body = await res.json();
    showAlert(body.message || 'Bulk job queued.', false);
    if (body.job_id) {
      startJobPolling(body.job_id);
    } else {
      setActionButtonsDisabled(false);
    }
  } catch (err) {
    showAlert(err.message || String(err));
    setActionButtonsDisabled(false);
    setBadge(el.bulkStatus, 'Bulk submission failed', 'rgba(248, 113, 113, 0.25)');
  }
}

async function fetchJob(jobId) {
  const res = await fetch(`/api/jobs/${jobId}`);
  if (!res.ok) throw new Error('Unable to retrieve job status');
  return res.json();
}

function stopJobPolling() {
  if (state.jobPoller) {
    clearInterval(state.jobPoller);
    state.jobPoller = null;
  }
}

function updateJobUI(job) {
  if (!job) return;
  const prevStatus = state.lastJobStatus;
  state.lastJobStatus = job.status;
  resetJobLog(job.logs.join('\n'));
  const newSnapshots = Array.isArray(job.details?.snapshots) ? job.details.snapshots.filter(Boolean) : [];
  if (newSnapshots.length && snapshotsChanged(newSnapshots)) {
    loadSnapshotMetadata(newSnapshots);
  }
  const epPlotFiles = Array.isArray(job.details?.epitope_plot_files) ? job.details.epitope_plot_files : [];
  const epPlotPaths = Array.isArray(job.details?.epitope_plots) ? job.details.epitope_plots : [];
  const plotDir = typeof job.details?.plot_dir === 'string' ? job.details.plot_dir : '';
  const epPlotNames = [];
  const seen = new Set();
  const makeSrc = (val) => {
    if (!val) return null;
    const text = String(val);
    // Never use `file://` for assets from an HTTP-served page; browsers typically block it.
    // Instead, always ask the backend to serve the basename from the bulk output directory.
    const baseName = text.split('/').pop();
    return `/api/bulk/file?name=${encodeURIComponent(baseName || text)}&t=${Date.now()}`;
  };
  epPlotFiles.forEach((name) => {
    const label = (name || '').split('/').pop();
    if (!label || seen.has(label)) return;
    seen.add(label);
    const candidate = plotDir ? `${plotDir.replace(/\/$/, '')}/${label}` : label;
    const src = makeSrc(candidate) || makeSrc(label);
    if (src) epPlotNames.push({ src, label });
  });
  epPlotPaths.forEach((path) => {
    const label = (path || '').split('/').pop();
    if (!label || seen.has(label)) return;
    seen.add(label);
    const src = makeSrc(path);
    if (src) epPlotNames.push({ src, label });
  });
  if (epPlotNames.length || state.epitopePlots.length) {
    renderEpitopePlots(epPlotNames);
  }
  if (el.jobMeta) {
    el.jobMeta.innerHTML = '';
    const parts = [];
    const addLink = (label, filename) => {
      const a = document.createElement('a');
      a.href = `/api/bulk/file?name=${encodeURIComponent(filename)}`;
      a.textContent = label;
      a.target = '_blank';
      a.rel = 'noopener noreferrer';
      parts.push(a);
    };
    if (job.details?.log_path) {
      const name = job.details.log_path.split('/').pop();
      if (name) addLink('Log', name);
    }
    if (job.details?.insights_filename) {
      addLink('Insights CSV', job.details.insights_filename);
    }
    if (job.details?.design_config_filename) {
      addLink('Design config CSV', job.details.design_config_filename);
    }
    if (job.details?.selection_filename) {
      addLink('Epitope input TSV', job.details.selection_filename);
    }
    if (job.details?.targets_table_filename) {
      addLink('Detected targets', job.details.targets_table_filename);
    }
    if (job.details?.snapshot_report_filename) {
      addLink('Snapshot report', job.details.snapshot_report_filename);
    }
    if (Array.isArray(job.details?.snapshots) && job.details.snapshots.length) {
      job.details.snapshots.forEach((snap) => {
        const name = (snap || '').trim();
        if (name) addLink('Snapshot', name);
      });
    }
    if (epPlotNames.length) {
      epPlotNames.forEach((name) => addLink('Epitope plot', name));
    }
    if (parts.length) {
      parts.forEach((node, idx) => {
        el.jobMeta.appendChild(node);
        if (idx < parts.length - 1) {
          const sep = document.createTextNode(' · ');
          el.jobMeta.appendChild(sep);
        }
      });
      el.jobMeta.hidden = false;
    } else {
      el.jobMeta.hidden = true;
    }
  }

  if (job.status === 'running' || job.status === 'pending') {
    setBadge(el.bulkStatus, job.status === 'running' ? 'Running…' : 'Queued…');
    setActionButtonsDisabled(true);
  } else if (job.status === 'success') {
    setBadge(el.bulkStatus, 'Bulk complete', 'rgba(134, 239, 172, 0.25)');
    showAlert(job.message || 'Bulk job finished.', false);
    setActionButtonsDisabled(false);
    if (job.details?.design_config_filename && prevStatus !== 'success') {
      const url = `/api/bulk/file?name=${encodeURIComponent(job.details.design_config_filename)}`;
      window.open(url, '_blank');
    }
    if (Array.isArray(job.details?.snapshots) && job.details.snapshots.length) {
      loadSnapshotMetadata(job.details.snapshots.filter(Boolean));
    } else {
      renderSnapshots([]);
    }
    renderEpitopePlots(epPlotNames);
    refreshDiversity({ silent: true });
    stopJobPolling();
  } else {
    setBadge(el.bulkStatus, 'Bulk failed', 'rgba(248, 113, 113, 0.25)');
    showAlert(job.message || 'Bulk job failed.');
    setActionButtonsDisabled(false);
    stopJobPolling();
  }
}

function startJobPolling(jobId) {
  state.currentJobId = jobId;
  state.lastJobStatus = null;
  stopJobPolling();
  renderSnapshots([]);
  renderEpitopePlots([]);
  renderDiversityPlots([]);
  state.diversityCsv = null;
  state.diversityHtml = null;
  const poll = async () => {
    try {
      const job = await fetchJob(jobId);
      updateJobUI(job);
    } catch (err) {
      appendLog(`Polling error: ${err.message || err}`);
      stopJobPolling();
      setActionButtonsDisabled(false);
    }
  };
  poll();
  state.jobPoller = setInterval(poll, 2000);
}

function scrollToSnapshots() {
  if (!el.snapshotSection) return;
  el.snapshotSection.hidden = false;
  el.snapshotSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
  if (!el.snapshotGrid || el.snapshotGrid.childElementCount === 0) {
    appendLog("No hotspot snapshots yet. Run a bulk job with PyMOL enabled to generate visualizations.");
  }
}

async function downloadSnapshotPdf() {
  try {
    if (!el.snapshotGrid || el.snapshotGrid.childElementCount === 0) {
      showAlert('No hotspot snapshots to export yet.');
      return;
    }

    const names = (state.snapshotNames || []).filter(Boolean);
    if (names.length) {
      const params = new URLSearchParams();
      params.append('names', names.join(','));
      const res = await fetch(`/api/pymol/snapshots/package?${params.toString()}`);
      if (res.ok) {
        const blob = await res.blob();
        const cd = res.headers.get('content-disposition') || '';
        const match = cd.match(/filename=\"?([^\";]+)\"?/i);
        const filename = match?.[1] || 'hotspot_snapshots.html';
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(url), 5000);
        return;
      }
    }

    // Fallback: open printable window
    const printable = el.snapshotGrid.cloneNode(true);
    printable.classList.remove('catalog-viewer');
    printable.style.maxHeight = 'none';
    printable.style.overflow = 'visible';
    printable.style.border = 'none';
    printable.style.background = 'transparent';
    const container = document.createElement('div');
    container.appendChild(printable);
    const html = `
      <!doctype html>
      <html>
        <head>
          <meta charset="utf-8">
          <title>Hotspot snapshots</title>
          <link rel="stylesheet" href="/static/css/app.css">
          <style>
            body { padding: 18px; background: #ffffff; }
            h1 { margin-bottom: 12px; }
            .catalog-viewer { max-height: none !important; overflow: visible !important; border: none; }
            .snapshot-card { page-break-inside: avoid; break-inside: avoid; }
            .snapshot-card + .snapshot-card { margin-top: 12px; }
            .snapshot-card img { max-width: 100%; height: auto; }
            @media print { body { padding: 0 12px; } }
          </style>
        </head>
        <body>
          <h1>Hotspot snapshots</h1>
          <p class="help-text" style="margin-top: 0;">Includes all hotspot images and annotations currently shown.</p>
          ${container.innerHTML}
        </body>
      </html>`;
    const popup = window.open('', '_blank', 'width=1000,height=1200');
    if (!popup || !popup.document) {
      showAlert('Pop-up blocked. Allow pop-ups to download the PDF.');
      return;
    }
    popup.document.open();
    popup.document.write(html);
    popup.document.close();
    popup.focus();
    const imgs = Array.from(popup.document.images || []);
    await Promise.all(imgs.map((img) => new Promise((resolve) => {
      if (img.complete) return resolve();
      img.onload = () => resolve();
      img.onerror = () => resolve();
    })));
    popup.print();
    setTimeout(() => popup.close(), 500);
  } catch (err) {
    showAlert(`Failed to prepare PDF: ${err.message || err}`);
  }
}

async function handleVisualizeEpitopes() {
  await startBulkRun({
    submitDesignsOverride: false,
    triggerEl: el.bulkVisualizeEpitopes,
    llmDelaySeconds: 5,
    decideScopeAttempts: 3,
    launchPymolOverride: false,
    renderSnapshotsOverride: true,
  });
  scrollToSnapshots();
}

function setupExampleTabs() {
  const containers = document.querySelectorAll('[data-example-tabs]');
  containers.forEach((container) => {
    const tabs = Array.from(container.querySelectorAll('[data-example-tab]'));
    const panels = Array.from(container.querySelectorAll('[data-example-panel]'));
    if (!tabs.length || !panels.length) return;

    const activate = (key) => {
      tabs.forEach((btn) => {
        const active = btn.dataset.exampleTab === key;
        btn.classList.toggle('is-active', active);
        btn.setAttribute('aria-selected', active ? 'true' : 'false');
      });
      panels.forEach((panel) => {
        panel.classList.toggle('is-active', panel.dataset.examplePanel === key);
      });
    };

    tabs.forEach((btn) => {
      btn.addEventListener('click', () => activate(btn.dataset.exampleTab));
    });

    const initial = tabs.find((btn) => btn.classList.contains('is-active')) || tabs[0];
    if (initial) activate(initial.dataset.exampleTab);
  });
}

function init() {
  loadClusterStatus();
  setupExampleTabs();
  if (el.bulkPreviewBtn) el.bulkPreviewBtn.addEventListener('click', () => previewBulkCsv({ silent: false }));
  if (el.bulkPreviewRefresh) el.bulkPreviewRefresh.addEventListener('click', () => previewBulkCsv({ silent: true }));
  if (el.bulkCsvInput) el.bulkCsvInput.addEventListener('input', scheduleBulkPreview);
  if ((el.bulkCsvInput?.value || '').trim()) {
    previewBulkCsv({ silent: true });
  }
  if (el.bulkVisualizeEpitopes) el.bulkVisualizeEpitopes.addEventListener('click', handleVisualizeEpitopes);
  if (el.bulkRunBtn) el.bulkRunBtn.addEventListener('click', () => startBulkRun());
  if (el.snapshotDownloadPdf) el.snapshotDownloadPdf.addEventListener('click', downloadSnapshotPdf);
  if (el.boltzRefresh) {
    el.boltzRefresh.addEventListener('click', () => {
      const hasPreview = Array.isArray(state.bulkPreviewRows) && state.bulkPreviewRows.length > 0;
      const hasCsv = Boolean((el.bulkCsvInput?.value || '').trim());
      if (!hasPreview && hasCsv) {
        previewBulkCsv({ silent: true });
      } else {
        loadBoltzConfigs({ silent: true });
      }
    });
  }
  if (el.boltzRegenerate) {
    el.boltzRegenerate.addEventListener('click', openRegenerateRangeModal);
  }
  if (el.diversityRefresh) el.diversityRefresh.addEventListener('click', () => refreshDiversity());
  if (el.diversityRebuild) el.diversityRebuild.addEventListener('click', forceRefreshDiversity);
  if (el.diversityDownloadCsv) el.diversityDownloadCsv.addEventListener('click', () => downloadDiversityFile('csv'));
  if (el.diversityDownloadHtml) el.diversityDownloadHtml.addEventListener('click', () => downloadDiversityFile('html'));
  if (el.boltzPlotDiversity) el.boltzPlotDiversity.addEventListener('click', () => plotBoltzAntigenDiversity());
  if (el.epitopeDiversityPlot) el.epitopeDiversityPlot.addEventListener('click', () => plotEpitopeDiversitySelection());
  if (el.boltzShowRunSelection) el.boltzShowRunSelection.addEventListener('click', showRunCommandSelection);
  if (el.boltzShowRunAll) el.boltzShowRunAll.addEventListener('click', showRunCommandAll);
  if (el.boltzRerunRange) {
    el.boltzRerunRange.addEventListener('click', openPipelineRerunRangeModal);
  }
  if (el.boltzShowRunRange) {
    el.boltzShowRunRange.addEventListener('click', openRunCommandRangeModal);
  }
  if (el.boltzRerunRangeConfirm) {
    el.boltzRerunRangeConfirm.addEventListener('click', (event) => submitPipelineRerunRangeModal(event.currentTarget));
  }
  if (el.boltzRerunRangeCancel) {
    el.boltzRerunRangeCancel.addEventListener('click', () => toggleModal(el.boltzRerunRangeModal, false));
  }
  if (el.boltzRerunRangeClose) {
    el.boltzRerunRangeClose.addEventListener('click', () => toggleModal(el.boltzRerunRangeModal, false));
  }
  if (el.boltzRunRangeConfirm) {
    el.boltzRunRangeConfirm.addEventListener('click', showRunCommandRange);
  }
  if (el.boltzRunRangeCancel) {
    el.boltzRunRangeCancel.addEventListener('click', () => toggleModal(el.boltzRunRangeModal, false));
  }
  if (el.boltzRunRangeClose) {
    el.boltzRunRangeClose.addEventListener('click', () => toggleModal(el.boltzRunRangeModal, false));
  }
  if (el.boltzRunRangeStart) {
    el.boltzRunRangeStart.addEventListener('input', updateRunCommandFilterNote);
    el.boltzRunRangeStart.addEventListener('change', updateRunCommandFilterNote);
  }
  if (el.boltzRunRangeEnd) {
    el.boltzRunRangeEnd.addEventListener('input', updateRunCommandFilterNote);
    el.boltzRunRangeEnd.addEventListener('change', updateRunCommandFilterNote);
  }
  if (el.boltzRunRangeMinAllowed) {
    el.boltzRunRangeMinAllowed.addEventListener('input', updateRunCommandFilterNote);
    el.boltzRunRangeMinAllowed.addEventListener('change', updateRunCommandFilterNote);
  }
  if (el.boltzRunRangeMinEpitopes) {
    el.boltzRunRangeMinEpitopes.addEventListener('input', updateRunCommandFilterNote);
    el.boltzRunRangeMinEpitopes.addEventListener('change', updateRunCommandFilterNote);
  }
  if (el.boltzRegenerateRangeConfirm) {
    el.boltzRegenerateRangeConfirm.addEventListener('click', (event) => {
      submitRegenerateRangeModal(event.currentTarget);
    });
  }
  if (el.boltzRegenerateRangeCancel) {
    el.boltzRegenerateRangeCancel.addEventListener('click', () => toggleModal(el.boltzRegenerateRangeModal, false));
  }
  if (el.boltzRegenerateRangeClose) {
    el.boltzRegenerateRangeClose.addEventListener('click', () => toggleModal(el.boltzRegenerateRangeModal, false));
  }
  if (el.boltzTable) el.boltzTable.addEventListener('click', handleBoltzTableClick);
  if (el.binderTable) el.binderTable.addEventListener('click', handleBinderTableClick);
  if (el.binderPagination) el.binderPagination.addEventListener('click', handleBinderPagination);
  if (el.binderRefresh) el.binderRefresh.addEventListener('click', () => refreshDiversity({ silent: false, page: state.binderPage || 1 }));
  if (el.binderExportBtn) el.binderExportBtn.addEventListener('click', () => toggleModal(el.binderExportModal, true));
  if (el.binderDownload) el.binderDownload.addEventListener('click', downloadBinderCsv);
  if (el.bulkAlgorithmBtn) el.bulkAlgorithmBtn.addEventListener('click', () => toggleModal(el.bulkAlgorithmModal, true));
  if (el.binderFilterPdb) {
    el.binderFilterPdb.addEventListener('input', () => scheduleBinderRefresh({ silent: true }));
    el.binderFilterPdb.addEventListener('change', () => scheduleBinderRefresh({ silent: true }));
  }
  if (el.binderFilterEpitope) {
    el.binderFilterEpitope.addEventListener('input', () => scheduleBinderRefresh({ silent: true }));
    el.binderFilterEpitope.addEventListener('change', () => scheduleBinderRefresh({ silent: true }));
  }
  if (el.binderFilterEngine) {
    el.binderFilterEngine.addEventListener('change', () => scheduleBinderRefresh({ silent: true }));
  }
  if (el.binderOrderBy) {
    el.binderOrderBy.addEventListener('change', () => {
      state.binderPage = 1;
      refreshDiversity({ silent: true, page: 1 });
    });
  }
  if (el.boltzConfigEngine) {
    el.boltzConfigEngine.addEventListener('change', () => {
      if (el.boltzConfigModal && !el.boltzConfigModal.hidden) {
        renderConfigModal();
      }
    });
  }
  if (el.boltzRunEngine) {
    if (el.boltzRunEngine.tagName === 'SELECT') {
      el.boltzRunEngine.addEventListener('change', () => {
        if (el.boltzRunModal && !el.boltzRunModal.hidden) {
          renderRunCommandModal();
        }
      });
    } else {
      el.boltzRunEngine.addEventListener('click', (evt) => {
        const btn = evt.target?.closest('button[data-engine]');
        if (!btn) return;
        setRunEngine(btn.dataset.engine || 'boltzgen');
        if (el.boltzRunModal && !el.boltzRunModal.hidden) {
          renderRunCommandModal();
        }
      });
    }
  }
  if (el.boltzConfigClose) el.boltzConfigClose.addEventListener('click', () => toggleModal(el.boltzConfigModal, false));
  if (el.boltzLogClose) el.boltzLogClose.addEventListener('click', () => toggleModal(el.boltzLogModal, false));
  if (el.boltzRunClose) el.boltzRunClose.addEventListener('click', () => toggleModal(el.boltzRunModal, false));
  if (el.bulkAlgorithmClose) el.bulkAlgorithmClose.addEventListener('click', () => toggleModal(el.bulkAlgorithmModal, false));
  if (el.binderExportClose) el.binderExportClose.addEventListener('click', () => toggleModal(el.binderExportModal, false));
  if (el.binderExportCancel) el.binderExportCancel.addEventListener('click', () => toggleModal(el.binderExportModal, false));
  if (el.binderExportConfirm) el.binderExportConfirm.addEventListener('click', exportSelectedBinders);
  if (el.pipelineRerunClose) el.pipelineRerunClose.addEventListener('click', () => toggleModal(el.pipelineRerunModal, false));
  if (el.pipelineRerunCancel) el.pipelineRerunCancel.addEventListener('click', () => toggleModal(el.pipelineRerunModal, false));
  if (el.pipelineRerunConfirm) el.pipelineRerunConfirm.addEventListener('click', () => submitPipelineRerun(el.pipelineRerunConfirm));
  if (el.epitopeReportDownload) el.epitopeReportDownload.addEventListener('click', downloadEpitopeReportHtml);
  if (el.epitopeHotspotCsvDownload) el.epitopeHotspotCsvDownload.addEventListener('click', downloadEpitopeHotspotCsv);
  if (el.epitopeCsvDownload) el.epitopeCsvDownload.addEventListener('click', downloadEpitopeCsv);
  document.addEventListener('click', (evt) => {
    const closeTarget = evt.target?.dataset?.close;
    if (closeTarget === 'boltz-config') toggleModal(el.boltzConfigModal, false);
    if (closeTarget === 'boltz-log') toggleModal(el.boltzLogModal, false);
    if (closeTarget === 'boltz-run') toggleModal(el.boltzRunModal, false);
    if (closeTarget === 'boltz-rerun-range') toggleModal(el.boltzRerunRangeModal, false);
    if (closeTarget === 'boltz-run-range') toggleModal(el.boltzRunRangeModal, false);
    if (closeTarget === 'boltz-regenerate-range') toggleModal(el.boltzRegenerateRangeModal, false);
    if (closeTarget === 'pipeline-rerun') toggleModal(el.pipelineRerunModal, false);
    if (closeTarget === 'bulk-algorithm') toggleModal(el.bulkAlgorithmModal, false);
    if (closeTarget === 'binder-export') toggleModal(el.binderExportModal, false);
  });

  renderBoltzConfigs();
  setRunEngine(getRunEngine());
  updateRunCommandFilterNote();
  refreshDiversity({ silent: true });
  hydrateSnapshotsFromCache({ silent: true });
}

document.addEventListener('DOMContentLoaded', init);
