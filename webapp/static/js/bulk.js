const state = {
  bulkPreviewRows: [],
  jobPoller: null,
  currentJobId: null,
  snapshotNames: [],
  boltzConfigs: [],
  epitopePlots: [],
  diversityPlots: [],
  lastJobStatus: null,
  commandDefaults: null,
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
  bulkPreviewMessage: null,
  binderFilterTimer: null,
  rfaConfigs: {},
  configContext: null,
  runContext: null,
  runMode: null,
  uiSettings: null,
  guiReadme: null,
  llmConversations: [],
  activeConversationId: null,
  llmHistory: [],
  llmMatchedRows: [],
  llmUnmatched: [],
  llmUnmatchedActions: {},
  llmPickedPdbIds: [],
  llmDeletedPdbIds: [],
  llmSuggestPending: false,
  llmUnmatchedPoller: null,
  llmViewMode: 'all',
  llmCatalogFilter: 'biotin',
  activeCatalogName: null,
};

const PIPELINE_RERUN_DELAY_MS = 3000; // 3 seconds
const LLM_LOCAL_STORAGE_KEY = 'initbinder.bulk.llm.state.v1';
const LLM_MAX_HISTORY_STORED = 200;
const LLM_HISTORY_WINDOW_FOR_API = 30;

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
  llmTargetPanel: document.querySelector('#llm-target-panel'),
  llmConversationSelect: document.querySelector('#llm-conversation-select'),
  llmNewConversation: document.querySelector('#llm-new-conversation'),
  llmViewMode: document.querySelector('#llm-view-mode'),
  llmCatalogFilter: document.querySelector('#llm-catalog-filter'),
  llmChatTranscript: document.querySelector('#llm-chat-transcript'),
  llmChatPrompt: document.querySelector('#llm-chat-prompt'),
  llmSuggestTargets: document.querySelector('#llm-suggest-targets'),
  llmSuggestLoading: document.querySelector('#llm-suggest-loading'),
  llmMatchedPanel: document.querySelector('#llm-matched-panel'),
  llmUnmatchedPanel: document.querySelector('#llm-unmatched-panel'),
  llmMatchedList: document.querySelector('#llm-matched-list'),
  llmUnmatchedList: document.querySelector('#llm-unmatched-list'),
  llmMatchedCount: document.querySelector('#llm-matched-count'),
  llmUnmatchedCount: document.querySelector('#llm-unmatched-count'),
  bulkReadmeBtn: document.querySelector('#bulk-readme-btn'),
  bulkSettingsBtn: document.querySelector('#bulk-settings-btn'),
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
  bulkSettingsModal: document.querySelector('#bulk-settings-modal'),
  bulkSettingsClose: document.querySelector('#bulk-settings-close'),
  bulkSettingsCancel: document.querySelector('#bulk-settings-cancel'),
  bulkSettingsSave: document.querySelector('#bulk-settings-save'),
  bulkReadmeModal: document.querySelector('#bulk-readme-modal'),
  bulkReadmeTitle: document.querySelector('#bulk-readme-title'),
  bulkReadmeBody: document.querySelector('#bulk-readme-body'),
  bulkReadmeClose: document.querySelector('#bulk-readme-close'),
  bulkSettingsMeta: document.querySelector('#bulk-settings-meta'),
  bulkSettingsSshAlias: document.querySelector('#bulk-settings-ssh-alias'),
  bulkSettingsRemoteRoot: document.querySelector('#bulk-settings-remote-root'),
  bulkSettingsTargetRoot: document.querySelector('#bulk-settings-target-root'),
  bulkSettingsConda: document.querySelector('#bulk-settings-conda'),
  bulkSettingsPymolPath: document.querySelector('#bulk-settings-pymol-path'),
  bulkSettingsPymolCondaEnv: document.querySelector('#bulk-settings-pymol-conda-env'),
  bulkSettingsBoltzPartition: document.querySelector('#bulk-settings-boltz-partition'),
  bulkSettingsBoltzAccount: document.querySelector('#bulk-settings-boltz-account'),
  bulkSettingsBoltzGpus: document.querySelector('#bulk-settings-boltz-gpus'),
  bulkSettingsBoltzCpus: document.querySelector('#bulk-settings-boltz-cpus'),
  bulkSettingsBoltzMem: document.querySelector('#bulk-settings-boltz-mem'),
  bulkSettingsBoltzTime: document.querySelector('#bulk-settings-boltz-time'),
  bulkSettingsBoltzDesigns: document.querySelector('#bulk-settings-boltz-designs'),
  bulkSettingsBoltzScaffolds: document.querySelector('#bulk-settings-boltz-scaffolds'),
  bulkSettingsOpenaiKey: document.querySelector('#bulk-settings-openai-key'),
  bulkSettingsOpenaiModel: document.querySelector('#bulk-settings-openai-model'),
  bulkSettingsInputPath: document.querySelector('#bulk-settings-input-path'),
  bulkSettingsAutoLoad: document.querySelector('#bulk-settings-auto-load'),
  bulkSettingsLoadDefault: document.querySelector('#bulk-settings-load-default'),
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

function normalizePdbId(value) {
  const text = String(value || '').trim().toUpperCase().replace(/[^0-9A-Z]/g, '');
  return text || '';
}

function normalizePdbIdForPymol(value) {
  const text = String(value || '').trim();
  if (!text) return '';
  const base = text.split(/[\\/]/).pop() || text;
  const stripped = base.replace(/\.(mmcif|cif|pdb)$/i, '');
  const cleaned = stripped.toUpperCase().replace(/[^0-9A-Z]/g, '');
  if (cleaned.length < 4) return '';
  return cleaned.slice(0, 4);
}

function catalogNameFromPath(value) {
  const text = String(value || '').trim();
  if (!text) return null;
  const base = text.split(/[\\/]/).pop() || '';
  if (!base || !/\.(csv|tsv)$/i.test(base)) return null;
  return base;
}

function getConfiguredLlmCatalogName() {
  const fromSettingsField = catalogNameFromPath(el.bulkSettingsInputPath?.value);
  if (fromSettingsField) return fromSettingsField;
  const fromUiSettings = catalogNameFromPath(state.uiSettings?.input?.default_input_path);
  if (fromUiSettings) return fromUiSettings;
  return null;
}

function hashUnmatchedKey(raw) {
  let hash = 0;
  const text = String(raw || '');
  for (let idx = 0; idx < text.length; idx += 1) {
    hash = (hash * 31 + text.charCodeAt(idx)) >>> 0;
  }
  return hash.toString(16);
}

function buildUnmatchedKey(entry = {}, fallbackIndex = 0) {
  const candidate = entry?.candidate || {};
  const basis = [
    candidate.target_name,
    candidate.protein_name,
    candidate.gene,
    candidate.uniprot,
    candidate.pdb_id,
    candidate.antigen_catalog,
    candidate.accession,
    entry?.reason,
    fallbackIndex,
  ]
    .map((value) => String(value || '').trim().toLowerCase())
    .join('|');
  return `unm_${hashUnmatchedKey(basis)}`;
}

function normalizeUnmatchedAction(raw = null) {
  const status = String(raw?.status || '').toLowerCase();
  const allowed = new Set(['idle', 'queued', 'running', 'success', 'failed']);
  return {
    status: allowed.has(status) ? status : 'idle',
    job_id: raw?.job_id ? String(raw.job_id) : null,
    message: raw?.message ? String(raw.message) : null,
    updated_at: Number(raw?.updated_at || llmNowTs()),
  };
}

function normalizeUnmatchedList(list = []) {
  const out = [];
  const seen = new Set();
  (Array.isArray(list) ? list : []).forEach((entry, idx) => {
    if (!entry || typeof entry !== 'object') return;
    const base = String(entry.unmatched_key || '').trim() || buildUnmatchedKey(entry, idx);
    let key = base;
    let suffix = 1;
    while (seen.has(key)) {
      key = `${base}_${suffix}`;
      suffix += 1;
    }
    seen.add(key);
    out.push({ ...entry, unmatched_key: key });
  });
  return out;
}

function llmNowTs() {
  return Date.now();
}

function llmConversationId() {
  return `conv_${llmNowTs()}_${Math.random().toString(36).slice(2, 8)}`;
}

function llmConversationTitleFromHistory(history = [], fallback = null) {
  const firstUser = (history || []).find((msg) => String(msg?.role || '').toLowerCase() === 'user');
  const text = String(firstUser?.content || '').trim();
  if (text) {
    return text.length > 48 ? `${text.slice(0, 48)}...` : text;
  }
  return fallback || 'Untitled conversation';
}

function buildLlmConversation(seed = {}) {
  const createdAt = Number(seed.created_at || llmNowTs());
  const history = Array.isArray(seed.history) ? seed.history : [];
  const matchedRows = Array.isArray(seed.matched_rows) ? seed.matched_rows : [];
  const unmatched = normalizeUnmatchedList(seed.unmatched);
  const unmatchedActionsRaw = (seed && typeof seed.unmatched_actions === 'object' && seed.unmatched_actions)
    ? seed.unmatched_actions
    : {};
  const unmatchedActions = {};
  Object.entries(unmatchedActionsRaw).forEach(([key, value]) => {
    if (!key) return;
    unmatchedActions[String(key)] = normalizeUnmatchedAction(value);
  });
  unmatched.forEach((entry) => {
    if (!entry?.unmatched_key) return;
    if (!unmatchedActions[entry.unmatched_key]) {
      unmatchedActions[entry.unmatched_key] = normalizeUnmatchedAction();
    }
  });
  const deleted = Array.isArray(seed.deleted_pdb_ids)
    ? seed.deleted_pdb_ids.map((id) => normalizePdbId(id)).filter(Boolean)
    : [];
  const picked = Array.isArray(seed.picked_pdb_ids)
    ? seed.picked_pdb_ids.map((id) => normalizePdbId(id)).filter(Boolean)
    : [];
  const title = String(seed.title || '').trim() || llmConversationTitleFromHistory(history, null);
  return {
    id: String(seed.id || llmConversationId()),
    title: title || 'Untitled conversation',
    created_at: createdAt,
    updated_at: Number(seed.updated_at || createdAt),
    history,
    matched_rows: matchedRows,
    unmatched,
    unmatched_actions: unmatchedActions,
    picked_pdb_ids: picked,
    deleted_pdb_ids: deleted,
    view_mode: seed.view_mode === 'picked' ? 'picked' : 'all',
    catalog_filter: seed.catalog_filter === 'all' ? 'all' : 'biotin',
    active_catalog_name: seed.active_catalog_name ? String(seed.active_catalog_name) : null,
  };
}

function getActiveConversation() {
  const list = Array.isArray(state.llmConversations) ? state.llmConversations : [];
  if (!state.activeConversationId) return null;
  return list.find((conv) => conv.id === state.activeConversationId) || null;
}

function syncStateFromConversation(conv) {
  if (!conv) return;
  state.llmSuggestPending = false;
  state.llmHistory = Array.isArray(conv.history) ? conv.history : [];
  state.llmMatchedRows = Array.isArray(conv.matched_rows) ? conv.matched_rows : [];
  state.llmUnmatched = normalizeUnmatchedList(conv.unmatched);
  state.llmUnmatchedActions = {};
  const actionMap = (conv && typeof conv.unmatched_actions === 'object' && conv.unmatched_actions)
    ? conv.unmatched_actions
    : {};
  state.llmUnmatched.forEach((entry) => {
    const key = String(entry?.unmatched_key || '').trim();
    if (!key) return;
    state.llmUnmatchedActions[key] = normalizeUnmatchedAction(actionMap[key]);
  });
  state.llmPickedPdbIds = Array.isArray(conv.picked_pdb_ids) ? conv.picked_pdb_ids : [];
  state.llmDeletedPdbIds = Array.isArray(conv.deleted_pdb_ids) ? conv.deleted_pdb_ids : [];
  state.llmViewMode = conv.view_mode === 'picked' ? 'picked' : 'all';
  state.llmCatalogFilter = conv.catalog_filter === 'all' ? 'all' : 'biotin';
  state.activeCatalogName = getConfiguredLlmCatalogName() || null;
}

function syncConversationFromState() {
  const conv = getActiveConversation();
  if (!conv) return;
  conv.history = Array.isArray(state.llmHistory) ? state.llmHistory : [];
  conv.matched_rows = Array.isArray(state.llmMatchedRows) ? state.llmMatchedRows : [];
  conv.unmatched = normalizeUnmatchedList(state.llmUnmatched);
  const normalizedActions = {};
  Object.entries(state.llmUnmatchedActions || {}).forEach(([key, value]) => {
    if (!key) return;
    normalizedActions[String(key)] = normalizeUnmatchedAction(value);
  });
  conv.unmatched_actions = normalizedActions;
  conv.picked_pdb_ids = Array.isArray(state.llmPickedPdbIds) ? state.llmPickedPdbIds : [];
  conv.deleted_pdb_ids = Array.isArray(state.llmDeletedPdbIds) ? state.llmDeletedPdbIds : [];
  conv.view_mode = state.llmViewMode === 'picked' ? 'picked' : 'all';
  conv.catalog_filter = state.llmCatalogFilter === 'all' ? 'all' : 'biotin';
  conv.active_catalog_name = getConfiguredLlmCatalogName() || null;
  conv.title = llmConversationTitleFromHistory(conv.history, conv.title || null);
  conv.updated_at = llmNowTs();
}

function renderConversationSelector() {
  if (!el.llmConversationSelect) return;
  const list = Array.isArray(state.llmConversations) ? state.llmConversations : [];
  const sorted = [...list].sort((a, b) => Number(b.updated_at || 0) - Number(a.updated_at || 0));
  el.llmConversationSelect.innerHTML = '';
  if (!sorted.length) {
    const opt = document.createElement('option');
    opt.value = '';
    opt.textContent = 'No conversations';
    el.llmConversationSelect.appendChild(opt);
    return;
  }
  sorted.forEach((conv, idx) => {
    const opt = document.createElement('option');
    opt.value = conv.id;
    const title = String(conv.title || '').trim() || `Conversation ${idx + 1}`;
    const ts = new Date(Number(conv.updated_at || conv.created_at || llmNowTs())).toLocaleString();
    const pickedCount = Array.isArray(conv.picked_pdb_ids) ? conv.picked_pdb_ids.length : 0;
    opt.textContent = `${title} · ${pickedCount} targets (${ts})`;
    el.llmConversationSelect.appendChild(opt);
  });
  el.llmConversationSelect.value = state.activeConversationId || sorted[0].id;
}

function ensureActiveConversation() {
  if (!Array.isArray(state.llmConversations)) {
    state.llmConversations = [];
  }
  let conv = getActiveConversation();
  if (conv) return conv;
  if (state.llmConversations.length) {
    state.activeConversationId = state.llmConversations[0].id;
    conv = getActiveConversation();
    if (conv) return conv;
  }
  const fresh = buildLlmConversation({
    title: 'New conversation',
  });
  state.llmConversations.push(fresh);
  state.activeConversationId = fresh.id;
  return fresh;
}

function recomputeLlmPickedPdbIds() {
  const deleted = new Set(
    (state.llmDeletedPdbIds || []).map((id) => normalizePdbId(id)).filter(Boolean),
  );
  const ids = Array.from(new Set(
    (state.llmMatchedRows || [])
      .map((row) => normalizePdbId(row?.resolved_pdb_id || row?.pdb_id))
      .filter((id) => id && !deleted.has(id)),
  ));
  state.llmPickedPdbIds = ids;
}

function setActiveConversation(conversationId, { persist = true, rerender = true } = {}) {
  const list = Array.isArray(state.llmConversations) ? state.llmConversations : [];
  const target = list.find((conv) => conv.id === conversationId);
  if (!target) return;
  syncConversationFromState();
  state.activeConversationId = target.id;
  syncStateFromConversation(target);
  recomputeLlmPickedPdbIds();
  setLlmSuggestPending(false, { persist: false });
  renderConversationSelector();
  setLlmViewMode(state.llmViewMode || 'all', { persist: false, rerender: false });
  setLlmCatalogFilter(state.llmCatalogFilter || 'biotin', { persist: false, rerender: false });
  if (persist) persistLlmState();
  ensureLlmUnmatchedPolling();
  if (rerender) {
    renderLlmTranscript();
    renderLlmMatchPanels();
    renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
    renderBoltzConfigs();
    updateRunCommandFilterNote();
  }
}

function createNewConversation() {
  syncConversationFromState();
  const fresh = buildLlmConversation({
    title: 'New conversation',
    view_mode: 'all',
  });
  state.llmConversations.push(fresh);
  state.activeConversationId = fresh.id;
  syncStateFromConversation(fresh);
  recomputeLlmPickedPdbIds();
  setLlmSuggestPending(false, { persist: false });
  setLlmViewMode('all', { persist: false, rerender: false });
  setLlmCatalogFilter('biotin', { persist: false, rerender: false });
  renderConversationSelector();
  renderLlmTranscript();
  renderLlmMatchPanels();
  persistLlmState();
  ensureLlmUnmatchedPolling();
  renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
  renderBoltzConfigs();
  updateRunCommandFilterNote();
}

function getLlmPickedPdbSet() {
  const ids = Array.isArray(state.llmPickedPdbIds) ? state.llmPickedPdbIds : [];
  return new Set(ids.map((id) => normalizePdbId(id)).filter(Boolean));
}

function normalizeLlmCatalogFilter(mode) {
  return String(mode || '').toLowerCase() === 'all' ? 'all' : 'biotin';
}

function parseBiotinFlag(value) {
  if (typeof value === 'boolean') return value;
  const raw = (value === null || value === undefined)
    ? ''
    : String(value).trim().toLowerCase();
  if (!raw) return null;
  if (['1', 'true', 't', 'yes', 'y'].includes(raw)) return true;
  if (['0', 'false', 'f', 'no', 'n'].includes(raw)) return false;
  return null;
}

function rowPassesCatalogFilter(row) {
  const mode = normalizeLlmCatalogFilter(state.llmCatalogFilter);
  if (mode === 'all') return true;
  const selection = String(row?.selection || '').trim().toLowerCase();
  if (selection === 'biotin') return true;
  return parseBiotinFlag(row?.biotinylated) === true;
}

function getVisibleBulkRows(rows = null) {
  const source = Array.isArray(rows) ? rows : (Array.isArray(state.bulkPreviewRows) ? state.bulkPreviewRows : []);
  const byCatalog = source.filter((row) => rowPassesCatalogFilter(row));
  if (state.llmViewMode !== 'picked') return byCatalog;
  const picked = getLlmPickedPdbSet();
  if (!picked.size) {
    const hasMatched = Array.isArray(state.llmMatchedRows) && state.llmMatchedRows.length > 0;
    return hasMatched ? [] : byCatalog;
  }
  return byCatalog.filter((row) => {
    const pdb = normalizePdbId(row?.resolved_pdb_id || row?.pdb_id);
    return pdb && picked.has(pdb);
  });
}

function getVisibleBoltzTargets(targets = null) {
  const source = Array.isArray(targets) ? targets : (Array.isArray(state.boltzConfigs) ? state.boltzConfigs : []);
  const byCatalog = source.filter((target) => rowPassesCatalogFilter(target));
  if (state.llmViewMode !== 'picked') return byCatalog;
  const picked = getLlmPickedPdbSet();
  if (!picked.size) {
    const hasMatched = Array.isArray(state.llmMatchedRows) && state.llmMatchedRows.length > 0;
    return hasMatched ? [] : byCatalog;
  }
  return byCatalog.filter((target) => {
    const pdb = normalizePdbId(target?.pdb_id);
    return pdb && picked.has(pdb);
  });
}

function persistLlmState() {
  try {
    syncConversationFromState();
    const list = Array.isArray(state.llmConversations) ? state.llmConversations : [];
    const payload = {
      version: 4,
      active_conversation_id: state.activeConversationId || null,
      conversations: list.map((conv) => ({
        id: conv.id,
        title: conv.title,
        created_at: conv.created_at,
        updated_at: conv.updated_at,
        history: Array.isArray(conv.history) ? conv.history : [],
        matched_rows: Array.isArray(conv.matched_rows) ? conv.matched_rows : [],
        unmatched: normalizeUnmatchedList(conv.unmatched),
        unmatched_actions: (conv && typeof conv.unmatched_actions === 'object' && conv.unmatched_actions)
          ? conv.unmatched_actions
          : {},
        picked_pdb_ids: Array.isArray(conv.picked_pdb_ids) ? conv.picked_pdb_ids : [],
        deleted_pdb_ids: Array.isArray(conv.deleted_pdb_ids) ? conv.deleted_pdb_ids : [],
        view_mode: conv.view_mode === 'picked' ? 'picked' : 'all',
        catalog_filter: conv.catalog_filter === 'all' ? 'all' : 'biotin',
        active_catalog_name: conv.active_catalog_name || null,
      })),
    };
    localStorage.setItem(LLM_LOCAL_STORAGE_KEY, JSON.stringify(payload));
  } catch (err) {
    // Ignore storage failures in private mode/quota limits.
  }
}

function restoreLlmState() {
  try {
    const raw = localStorage.getItem(LLM_LOCAL_STORAGE_KEY);
    if (!raw) {
      const fresh = buildLlmConversation({ title: 'New conversation' });
      state.llmConversations = [fresh];
      state.activeConversationId = fresh.id;
      syncStateFromConversation(fresh);
      return;
    }
    const data = JSON.parse(raw);
    if (!data || typeof data !== 'object') {
      const fresh = buildLlmConversation({ title: 'New conversation' });
      state.llmConversations = [fresh];
      state.activeConversationId = fresh.id;
      syncStateFromConversation(fresh);
      return;
    }
    if (Array.isArray(data.conversations)) {
      state.llmConversations = data.conversations.map((conv) => buildLlmConversation(conv));
      state.activeConversationId = String(data.active_conversation_id || '');
    } else {
      // Legacy single-conversation payload migration.
      const migrated = buildLlmConversation({
        title: llmConversationTitleFromHistory(data.history || [], 'Migrated conversation'),
        history: Array.isArray(data.history) ? data.history : [],
        matched_rows: Array.isArray(data.matched_rows) ? data.matched_rows : [],
        unmatched: Array.isArray(data.unmatched) ? data.unmatched : [],
        picked_pdb_ids: Array.isArray(data.picked_pdb_ids) ? data.picked_pdb_ids : [],
        view_mode: data.view_mode === 'picked' ? 'picked' : 'all',
        catalog_filter: data.catalog_filter === 'all' ? 'all' : 'biotin',
        active_catalog_name: data.active_catalog_name ? String(data.active_catalog_name) : null,
      });
      state.llmConversations = [migrated];
      state.activeConversationId = migrated.id;
    }
    ensureActiveConversation();
    syncStateFromConversation(getActiveConversation());
    recomputeLlmPickedPdbIds();
  } catch (err) {
    const fresh = buildLlmConversation({ title: 'New conversation' });
    state.llmConversations = [fresh];
    state.activeConversationId = fresh.id;
    syncStateFromConversation(fresh);
  }
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
    wrapper.className = 'bulk-command-block';

    const header = document.createElement('div');
    header.className = 'bulk-command-header';

    const heading = document.createElement('div');
    heading.className = 'bulk-command-title';
    heading.textContent = 'Commands';
    header.appendChild(heading);

    const copyBtn = document.createElement('button');
    copyBtn.type = 'button';
    copyBtn.className = 'ghost';
    copyBtn.textContent = 'Copy';
    copyBtn.addEventListener('click', () => copyTextToClipboard(String(text), copyBtn));
    header.appendChild(copyBtn);

    const block = document.createElement('pre');
    block.className = 'bulk-command-pre';
    block.textContent = String(text).trim() || '—';

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
    wrapper.className = 'bulk-command-block';

    const header = document.createElement('div');
    header.className = 'bulk-command-header';

    const heading = document.createElement('div');
    heading.className = 'bulk-command-title';
    heading.textContent = title;
    header.appendChild(heading);

    const copyBtn = document.createElement('button');
    copyBtn.type = 'button';
    copyBtn.className = 'ghost';
    copyBtn.textContent = 'Copy';
    copyBtn.addEventListener('click', () => copyTextToClipboard(commandText, copyBtn));
    header.appendChild(copyBtn);

    const block = document.createElement('pre');
    block.className = 'bulk-command-pre';
    block.textContent = commandText || '—';

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
    btn.classList.remove('pymol-unavailable');
    btn.title = '';
    return;
  }
  btn.disabled = true;
  btn.classList.add('pymol-unavailable');
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
    img.className = 'bulk-image';
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
      const swatch = document.createElement('span');
      swatch.className = 'bulk-swatch';
      swatch.style.backgroundColor = color;
      li.appendChild(swatch);
      const nameEl = document.createElement('strong');
      nameEl.textContent = name;
      li.appendChild(nameEl);
      li.appendChild(document.createTextNode(`: ${residueText}`));
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
    img.className = 'bulk-image';
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
    dirNote.className = 'help-text bulk-diversity-note';
    dirNote.textContent = `Output directory: ${outputDir}`;
    el.diversityGrid.appendChild(dirNote);
  }
  const files = Array.isArray(state.diversityFiles) ? state.diversityFiles : [];
  if (files.length) {
    const details = document.createElement('details');
    details.className = 'bulk-diversity-files';
    const summary = document.createElement('summary');
    summary.className = 'bulk-diversity-summary';
    summary.textContent = `Detected ${files.length} all_designs_metrics.csv file${files.length === 1 ? '' : 's'}`;
    details.appendChild(summary);

    const pre = document.createElement('pre');
    pre.className = 'bulk-diversity-pre';
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
    empty.className = 'help-text bulk-doc-help-zero';
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
      img.className = 'bulk-image';
      imgBox.appendChild(img);
      card.appendChild(imgBox);
    }

    const legend = document.createElement('div');
    legend.className = 'snapshot-note';
    const colors = plot.epitope_colors || {};
    if (Object.keys(colors).length) {
      Object.entries(colors).forEach(([name, color]) => {
        const item = document.createElement('span');
        item.className = 'bulk-legend-item';
        const dot = document.createElement('span');
        dot.className = 'bulk-legend-dot';
        dot.style.backgroundColor = color;
        item.appendChild(dot);
        item.appendChild(document.createTextNode(name));
        legend.appendChild(item);
      });
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
    td.className = 'bulk-empty-row';
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
      pymolBtn.className = 'bulk-pymol-btn';
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
    const plotExports = Array.isArray(body.plot_exports) ? body.plot_exports : [];
    for (const item of plotExports) {
      if (!item || typeof item !== 'object') continue;
      const files = [item.png_name, item.svg_name, item.map_csv_name];
      for (const name of files) {
        if (!name) continue;
        await downloadNamedFile(String(name));
      }
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
  const normalizedPdbId = normalizePdbIdForPymol(row.pdb_id);
  if (!normalizedPdbId) {
    showAlert(`Invalid PDB ID for PyMOL: ${row.pdb_id || 'missing value'}.`);
    return;
  }
  try {
    const epLabel = row.epitope_id || row.epitope || null;
    const payload = {
      pdb_id: normalizedPdbId,
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
  const normalizedPdbId = normalizePdbIdForPymol(pdbId);
  if (!normalizedPdbId) {
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
    const res = await fetch(`/api/targets/${encodeURIComponent(normalizedPdbId)}/pymol/hotspots`, {
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
    showAlert(`PyMOL hotspot bundle ready (${normalizedPdbId} · ${scopeText})`, false);
    appendLog(body.message || `PyMOL bundle ready for ${normalizedPdbId}`);
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
  const targets = getVisibleBoltzTargets();
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
  const targets = getVisibleBoltzTargets();
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
  const targets = getVisibleBoltzTargets();
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
  const targets = getVisibleBoltzTargets();
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
  await loadCommandDefaults();
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
  const allTargets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  const targets = getVisibleBoltzTargets(allTargets);
  if (!allTargets.length) {
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
  if (!targets.length) {
    const emptyRow = document.createElement('tr');
    const emptyCell = document.createElement('td');
    emptyCell.colSpan = 8;
    emptyCell.className = 'empty-note';
    const viewNote = state.llmViewMode === 'picked'
      ? 'No LLM-picked targets are currently visible.'
      : 'No targets are currently visible.';
    const filterNote = state.llmCatalogFilter === 'biotin'
      ? ' Try Catalog filter -> All.'
      : '';
    emptyCell.textContent = `${viewNote}${filterNote} Switch View mode to All targets if needed.`;
    emptyRow.appendChild(emptyCell);
    tbody.appendChild(emptyRow);
    if (el.boltzSummary) {
      const modeLabel = state.llmViewMode === 'picked' ? 'LLM-picked' : 'all-target';
      const catalogLabel = state.llmCatalogFilter === 'biotin' ? 'biotin' : 'all';
      el.boltzSummary.textContent = `Showing 0 of ${allTargets.length} target${allTargets.length === 1 ? '' : 's'} (${modeLabel}, catalog=${catalogLabel})`;
      el.boltzSummary.hidden = false;
    }
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
    pymolBtn.className = 'bulk-pymol-btn';
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
        if (cellIdx === 5) td.classList.add('epitope-row-indent');
        epRow.appendChild(td);
      });
      if (cfg?.hotspot_surface_ok === false) {
        epRow.classList.add('is-filtered');
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
      epPymol.className = 'bulk-pymol-btn';
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
    const baseSummary = `${configCount} config${configCount === 1 ? '' : 's'} across ${targets.length} target${targets.length === 1 ? '' : 's'}`;
    const filtered = targets.length !== allTargets.length;
    const modeLabel = state.llmViewMode === 'picked' ? 'LLM-picked' : 'all-target';
    const catalogLabel = state.llmCatalogFilter === 'biotin' ? 'biotin' : 'all';
    el.boltzSummary.textContent = filtered
      ? `${baseSummary} (filtered from ${allTargets.length}; mode=${modeLabel}; catalog=${catalogLabel})`
      : `${baseSummary} (mode=${modeLabel}; catalog=${catalogLabel})`;
    el.boltzSummary.hidden = false;
  }
  el.boltzPanel.hidden = false;
}

async function loadCommandDefaults() {
  if (state.commandDefaults) return state.commandDefaults;
  try {
    const res = await fetch('/api/bulk/command-defaults');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    state.commandDefaults = await res.json();
  } catch (err) {
    state.commandDefaults = {};
  }
  return state.commandDefaults;
}

function parseOptionalInteger(raw, { min = null, max = null } = {}) {
  const text = String(raw ?? '').trim();
  if (!text) return null;
  const num = Number(text);
  if (!Number.isFinite(num)) return null;
  const rounded = Math.round(num);
  if (min !== null && rounded < min) return null;
  if (max !== null && rounded > max) return null;
  return rounded;
}

function normalizeOptionalText(raw) {
  const text = String(raw ?? '').trim();
  return text || null;
}

function setTextInputValue(inputEl, value) {
  if (!inputEl) return;
  inputEl.value = value ?? '';
}

function normalizeTextList(rawList) {
  if (!Array.isArray(rawList)) return [];
  return rawList
    .map((item) => String(item ?? '').trim())
    .filter((item) => item.length > 0);
}

function parseMultilineList(raw) {
  const text = String(raw ?? '').trim();
  if (!text) return [];
  return text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0);
}

function applyBulkUiSettingsToForm(payload = {}) {
  const cluster = payload.cluster || {};
  const boltzgen = payload.boltzgen || {};
  const llm = payload.llm || {};
  const input = payload.input || {};
  setTextInputValue(el.bulkSettingsSshAlias, cluster.ssh_config_alias);
  setTextInputValue(el.bulkSettingsRemoteRoot, cluster.remote_root);
  setTextInputValue(el.bulkSettingsTargetRoot, cluster.target_root);
  setTextInputValue(el.bulkSettingsConda, cluster.conda_activate);
  setTextInputValue(el.bulkSettingsPymolPath, cluster.pymol_path);
  setTextInputValue(el.bulkSettingsPymolCondaEnv, cluster.pymol_conda_env);

  setTextInputValue(el.bulkSettingsBoltzPartition, boltzgen.partition);
  setTextInputValue(el.bulkSettingsBoltzAccount, boltzgen.account);
  setTextInputValue(el.bulkSettingsBoltzGpus, boltzgen.gpus);
  setTextInputValue(el.bulkSettingsBoltzCpus, boltzgen.cpus);
  setTextInputValue(el.bulkSettingsBoltzMem, boltzgen.mem_gb);
  setTextInputValue(el.bulkSettingsBoltzTime, boltzgen.time_hours);
  setTextInputValue(el.bulkSettingsBoltzDesigns, boltzgen.default_num_designs);
  setTextInputValue(el.bulkSettingsBoltzScaffolds, normalizeTextList(boltzgen.nanobody_scaffolds).join('\n'));

  setTextInputValue(el.bulkSettingsOpenaiKey, llm.openai_api_key);
  setTextInputValue(el.bulkSettingsOpenaiModel, llm.openai_model);
  setTextInputValue(el.bulkSettingsInputPath, input.default_input_path);
  state.activeCatalogName = getConfiguredLlmCatalogName();
  if (el.bulkSettingsAutoLoad) el.bulkSettingsAutoLoad.checked = Boolean(input.auto_load_default_input);

  if (el.bulkSettingsMeta) {
    const path = payload.local_config_path || 'cfg/webapp.local.yaml';
    el.bulkSettingsMeta.innerHTML = `Edits are saved to <code>${path}</code>.`;
  }

  const defaultDesignCount = parseOptionalInteger(boltzgen.default_num_designs, { min: 1, max: 50000 });
  if (defaultDesignCount && el.boltzDesignCount) {
    el.boltzDesignCount.value = String(defaultDesignCount);
  }
  const defaultTimeHours = parseOptionalInteger(boltzgen.time_hours, { min: 1, max: 240 });
  if (defaultTimeHours && el.boltzTimeHours) {
    el.boltzTimeHours.value = String(defaultTimeHours);
  }
}

function buildBulkUiSettingsPayload() {
  return {
    cluster: {
      mock: false,
      ssh_config_alias: normalizeOptionalText(el.bulkSettingsSshAlias?.value),
      remote_root: normalizeOptionalText(el.bulkSettingsRemoteRoot?.value),
      target_root: normalizeOptionalText(el.bulkSettingsTargetRoot?.value),
      conda_activate: normalizeOptionalText(el.bulkSettingsConda?.value),
      pymol_path: normalizeOptionalText(el.bulkSettingsPymolPath?.value),
      pymol_conda_env: normalizeOptionalText(el.bulkSettingsPymolCondaEnv?.value),
    },
    boltzgen: {
      partition: normalizeOptionalText(el.bulkSettingsBoltzPartition?.value),
      account: normalizeOptionalText(el.bulkSettingsBoltzAccount?.value),
      gpus: normalizeOptionalText(el.bulkSettingsBoltzGpus?.value),
      cpus: parseOptionalInteger(el.bulkSettingsBoltzCpus?.value, { min: 1, max: 256 }),
      mem_gb: parseOptionalInteger(el.bulkSettingsBoltzMem?.value, { min: 1, max: 2048 }),
      time_hours: parseOptionalInteger(el.bulkSettingsBoltzTime?.value, { min: 1, max: 240 }),
      default_num_designs: parseOptionalInteger(el.bulkSettingsBoltzDesigns?.value, { min: 1, max: 50000 }),
      nanobody_scaffolds: parseMultilineList(el.bulkSettingsBoltzScaffolds?.value),
    },
    llm: {
      openai_api_key: normalizeOptionalText(el.bulkSettingsOpenaiKey?.value),
      openai_model: normalizeOptionalText(el.bulkSettingsOpenaiModel?.value),
    },
    input: {
      default_input_path: normalizeOptionalText(el.bulkSettingsInputPath?.value),
      auto_load_default_input: Boolean(el.bulkSettingsAutoLoad?.checked),
    },
  };
}

async function loadBulkUiSettings({ autoLoadInput = true } = {}) {
  try {
    const res = await fetch('/api/bulk/ui-config');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const payload = await res.json();
    state.uiSettings = payload;
    applyBulkUiSettingsToForm(payload);
    if (autoLoadInput) {
      const configuredInput = normalizeOptionalText(payload?.input?.default_input_path);
      if (configuredInput) {
        await loadDefaultInputFile({ force: true, silent: false });
      }
    }
  } catch (err) {
    state.uiSettings = null;
  }
}

async function saveBulkUiSettings(triggerEl = null) {
  const payload = buildBulkUiSettingsPayload();
  if (triggerEl) triggerEl.disabled = true;
  try {
    const res = await fetch('/api/bulk/ui-config', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const body = await res.json().catch(() => ({}));
    if (!res.ok) {
      const message = body?.detail || body?.message || `HTTP ${res.status}`;
      throw new Error(message);
    }
    state.uiSettings = body;
    applyBulkUiSettingsToForm(body);
    state.commandDefaults = null;
    await loadCommandDefaults();
    const configuredInput = normalizeOptionalText(body?.input?.default_input_path);
    if (configuredInput) {
      await loadDefaultInputFile({ force: true, silent: false });
    }
    showAlert('Bulk settings saved.', false);
    toggleModal(el.bulkSettingsModal, false);
  } catch (err) {
    showAlert(`Failed to save settings: ${err.message || err}`);
  } finally {
    if (triggerEl) triggerEl.disabled = false;
  }
}

async function loadDefaultInputFile({ force = false, silent = false } = {}) {
  const hasInput = Boolean((el.bulkCsvInput?.value || '').trim());
  if (!force && hasInput) return;
  try {
    const res = await fetch('/api/bulk/default-input');
    const body = await res.json().catch(() => ({}));
    if (!res.ok) {
      const message = body?.detail || body?.message || `HTTP ${res.status}`;
      throw new Error(message);
    }
    if (el.bulkCsvInput) el.bulkCsvInput.value = body.text || '';
    state.bulkPreviewSig = null;
    await previewBulkCsv({ silent: true });
    if (!silent) {
      const shownPath = body.path || 'configured default input';
      showAlert(`Loaded default input from ${shownPath}.`, false);
    }
  } catch (err) {
    if (!silent) showAlert(`Failed to load default input: ${err.message || err}`);
  }
}

async function openBulkReadmeModal() {
  if (el.bulkReadmeBody) {
    el.bulkReadmeBody.textContent = 'Loading README...';
  }
  toggleModal(el.bulkReadmeModal, true);
  try {
    if (!state.guiReadme) {
      const res = await fetch('/api/bulk/readme-gui');
      const body = await res.json().catch(() => ({}));
      if (!res.ok) {
        throw new Error(body?.detail || body?.message || `HTTP ${res.status}`);
      }
      state.guiReadme = body;
    }
    const path = state.guiReadme?.path || 'README_GUI.md';
    const text = state.guiReadme?.text || 'README_GUI.md is empty.';
    if (el.bulkReadmeTitle) {
      el.bulkReadmeTitle.textContent = `GUI README · ${path.split('/').pop() || path}`;
    }
    if (el.bulkReadmeBody) {
      el.bulkReadmeBody.textContent = text;
    }
  } catch (err) {
    if (el.bulkReadmeTitle) {
      el.bulkReadmeTitle.textContent = 'GUI README';
    }
    if (el.bulkReadmeBody) {
      el.bulkReadmeBody.textContent = `Failed to load README: ${err.message || err}`;
    }
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
    const previewMap = new Map();
    (state.bulkPreviewRows || []).forEach((row) => {
      const key = (row.resolved_pdb_id || row.pdb_id || '').toUpperCase();
      if (key) previewMap.set(key, row);
    });
    targets.forEach((target) => {
      const key = (target.pdb_id || '').toUpperCase();
      const previewRow = previewMap.get(key);
      if (previewRow) {
        target.preset_name = previewRow.preset_name || previewRow.gene || '';
        target.selection = previewRow.selection ?? null;
        target.biotinylated = previewRow.biotinylated ?? null;
        target.tags = previewRow.tags ?? null;
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
  const visibleTargets = getVisibleBoltzTargets();
  const candidates = Array.isArray(visibleTargets) && visibleTargets.length
    ? visibleTargets
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
  const targets = getVisibleBoltzTargets();
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
  const info = state.commandDefaults || {};
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
  await loadCommandDefaults();
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
  await loadCommandDefaults();
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
  const targets = getVisibleBoltzTargets();
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
  await loadCommandDefaults();
  const engine = getRunEngine();
  if (engine === 'rfantibody') {
    await showRfaRunCommandAll();
    return;
  }
  const targets = getVisibleBoltzTargets();
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
  await loadCommandDefaults();

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
  const allRows = Array.isArray(rows) ? rows : [];
  const list = getVisibleBulkRows(allRows);
  if (!allRows.length) {
    el.bulkPreviewPanel.hidden = true;
    return;
  }
  if (!list.length) {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = 10;
    td.className = 'empty-note';
    const catalogSuffix = state.llmCatalogFilter === 'biotin'
      ? ' Try Catalog filter -> All.'
      : '';
    td.textContent = state.llmViewMode === 'picked'
      ? `No LLM-picked targets found in current preview.${catalogSuffix} Switch View mode to All targets if needed.`
      : `No targets detected.${catalogSuffix}`;
    tr.appendChild(td);
    tbody.appendChild(tr);
  }
  list.forEach((row) => {
    const tr = document.createElement('tr');
    const biotinValue = parseBiotinFlag(row?.biotinylated);
    const biotinText = biotinValue === true ? 'Yes' : (biotinValue === false ? 'No' : '—');
    const values = [
      row.raw_index ?? '',
      row.preset_name || '',
      row.protein_name || '',
      row.antigen_url || '',
      row.accession || '',
      row.selection || '—',
      biotinText,
      row.tags || '—',
      row.resolved_pdb_id || row.pdb_id || '—',
      row.expression_host || '—',
    ];
    values.forEach((value) => {
      const td = document.createElement('td');
      td.textContent = value;
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
  if (el.bulkPreviewSummary) {
    const baseSummary = summary || `${allRows.length} rows parsed`;
    const viewSummary = state.llmViewMode === 'picked'
      ? `showing ${list.length} LLM-picked`
      : `showing ${list.length} targets`;
    const catalogSummary = state.llmCatalogFilter === 'biotin'
      ? 'Catalog filter: biotin only'
      : 'Catalog filter: all';
    el.bulkPreviewSummary.textContent = `${baseSummary} · ${viewSummary} · ${catalogSummary}`;
    el.bulkPreviewSummary.hidden = false;
  }
  el.bulkPreviewPanel.hidden = false;
}

function renderLlmSuggestLoading() {
  if (el.llmSuggestLoading) {
    el.llmSuggestLoading.hidden = !state.llmSuggestPending;
  }
}

function setLlmSuggestPending(next, { persist = true } = {}) {
  state.llmSuggestPending = Boolean(next);
  if (el.llmSuggestTargets) {
    el.llmSuggestTargets.disabled = state.llmSuggestPending;
    el.llmSuggestTargets.setAttribute('aria-busy', state.llmSuggestPending ? 'true' : 'false');
  }
  if (el.llmTargetPanel) {
    el.llmTargetPanel.setAttribute('aria-busy', state.llmSuggestPending ? 'true' : 'false');
  }
  if (el.llmChatPrompt) {
    el.llmChatPrompt.setAttribute('aria-busy', state.llmSuggestPending ? 'true' : 'false');
  }
  renderLlmSuggestLoading();
  if (persist) persistLlmState();
}

function ensureLlmUnmatchedActions() {
  const next = {};
  const current = (state && typeof state.llmUnmatchedActions === 'object' && state.llmUnmatchedActions)
    ? state.llmUnmatchedActions
    : {};
  const normalized = normalizeUnmatchedList(state.llmUnmatched);
  state.llmUnmatched = normalized;
  normalized.forEach((entry) => {
    const key = String(entry?.unmatched_key || '').trim();
    if (!key) return;
    next[key] = normalizeUnmatchedAction(current[key]);
  });
  state.llmUnmatchedActions = next;
}

function getLlmUnmatchedAction(unmatchedKey) {
  const key = String(unmatchedKey || '').trim();
  if (!key) return normalizeUnmatchedAction();
  ensureLlmUnmatchedActions();
  return state.llmUnmatchedActions[key] || normalizeUnmatchedAction();
}

function setLlmUnmatchedAction(unmatchedKey, patch = {}) {
  const key = String(unmatchedKey || '').trim();
  if (!key) return;
  ensureLlmUnmatchedActions();
  const prev = state.llmUnmatchedActions[key] || normalizeUnmatchedAction();
  state.llmUnmatchedActions[key] = normalizeUnmatchedAction({
    ...prev,
    ...patch,
    updated_at: llmNowTs(),
  });
}

function mergeMatchedRowsInState(rows = []) {
  const mergedRows = new Map();
  (Array.isArray(state.llmMatchedRows) ? state.llmMatchedRows : []).forEach((row) => {
    const key = normalizePdbId(row?.resolved_pdb_id || row?.pdb_id)
      || `row_${row?.raw_index || ''}_${String(row?.preset_name || '')}`;
    mergedRows.set(key, row);
  });
  (Array.isArray(rows) ? rows : []).forEach((row) => {
    const key = normalizePdbId(row?.resolved_pdb_id || row?.pdb_id)
      || `row_${row?.raw_index || ''}_${String(row?.preset_name || '')}`;
    mergedRows.set(key, row);
  });
  state.llmMatchedRows = Array.from(mergedRows.values());
}

function mergeMatchedRowsIntoPreview(rows = []) {
  if (!Array.isArray(state.bulkPreviewRows) || !state.bulkPreviewRows.length) {
    state.bulkPreviewRows = Array.isArray(rows) ? [...rows] : [];
    return;
  }
  const existing = Array.isArray(state.bulkPreviewRows) ? state.bulkPreviewRows : [];
  const merged = [...existing];
  const existingIds = new Set(
    existing
      .map((row) => normalizePdbId(row?.resolved_pdb_id || row?.pdb_id))
      .filter(Boolean),
  );
  (Array.isArray(rows) ? rows : []).forEach((row) => {
    const pdb = normalizePdbId(row?.resolved_pdb_id || row?.pdb_id);
    if (pdb && !existingIds.has(pdb)) {
      merged.push(row);
      existingIds.add(pdb);
    }
  });
  state.bulkPreviewRows = merged;
}

function removeUnmatchedEntry(unmatchedKey) {
  const key = String(unmatchedKey || '').trim();
  if (!key) return;
  state.llmUnmatched = (Array.isArray(state.llmUnmatched) ? state.llmUnmatched : [])
    .filter((entry) => String(entry?.unmatched_key || '').trim() !== key);
  if (state.llmUnmatchedActions && typeof state.llmUnmatchedActions === 'object') {
    delete state.llmUnmatchedActions[key];
  }
}

function stopLlmUnmatchedPolling() {
  if (state.llmUnmatchedPoller) {
    clearInterval(state.llmUnmatchedPoller);
    state.llmUnmatchedPoller = null;
  }
}

function hasActiveLlmUnmatchedJobs() {
  ensureLlmUnmatchedActions();
  return Object.values(state.llmUnmatchedActions || {}).some((entry) => {
    const status = String(entry?.status || '').toLowerCase();
    return status === 'queued' || status === 'running';
  });
}

async function pollLlmUnmatchedJobs() {
  ensureLlmUnmatchedActions();
  const actionEntries = Object.entries(state.llmUnmatchedActions || {})
    .filter(([, entry]) => {
      const status = String(entry?.status || '').toLowerCase();
      return (status === 'queued' || status === 'running') && entry?.job_id;
    });
  if (!actionEntries.length) {
    stopLlmUnmatchedPolling();
    return;
  }
  const activeLogJobId = String(actionEntries[0]?.[1]?.job_id || '').trim();
  if (activeLogJobId) {
    await refreshDiscoveryJobLog(activeLogJobId);
  }

  let didUpdate = false;
  for (const [key, action] of actionEntries) {
    try {
      const res = await fetch(`/api/bulk/llm-targets/unmatched/discover/${encodeURIComponent(action.job_id)}`);
      const body = await res.json().catch(() => ({}));
      if (!res.ok) {
        throw new Error(body?.detail || body?.message || `Discovery status failed (${res.status})`);
      }
      const status = String(body?.status || '').toLowerCase();
      if (status === 'pending' || status === 'running') {
        const phase = String(body?.phase || '').trim();
        const attempts = Array.isArray(body?.attempts) ? body.attempts.length : 0;
        const queryPreview = Array.isArray(body?.planned_queries) && body.planned_queries.length
          ? ` · ${body.planned_queries.slice(0, 2).join(', ')}`
          : '';
        const phasePrefix = phase ? `[${phase}] ` : '';
        const attemptSuffix = attempts ? ` · attempts: ${attempts}` : '';
        setLlmUnmatchedAction(key, {
          status: status === 'pending' ? 'queued' : 'running',
          message: `${phasePrefix}${body?.message || 'Discovery in progress...'}${queryPreview}${attemptSuffix}`,
          job_id: action.job_id,
        });
        didUpdate = true;
        continue;
      }
      if (status === 'success') {
        await refreshDiscoveryJobLog(action.job_id);
        if (body?.matched_row && typeof body.matched_row === 'object') {
          mergeMatchedRowsInState([body.matched_row]);
          mergeMatchedRowsIntoPreview([body.matched_row]);
          recomputeLlmPickedPdbIds();
        }
        removeUnmatchedEntry(key);
        didUpdate = true;
        continue;
      }
      setLlmUnmatchedAction(key, {
        status: 'failed',
        message: [
          body?.phase ? `[${String(body.phase).trim()}]` : '',
          body?.failure_reason || body?.message || 'Discovery failed.',
        ].filter(Boolean).join(' '),
        job_id: action.job_id,
      });
      await refreshDiscoveryJobLog(action.job_id);
      didUpdate = true;
    } catch (err) {
      setLlmUnmatchedAction(key, {
        status: 'failed',
        message: err?.message || 'Discovery polling failed.',
        job_id: action.job_id,
      });
      await refreshDiscoveryJobLog(action.job_id);
      didUpdate = true;
    }
  }

  if (didUpdate) {
    ensureLlmUnmatchedActions();
    renderLlmMatchPanels();
    renderConversationSelector();
    persistLlmState();
    renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
    loadBoltzConfigs({ silent: true });
    updateRunCommandFilterNote();
  }
  if (!hasActiveLlmUnmatchedJobs()) {
    stopLlmUnmatchedPolling();
  }
}

function ensureLlmUnmatchedPolling() {
  if (!hasActiveLlmUnmatchedJobs()) {
    stopLlmUnmatchedPolling();
    return;
  }
  if (state.llmUnmatchedPoller) return;
  pollLlmUnmatchedJobs();
  state.llmUnmatchedPoller = setInterval(pollLlmUnmatchedJobs, 2500);
}

function buildLlmMessageNode(role, text, { pending = false } = {}) {
  const safeRole = role === 'assistant' ? 'assistant' : 'user';
  const box = document.createElement('div');
  box.className = `llm-msg llm-msg-${safeRole}${pending ? ' llm-msg-pending' : ''}`;

  const head = document.createElement('div');
  head.className = 'llm-msg-head';
  const badge = document.createElement('span');
  badge.className = `llm-role-pill ${safeRole}`;
  badge.textContent = safeRole === 'assistant' ? 'Assistant' : 'You';
  head.appendChild(badge);
  box.appendChild(head);

  const body = document.createElement('div');
  body.className = 'llm-msg-body';
  body.textContent = String(text || '').trim();
  box.appendChild(body);
  return box;
}

function llmActionStatusLabel(status) {
  const normalized = String(status || 'idle').toLowerCase();
  if (normalized === 'queued') return 'Queued';
  if (normalized === 'running') return 'Running';
  if (normalized === 'failed') return 'Failed';
  if (normalized === 'success') return 'Added';
  return 'Ready';
}

function hasLlmAssistantResponse() {
  const history = Array.isArray(state.llmHistory) ? state.llmHistory : [];
  const hasAssistantInHistory = history.some((msg) => (
    String(msg?.role || '').trim().toLowerCase() === 'assistant'
    && String(msg?.content || '').trim()
  ));
  if (hasAssistantInHistory) return true;
  const hasMatched = Array.isArray(state.llmMatchedRows) && state.llmMatchedRows.length > 0;
  if (hasMatched) return true;
  return Array.isArray(state.llmUnmatched) && state.llmUnmatched.length > 0;
}

function renderLlmTranscript() {
  if (!el.llmChatTranscript) return;
  const history = Array.isArray(state.llmHistory) ? state.llmHistory : [];
  el.llmChatTranscript.innerHTML = '';
  if (!history.length) {
    const empty = document.createElement('div');
    empty.className = 'llm-empty';
    empty.textContent = 'No conversation yet.';
    el.llmChatTranscript.appendChild(empty);
    return;
  }
  history.forEach((msg) => {
    const role = String(msg?.role || '').trim().toLowerCase() === 'assistant' ? 'assistant' : 'user';
    const text = String(msg?.content || '').trim();
    if (!text) return;
    el.llmChatTranscript.appendChild(buildLlmMessageNode(role, text));
  });
  if (state.llmSuggestPending) {
    el.llmChatTranscript.appendChild(
      buildLlmMessageNode('assistant', 'Preparing suggestions...', { pending: true }),
    );
  }
  el.llmChatTranscript.scrollTop = el.llmChatTranscript.scrollHeight;
}

function renderLlmMatchPanels() {
  const showPanels = hasLlmAssistantResponse();
  const matchedPanel = el.llmMatchedPanel || el.llmMatchedList?.closest('.llm-list-panel');
  const unmatchedPanel = el.llmUnmatchedPanel || el.llmUnmatchedList?.closest('.llm-list-panel');
  if (matchedPanel) matchedPanel.hidden = !showPanels;
  if (unmatchedPanel) unmatchedPanel.hidden = !showPanels;

  const matchedList = Array.isArray(state.llmMatchedRows) ? state.llmMatchedRows : [];
  const deletedSet = new Set(
    (state.llmDeletedPdbIds || []).map((id) => normalizePdbId(id)).filter(Boolean),
  );
  ensureLlmUnmatchedActions();
  const activeCount = Array.isArray(state.llmPickedPdbIds) ? state.llmPickedPdbIds.length : 0;
  const totalCount = matchedList.length;
  if (el.llmMatchedCount) {
    el.llmMatchedCount.textContent = `${activeCount}/${totalCount}`;
  }
  if (el.llmUnmatchedCount) {
    const count = Array.isArray(state.llmUnmatched) ? state.llmUnmatched.length : 0;
    el.llmUnmatchedCount.textContent = String(count);
  }

  if (el.llmMatchedList) {
    el.llmMatchedList.innerHTML = '';
    if (!matchedList.length) {
      const empty = document.createElement('div');
      empty.className = 'llm-empty';
      empty.textContent = 'No matched targets yet.';
      el.llmMatchedList.appendChild(empty);
    } else {
      matchedList.forEach((row) => {
        const item = document.createElement('div');
        item.className = 'item';
        const pdb = normalizePdbId(row?.resolved_pdb_id || row?.pdb_id) || '—';
        const name = row?.protein_name || row?.preset_name || 'Unnamed target';
        const acc = row?.accession ? `Accession: ${row.accession}` : '';
        const preset = row?.preset_name ? `Preset: ${row.preset_name}` : '';
        const isDeleted = deletedSet.has(pdb);

        const rowWrap = document.createElement('div');
        rowWrap.className = 'llm-item-row';

        const main = document.createElement('div');
        main.className = 'llm-item-main';
        const label = document.createElement('span');
        label.className = `llm-item-title label${isDeleted ? ' deleted' : ''}`;
        label.textContent = `${pdb} · ${name}`;
        main.appendChild(label);
        const metaParts = [acc, preset].filter(Boolean);
        if (metaParts.length) {
          const meta = document.createElement('span');
          meta.className = 'llm-item-meta';
          meta.textContent = metaParts.join(' · ');
          main.appendChild(meta);
        }
        rowWrap.appendChild(main);

        const toggleBtn = document.createElement('button');
        toggleBtn.type = 'button';
        toggleBtn.className = 'mini-btn';
        toggleBtn.dataset.action = 'toggle-picked-target';
        toggleBtn.dataset.pdbId = pdb;
        toggleBtn.textContent = isDeleted ? 'Undelete' : 'Delete';
        const actions = document.createElement('div');
        actions.className = 'llm-item-actions';
        actions.appendChild(toggleBtn);
        rowWrap.appendChild(actions);

        item.appendChild(rowWrap);
        el.llmMatchedList.appendChild(item);
      });
    }
  }

  if (el.llmUnmatchedList) {
    const list = Array.isArray(state.llmUnmatched) ? state.llmUnmatched : [];
    el.llmUnmatchedList.innerHTML = '';
    if (!list.length) {
      const empty = document.createElement('div');
      empty.className = 'llm-empty';
      empty.textContent = 'No unmatched suggestions.';
      el.llmUnmatchedList.appendChild(empty);
    } else {
      list.forEach((entry) => {
        const item = document.createElement('div');
        item.className = 'item';
        const key = String(entry?.unmatched_key || '').trim() || buildUnmatchedKey(entry);
        const action = getLlmUnmatchedAction(key);
        const candidate = entry?.candidate || {};
        const label = candidate.target_name
          || candidate.protein_name
          || candidate.gene
          || candidate.uniprot
          || candidate.pdb_id
          || 'Unknown target';
        const nearest = Array.isArray(entry?.nearest) ? entry.nearest : [];
        const hints = nearest
          .slice(0, 2)
          .map((hit) => {
            const row = hit?.row || {};
            return row?.resolved_pdb_id || row?.pdb_id || row?.preset_name || '';
          })
          .filter(Boolean)
          .join(', ');
        const rowWrap = document.createElement('div');
        rowWrap.className = 'llm-item-row';
        const actionStatus = String(action?.status || 'idle').toLowerCase();
        const statusClass = ['queued', 'running', 'failed', 'success'].includes(actionStatus) ? actionStatus : 'idle';

        const main = document.createElement('div');
        main.className = 'llm-item-main';
        const textWrap = document.createElement('span');
        textWrap.className = 'llm-item-title label';
        textWrap.textContent = label;
        main.appendChild(textWrap);
        const meta = document.createElement('span');
        meta.className = 'llm-item-meta';
        const reason = entry?.reason ? `Reason: ${String(entry.reason).trim()}` : '';
        const matchHint = hints ? `Nearest: ${hints}` : 'No catalog match';
        meta.textContent = [matchHint, reason].filter(Boolean).join(' · ');
        main.appendChild(meta);
        rowWrap.appendChild(main);

        const actions = document.createElement('div');
        actions.className = 'llm-item-actions';
        const statusBadge = document.createElement('span');
        statusBadge.className = `llm-status-badge ${statusClass}`;
        statusBadge.textContent = llmActionStatusLabel(statusClass);
        actions.appendChild(statusBadge);

        const actionBtn = document.createElement('button');
        actionBtn.type = 'button';
        actionBtn.className = 'mini-btn';
        actionBtn.dataset.action = 'discover-unmatched-target';
        actionBtn.dataset.unmatchedKey = key;
        if (actionStatus === 'running' || actionStatus === 'queued') {
          actionBtn.textContent = 'Running...';
          actionBtn.disabled = true;
        } else if (actionStatus === 'failed') {
          actionBtn.textContent = 'Retry';
          actionBtn.disabled = false;
        } else if (actionStatus === 'success') {
          actionBtn.textContent = 'Added';
          actionBtn.disabled = true;
        } else {
          actionBtn.textContent = 'Discover';
          actionBtn.disabled = false;
        }
        actions.appendChild(actionBtn);
        rowWrap.appendChild(actions);
        item.appendChild(rowWrap);

        if (action?.message) {
          const meta = document.createElement('div');
          meta.className = 'llm-item-status';
          meta.textContent = action.message;
          item.appendChild(meta);
        }
        el.llmUnmatchedList.appendChild(item);
      });
    }
  }
}

function toggleMatchedTargetDeleted(pdbId) {
  const key = normalizePdbId(pdbId);
  if (!key) return;
  const next = new Set(
    (state.llmDeletedPdbIds || []).map((id) => normalizePdbId(id)).filter(Boolean),
  );
  if (next.has(key)) {
    next.delete(key);
  } else {
    next.add(key);
  }
  state.llmDeletedPdbIds = Array.from(next);
  recomputeLlmPickedPdbIds();
  renderLlmMatchPanels();
  renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
  renderBoltzConfigs();
  updateRunCommandFilterNote();
  persistLlmState();
}

function handleLlmMatchedListClick(event) {
  const btn = event.target.closest('button[data-action="toggle-picked-target"]');
  if (!btn) return;
  toggleMatchedTargetDeleted(btn.dataset.pdbId || '');
}

function findUnmatchedEntry(unmatchedKey) {
  const key = String(unmatchedKey || '').trim();
  if (!key) return null;
  const list = Array.isArray(state.llmUnmatched) ? state.llmUnmatched : [];
  return list.find((entry) => String(entry?.unmatched_key || '').trim() === key) || null;
}

async function submitUnmatchedDiscovery(unmatchedKey) {
  const key = String(unmatchedKey || '').trim();
  const entry = findUnmatchedEntry(key);
  if (!entry) return;
  const catalogName = String(getConfiguredLlmCatalogName() || '').trim();
  if (!catalogName) {
    showAlert('Set Config -> Input TSV defaults -> Default input file path to a catalog .tsv/.csv file in targets_catalog.');
    return;
  }
  setLlmUnmatchedAction(key, {
    status: 'queued',
    message: 'Queued discovery...',
  });
  renderLlmMatchPanels();
  persistLlmState();
  try {
    const historyForDiscover = Array.isArray(state.llmHistory)
      ? state.llmHistory.slice(-6)
      : [];
    const res = await fetch('/api/bulk/llm-targets/unmatched/discover', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        catalog_name: catalogName,
        unmatched_key: key,
        candidate: entry.candidate || {},
        history: historyForDiscover,
        vendor_scope: 'both',
        planning_mode: 'balanced',
        max_targets: 3,
        launch_browser: true,
      }),
    });
    const body = await res.json().catch(() => ({}));
    if (!res.ok) {
      throw new Error(body?.detail || body?.message || `Discover request failed (${res.status})`);
    }
    setLlmUnmatchedAction(key, {
      status: 'queued',
      job_id: body?.job_id ? String(body.job_id) : null,
      message: body?.message || 'Discovery queued.',
    });
    if (body?.job_id) {
      appendLog(`[discover] queued unmatched_key=${key} job_id=${body.job_id}`);
      await refreshDiscoveryJobLog(body.job_id);
    }
    ensureLlmUnmatchedPolling();
    renderLlmMatchPanels();
    renderConversationSelector();
    persistLlmState();
  } catch (err) {
    setLlmUnmatchedAction(key, {
      status: 'failed',
      message: err?.message || 'Failed to queue discovery.',
    });
    renderLlmMatchPanels();
    persistLlmState();
    showAlert(err?.message || 'Failed to queue discovery.');
  }
}

function handleLlmUnmatchedListClick(event) {
  const btn = event.target.closest('button[data-action="discover-unmatched-target"]');
  if (!btn) return;
  submitUnmatchedDiscovery(btn.dataset.unmatchedKey || '');
}

function setLlmViewMode(mode, { persist = true, rerender = true } = {}) {
  const next = mode === 'picked' ? 'picked' : 'all';
  state.llmViewMode = next;
  if (el.llmViewMode) {
    el.llmViewMode.value = next;
  }
  if (persist) persistLlmState();
  if (rerender) {
    renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
    renderBoltzConfigs();
    updateRunCommandFilterNote();
  }
}

function setLlmCatalogFilter(mode, { persist = true, rerender = true } = {}) {
  const next = normalizeLlmCatalogFilter(mode);
  state.llmCatalogFilter = next;
  if (el.llmCatalogFilter) {
    el.llmCatalogFilter.value = next;
  }
  if (persist) persistLlmState();
  if (rerender) {
    renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
    renderBoltzConfigs();
    updateRunCommandFilterNote();
  }
}

function applyLlmSuggestResponse(body = {}, { promptText = '', appendUserMessage = true } = {}) {
  const assistantMessage = String(body?.assistant_message || body?.message || '').trim();
  const matchedRows = Array.isArray(body?.matched_rows) ? body.matched_rows : [];
  const unmatched = normalizeUnmatchedList(body?.unmatched);
  const userText = String(promptText || '').trim();
  if (appendUserMessage && userText) {
    state.llmHistory.push({ role: 'user', content: userText });
  }
  if (assistantMessage) {
    state.llmHistory.push({ role: 'assistant', content: assistantMessage });
  }
  if (state.llmHistory.length > LLM_MAX_HISTORY_STORED) {
    state.llmHistory = state.llmHistory.slice(-LLM_MAX_HISTORY_STORED);
  }

  mergeMatchedRowsInState(matchedRows);
  state.llmUnmatched = unmatched;
  ensureLlmUnmatchedActions();
  recomputeLlmPickedPdbIds();

  if (!Array.isArray(state.bulkPreviewRows) || !state.bulkPreviewRows.length) {
    state.bulkPreviewRows = Array.isArray(state.llmMatchedRows) ? [...state.llmMatchedRows] : [];
    state.bulkPreviewMessage = body?.message || '';
  } else {
    mergeMatchedRowsIntoPreview(matchedRows);
  }
  if (body?.catalog_name) {
    state.activeCatalogName = String(body.catalog_name);
  }
  setLlmViewMode('picked', { persist: false, rerender: false });
  setLlmSuggestPending(false, { persist: false });
  renderLlmTranscript();
  renderLlmMatchPanels();
  renderConversationSelector();
  persistLlmState();
  ensureLlmUnmatchedPolling();
  renderBulkPreview(state.bulkPreviewRows || [], state.bulkPreviewMessage || '');
  loadBoltzConfigs({ silent: true });
  updateRunCommandFilterNote();
}

async function submitLlmTargetSuggest(triggerEl = null) {
  const prompt = String(el.llmChatPrompt?.value || '').trim();
  const catalogName = String(getConfiguredLlmCatalogName() || '').trim();
  if (!catalogName) {
    showAlert('Set Config -> Input TSV defaults -> Default input file path to a catalog .tsv/.csv file in targets_catalog.');
    return;
  }
  if (!prompt) {
    showAlert('Enter a prompt before suggesting targets.');
    return;
  }
  if (state.llmSuggestPending) return;
  const historyForApi = Array.isArray(state.llmHistory)
    ? state.llmHistory.slice(-LLM_HISTORY_WINDOW_FOR_API)
    : [];
  state.llmHistory.push({ role: 'user', content: prompt });
  if (state.llmHistory.length > LLM_MAX_HISTORY_STORED) {
    state.llmHistory = state.llmHistory.slice(-LLM_MAX_HISTORY_STORED);
  }
  if (el.llmChatPrompt) el.llmChatPrompt.value = '';
  setLlmSuggestPending(true, { persist: false });
  renderLlmTranscript();
  renderConversationSelector();
  persistLlmState();
  if (triggerEl) triggerEl.disabled = true;
  try {
    const payload = {
      catalog_name: catalogName,
      prompt,
      history: historyForApi,
      max_candidates: 24,
    };
    const res = await fetch('/api/bulk/llm-targets/suggest', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    const body = await res.json().catch(() => ({}));
    if (!res.ok) {
      throw new Error(body?.detail || body?.message || `LLM suggest failed (${res.status})`);
    }
    applyLlmSuggestResponse(body, { promptText: prompt, appendUserMessage: false });
    showAlert(body?.message || 'LLM target suggestions ready.', false);
  } catch (err) {
    setLlmSuggestPending(false, { persist: false });
    renderLlmTranscript();
    persistLlmState();
    showAlert(err.message || 'Failed to suggest targets.');
  } finally {
    if (triggerEl) triggerEl.disabled = false;
  }
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
    state.bulkPreviewMessage = body.message || '';
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

async function refreshDiscoveryJobLog(jobId) {
  const id = String(jobId || '').trim();
  if (!id) return;
  // Avoid clobbering the primary bulk run log stream while bulk polling is active.
  if (state.jobPoller && state.currentJobId && state.currentJobId !== id) return;
  try {
    const job = await fetchJob(id);
    const logs = Array.isArray(job?.logs) ? job.logs : [];
    if (logs.length) {
      resetJobLog(logs.join('\n'));
    } else {
      resetJobLog(`[discover] Job ${id} has started. Waiting for log lines...`);
    }
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.bulkStatus, 'Discovery running…');
    } else if (job.status === 'success') {
      setBadge(el.bulkStatus, 'Discovery complete', 'rgba(134, 239, 172, 0.25)');
    } else if (job.status === 'failed' || job.status === 'canceled') {
      setBadge(el.bulkStatus, 'Discovery failed', 'rgba(248, 113, 113, 0.25)');
    }
  } catch (err) {
    appendLog(`[discover] Unable to retrieve log stream for ${id}: ${err?.message || err}`);
  }
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
      addLink('Selected targets', job.details.targets_table_filename);
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
    printable.classList.add('bulk-printable-grid');
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
          <p class="help-text bulk-print-note">Includes all hotspot images and annotations currently shown.</p>
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
  restoreLlmState();
  ensureActiveConversation();
  syncStateFromConversation(getActiveConversation());
  recomputeLlmPickedPdbIds();
  setLlmSuggestPending(false, { persist: false });
  renderConversationSelector();
  renderLlmSuggestLoading();
  renderLlmTranscript();
  renderLlmMatchPanels();
  setLlmViewMode(state.llmViewMode || 'all', { persist: false, rerender: false });
  setLlmCatalogFilter(state.llmCatalogFilter || 'biotin', { persist: false, rerender: false });
  ensureLlmUnmatchedPolling();
  loadCommandDefaults();
  loadBulkUiSettings({ autoLoadInput: true });
  setupExampleTabs();
  if (el.llmConversationSelect) {
    el.llmConversationSelect.addEventListener('change', () => {
      const id = String(el.llmConversationSelect.value || '').trim();
      if (!id) return;
      setActiveConversation(id, { persist: true, rerender: true });
    });
  }
  if (el.llmNewConversation) {
    el.llmNewConversation.addEventListener('click', () => {
      createNewConversation();
      if (el.llmChatPrompt) el.llmChatPrompt.value = '';
    });
  }
  if (el.llmSuggestTargets) {
    el.llmSuggestTargets.addEventListener('click', () => submitLlmTargetSuggest(el.llmSuggestTargets));
  }
  if (el.llmViewMode) {
    el.llmViewMode.addEventListener('change', () => setLlmViewMode(el.llmViewMode.value));
  }
  if (el.llmCatalogFilter) {
    el.llmCatalogFilter.addEventListener('change', () => setLlmCatalogFilter(el.llmCatalogFilter.value));
  }
  if (el.llmMatchedList) {
    el.llmMatchedList.addEventListener('click', handleLlmMatchedListClick);
  }
  if (el.llmUnmatchedList) {
    el.llmUnmatchedList.addEventListener('click', handleLlmUnmatchedListClick);
  }
  if (el.llmChatPrompt) {
    el.llmChatPrompt.addEventListener('keydown', (evt) => {
      if (evt.key === 'Enter' && (evt.metaKey || evt.ctrlKey)) {
        evt.preventDefault();
        submitLlmTargetSuggest(el.llmSuggestTargets);
      }
    });
  }
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
  if (el.bulkReadmeBtn) el.bulkReadmeBtn.addEventListener('click', openBulkReadmeModal);
  if (el.bulkAlgorithmBtn) el.bulkAlgorithmBtn.addEventListener('click', () => toggleModal(el.bulkAlgorithmModal, true));
  if (el.bulkSettingsBtn) {
    el.bulkSettingsBtn.addEventListener('click', async () => {
      await loadBulkUiSettings({ autoLoadInput: false });
      toggleModal(el.bulkSettingsModal, true);
    });
  }
  if (el.bulkSettingsLoadDefault) {
    el.bulkSettingsLoadDefault.addEventListener('click', async () => {
      await loadDefaultInputFile({ force: true, silent: false });
    });
  }
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
  if (el.bulkReadmeClose) el.bulkReadmeClose.addEventListener('click', () => toggleModal(el.bulkReadmeModal, false));
  if (el.bulkSettingsClose) {
    el.bulkSettingsClose.addEventListener('click', () => {
      applyBulkUiSettingsToForm(state.uiSettings || {});
      toggleModal(el.bulkSettingsModal, false);
    });
  }
  if (el.bulkSettingsCancel) {
    el.bulkSettingsCancel.addEventListener('click', () => {
      applyBulkUiSettingsToForm(state.uiSettings || {});
      toggleModal(el.bulkSettingsModal, false);
    });
  }
  if (el.bulkSettingsSave) {
    el.bulkSettingsSave.addEventListener('click', () => saveBulkUiSettings(el.bulkSettingsSave));
  }
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
    if (closeTarget === 'bulk-readme') toggleModal(el.bulkReadmeModal, false);
    if (closeTarget === 'bulk-settings') {
      applyBulkUiSettingsToForm(state.uiSettings || {});
      toggleModal(el.bulkSettingsModal, false);
    }
    if (closeTarget === 'binder-export') toggleModal(el.binderExportModal, false);
  });

  renderBoltzConfigs();
  setRunEngine(getRunEngine());
  updateRunCommandFilterNote();
  refreshDiversity({ silent: true });
  hydrateSnapshotsFromCache({ silent: true });
}

document.addEventListener('DOMContentLoaded', init);
