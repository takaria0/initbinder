const DEFAULT_TABLE_SORT = { key: 'iptm', dir: 'desc' };
const DEFAULT_RFDIFF_CROP_RADIUS = 14;
const DEFAULT_BOLTZGEN_CROP_RADIUS = 14;
const SCATTER_COLOR_PALETTE = [
  '#2563eb',
  '#f97316',
  '#22c55e',
  '#a855f7',
  '#ef4444',
  '#0ea5e9',
  '#14b8a6',
  '#facc15',
];

const SCATTER_FIELDS = {
  rmsd_diego: {
    label: 'Binder RMSD (Å)',
    shortLabel: 'RMSD',
    defaultMin: 0,
    defaultMax: 5,
    tickDecimals: 1,
    valueDecimals: 3,
    threshold: { direction: 'max', default: 3.7, min: 0, max: null, step: 0.1, decimals: 2 },
  },
  iptm: {
    label: 'ipTM',
    shortLabel: 'ipTM',
    defaultMin: 0,
    defaultMax: 1,
    tickDecimals: 2,
    valueDecimals: 3,
    threshold: { direction: 'min', default: 0.8, min: 0, max: 1, step: 0.01, decimals: 2 },
  },
  ipsae_min: {
    label: 'ipSAE_min',
    shortLabel: 'ipSAE_min',
    defaultMin: 0,
    defaultMax: 1,
    tickDecimals: 2,
    valueDecimals: 3,
    threshold: { direction: 'min', default: null, min: 0, max: 1, step: 0.01, decimals: 2 },
  },
};

const DEFAULT_SCATTER_THRESHOLDS = Object.fromEntries(
  Object.entries(SCATTER_FIELDS).map(([key, config]) => {
    const defaultValue = config.threshold && Number.isFinite(config.threshold.default)
      ? config.threshold.default
      : null;
    return [key, Number.isFinite(defaultValue) ? defaultValue : null];
  }),
);

const RESULTS_TABLE_DISPLAY_LIMIT = 2000;

const ENGINE_FIELD_BLUEPRINT = [
  {
    field_id: 'total_designs',
    label: 'Total designs',
    description: 'Split evenly across detected epitope arms.',
    visible: true,
    debug_only: false,
  },
  {
    field_id: 'num_sequences',
    label: 'Sequences per backbone',
    description: 'ProteinMPNN sequences to sample per RFdiffusion backbone.',
    visible: true,
    debug_only: false,
  },
  {
    field_id: 'temperature',
    label: 'RFdiffusion temperature',
    description: 'Controls sampling diversity for RFdiffusion (0.0–1.0).',
    visible: true,
    debug_only: false,
  },
  {
    field_id: 'binder_chain_id',
    label: 'Binder chain ID',
    description: 'Override binder chain for downstream assessment (debug only).',
    visible: true,
    debug_only: true,
  },
  {
    field_id: 'af3_seed',
    label: 'AlphaFold 3 seed',
    description: 'Seed forwarded to AlphaFold 3 stage for reproducibility.',
    visible: true,
    debug_only: false,
  },
  {
    field_id: 'rfdiff_crop_radius',
    label: 'RFdiffusion crop radius',
    description: 'Enable to crop the prepared target around hotspots for RFdiffusion.',
    visible: true,
    debug_only: true,
  },
  {
    field_id: 'boltzgen_crop_radius',
    label: 'BoltzGen crop radius (Å)',
    description: 'Crop target around hotspots for BoltzGen configs/specs.',
    visible: false,
    debug_only: false,
  },
];

const DEFAULT_ENGINE_OPTIONS = [
  {
    engine_id: 'rfantibody',
    label: 'RFantibody (RFdiffusion → MPNN → AF3)',
    description: 'Legacy antibody pipeline using RFdiffusion, ProteinMPNN, and AlphaFold3.',
    is_default: true,
    fields: ENGINE_FIELD_BLUEPRINT.map((field) => ({ ...field })),
  },
  {
    engine_id: 'boltzgen',
    label: 'BoltzGen Diffusion',
    description: 'BoltzGen generative diffusion pipeline for protein binder design.',
    is_default: false,
    fields: ENGINE_FIELD_BLUEPRINT.map((field) => {
      if (field.field_id === 'total_designs') {
        return {
          ...field,
          label: 'Total designs',
          description: 'Total BoltzGen designs (single spec; no epitope splitting).',
          visible: true,
          debug_only: false,
        };
      }
      if (field.field_id === 'boltzgen_crop_radius') {
        return {
          ...field,
          label: 'BoltzGen crop radius (Å)',
          description: 'Crop target around hotspots for BoltzGen configs/specs.',
          visible: true,
          debug_only: false,
        };
      }
      return {
        ...field,
        description: '',
        visible: false,
      };
    }),
  },
];

const COUNTER_SELECTION_LIBRARY = {
  macs_streptavidin: {
    label: 'MACS streptavidin microbeads + biotinylated antigen',
    description:
      'Magnetic bead enrichment using streptavidin microbeads loaded with biotinylated antigen.',
    recommendations: [
      'Run streptavidin microbeads without antigen to deplete bead- and biotin-binding clones before the positive pull-down.',
      'Perform a mock pull-down with beads pre-blocked using excess free biotin to measure anti-biotin binders prior to enrichment.',
      'Include an OrthoRep host-only passage to monitor background binding after the MACS step.',
    ],
    notes: [
      'Load antigen at saturating levels but remove unbound biotinylated antigen so it does not leach off the beads during selection.',
      'Keep bead-to-cell ratios consistent between maturation cycles to maintain selection stringency.',
    ],
  },
  facs_anti_biotin: {
    label: 'FACS anti-biotin AF647 sorting with biotinylated antigen',
    description:
      'Two-color FACS sorting with anti-biotin AF647 secondary detecting nanobodies captured on biotinylated antigen.',
    recommendations: [
      'Sort a control sample stained with anti-biotin AF647 but without antigen to gate out clones that bind the detection reagent.',
      'Pre-clear the library with anti-biotin AF647-coated beads or antibody-only mixes before the positive FACS sort.',
      'Carry a secondary-only control (no antigen, no primary) every round to normalize OrthoRep background fluorescence.',
    ],
    notes: [
      'Titrate the AF647 conjugate to stay in the linear range and avoid fluorophore aggregation artifacts.',
      'Record the gating template (FSC/SSC, viability, negative control windows) and reuse it across maturation cycles.',
    ],
  },
  custom: {
    label: 'Custom / other selection workflow',
    description:
      'Capture bespoke counter-selection reminders for non-standard enrichment schemes.',
    recommendations: [
      'Design a negative selection that mimics your dominant off-target binder (tags, carrier proteins, scaffold, bead chemistry).',
      'Include a no-antigen passage on naïve OrthoRep host cells to track spontaneous display or autofluorescence enrichment.',
    ],
    notes: [
      'Document buffer additives (detergents, blockers, competitors) together with selection history so future rounds can replicate success.',
    ],
  },
};

const GLOBAL_COUNTER_NOTES = [
  'Confirm Golden Gate fragments remain free of BsaI/BbsI/BsmBI motifs before ordering DNA; recode any variants that reintroduce these sites.',
  'Log counter-selection timing and reagent lots after each maturation cycle so subsequent OrthoRep rounds stay reproducible.',
];

const COUNTER_TAG_RULES = [
  {
    id: 'trxa',
    chip: 'TrxA fusion detected',
    matches: (ctx) =>
      ctx.tagsLower.some((tag) => /trx|thioredoxin/.test(tag)) || /trx|thioredoxin/.test(ctx.text),
    recommendations: [
      'Run a counter-selection against thioredoxin (TrxA) carrier protein without the antigen insert to remove tag binders.',
    ],
    notes: [
      'Thioredoxin carriers are sticky; keep a TrxA-only control in each maturation round to spot false positives.',
    ],
    context: ['Vendor metadata indicates a thioredoxin (TrxA) fusion tag on the antigen.'],
  },
  {
    id: 'his',
    chip: 'His tag present',
    matches: (ctx) =>
      ctx.tagsLower.some((tag) => tag.includes('his')) || /his[-\s]?tag|6xhis|his6/.test(ctx.text),
    recommendations: [
      'Pre-clear the library on Ni-NTA or include soluble poly-His peptide/imidazole to suppress histidine-tag binders.',
    ],
    notes: [
      'Histidine tags quickly enrich sticky clones; track Ni-NTA depletion efficiency after each cycle.',
    ],
    context: ['Histidine affinity tag detected in antigen metadata.'],
  },
  {
    id: 'fc',
    chip: 'Fc fusion partner',
    matches: (ctx) =>
      ctx.tagsLower.some((tag) => /(^|[^a-z])fc($|[^a-z])/.test(tag)) || /(igg|fc\b)/.test(ctx.text),
    recommendations: [
      'Counter-select with human IgG Fc or Fc-blocking reagent to remove clones that prefer the Fc fusion.',
    ],
    notes: [
      'Fc fusions attract Protein A/G-like binders; include an Fc-only control during enrichment.',
    ],
    context: ['Antigen appears to include an Fc fusion partner.'],
  },
  {
    id: 'gst',
    chip: 'GST fusion partner',
    matches: (ctx) =>
      ctx.tagsLower.some((tag) => tag.includes('gst')) || /glutathione|gst/.test(ctx.text),
    recommendations: [
      'Include glutathione resin or GST-only protein counter selections to deplete GST binders.',
    ],
    notes: [
      'GST is highly sticky; verify clones do not track with glutathione resin over time.',
    ],
    context: ['Detected GST/glutathione fusion tag in antigen metadata.'],
  },
  {
    id: 'biotin',
    chip: 'Biotinylated tag',
    matches: (ctx) =>
      ctx.tagsLower.some((tag) => tag.includes('biotin') || tag.includes('avi')) || ctx.text.includes('biotin'),
    recommendations: [
      'Include a streptavidin-only or free biotin counter-selection to identify anti-biotin binders early.',
    ],
    notes: [
      'Track free biotin concentration in buffers; insufficient blocking increases bead-binding artifacts.',
    ],
    context: ['Biotin/avi-tag annotation detected; monitor for anti-biotin clones.'],
  },
  {
    id: 'membrane_scaffold',
    chip: 'Membrane scaffold',
    matches: (ctx) =>
      /nanodisc|liposome|amphipol|detergent|micelle|membrane/.test(ctx.text),
    recommendations: [
      'Run a counter-selection with the membrane scaffold or detergent micelle alone (no target protein) to eliminate scaffold binders.',
    ],
    notes: [
      'Carrier scaffolds can dominate affinity signals; document detergent composition and concentration for reproducibility.',
    ],
    context: ['Vendor notes mention a membrane scaffold or detergent carrier.'],
  },
  {
    id: 'albumin',
    chip: 'Carrier protein (BSA/albumin)',
    matches: (ctx) => ctx.tagsLower.some((tag) => tag.includes('bsa') || tag.includes('albumin')) || /bsa|albumin/.test(ctx.text),
    recommendations: [
      'Include BSA/albumin-only counter selections or add soluble carrier protein during positive selections to soak up carrier binders.',
    ],
    notes: [
      'Albumin carriers frequently drive false positives; watch for clones that persist even without antigen present.',
    ],
    context: ['Carrier protein (BSA/albumin) annotated in antigen metadata.'],
  },
];

const state = {
  currentPdb: null,
  designEngines: [],
  selectedDesignEngine: 'rfantibody',
  activeDesignFields: new Set(),
  jobPoller: null,
  jobContext: null,
  rankings: [],
  rankingsResponse: null,
  rankingsTotal: 0,
  rankingsShowAll: false,
  rankingsDisplayLimit: RESULTS_TABLE_DISPLAY_LIMIT,
  resultsSource: 'af3',
  tableSort: { ...DEFAULT_TABLE_SORT },
  scatterLayout: [],
  scatterPoints: [],
  scatterAxes: { x: 'rmsd_diego', y: 'iptm' },
  scatterThresholds: { ...DEFAULT_SCATTER_THRESHOLDS },
  scatterStats: null,
  scatterColorMap: new Map(),
  selectedDesign: null,
  jobs: [],
  selectedJobId: null,
  runHistory: {},
  boltzRunHistory: {},
  boltzActiveSpec: null,
  jobRefreshTimer: null,
  presets: [],
  activePresetId: null,
  debugMode: false,
  rankingsPollTimer: null,
  rankingsFetching: false,
  activeRunLabel: '',
  galleryAvailable: false,
  targetStatus: null,
  targetDetails: null,
  lastClusterSyncAt: 0,
  pendingClusterSync: null,
  pendingClusterSyncTimer: null,
  clusterRetryTimer: null,
  analysisResults: null,
  analysisRunning: false,
  analysisContext: null,
  librarySummary: null,
  libraryAlignmentMode: 'aa',
  libraryJobId: null,
  counterSelectionMethod: '',
  counterSelectionPlan: { recommendations: [], notes: [], context: [], chips: [] },
  catalogFiles: [],
  catalogSelectedFile: null,
  catalogHeaders: [],
  catalogRows: [],
  catalogSelectedRowIndex: null,
};

state.activeDesignFields = new Set(
  ENGINE_FIELD_BLUEPRINT.filter((field) => field.visible !== false).map((field) => field.field_id),
);

const SYNC_RESULTS_MIN_INTERVAL_MS = 600000;

function escapeHtml(value) {
  return String(value ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function getResultsSource() {
  if (el.resultsSource) {
    const value = (el.resultsSource.value || '').trim().toLowerCase();
    if (value === 'boltzgen' || value === 'af3') {
      state.resultsSource = value;
      return value;
    }
  }
  return state.resultsSource || 'af3';
}

function setResultsSource(source) {
  const normalized = source === 'boltzgen' ? 'boltzgen' : 'af3';
  state.resultsSource = normalized;
  if (el.resultsSource) {
    el.resultsSource.value = normalized;
  }
  const isBoltz = normalized === 'boltzgen';
  if (el.boltzSpecWrapper) {
    el.boltzSpecWrapper.hidden = !isBoltz;
  }
  if (el.syncResultsBtn) {
    el.syncResultsBtn.hidden = isBoltz;
  }
  if (el.syncBoltzResultsBtn) {
    el.syncBoltzResultsBtn.hidden = !isBoltz;
  }
  if (el.boltzRunList) {
    el.boltzRunList.hidden = !isBoltz;
  }
  if (isBoltz && state.currentPdb) {
    fetchBoltzRunHistory(state.currentPdb);
  }
}

const el = {
  targetForm: document.querySelector('#target-form'),
  clusterStatus: document.querySelector('#cluster-status'),
  debugToggle: document.querySelector('#debug-toggle'),
  pdbInput: document.querySelector('#pdb-id'),
  antigenInput: document.querySelector('#antigen-url'),
  targetName: document.querySelector('#target-name'),
  targetEpitopes: document.querySelector('#target-epitopes'),
  targetDecidePrompt: document.querySelector('#decide-scope-prompt'),
  presetIdInput: document.querySelector('#preset-id'),
  presetList: document.querySelector('#preset-list'),
  presetRefresh: document.querySelector('#presets-refresh'),
  forceCheckbox: document.querySelector('#force-refresh'),
  decideCheckbox: document.querySelector('#skip-decide'),
  prepCheckbox: document.querySelector('#skip-prep'),
  targetSubmit: document.querySelector('#target-submit'),
  targetBadge: document.querySelector('#target-status-badge'),
  pymolHotspots: document.querySelector('#pymol-hotspots'),
  jobLog: document.querySelector('#job-log'),
  jobAlert: document.querySelector('#job-alert'),
  refreshAlignmentBtn: document.querySelector('#refresh-alignment'),
  alignmentMeta: document.querySelector('#alignment-meta'),
  alignmentView: document.querySelector('#alignment-view'),
  designSubmit: document.querySelector('#design-submit'),
  designStatus: document.querySelector('#design-status'),
  designEngineSelect: document.querySelector('#design-engine'),
  designEngineHint: document.querySelector('#design-engine-hint'),
  designTotal: document.querySelector('#design-total'),
  designNumSeq: document.querySelector('#design-num-seq'),
  designTemp: document.querySelector('#design-temp'),
  designRunLabel: document.querySelector('#design-run-label'),
  designBinderChain: document.querySelector('#design-binder-chain'),
  designAf3Seed: document.querySelector('#design-af3-seed'),
  designEnableRfdiffCrop: document.querySelector('#design-enable-rfdiff-crop'),
  designBoltzBinding: document.querySelector('#design-boltz-binding'),
  designBoltzCropRadius: document.querySelector('#design-boltz-crop-radius'),
  designBoltzBindingRow: document.querySelector('#design-boltz-binding-row'),
  refreshResultsBtn: document.querySelector('#refresh-results'),
  resultsMeta: document.querySelector('#results-meta'),
  binderDetail: document.querySelector('#binder-detail'),
  resultsTableWrapper: document.querySelector('#results-table-wrapper'),
  resultsTable: document.querySelector('#results-table'),
  resultsTruncate: document.querySelector('#results-truncate'),
  resultsTruncateText: document.querySelector('#results-truncate-text'),
  resultsTruncateToggle: document.querySelector('#results-truncate-toggle'),
  scatterCanvas: document.querySelector('#scatter-plot'),
  scatterHistX: document.querySelector('#scatter-hist-x'),
  scatterHistY: document.querySelector('#scatter-hist-y'),
  scatterXAxis: document.querySelector('#scatter-axis-x'),
  scatterYAxis: document.querySelector('#scatter-axis-y'),
  scatterThresholdX: document.querySelector('#scatter-threshold-x'),
  scatterThresholdY: document.querySelector('#scatter-threshold-y'),
  scatterThresholdXLabel: document.querySelector('#scatter-threshold-x-label'),
  scatterThresholdYLabel: document.querySelector('#scatter-threshold-y-label'),
  scatterThresholdCount: document.querySelector('#scatter-threshold-count'),
  scatterLegend: document.querySelector('#scatter-legend'),
  scatterRefresh: document.querySelector('#scatter-refresh'),
  plotContainer: document.querySelector('#plot-container'),
  exportOpen: document.querySelector('#export-open'),
  exportModal: document.querySelector('#export-modal'),
  exportClose: document.querySelector('#export-close'),
  exportBackdrop: document.querySelector('#export-modal .modal-backdrop'),
  exportButton: document.querySelector('#export-run'),
  exportTopN: document.querySelector('#export-topn'),
  exportCodonHost: document.querySelector('#export-codon-host'),
  exportGc: document.querySelector('#export-gc'),
  exportGcWindow: document.querySelector('#export-gc-window'),
  exportPrefix: document.querySelector('#export-prefix'),
  exportSuffix: document.querySelector('#export-suffix'),
  exportDnac: document.querySelector('#export-dnachisel'),
  resultsRunLabel: document.querySelector('#results-run-label'),
  resultsSource: document.querySelector('#results-source'),
  boltzSpec: document.querySelector('#boltz-spec'),
  boltzSpecWrapper: document.querySelector('#boltz-spec-wrapper'),
  resultsLimit: document.querySelector('#results-limit'),
  pymolTop: document.querySelector('#pymol-top'),
  pymolMovie: document.querySelector('#pymol-movie'),
  syncResultsBtn: document.querySelector('#sync-results'),
  syncBoltzResultsBtn: document.querySelector('#sync-boltz-results'),
  assessSubmit: document.querySelector('#assess-submit'),
  jobList: document.querySelector('#job-list'),
  refreshJobsBtn: document.querySelector('#refresh-jobs'),
  runHistory: document.querySelector('#run-history'),
  boltzRunList: document.querySelector('#boltz-run-list'),
  runLabelOptions: document.querySelector('#run-label-options'),
  activeRunName: document.querySelector('#active-run-name'),
  selectionCard: document.querySelector('#current-selection-card'),
  selectionPdb: document.querySelector('#current-selection-pdb'),
  selectionName: document.querySelector('#current-selection-name'),
  selectionAntigen: document.querySelector('#current-selection-antigen'),
  epitopeList: document.querySelector('#epitope-list'),
  chainList: document.querySelector('#chain-list'),
  antigenDetails: document.querySelector('#antigen-details'),
  analysisPanel: document.querySelector('#analysis-panel'),
  analysisStatus: document.querySelector('#analysis-status'),
  analysisPlots: document.querySelector('#analysis-plots'),
  analysisMatrix: document.querySelector('#analysis-matrix'),
  analysisRun: document.querySelector('#analysis-run'),
  analysisLogs: document.querySelector('#analysis-logs'),
  libraryRun: document.querySelector('#library-run'),
  libraryTopN: document.querySelector('#library-topn'),
  libraryCodonHost: document.querySelector('#library-codon-host'),
  librarySequenceColumn: document.querySelector('#library-sequence-column'),
  libraryUpstream: document.querySelector('#library-upstream'),
  libraryDownstream: document.querySelector('#library-downstream'),
  libraryStatus: document.querySelector('#library-status'),
  librarySummary: document.querySelector('#library-summary'),
  libraryDownloads: document.querySelector('#library-downloads'),
  libraryAlignment: document.querySelector('#library-alignment'),
  libraryAlignmentToggle: document.querySelector('#library-alignment-toggle'),
  libraryAlignmentLegend: document.querySelector('#library-alignment-legend'),
  libraryAlignmentTable: document.querySelector('#library-alignment-table'),
  libraryAssembly: document.querySelector('#library-assembly'),
  libraryPreview: document.querySelector('#library-preview'),
  counterSelectionSection: document.querySelector('#counter-selection-section'),
  counterSelectionMethod: document.querySelector('#counter-selection-method'),
  counterSelectionList: document.querySelector('#counter-selection-list'),
  counterSelectionNotes: document.querySelector('#counter-selection-notes'),
  counterSelectionContext: document.querySelector('#counter-selection-context'),
  counterSelectionContextList: document.querySelector('#counter-selection-context-list'),
  counterSelectionChips: document.querySelector('#counter-selection-chips'),
  counterSelectionDescription: document.querySelector('#counter-selection-description'),
  catalogOpen: document.querySelector('#open-catalog-browser'),
  catalogModal: document.querySelector('#catalog-modal'),
  catalogClose: document.querySelector('#catalog-close'),
  catalogBackdrop: document.querySelector('#catalog-modal .modal-backdrop'),
  catalogList: document.querySelector('#catalog-browser-list'),
  catalogViewer: document.querySelector('#catalog-browser-viewer'),
  catalogSelectedName: document.querySelector('#catalog-browser-selected-name'),
  catalogSelectedInfo: document.querySelector('#catalog-browser-selected-info'),
  catalogStatus: document.querySelector('#catalog-browser-status'),
  catalogRefresh: document.querySelector('#catalog-browser-refresh'),
  catalogLimit: document.querySelector('#catalog-browser-limit'),
  catalogDownload: document.querySelector('#catalog-browser-download'),
};

const DESIGN_FIELD_ELEMENTS = (() => {
  const blueprintMap = Object.fromEntries(ENGINE_FIELD_BLUEPRINT.map((field) => [field.field_id, field]));
  const base = {
    total_designs: { input: el.designTotal },
    num_sequences: { input: el.designNumSeq },
    temperature: { input: el.designTemp },
    binder_chain_id: { input: el.designBinderChain },
    af3_seed: { input: el.designAf3Seed },
    rfdiff_crop_radius: { input: el.designEnableRfdiffCrop },
    boltzgen_crop_radius: { input: el.designBoltzCropRadius },
  };

  Object.entries(base).forEach(([fieldId, entry]) => {
    const wrapper = document.querySelector(`[data-engine-field="${fieldId}"]`);
    const container = document.querySelector(`[data-engine-field-container="${fieldId}"]`)
      || (wrapper ? wrapper.closest('[data-engine-field-container]') : null);
    const label = wrapper?.querySelector(`[data-engine-label="${fieldId}"]`) || wrapper?.querySelector('.field-label') || null;
    const help = wrapper?.querySelector(`[data-engine-help="${fieldId}"]`) || null;

    entry.wrapper = wrapper || null;
    entry.container = container || null;
    entry.label = label;
    entry.help = help;
    entry.defaultLabel = label ? label.textContent.trim() : (blueprintMap[fieldId]?.label || '');
    entry.defaultDescription = help ? help.textContent.trim() : (blueprintMap[fieldId]?.description || '');
    entry.defaultVisible = blueprintMap[fieldId]?.visible !== false;
  });

  return base;
})();

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

function getEpitopeLabel(row) {
  if (!row || typeof row !== 'object') return '';
  const metadata = row.metadata && typeof row.metadata === 'object' ? row.metadata : {};
  const value =
    metadata.epitope ||
    metadata.epitope_name ||
    metadata.epitopeName ||
    metadata.arm ||
    metadata.target_epitope ||
    '';
  if (typeof value === 'string') return value.trim();
  if (typeof value === 'number' || typeof value === 'boolean') return String(value);
  return '';
}

function ensureValidTableSort(rows) {
  if (!state.tableSort || !state.tableSort.key) {
    state.tableSort = { ...DEFAULT_TABLE_SORT };
    return;
  }
  if (!Array.isArray(rows) || rows.length === 0) {
    return;
  }
  const firstRow = rows[0];
  if (firstRow && !Object.prototype.hasOwnProperty.call(firstRow, state.tableSort.key)) {
    state.tableSort = { ...DEFAULT_TABLE_SORT };
  }
}

function getScatterFieldConfig(key) {
  if (!key) return null;
  return Object.prototype.hasOwnProperty.call(SCATTER_FIELDS, key) ? SCATTER_FIELDS[key] : null;
}

function getActiveScatterAxes() {
  const axes = state.scatterAxes || {};
  const x = getScatterFieldConfig(axes.x) ? axes.x : 'rmsd_diego';
  const y = getScatterFieldConfig(axes.y) ? axes.y : 'iptm';
  return { x, y };
}

function coerceScatterNumber(value) {
  if (typeof value === 'number') {
    return Number.isFinite(value) ? value : Number.NaN;
  }
  if (typeof value === 'string') {
    const trimmed = value.trim();
    if (!trimmed) return Number.NaN;
    const parsed = Number.parseFloat(trimmed);
    return Number.isFinite(parsed) ? parsed : Number.NaN;
  }
  return Number.NaN;
}

function buildScatterPointsFromRows(rows) {
  if (!Array.isArray(rows)) return [];
  return rows
    .map((row) => {
      if (!row || typeof row !== 'object') return null;
      const iptm = coerceScatterNumber(row.iptm);
      const rmsd = coerceScatterNumber(row.rmsd_diego);
      const ipsae = coerceScatterNumber(row.ipsae_min);
      if (!Number.isFinite(iptm) && !Number.isFinite(rmsd) && !Number.isFinite(ipsae)) {
        return null;
      }
      const hotspot = coerceScatterNumber(row.hotspot_min_distance);
      const metadata = row.metadata && typeof row.metadata === 'object' ? row.metadata : {};
      return {
        design_name: row.design_name || metadata.design_name || '',
        iptm: Number.isFinite(iptm) ? iptm : null,
        rmsd_diego: Number.isFinite(rmsd) ? rmsd : null,
        ipsae_min: Number.isFinite(ipsae) ? ipsae : null,
        hotspot_min_distance: Number.isFinite(hotspot) ? hotspot : null,
        metadata,
      };
    })
    .filter(Boolean);
}

function resolveScatterPoints() {
  const scatter = state.rankingsResponse?.scatter;
  if (Array.isArray(scatter) && scatter.length > 0) {
    return scatter;
  }
  return buildScatterPointsFromRows(state.rankings);
}

function inferBoltzSpecName(payload) {
  if (!payload) return null;
  const rows = Array.isArray(payload.rows) ? payload.rows : [];
  for (let i = 0; i < rows.length; i += 1) {
    const meta = rows[i]?.metadata;
    if (meta && typeof meta.spec === 'string' && meta.spec.trim()) {
      return meta.spec.trim();
    }
  }
  const sourcePath = payload.source_path || payload.sourcePath;
  if (typeof sourcePath === 'string' && sourcePath) {
    const parts = sourcePath.split(/[\\/]+/).filter(Boolean);
    const idx = parts.lastIndexOf('final_ranked_designs');
    if (idx > 0 && idx - 1 < parts.length) {
      return parts[idx - 1];
    }
    if (parts.length >= 3) {
      return parts[parts.length - 3];
    }
  }
  return null;
}

function getScatterFieldValue(row, key) {
  if (!row || !key) return Number.NaN;
  const value = row[key];
  if (typeof value === 'number') return value;
  if (typeof value === 'string' && value.trim() !== '') {
    const parsed = Number.parseFloat(value);
    return Number.isFinite(parsed) ? parsed : Number.NaN;
  }
  return Number.NaN;
}

function formatScatterValue(fieldKey, value, fallbackDecimals = 2) {
  if (!Number.isFinite(value)) return 'NA';
  const config = getScatterFieldConfig(fieldKey);
  if (config && Number.isInteger(config.valueDecimals)) {
    return value.toFixed(config.valueDecimals);
  }
  return value.toFixed(fallbackDecimals);
}

function formatScatterThreshold(fieldKey, value) {
  if (!Number.isFinite(value)) return 'NA';
  const config = getScatterFieldConfig(fieldKey);
  if (config && config.threshold && Number.isInteger(config.threshold.decimals)) {
    return value.toFixed(config.threshold.decimals);
  }
  if (config && Number.isInteger(config.tickDecimals)) {
    return value.toFixed(config.tickDecimals);
  }
  return value.toFixed(2);
}

function clampScatterThreshold(fieldKey, value) {
  const config = getScatterFieldConfig(fieldKey);
  if (!config || !config.threshold) return value;
  let result = value;
  const { min, max, decimals } = config.threshold;
  if (Number.isFinite(min)) result = Math.max(min, result);
  if (Number.isFinite(max)) result = Math.min(max, result);
  if (Number.isInteger(decimals)) {
    result = Number.parseFloat(result.toFixed(decimals));
  }
  return result;
}

function updateScatterControls() {
  const axes = getActiveScatterAxes();
  state.scatterAxes = axes;
  if (el.scatterXAxis && el.scatterXAxis.value !== axes.x) {
    el.scatterXAxis.value = axes.x;
  }
  if (el.scatterYAxis && el.scatterYAxis.value !== axes.y) {
    el.scatterYAxis.value = axes.y;
  }
  updateScatterThresholdInputs();
}

function updateScatterThresholdInputs() {
  updateSingleScatterThresholdInput('x');
  updateSingleScatterThresholdInput('y');
}

function updateSingleScatterThresholdInput(axis) {
  const axes = state.scatterAxes || {};
  const fieldKey = axis === 'y' ? axes.y : axes.x;
  const config = getScatterFieldConfig(fieldKey);
  const input = axis === 'y' ? el.scatterThresholdY : el.scatterThresholdX;
  const label = axis === 'y' ? el.scatterThresholdYLabel : el.scatterThresholdXLabel;
  if (!input || !label) return;

  if (!config) {
    label.textContent = axis === 'y' ? 'Y threshold' : 'X threshold';
    input.value = '';
    input.disabled = true;
    input.placeholder = 'N/A';
    input.removeAttribute('min');
    input.removeAttribute('max');
    input.removeAttribute('step');
    return;
  }

  const thresholdConfig = config.threshold;
  if (thresholdConfig) {
    const directionSymbol = thresholdConfig.direction === 'max' ? '≤'
      : thresholdConfig.direction === 'min' ? '≥'
        : '';
    label.textContent = `${config.label} threshold${directionSymbol ? ` (${directionSymbol})` : ''}`;
    input.disabled = false;
    if (Number.isFinite(thresholdConfig.min)) input.min = thresholdConfig.min;
    else input.removeAttribute('min');
    if (Number.isFinite(thresholdConfig.max)) input.max = thresholdConfig.max;
    else input.removeAttribute('max');
    if (Number.isFinite(thresholdConfig.step)) input.step = thresholdConfig.step;
    else input.removeAttribute('step');
    const current = state.scatterThresholds?.[fieldKey];
    if (Number.isFinite(current)) {
      input.value = current;
    } else {
      input.value = '';
    }
    input.placeholder = '';
  } else {
    label.textContent = `${config.label} threshold`;
    input.disabled = true;
    input.value = '';
    input.placeholder = 'N/A';
    input.removeAttribute('min');
    input.removeAttribute('max');
    input.removeAttribute('step');
  }
}

function handleScatterAxisChange(axis, fieldKey) {
  if (!axis) return;
  if (!getScatterFieldConfig(fieldKey)) {
    updateScatterControls();
    return;
  }
  const nextAxes = { ...getActiveScatterAxes(), [axis]: fieldKey };
  state.scatterAxes = nextAxes;
  updateScatterControls();
  computeScatterLayout(state.scatterPoints || []);
  renderScatter();
}

function timestampString() {
  const now = new Date();
  const pad = (n) => String(n).padStart(2, '0');
  return `${now.getFullYear()}${pad(now.getMonth() + 1)}${pad(now.getDate())}_${pad(now.getHours())}${pad(now.getMinutes())}${pad(now.getSeconds())}`;
}

function setDebugMode(enabled) {
  state.debugMode = Boolean(enabled);
  document.body.classList.toggle('debug-mode', state.debugMode);
}

function updateActiveRunDisplay() {
  if (!el.activeRunName) return;
  el.activeRunName.textContent = state.activeRunLabel || 'latest';
}

function dedupeList(items = []) {
  const seen = new Set();
  const output = [];
  items.forEach((item) => {
    const text = typeof item === 'string' ? item.trim() : '';
    if (!text) return;
    const key = text.toLowerCase();
    if (seen.has(key)) return;
    seen.add(key);
    output.push(text);
  });
  return output;
}

function extractAntigenCounterInsights(target = state.targetDetails) {
  const plan = { recommendations: [], notes: [], context: [], chips: [] };
  if (!target || !target.antigen) {
    return plan;
  }
  const antigen = target.antigen || {};
  const rawTags = Array.isArray(antigen.tags)
    ? antigen.tags
    : antigen.tags
      ? [antigen.tags]
      : [];
  const tags = rawTags
    .map((tag) => String(tag).trim())
    .filter(Boolean);
  const tagsLower = tags.map((tag) => tag.toLowerCase());
  const textParts = [
    antigen.notes,
    antigen.product_form,
    antigen.catalog,
    antigen.expression_host,
  ]
    .map((part) => (part == null ? '' : String(part).toLowerCase()))
    .filter(Boolean);
  const text = textParts.join(' ');

  if (tags.length) {
    plan.context.push(`Vendor tags: ${tags.join(', ')}`);
    tags.forEach((tag) => {
      plan.chips.push(`Tag: ${tag}`);
    });
  }
  if (antigen.expression_host) {
    plan.context.push(`Expression host: ${antigen.expression_host}`);
  }
  if (antigen.product_form) {
    plan.context.push(`Product form: ${antigen.product_form}`);
  }

  const ctx = { tagsLower, text, antigen, target, tags };
  COUNTER_TAG_RULES.forEach((rule) => {
    try {
      if (!rule || typeof rule.matches !== 'function') return;
      if (!rule.matches(ctx)) return;
      if (rule.chip) plan.chips.push(rule.chip);
      if (Array.isArray(rule.recommendations)) {
        plan.recommendations.push(...rule.recommendations);
      }
      if (Array.isArray(rule.notes)) {
        plan.notes.push(...rule.notes);
      }
      if (Array.isArray(rule.context)) {
        plan.context.push(...rule.context);
      }
    } catch (err) {
      console.warn('Failed to evaluate counter-selection rule', rule?.id, err);
    }
  });

  plan.recommendations = dedupeList(plan.recommendations);
  plan.notes = dedupeList(plan.notes);
  plan.context = dedupeList(plan.context);
  plan.chips = dedupeList(plan.chips);
  return plan;
}

function computeCounterSelectionPlan(selectionKey, target = state.targetDetails) {
  const plan = { recommendations: [], notes: [], context: [], chips: [] };
  const config = selectionKey ? COUNTER_SELECTION_LIBRARY[selectionKey] : null;
  if (config) {
    if (Array.isArray(config.recommendations)) {
      plan.recommendations.push(...config.recommendations);
    }
    if (Array.isArray(config.notes)) {
      plan.notes.push(...config.notes);
    }
  }
  plan.notes.push(...GLOBAL_COUNTER_NOTES);

  const antigenInsights = extractAntigenCounterInsights(target);
  plan.recommendations.push(...antigenInsights.recommendations);
  plan.notes.push(...antigenInsights.notes);
  plan.context.push(...antigenInsights.context);
  plan.chips.push(...antigenInsights.chips);

  plan.recommendations = dedupeList(plan.recommendations);
  plan.notes = dedupeList(plan.notes);
  plan.context = dedupeList(plan.context);
  plan.chips = dedupeList(plan.chips);
  return plan;
}

function renderCounterList(listEl, items, emptyText) {
  if (!listEl) return;
  listEl.innerHTML = '';
  if (!Array.isArray(items) || items.length === 0) {
    const li = document.createElement('li');
    li.className = 'counter-empty';
    li.textContent = emptyText;
    listEl.appendChild(li);
    return;
  }
  items.forEach((item) => {
    const li = document.createElement('li');
    li.textContent = item;
    listEl.appendChild(li);
  });
}

function renderCounterSelectionPanel() {
  if (!el.counterSelectionSection) return;
  const method = state.counterSelectionMethod || '';
  const config = method ? COUNTER_SELECTION_LIBRARY[method] : null;
  const plan = state.counterSelectionPlan || { recommendations: [], notes: [], context: [], chips: [] };

  if (el.counterSelectionMethod && el.counterSelectionMethod.value !== method) {
    el.counterSelectionMethod.value = method;
  }

  if (el.counterSelectionDescription) {
    let description = 'Select a primary selection workflow to surface counter selections and pitfalls to keep OrthoRep affinity maturation on track.';
    if (config && config.description) {
      description = config.description;
    } else if (config && config.label) {
      description = config.label;
    }
    el.counterSelectionDescription.textContent = description;
  }

  const recEmptyText = method
    ? 'No counter selections recorded yet for this workflow.'
    : 'Choose a primary selection to see counter-selection suggestions.';
  renderCounterList(el.counterSelectionList, plan.recommendations, recEmptyText);

  const noteEmptyText = method
    ? 'No additional pitfalls detected yet.'
    : 'Select a workflow to surface technical reminders.';
  renderCounterList(el.counterSelectionNotes, plan.notes, noteEmptyText);

  if (el.counterSelectionChips) {
    el.counterSelectionChips.innerHTML = '';
    if (Array.isArray(plan.chips) && plan.chips.length) {
      plan.chips.forEach((chip) => {
        const span = document.createElement('span');
        span.className = 'counter-chip';
        span.textContent = chip;
        el.counterSelectionChips.appendChild(span);
      });
    }
  }

  if (el.counterSelectionContextList) {
    const contextItems = Array.isArray(plan.context) ? plan.context : [];
    el.counterSelectionContextList.innerHTML = '';
    if (contextItems.length === 0) {
      const li = document.createElement('li');
      li.className = 'counter-empty';
      li.textContent = 'No antigen-specific risk factors detected yet.';
      el.counterSelectionContextList.appendChild(li);
    } else {
      contextItems.forEach((item) => {
        const li = document.createElement('li');
        li.textContent = item;
        el.counterSelectionContextList.appendChild(li);
      });
    }
  }

  if (el.counterSelectionContext) {
    const hasContext = (plan.context && plan.context.length > 0) || (plan.chips && plan.chips.length > 0);
    el.counterSelectionContext.hidden = !hasContext;
  }
}

function updateCounterSelectionPlan() {
  state.counterSelectionPlan = computeCounterSelectionPlan(state.counterSelectionMethod || '', state.targetDetails);
  renderCounterSelectionPanel();
}

function renderTargetInsights(status = state.targetStatus) {
  if (!el.epitopeList || !el.chainList || !el.antigenDetails) return;

  const target = status || {};
  state.targetDetails = target;

  const epitopes = Array.isArray(target.epitopes) ? target.epitopes : [];
  el.epitopeList.innerHTML = '';
  if (epitopes.length === 0) {
    const li = document.createElement('li');
    li.className = 'insight-empty';
    li.textContent = 'No epitope metadata available yet.';
    el.epitopeList.appendChild(li);
  } else {
    epitopes.forEach((ep) => {
      const li = document.createElement('li');
      const title = document.createElement('span');
      title.className = 'insight-title';
      title.textContent = ep.name || 'Unnamed epitope';
      li.appendChild(title);
      const residues = Array.isArray(ep.residues) ? ep.residues.filter(Boolean) : [];
      const residueLine = document.createElement('div');
      residueLine.textContent = residues.length ? `Residues: ${residues.join(', ')}` : 'Residues: —';
      li.appendChild(residueLine);
      if (typeof ep.surface_exposed === 'boolean') {
        const exposedLine = document.createElement('div');
        exposedLine.textContent = ep.surface_exposed ? 'Surface exposed: yes' : 'Surface exposed: no';
        li.appendChild(exposedLine);
      }
      if (typeof ep.hotspot_count === 'number') {
        const hotspotLine = document.createElement('div');
        hotspotLine.textContent = `Hotspots: ${ep.hotspot_count}`;
        li.appendChild(hotspotLine);
      }
      if (ep.rationale) {
        const rationaleLine = document.createElement('div');
        rationaleLine.className = 'insight-subtext';
        rationaleLine.textContent = `Rationale: ${ep.rationale}`;
        li.appendChild(rationaleLine);
      }
      el.epitopeList.appendChild(li);
    });
  }

  const chains = Array.isArray(target.chains) ? target.chains : [];
  el.chainList.innerHTML = '';
  if (chains.length === 0) {
    const li = document.createElement('li');
    li.className = 'insight-empty';
    li.textContent = 'No chain annotations available.';
    el.chainList.appendChild(li);
  } else {
    chains.forEach((chain) => {
      const li = document.createElement('li');
      if (chain.is_primary) {
        li.classList.add('chain-primary');
      }
      const title = document.createElement('span');
      title.className = 'insight-title';
      const chainId = chain.id || chain.chain_id || chain.chain || '—';
      title.textContent = `Chain ${chainId}`;
      li.appendChild(title);
      const summaryText = chain.summary || chain.description || '';
      if (summaryText) {
        const summary = document.createElement('div');
        summary.className = 'chain-summary';
        summary.textContent = summaryText;
        li.appendChild(summary);
      }
      const metaBits = [];
      if (chain.role) metaBits.push(chain.role);
      if (chain.length) metaBits.push(`${chain.length} aa`);
      if (chain.organism) metaBits.push(chain.organism);
      if (chain.polymer_type) metaBits.push(chain.polymer_type);
      if (metaBits.length) {
        const meta = document.createElement('div');
        meta.className = 'chain-meta';
        meta.textContent = metaBits.join(' · ');
        li.appendChild(meta);
      }
      if (chain.context && chain.context !== summaryText) {
        const context = document.createElement('div');
        context.className = 'insight-subtext';
        context.textContent = chain.context;
        li.appendChild(context);
      }
      if (chain.function) {
        const functionLine = document.createElement('div');
        functionLine.className = 'insight-subtext';
        functionLine.textContent = `Function: ${chain.function}`;
        li.appendChild(functionLine);
      }
      if (Array.isArray(chain.synonyms) && chain.synonyms.length) {
        const synonymsLine = document.createElement('div');
        synonymsLine.className = 'insight-subtext';
        synonymsLine.textContent = `Also known as: ${chain.synonyms.join(', ')}`;
        li.appendChild(synonymsLine);
      }
      el.chainList.appendChild(li);
    });
  }

  const antigen = target.antigen || {};
  el.antigenDetails.innerHTML = '';
  const detailEntries = [];
  if (antigen.accession) detailEntries.push({ key: 'Accession', value: antigen.accession });
  if (antigen.expressed_range) detailEntries.push({ key: 'Expressed range', value: antigen.expressed_range });
  if (antigen.expressed_length) detailEntries.push({ key: 'Expressed length', value: `${antigen.expressed_length} aa` });
  if (antigen.allowed_range) detailEntries.push({ key: 'Allowed epitope range', value: antigen.allowed_range });
  if (antigen.length) detailEntries.push({ key: 'Full length', value: `${antigen.length} aa` });
  if (antigen.molecular_weight_kda) {
    const massValue = typeof antigen.molecular_weight_kda === 'number'
      ? `${antigen.molecular_weight_kda} kDa`
      : String(antigen.molecular_weight_kda);
    detailEntries.push({ key: 'Molecular mass', value: massValue });
  }
  if (antigen.expression_host) detailEntries.push({ key: 'Expression host', value: antigen.expression_host });
  if (antigen.product_form) detailEntries.push({ key: 'Product form', value: antigen.product_form });
  if (antigen.tags) {
    const tags = Array.isArray(antigen.tags) ? antigen.tags.join(', ') : antigen.tags;
    if (tags) detailEntries.push({ key: 'Tags', value: tags });
  }
  if (antigen.catalog) detailEntries.push({ key: 'Catalog', value: antigen.catalog });
  if (antigen.notes) detailEntries.push({ key: 'Notes', value: antigen.notes });

  if (detailEntries.length === 0) {
    const empty = document.createElement('div');
    empty.className = 'empty-note';
    empty.textContent = 'No antigen metadata captured yet.';
    el.antigenDetails.appendChild(empty);
  } else {
    detailEntries.forEach((entry) => {
      const dt = document.createElement('dt');
      dt.textContent = entry.key;
      const dd = document.createElement('dd');
      if (entry.type === 'link' && typeof entry.value === 'string' && /^https?:/i.test(entry.value)) {
        const link = document.createElement('a');
        link.href = entry.value;
        link.target = '_blank';
        link.rel = 'noopener noreferrer';
        link.textContent = entry.value;
        dd.appendChild(link);
      } else {
        dd.textContent = entry.value;
      }
      el.antigenDetails.appendChild(dt);
      el.antigenDetails.appendChild(dd);
    });
  }

  updateCounterSelectionPlan();
}

function renderCurrentSelection() {
  if (!el.selectionCard) return;
  const pdb = state.currentPdb || '';
  if (el.selectionPdb) {
    el.selectionPdb.textContent = pdb || '—';
  }
  let presetName = '';
  if (state.activePresetId) {
    const preset = state.presets.find((p) => p.id === state.activePresetId);
    if (preset && preset.name) presetName = preset.name;
  }
  if (!presetName && state.targetDetails && state.targetDetails.target_name) {
    presetName = state.targetDetails.target_name;
  }
  if (!presetName && el.targetName && el.targetName.value) {
    presetName = el.targetName.value.trim();
  }
  if (el.selectionName) {
    el.selectionName.textContent = `Preset: ${presetName || '—'}`;
  }
  const antigenUrl = (() => {
    if (state.targetDetails && state.targetDetails.antigen && state.targetDetails.antigen.url) {
      return state.targetDetails.antigen.url;
    }
    if (el.antigenInput && el.antigenInput.value) {
      return el.antigenInput.value.trim();
    }
    return '';
  })();
  if (el.selectionAntigen) {
    el.selectionAntigen.textContent = 'Antigen URL: ';
    if (antigenUrl && /^https?:/i.test(antigenUrl)) {
      const link = document.createElement('a');
      link.href = antigenUrl;
      link.target = '_blank';
      link.rel = 'noopener noreferrer';
      link.textContent = antigenUrl;
      el.selectionAntigen.appendChild(link);
    } else if (antigenUrl) {
      const span = document.createElement('span');
      span.textContent = antigenUrl;
      el.selectionAntigen.appendChild(span);
    } else {
      const span = document.createElement('span');
      span.textContent = '—';
      el.selectionAntigen.appendChild(span);
    }
  }
}

function renderPresets() {
  if (!el.presetList) return;
  el.presetList.innerHTML = '';
  if (!state.presets || state.presets.length === 0) {
    const span = document.createElement('span');
    span.className = 'empty-note';
    span.textContent = 'No saved pairs yet.';
    el.presetList.appendChild(span);
    return;
  }
  state.presets.forEach((preset) => {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.textContent = preset.name || preset.pdb_id;
    btn.classList.toggle('active', preset.id === state.activePresetId);
    const metaParts = [preset.pdb_id];
    if (preset.antigen_url) metaParts.push(preset.antigen_url);
    if (preset.num_epitopes) metaParts.push(`${preset.num_epitopes} epitopes`);
    btn.title = metaParts.join(' · ');
    btn.addEventListener('click', () => selectPreset(preset));
    el.presetList.appendChild(btn);
  });
}

function getAvailableDesignEngines() {
  return state.designEngines.length > 0 ? state.designEngines : DEFAULT_ENGINE_OPTIONS;
}

function getDesignEngineMeta(engineId) {
  const engines = getAvailableDesignEngines();
  return engines.find((engine) => engine.engine_id === engineId) || null;
}

function normalizeEngineInfo(raw) {
  const fallback = DEFAULT_ENGINE_OPTIONS.find((engine) => engine.engine_id === (raw?.engine_id || ''))
    || DEFAULT_ENGINE_OPTIONS[0];
  const fallbackFields = fallback?.fields || [];
  const rawFields = Array.isArray(raw?.fields) ? raw.fields : [];
  const overrideMap = new Map(rawFields.map((field) => [field.field_id, field]));
  const mergedFields = fallbackFields.map((fallbackField) => {
    const override = overrideMap.get(fallbackField.field_id);
    if (override) overrideMap.delete(fallbackField.field_id);
    return {
      field_id: fallbackField.field_id,
      label: override?.label || fallbackField.label,
      description: override?.description ?? fallbackField.description ?? '',
      visible: override?.visible !== undefined ? Boolean(override.visible) : fallbackField.visible !== false,
      debug_only: override?.debug_only !== undefined ? Boolean(override.debug_only) : Boolean(fallbackField.debug_only),
    };
  });
  overrideMap.forEach((override, fieldId) => {
    mergedFields.push({
      field_id: fieldId,
      label: override.label || fieldId,
      description: override.description || '',
      visible: override.visible !== false,
      debug_only: Boolean(override.debug_only),
    });
  });
  return {
    engine_id: raw?.engine_id || fallback.engine_id,
    label: raw?.label || fallback.label,
    description: raw?.description || fallback.description || '',
    is_default: raw?.is_default !== undefined ? Boolean(raw.is_default) : Boolean(fallback.is_default),
    fields: mergedFields,
  };
}

function applyDesignEngineFieldConfig(engineMeta) {
  const meta = engineMeta || getDesignEngineMeta(state.selectedDesignEngine) || DEFAULT_ENGINE_OPTIONS[0];
  const fieldMap = new Map((Array.isArray(meta?.fields) ? meta.fields : []).map((field) => [field.field_id, field]));
  const activeFields = new Set();

  Object.entries(DESIGN_FIELD_ELEMENTS).forEach(([fieldId, refs]) => {
    if (!refs) return;
    const config = fieldMap.get(fieldId);
    const visible = config ? config.visible !== false : refs.defaultVisible;
    const labelText = config?.label || refs.defaultLabel || '';
    const description = config?.description || '';

    if (refs.label && labelText) {
      refs.label.textContent = labelText;
    }
    if (refs.help) {
      if (description) {
        refs.help.textContent = description;
        refs.help.hidden = false;
      } else {
        refs.help.textContent = '';
        refs.help.hidden = true;
      }
    }

    const target = refs.container || refs.wrapper;
    if (target) {
      if (visible) {
        target.classList.remove('engine-hidden');
        target.hidden = false;
      } else {
        target.classList.add('engine-hidden');
        target.hidden = true;
      }
    }

    if (refs.input) {
      refs.input.disabled = !visible;
      if (!visible && refs.input.type === 'checkbox') {
        refs.input.checked = false;
      }
    }

    if (visible) {
      activeFields.add(fieldId);
    }
  });

  state.activeDesignFields = activeFields;
}

function isDesignFieldActive(fieldId) {
  if (!(state.activeDesignFields instanceof Set)) return true;
  return state.activeDesignFields.has(fieldId);
}

function updateDesignEngineHint() {
  if (!el.designEngineHint) return;
  const engines = getAvailableDesignEngines();
  const selected = getDesignEngineMeta(state.selectedDesignEngine);
  const fallback = selected || engines.find((engine) => engine.is_default) || engines[0];
  const description = (selected && selected.description) || (fallback && fallback.description) || '';
  if (description) {
    el.designEngineHint.textContent = description;
    el.designEngineHint.hidden = false;
  } else {
    el.designEngineHint.textContent = '';
    el.designEngineHint.hidden = true;
  }
}

function renderDesignEngineOptions() {
  if (!el.designEngineSelect) return;
  const select = el.designEngineSelect;
  select.innerHTML = '';
  const engines = getAvailableDesignEngines();
  const defaultEngine = engines.find((engine) => engine.is_default) || engines[0] || null;
  if (!state.selectedDesignEngine || !engines.some((engine) => engine.engine_id === state.selectedDesignEngine)) {
    state.selectedDesignEngine = defaultEngine ? defaultEngine.engine_id : 'rfantibody';
  }
  engines.forEach((engine) => {
    const option = document.createElement('option');
    option.value = engine.engine_id;
    option.textContent = engine.label || engine.engine_id;
    if (engine.description) option.title = engine.description;
    select.appendChild(option);
  });
  select.value = state.selectedDesignEngine;
  select.disabled = engines.length <= 1;
  updateDesignEngineHint();
  applyDesignEngineFieldConfig(getDesignEngineMeta(state.selectedDesignEngine));
  toggleBoltzBinding(state.selectedDesignEngine === 'boltzgen');
}

function toggleBoltzBinding(show) {
  if (!el.designBoltzBindingRow) return;
  el.designBoltzBindingRow.hidden = !show;
}

async function loadDesignEngines() {
  if (!el.designEngineSelect) return;
  let previous = state.selectedDesignEngine;
  try {
    const res = await fetch('/api/designs/engines', { headers: { 'Cache-Control': 'no-cache' } });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const payload = await res.json();
    const engines = Array.isArray(payload?.engines) ? payload.engines : [];
    if (engines.length) {
      const normalized = engines.map((engine) => normalizeEngineInfo(engine));
      state.designEngines = normalized;
      if (!normalized.some((engine) => engine.engine_id === previous)) {
        const fallback = normalized.find((engine) => engine.is_default) || normalized[0];
        state.selectedDesignEngine = fallback ? fallback.engine_id : previous;
      }
    } else {
      state.designEngines = [];
    }
  } catch (err) {
    console.warn('Failed to load design engines', err);
    state.designEngines = [];
    state.selectedDesignEngine = previous;
  }
  renderDesignEngineOptions();
}

function applyTargetStatus(status = state.targetStatus) {
  const effectiveStatus = status || state.targetStatus;
  const isTargetJobActive = state.jobContext === 'target' && Boolean(state.jobPoller);
  if (effectiveStatus && el.pymolHotspots) {
    const readyForHotspots = Boolean(effectiveStatus.has_prep);
    el.pymolHotspots.disabled = isTargetJobActive || !readyForHotspots;
  }
  if (el.designSubmit && el.designStatus && !isTargetJobActive) {
    const readyForDesigns = Boolean(effectiveStatus && effectiveStatus.has_prep);
    if (readyForDesigns) {
      el.designSubmit.disabled = false;
      setBadge(el.designStatus, 'Ready for submission');
    } else {
      el.designSubmit.disabled = true;
      setBadge(el.designStatus, 'Awaiting target prep');
    }
  }
}

async function fetchTargetStatus(pdbId = state.currentPdb, { silent = true } = {}) {
  const target = (pdbId || '').trim();
  if (!target) return null;
  try {
    const res = await fetch(`/api/targets/${target}/status`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const payload = await res.json();
    state.targetStatus = payload;
    applyTargetStatus(payload);
    renderTargetInsights(payload);
    renderCurrentSelection();
    return payload;
  } catch (err) {
    if (!silent) {
      console.warn('Failed to retrieve target status', err);
    }
    return null;
  }
}

function syncActivePreset() {
  if (!state.presets || state.presets.length === 0) {
    state.activePresetId = null;
    renderPresets();
    renderCurrentSelection();
    return;
  }
  const pdb = state.currentPdb || (el.pdbInput?.value || '').trim().toUpperCase();
  const antigen = (el.antigenInput?.value || '').trim();
  const match = state.presets.find((preset) => {
    if (!pdb || preset.pdb_id !== pdb) return false;
    const presetUrl = (preset.antigen_url || '').trim();
    return presetUrl === antigen;
  });
  if (match) {
    state.activePresetId = match.id;
    if (el.presetIdInput) el.presetIdInput.value = match.id;
    if (el.targetName && !el.targetName.value) el.targetName.value = match.name || '';
    if (el.targetEpitopes && !el.targetEpitopes.value && match.num_epitopes) {
      el.targetEpitopes.value = match.num_epitopes;
    }
  } else {
    state.activePresetId = null;
    if (el.presetIdInput) el.presetIdInput.value = '';
  }
  renderPresets();
  renderCurrentSelection();
}

async function touchPreset(preset) {
  try {
    const payload = {
      preset_id: preset.id,
      name: preset.name,
      pdb_id: preset.pdb_id,
      antigen_url: preset.antigen_url,
      num_epitopes: preset.num_epitopes,
    };
    await fetch('/api/targets/presets', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
  } catch (err) {
    console.warn('Failed to update preset usage', err);
  }
}

function selectPreset(preset) {
  state.activePresetId = preset.id;
  if (el.presetIdInput) el.presetIdInput.value = preset.id;
  if (el.targetName) el.targetName.value = preset.name || '';
  if (el.antigenInput) el.antigenInput.value = preset.antigen_url || '';
  if (el.targetEpitopes) {
    if (preset.num_epitopes) {
      el.targetEpitopes.value = preset.num_epitopes;
    } else {
      el.targetEpitopes.value = '';
    }
  }
  if (el.pdbInput) el.pdbInput.value = preset.pdb_id;
  setCurrentPdb(preset.pdb_id, { loadAlignment: true });
  touchPreset(preset);
  renderPresets();
}

async function loadPresets() {
  try {
    const res = await fetch('/api/targets/presets');
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const payload = await res.json();
    const list = Array.isArray(payload)
      ? payload
      : Array.isArray(payload.presets)
        ? payload.presets
        : [];
    state.presets = list.map((entry) => ({
      id: entry.id,
      name: entry.name,
      pdb_id: (entry.pdb_id || '').toUpperCase(),
      antigen_url: entry.antigen_url || '',
      num_epitopes: entry.num_epitopes || null,
      last_used: entry.last_used || 0,
    }));
    state.presets.sort((a, b) => (b.last_used || 0) - (a.last_used || 0));
    syncActivePreset();
  } catch (err) {
    console.warn('Failed to load presets', err);
    renderPresets();
  }
}

function stopRankingsPolling() {
  if (state.rankingsPollTimer) {
    clearInterval(state.rankingsPollTimer);
    state.rankingsPollTimer = null;
  }
}

function scheduleRankingsPolling() {
  stopRankingsPolling();
}

function updateCatalogStatus(message, color = null) {
  setBadge(el.catalogStatus, message, color);
}

function resetCatalogPreview() {
  state.catalogHeaders = [];
  state.catalogRows = [];
  state.catalogSelectedRowIndex = null;
  if (el.catalogSelectedName) {
    el.catalogSelectedName.textContent = 'Select a file to preview';
  }
  if (el.catalogSelectedInfo) {
    el.catalogSelectedInfo.hidden = true;
  }
  if (el.catalogViewer) {
    el.catalogViewer.innerHTML = '';
  }
  if (el.catalogDownload) {
    el.catalogDownload.disabled = true;
  }
  updateCatalogStatus('', null);
}

function formatBytes(bytes) {
  if (!Number.isFinite(bytes) || bytes <= 0) {
    return '0 B';
  }
  const units = ['B', 'KB', 'MB', 'GB'];
  const idx = Math.min(units.length - 1, Math.floor(Math.log(bytes) / Math.log(1024)));
  const value = bytes / Math.pow(1024, idx);
  return `${value.toFixed(value >= 10 || idx === 0 ? 0 : 1)} ${units[idx]}`;
}

function formatTimestamp(ts) {
  if (!Number.isFinite(ts)) {
    return '';
  }
  return new Date(ts * 1000).toLocaleString();
}

function renderCatalogFiles() {
  if (!el.catalogList) return;
  el.catalogList.innerHTML = '';
  if (!state.catalogFiles.length) {
    const empty = document.createElement('li');
    empty.className = 'catalog-empty';
    empty.textContent = 'No catalog files found yet.';
    el.catalogList.appendChild(empty);
    return;
  }
  state.catalogFiles.forEach((file) => {
    const item = document.createElement('li');
    item.className = 'catalog-item';
    const button = document.createElement('button');
    button.type = 'button';
    button.textContent = file.name;
    button.dataset.filename = file.name;
    button.className = 'catalog-button';
    if (state.catalogSelectedFile === file.name) {
      button.classList.add('active');
    }
    button.addEventListener('click', () => {
      if (state.catalogSelectedFile !== file.name) {
        state.catalogSelectedFile = file.name;
      }
      selectCatalogFile(file.name);
    });
    const meta = document.createElement('span');
    meta.className = 'catalog-item-meta';
    const size = Number.isFinite(file.size_bytes) ? formatBytes(file.size_bytes) : '';
    const modified = Number.isFinite(file.modified_at) ? formatTimestamp(file.modified_at) : '';
    meta.textContent = [size, modified].filter(Boolean).join(' • ');
    item.appendChild(button);
    item.appendChild(meta);
    el.catalogList.appendChild(item);
  });
}

async function loadCatalogFiles() {
  if (!el.catalogRefresh) return;
  try {
    el.catalogRefresh.disabled = true;
    const response = await fetch('/api/target-generation/catalog');
    if (!response.ok) {
      throw new Error(`Failed to load catalog (HTTP ${response.status})`);
    }
    const payload = await response.json();
    state.catalogFiles = Array.isArray(payload.files) ? payload.files : [];
    renderCatalogFiles();
    if (
      state.catalogSelectedFile &&
      !state.catalogFiles.some((file) => file && file.name === state.catalogSelectedFile)
    ) {
      state.catalogSelectedFile = null;
      resetCatalogPreview();
    }
    if (!state.catalogSelectedFile && state.catalogFiles.length) {
      state.catalogSelectedFile = state.catalogFiles[0].name;
      selectCatalogFile(state.catalogSelectedFile);
    }
  } catch (err) {
    console.error('Failed to load catalog list', err);
    resetCatalogPreview();
    updateCatalogStatus(err.message || 'Unable to load catalog files.', 'rgba(248, 113, 113, 0.2)');
  } finally {
    el.catalogRefresh.disabled = false;
  }
}

function renderCatalogPreview(headers, rows, truncated) {
  if (!el.catalogViewer) return;
  el.catalogViewer.innerHTML = '';
  const table = document.createElement('table');
  table.className = 'catalog-table';
  if (headers && headers.length) {
    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    headers.forEach((header) => {
      const th = document.createElement('th');
      th.textContent = header;
      headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);
  }
  const tbody = document.createElement('tbody');
  rows.forEach((row, index) => {
    const tr = document.createElement('tr');
    tr.tabIndex = 0;
    if (state.catalogSelectedRowIndex === index) {
      tr.classList.add('selected');
    }
    row.forEach((cell) => {
      const td = document.createElement('td');
      td.textContent = cell;
      tr.appendChild(td);
    });
    const handleSelect = () => handleCatalogRowClick(index);
    tr.addEventListener('click', handleSelect);
    tr.addEventListener('keydown', (event) => {
      if (event.key === 'Enter' || event.key === ' ') {
        event.preventDefault();
        handleSelect();
      }
    });
    tbody.appendChild(tr);
  });
  table.appendChild(tbody);
  if (truncated) {
    const caption = document.createElement('caption');
    caption.textContent = 'Preview truncated. Increase the row limit to view more records.';
    table.appendChild(caption);
  }
  el.catalogViewer.appendChild(table);
}

function extractCatalogValue(headerMap, row, keys) {
  for (const key of keys) {
    const idx = headerMap.get(key);
    if (idx == null) continue;
    const value = row[idx];
    if (value == null) continue;
    const trimmed = String(value).trim();
    if (trimmed) {
      return trimmed;
    }
  }
  return '';
}

function applyCatalogRow(row) {
  if (!Array.isArray(row) || !state.catalogHeaders.length) {
    return { applied: false };
  }
  const headerMap = new Map();
  state.catalogHeaders.forEach((header, index) => {
    if (header == null) return;
    headerMap.set(String(header).trim().toLowerCase(), index);
  });
  const pdbRaw = extractCatalogValue(headerMap, row, [
    'chosen_pdb',
    'pdb_id',
    'pdb',
    'pdbid',
  ]);
  const antigenRaw = extractCatalogValue(headerMap, row, [
    'antigen_url',
    'antigen url',
    'vendor_url',
    'product_url',
    'url',
  ]);
  const gene = extractCatalogValue(headerMap, row, ['gene', 'gene_symbol', 'gene name']);
  const protein = extractCatalogValue(headerMap, row, ['protein_name', 'protein', 'target']);
  let appliedPdb = '';
  let appliedAntigen = '';
  if (antigenRaw && el.antigenInput) {
    el.antigenInput.value = antigenRaw;
    appliedAntigen = antigenRaw;
  }
  if (pdbRaw) {
    const normalized = pdbRaw.trim().toUpperCase();
    const condensed = normalized.replace(/[^0-9A-Z]/g, '');
    const truncated = condensed.length >= 4 ? condensed.slice(0, 4) : normalized;
    setCurrentPdb(truncated, { updateInput: true });
    appliedPdb = truncated;
  }
  if (!appliedPdb && !appliedAntigen) {
    return { applied: false };
  }
  const parts = [];
  if (gene) parts.push(gene);
  if (protein && protein !== gene) parts.push(protein);
  if (appliedPdb) parts.push(`PDB ${appliedPdb.slice(0, 4)}`);
  if (appliedAntigen) parts.push('Antigen URL copied');
  const message = parts.length ? `${parts.join(' · ')} ready in Target Setup.` : 'Target form updated.';
  showAlert(message, false);
  return { applied: true, appliedPdb, appliedAntigen };
}

function handleCatalogRowClick(index) {
  if (!state.catalogRows.length || index < 0 || index >= state.catalogRows.length) {
    return;
  }
  state.catalogSelectedRowIndex = index;
  if (el.catalogViewer) {
    const rows = el.catalogViewer.querySelectorAll('tbody tr');
    rows.forEach((tr, idx) => {
      if (idx === index) tr.classList.add('selected');
      else tr.classList.remove('selected');
    });
  }
  const result = applyCatalogRow(state.catalogRows[index]);
  if (result.applied) {
    updateCatalogStatus('Target form updated from catalog row.', 'rgba(34, 197, 94, 0.22)');
  } else {
    updateCatalogStatus('Row is missing a PDB ID or antigen URL.', 'rgba(248, 113, 113, 0.2)');
  }
}

async function selectCatalogFile(name) {
  if (!name) return;
  const limit = Number.parseInt(el.catalogLimit?.value ?? '', 10);
  state.catalogSelectedFile = name;
  state.catalogSelectedRowIndex = null;
  updateCatalogStatus('', null);
  if (el.catalogSelectedName) {
    el.catalogSelectedName.textContent = name;
  }
  if (el.catalogSelectedInfo) {
    el.catalogSelectedInfo.hidden = true;
  }
  if (el.catalogDownload) {
    el.catalogDownload.disabled = false;
  }
  if (el.catalogViewer) {
    el.catalogViewer.innerHTML = '';
    const spinner = document.createElement('div');
    spinner.className = 'catalog-spinner';
    spinner.textContent = 'Loading preview…';
    el.catalogViewer.appendChild(spinner);
  }
  try {
    const queryLimit = Number.isFinite(limit) ? limit : 200;
    const response = await fetch(
      `/api/target-generation/catalog/${encodeURIComponent(name)}?limit=${queryLimit}`,
    );
    if (!response.ok) {
      throw new Error(`Preview failed (HTTP ${response.status})`);
    }
    const payload = await response.json();
    state.catalogHeaders = Array.isArray(payload.headers) ? payload.headers : [];
    state.catalogRows = Array.isArray(payload.rows) ? payload.rows : [];
    renderCatalogPreview(state.catalogHeaders, state.catalogRows, Boolean(payload.truncated));
    const extra = [];
    if (Number.isFinite(payload.total_rows)) {
      extra.push(`${payload.total_rows} rows`);
    }
    if (Number.isFinite(payload.displayed_rows)) {
      extra.push(`${payload.displayed_rows} shown`);
    }
    if (el.catalogSelectedInfo) {
      if (extra.length) {
        el.catalogSelectedInfo.textContent = extra.join(' • ');
        el.catalogSelectedInfo.hidden = false;
      } else {
        el.catalogSelectedInfo.hidden = true;
      }
    }
    if (!state.catalogRows.length) {
      updateCatalogStatus('No rows found in this file.', 'rgba(248, 113, 113, 0.2)');
    }
  } catch (err) {
    console.error('Failed to load catalog preview', err);
    resetCatalogPreview();
    if (el.catalogViewer) {
      const errorBox = document.createElement('div');
      errorBox.className = 'catalog-error';
      errorBox.textContent = err.message || 'Failed to load preview.';
      el.catalogViewer.appendChild(errorBox);
    }
    updateCatalogStatus(err.message || 'Unable to load preview.', 'rgba(248, 113, 113, 0.2)');
  }
}

function downloadCatalogFile() {
  if (!state.catalogSelectedFile) return;
  window.open(
    `/api/target-generation/catalog/${encodeURIComponent(state.catalogSelectedFile)}/file`,
    '_blank',
  );
}

function openCatalogModal() {
  if (!el.catalogModal) return;
  el.catalogModal.hidden = false;
  document.body.classList.add('modal-open');
  if (!state.catalogFiles.length) {
    loadCatalogFiles();
  } else {
    renderCatalogFiles();
    if (state.catalogSelectedFile) {
      selectCatalogFile(state.catalogSelectedFile);
    } else if (state.catalogFiles.length) {
      state.catalogSelectedFile = state.catalogFiles[0].name;
      selectCatalogFile(state.catalogSelectedFile);
    }
  }
}

function closeCatalogModal() {
  if (!el.catalogModal) return;
  el.catalogModal.hidden = true;
  document.body.classList.remove('modal-open');
}

function openExportModal() {
  if (!el.exportModal) return;
  el.exportModal.hidden = false;
  document.body.classList.add('modal-open');
}

function closeExportModal() {
  if (!el.exportModal) return;
  el.exportModal.hidden = true;
  document.body.classList.remove('modal-open');
}

function setCurrentPdb(pdbId, options = {}) {
  if (!pdbId) return;
  const upper = pdbId.toUpperCase();
  const changed = state.currentPdb !== upper;
  state.currentPdb = upper;
  state.targetStatus = null;
  if (options.updateInput !== false && el.pdbInput) {
    el.pdbInput.value = upper;
  }
  if (el.pymolHotspots) {
    el.pymolHotspots.disabled = true;
  }
  state.activeRunLabel = '';
  if (changed || options.forceReset) {
    resetResultsPanel();
  }
  resetAnalysisPanel({ disableButton: true });
  updateActiveRunDisplay();
  syncActivePreset();
  renderCurrentSelection();
  renderRunHistory(upper);
  renderBoltzRunHistory(upper);
  fetchRunHistory(upper);
  if (getResultsSource() === 'boltzgen') {
    fetchBoltzRunHistory(upper);
  }
  if (options.loadAlignment) {
    loadAlignment(upper);
  }
  scheduleRankingsPolling();
  fetchTargetStatus(upper, { silent: true });
}

function refreshRunLabel(force = false) {
  if (!el.designRunLabel) return;
  if (force || !el.designRunLabel.value) {
    el.designRunLabel.value = timestampString();
  }
}

function setClusterStatus(mode, text) {
  if (!el.clusterStatus) return;
  el.clusterStatus.classList.remove('ok', 'warn', 'error');
  if (mode) {
    el.clusterStatus.classList.add(mode);
  }
  el.clusterStatus.textContent = text;
}

async function updateClusterStatus() {
  if (!el.clusterStatus) return;
  try {
    const res = await fetch('/api/cluster/status');
    if (!res.ok) {
      throw new Error(`HTTP ${res.status}`);
    }
    const payload = await res.json();
    if (payload.mock) {
      setClusterStatus('ok', 'Cluster: mock mode (no remote execution)');
      return;
    }
    if (payload.control_master) {
      const remoteRoot = payload.remote_root || 'remote';
      if (payload.remote_root_exists) {
        setClusterStatus('ok', `Cluster: connected - ${remoteRoot}`);
      } else {
        setClusterStatus('warn', `Cluster: connected - missing ${remoteRoot}`);
      }
    } else {
      const msg = payload.message || 'control master not active';
      const controlPath = payload.control_path || '';
      const controlPersist = payload.control_persist || '';
      const sshTarget = payload.ssh_target || '';
      const args = [];
      if (controlPath) args.push(`-o ControlPath=${controlPath}`);
      args.push('-o ControlMaster=auto');
      if (controlPersist) args.push(`-o ControlPersist=${controlPersist}`);
      const command = sshTarget
        ? `ssh ${args.join(' ')} -MNf ${sshTarget}`
        : 'ssh -o ControlMaster=auto -MNf <cluster>';
      setClusterStatus('error', `Cluster: not connected - ${msg}. Run ${command}`);
    }
  } catch (err) {
    setClusterStatus('error', `Cluster: status unavailable (${err.message || err})`);
  }
}

async function fetchJobList() {
  if (!el.jobList) return;
  try {
    const res = await fetch('/api/jobs');
    if (!res.ok) {
      throw new Error(`HTTP ${res.status}`);
    }
    state.jobs = await res.json();
    renderJobList();
    if (state.currentPdb) {
      fetchRunHistory(state.currentPdb);
      if (getResultsSource() === 'boltzgen') {
        fetchBoltzRunHistory(state.currentPdb);
      }
    }
  } catch (err) {
    console.error('Failed to fetch jobs', err);
  }
}

async function fetchRunHistory(pdbId) {
  if (!pdbId || !el.runHistory) return;
  const upper = pdbId.toUpperCase();
  try {
    const res = await fetch(`/api/targets/${upper}/runs`);
    if (!res.ok) {
      throw new Error(`HTTP ${res.status}`);
    }
    const payload = await res.json();
    state.runHistory[upper] = Array.isArray(payload) ? payload : [];
    if (state.currentPdb === upper) {
      renderRunHistory(upper);
    }
  } catch (err) {
    console.error('Failed to fetch run history', err);
    if (!state.runHistory[upper]) {
      state.runHistory[upper] = [];
    }
    if (state.currentPdb === upper) {
      renderRunHistory(upper);
    }
  }
}

async function fetchBoltzRunHistory(pdbId) {
  if (!pdbId || !el.boltzRunList) return;
  const upper = pdbId.toUpperCase();
  try {
    const res = await fetch(`/api/targets/${upper}/boltzgen/runs`);
    if (!res.ok) {
      throw new Error(`HTTP ${res.status}`);
    }
    const payload = await res.json();
    state.boltzRunHistory[upper] = Array.isArray(payload?.runs) ? payload.runs : [];
    if (state.currentPdb === upper) {
      renderBoltzRunHistory(upper);
    }
  } catch (err) {
    console.warn('Failed to fetch BoltzGen runs', err);
    if (!state.boltzRunHistory[upper]) {
      state.boltzRunHistory[upper] = [];
    }
    if (state.currentPdb === upper) {
      renderBoltzRunHistory(upper);
    }
  }
}

function highlightRunChip(selectedLabel = state.activeRunLabel) {
  if (!el.runHistory) return;
  const target = (selectedLabel || '').trim();
  el.runHistory.querySelectorAll('button[data-run-label]').forEach((btn) => {
    btn.classList.toggle('active', btn.dataset.runLabel === target && target !== '');
  });
}

function renderRunHistory(pdbId) {
  if (!el.runHistory) return;
  const upper = pdbId ? pdbId.toUpperCase() : state.currentPdb;
  const runs = (upper && state.runHistory[upper]) || [];

  el.runHistory.innerHTML = '';
  el.runHistory.classList.toggle('empty', runs.length === 0);

  if (!runs.length) {
    el.runHistory.textContent = 'No assessments found locally or on the cluster.';
    if (el.refreshResultsBtn) el.refreshResultsBtn.disabled = true;
    if (el.syncResultsBtn) el.syncResultsBtn.disabled = !state.currentPdb;
  } else {
    const frag = document.createDocumentFragment();
    let selected = state.activeRunLabel || '';
    if (!state.activeRunLabel) {
      const firstRunLabel = runs[0]?.run_label || '';
      if (firstRunLabel) {
        state.activeRunLabel = firstRunLabel;
        if (el.resultsRunLabel) el.resultsRunLabel.value = firstRunLabel;
        updateActiveRunDisplay();
        selected = state.activeRunLabel || '';
      }
    }
    runs.forEach((run) => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.dataset.runLabel = run.run_label;
      btn.dataset.origin = run.origin || 'local';
      btn.classList.toggle('has-local', !!run.available_local);
      btn.classList.toggle('has-remote', !!run.available_remote);
      btn.classList.toggle('active', run.run_label === selected && selected !== '');

      const labelSpan = document.createElement('span');
      labelSpan.className = 'label';
      labelSpan.textContent = run.run_label;

      const originSpan = document.createElement('span');
      originSpan.className = 'origin';
      if (run.available_local && run.available_remote) {
        originSpan.textContent = 'local + cluster';
      } else if (run.available_local) {
        originSpan.textContent = 'local';
      } else if (run.available_remote) {
        originSpan.textContent = 'cluster';
      } else {
        originSpan.textContent = 'pending';
      }

      const infoParts = [`Updated ${new Date(run.updated_at * 1000).toLocaleString()}`];
      if (typeof run.total_rows === 'number') {
        infoParts.push(`${run.total_rows} rows`);
      }
      if (run.local_path) {
        infoParts.push(`local: ${run.local_path}`);
      }
      if (run.remote_path) {
        infoParts.push(`cluster: ${run.remote_path}`);
      }
      btn.title = infoParts.join(' · ');

      btn.append(labelSpan, originSpan);
      btn.addEventListener('click', () => {
        state.activeRunLabel = run.run_label || '';
        updateActiveRunDisplay();
        if (el.resultsRunLabel) {
          el.resultsRunLabel.value = run.run_label || '';
        }
        highlightRunChip(run.run_label);
        if (!run.available_local) {
          const originHint = run.available_remote
            ? 'Use “Download assessments” to copy this run locally before reloading.'
            : 'No rankings available for this run yet.';
          showAlert(originHint);
        }
      });
      frag.appendChild(btn);
    });
    el.runHistory.appendChild(frag);
    if (el.refreshResultsBtn) el.refreshResultsBtn.disabled = false;
    if (el.syncResultsBtn) el.syncResultsBtn.disabled = false;
  }

  highlightRunChip();

  if (el.runLabelOptions) {
    el.runLabelOptions.innerHTML = '';
    runs.forEach((run) => {
      const option = document.createElement('option');
      option.value = run.run_label;
      el.runLabelOptions.appendChild(option);
    });
  }

  if (el.resultsRunLabel) {
    const firstLocal = runs.find((run) => run.available_local);
    if (firstLocal && !el.resultsRunLabel.value) {
      el.resultsRunLabel.placeholder = firstLocal.run_label;
    } else if (runs.length === 0) {
      el.resultsRunLabel.placeholder = 'latest';
    } else if (!el.resultsRunLabel.value && runs[0]) {
      el.resultsRunLabel.placeholder = runs[0].run_label;
    }
  }

  if (el.assessSubmit) {
    const hasCandidates = runs.length > 0 || !!(el.resultsRunLabel && el.resultsRunLabel.value.trim());
    const assessRunning = state.jobContext === 'assess' && state.jobPoller;
    if (!hasCandidates && !assessRunning) {
      el.assessSubmit.disabled = true;
    } else if (hasCandidates && !assessRunning) {
      el.assessSubmit.disabled = false;
    }
  }
}

function renderBoltzRunHistory(pdbId) {
  if (!el.boltzRunList) return;
  const isBoltz = getResultsSource() === 'boltzgen';
  el.boltzRunList.hidden = !isBoltz;
  if (!isBoltz) return;
  const upper = pdbId ? pdbId.toUpperCase() : state.currentPdb;
  const runs = (upper && state.boltzRunHistory[upper]) || [];
  el.boltzRunList.innerHTML = '';
  if (!runs.length) {
    el.boltzRunList.textContent = 'No BoltzGen runs found locally. Sync results to populate this list.';
    return;
  }
  const frag = document.createDocumentFragment();
  runs.forEach((run) => {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'run-chip';
    btn.dataset.runLabel = run.run_label;
    const label = document.createElement('span');
    label.className = 'label';
    label.textContent = run.run_label;
    const specs = document.createElement('span');
    specs.className = 'origin';
    const specText = (run.specs || [])
      .map((spec) => (spec.has_metrics ? `✔ ${spec.name}` : `… ${spec.name}`))
      .join(', ');
    specs.textContent = specText || 'no specs detected';
    btn.title = `Updated ${new Date(run.updated_at * 1000).toLocaleString()}`;
    btn.append(label, specs);
    btn.addEventListener('click', () => {
      if (el.resultsRunLabel) {
        el.resultsRunLabel.value = run.run_label;
      }
      if (el.boltzSpec && Array.isArray(run.specs) && run.specs.length === 1) {
        el.boltzSpec.value = run.specs[0].name;
      }
      highlightRunChip(run.run_label);
    });
    frag.appendChild(btn);
  });
  el.boltzRunList.appendChild(frag);
}

function renderJobList() {
  if (!el.jobList) return;
  el.jobList.innerHTML = '';
  state.jobs.forEach((job) => {
    const li = document.createElement('li');
    li.dataset.jobId = job.job_id;
    if (state.selectedJobId === job.job_id) {
      li.classList.add('selected');
    }
    const label = document.createElement('div');
    label.className = 'label';
    label.textContent = job.label || job.job_id;
    const meta = document.createElement('div');
    meta.className = 'meta';
    const parts = [];
    if (job.pdb_id) parts.push(job.pdb_id.toUpperCase());
    parts.push(job.kind);
    parts.push(job.status);
    if (job.run_label) parts.push(`run ${job.run_label}`);
    const created = new Date(job.created_at * 1000).toLocaleString();
    parts.push(created);
    meta.textContent = parts.join(' · ');
    li.appendChild(label);
    li.appendChild(meta);
    li.addEventListener('click', () => {
      state.selectedJobId = job.job_id;
      renderJobList();
      if (job.pdb_id) {
        setCurrentPdb(job.pdb_id);
      }
      if (job.run_label && el.resultsRunLabel) {
        el.resultsRunLabel.value = job.run_label;
      }
      startJobPolling(job.job_id, null, { manual: true });
      highlightRunChip(job.run_label || '');
    });
    el.jobList.appendChild(li);
  });

  if (!state.currentPdb && state.jobs.length > 0) {
    const first = state.jobs[0];
    if (first.pdb_id) {
      setCurrentPdb(first.pdb_id);
      if (first.run_label && el.resultsRunLabel) {
        el.resultsRunLabel.value = first.run_label;
        highlightRunChip(first.run_label);
      }
    }
  }
}

function showAlert(message, isError = true) {
  if (!el.jobAlert) return;
  el.jobAlert.textContent = message;
  el.jobAlert.hidden = false;
  el.jobAlert.classList.toggle('success', !isError);
}

function resetJobLog(message = 'Waiting for commands...') {
  el.jobLog.textContent = message;
}

function appendLog(line) {
  if (!line) return;
  const atBottom = Math.abs(el.jobLog.scrollHeight - el.jobLog.scrollTop - el.jobLog.clientHeight) < 8;
  el.jobLog.textContent += `\n${line}`;
  if (atBottom) {
    el.jobLog.scrollTop = el.jobLog.scrollHeight;
  }
}

async function queueTargetInit(event) {
  event.preventDefault();
  const pdbId = el.pdbInput.value.trim().toUpperCase();
  if (!pdbId) return;
  setCurrentPdb(pdbId);
  setBadge(el.targetBadge, 'Queued…');
  el.targetSubmit.disabled = true;
  el.pymolHotspots.disabled = true;
  state.targetStatus = null;
  el.jobAlert.hidden = true;
  resetJobLog('Submitting init-target job…');

  const payload = {
    pdb_id: pdbId,
    antigen_url: el.antigenInput.value.trim() || null,
    preset_name: el.targetName?.value.trim() || null,
    num_epitopes: (() => {
      const raw = el.targetEpitopes?.value || '';
      if (!raw) return null;
      const parsed = Number(raw);
      if (!Number.isFinite(parsed) || parsed <= 0) return null;
      return Math.round(parsed);
    })(),
    decide_scope_prompt: (() => {
      const value = el.targetDecidePrompt?.value || '';
      const trimmed = value.trim();
      return trimmed ? trimmed : null;
    })(),
    force_refresh: el.forceCheckbox.checked,
    run_decide_scope: el.decideCheckbox.checked,
    run_prep: el.prepCheckbox.checked,
  };

  try {
    const res = await fetch('/api/targets/init', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Request failed with ${res.status}`);
    }
    const body = await res.json();
    state.selectedJobId = body.job_id;
    startJobPolling(body.job_id, 'target');
    fetchJobList();
  } catch (err) {
    setBadge(el.targetBadge, 'Error', 'rgba(248, 113, 113, 0.2)');
    showAlert(err.message || String(err));
    el.targetSubmit.disabled = false;
  }
}

async function queueDesignRun() {
  if (!state.currentPdb) {
    showAlert('Initialize a target first.');
    return;
  }
  const payload = {
    pdb_id: state.currentPdb,
    model_engine: state.selectedDesignEngine || 'rfantibody',
  };

  const totalDesignsRaw = Number(el.designTotal.value) || 90;
  payload.total_designs = Math.max(1, Math.round(totalDesignsRaw));

  if (isDesignFieldActive('num_sequences')) {
    const sequences = Number(el.designNumSeq.value) || 1;
    payload.num_sequences = Math.max(1, Math.round(sequences));
  }

  if (isDesignFieldActive('temperature')) {
    const tempRaw = Number(el.designTemp.value);
    payload.temperature = Number.isFinite(tempRaw) ? tempRaw : 0.1;
  }

  if (isDesignFieldActive('binder_chain_id')) {
    const binder = el.designBinderChain?.value.trim() || null;
    payload.binder_chain_id = binder || null;
  }

  if (isDesignFieldActive('af3_seed')) {
    const rawSeed = Number(el.designAf3Seed?.value ?? 1);
    const fallbackSeed = 1;
    payload.af3_seed = Number.isFinite(rawSeed) ? Math.max(0, Math.floor(rawSeed)) : fallbackSeed;
  }

  payload.run_label = el.designRunLabel.value.trim() || null;
  payload.run_assess = true;

  if (isDesignFieldActive('rfdiff_crop_radius')) {
    payload.rfdiff_crop_radius = el.designEnableRfdiffCrop && el.designEnableRfdiffCrop.checked
      ? DEFAULT_RFDIFF_CROP_RADIUS
      : null;
  }
  if (isDesignFieldActive('boltzgen_crop_radius')) {
    const raw = Number(el.designBoltzCropRadius?.value ?? DEFAULT_BOLTZGEN_CROP_RADIUS);
    payload.boltzgen_crop_radius = Number.isFinite(raw) && raw > 0 ? raw : null;
  }
  if (state.selectedDesignEngine === 'boltzgen' && el.designBoltzBinding) {
    const binding = el.designBoltzBinding.value.trim();
    if (binding) {
      payload.boltz_binding = binding;
    }
  }

  el.designSubmit.disabled = true;
  setBadge(el.designStatus, 'Submitting…');
  resetJobLog('Submitting design pipeline…');
  el.jobAlert.hidden = true;
  if (el.targetSubmit) el.targetSubmit.disabled = true;
  if (el.assessSubmit) el.assessSubmit.disabled = true;

  try {
    const res = await fetch('/api/designs/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Request failed with ${res.status}`);
    }
    const body = await res.json();
    state.selectedJobId = body.job_id;
    startJobPolling(body.job_id, 'design');
    fetchJobList();
  } catch (err) {
    el.designSubmit.disabled = false;
    setBadge(el.designStatus, 'Error', 'rgba(248, 113, 113, 0.2)');
    showAlert(err.message || String(err));
    if (el.targetSubmit) el.targetSubmit.disabled = false;
    if (el.assessSubmit) el.assessSubmit.disabled = false;
  }
}

async function queueAssessmentRun() {
  if (!state.currentPdb) {
    showAlert('Initialize a target first.');
    return;
  }
  let runLabel = el.resultsRunLabel?.value.trim() || '';
  if (!runLabel && state.activeRunLabel) {
    runLabel = state.activeRunLabel;
  }
  if (!runLabel) {
    showAlert('Enter a run label to assess.');
    return;
  }
  const binderChain = el.designBinderChain?.value.trim() || 'H';
  const payload = {
    pdb_id: state.currentPdb,
    run_label: runLabel,
    binder_chain_id: binderChain,
    include_keyword: runLabel,
  };

  if (el.assessSubmit) el.assessSubmit.disabled = true;
  setBadge(el.resultsMeta, 'Submitting assessment…');
  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/assess`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Request failed with ${res.status}`);
    }
    const body = await res.json();
    state.selectedJobId = body.job_id;
    startJobPolling(body.job_id, 'assess');
    fetchJobList();
  } catch (err) {
    if (el.assessSubmit) el.assessSubmit.disabled = false;
    setBadge(el.resultsMeta, 'Assessment submission failed', 'rgba(248, 113, 113, 0.25)');
    showAlert(err.message || String(err));
  }
}

async function runExport() {
  if (!state.currentPdb || !state.rankingsResponse) {
    showAlert('Load rankings before exporting.');
    return;
  }
  if (!ensureAf3Results('Exporting sequences')) {
    return;
  }
  const payload = {
    pdb_id: state.currentPdb,
    rankings_path: state.rankingsResponse.source_path,
    top_n: Number(el.exportTopN.value) || 48,
    codon_host: el.exportCodonHost.value,
    gc_target: el.exportGc.value ? Number(el.exportGc.value) : null,
    gc_window: el.exportGcWindow.value ? Number(el.exportGcWindow.value) : null,
    prefix_raw: el.exportPrefix.value || null,
    suffix_raw: el.exportSuffix.value || null,
    use_dnachisel: el.exportDnac.checked,
  };

  el.exportButton.disabled = true;
  if (el.exportOpen) el.exportOpen.disabled = true;
  resetJobLog('Submitting export job…');
  el.jobAlert.hidden = true;

  try {
    const res = await fetch('/api/exports', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Request failed with ${res.status}`);
    }
    const body = await res.json();
    state.selectedJobId = body.job_id;
    startJobPolling(body.job_id, 'export');
    fetchJobList();
    closeExportModal();
  } catch (err) {
    showAlert(err.message || String(err));
    el.exportButton.disabled = false;
    if (el.exportOpen) el.exportOpen.disabled = false;
  }
}

async function runGoldenGatePlanner() {
  if (!state.currentPdb || !state.rankingsResponse) {
    showAlert('Load rankings before building the Golden Gate plan.');
    return;
  }
  if (!ensureAf3Results('Golden Gate planning')) {
    return;
  }

  const sequenceColumn = el.librarySequenceColumn?.value.trim();
  const payload = {
    pdb_id: state.currentPdb,
    rankings_path: state.rankingsResponse.source_path,
    top_n: Number(el.libraryTopN?.value || 0) || 48,
    codon_host: el.libraryCodonHost?.value || 'yeast',
    sequence_column: sequenceColumn ? sequenceColumn : null,
    upstream_flank: (el.libraryUpstream?.value || 'GGAG').trim().toUpperCase(),
    downstream_flank: (el.libraryDownstream?.value || 'CGCT').trim().toUpperCase(),
  };

  if (payload.upstream_flank.length < 2 || payload.downstream_flank.length < 2) {
    showAlert('Flanking sequences should be at least 2 bases long.');
    return;
  }

  state.librarySummary = null;
  renderLibrarySummary(null);
  if (el.libraryRun) el.libraryRun.disabled = true;
  if (el.libraryStatus) setBadge(el.libraryStatus, 'Submitting…');
  resetJobLog('Submitting Golden Gate planner…');
  el.jobAlert.hidden = true;

  try {
    const res = await fetch('/api/golden-gate', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Request failed with ${res.status}`);
    }
    const body = await res.json();
    state.selectedJobId = body.job_id;
    startJobPolling(body.job_id, 'library');
    fetchJobList();
  } catch (err) {
    const message = err && err.message ? err.message : String(err);
    showAlert(message);
    if (el.libraryRun) el.libraryRun.disabled = false;
    if (el.libraryStatus) setBadge(el.libraryStatus, `Failed — ${message}`, 'rgba(248, 113, 113, 0.25)');
  }
}

async function fetchRankings(options = {}) {
  const { silent = false } = options;
  if (!state.currentPdb) {
    if (!silent) showAlert('Initialize target first.');
    return;
  }
  if (state.clusterRetryTimer) {
    clearTimeout(state.clusterRetryTimer);
    state.clusterRetryTimer = null;
  }
  if (state.rankingsFetching) return;
  state.rankingsFetching = true;

  let runLabel = el.resultsRunLabel?.value.trim() || '';
  if (!runLabel && state.activeRunLabel) {
    runLabel = state.activeRunLabel;
  }
  const source = getResultsSource();
  const limitVal = Number(el.resultsLimit?.value) || null;
  const params = new URLSearchParams();
  if (runLabel) params.append('run_label', runLabel);
  if (limitVal) params.append('limit', String(limitVal));
  let requestedSpec = null;
  if (source === 'boltzgen') {
    const specValue = el.boltzSpec?.value.trim();
    if (specValue) {
      params.append('spec', specValue);
      requestedSpec = specValue;
    }
  }

  if (!silent && el.refreshResultsBtn) {
    el.refreshResultsBtn.disabled = true;
  }
  if (!silent) {
    const label = source === 'boltzgen' ? 'Loading BoltzGen metrics' : 'Loading rankings';
    setBadge(el.resultsMeta, `${label} for ${state.currentPdb}…`);
    el.resultsMeta.hidden = false;
  }

  const requestedSignature = analysisSignature(
    state.currentPdb,
    runLabel || '',
    state.rankingsResponse?.source_path || '',
    state.rankingsResponse?.engine_id || 'rfantibody',
  );
  const shouldPreserveAnalysis =
    state.analysisContext && state.analysisContext.signature === requestedSignature;

  if (!shouldPreserveAnalysis) {
    resetAnalysisPanel({ disableButton: true });
  } else if (el.analysisRun) {
    el.analysisRun.disabled = true;
  }

  try {
    const endpoint =
      source === 'boltzgen'
        ? `/api/targets/${state.currentPdb}/boltzgen/results?${params.toString()}`
        : `/api/targets/${state.currentPdb}/rankings?${params.toString()}`;
    const res = await fetch(endpoint);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    const newSignature = analysisSignature(
      state.currentPdb,
      payload.run_label || '',
      payload.source_path || '',
      payload.engine_id || 'rfantibody',
    );
    const contextMatches =
      state.analysisContext && state.analysisContext.signature === newSignature;
    const previousSummary = state.librarySummary;
    const newSource = payload.source_path ? String(payload.source_path) : '';
    const keepSummary =
      Boolean(previousSummary && previousSummary.rankings_source && newSource) &&
      String(previousSummary.rankings_source) === newSource;
    state.rankings = payload.rows || [];
    state.rankingsTotal = state.rankings.length;
    state.rankingsDisplayLimit = RESULTS_TABLE_DISPLAY_LIMIT;
    state.rankingsShowAll = state.rankings.length <= state.rankingsDisplayLimit;
    state.rankingsResponse = payload;
    if (!keepSummary) {
      state.librarySummary = null;
    }
    ensureValidTableSort(state.rankings);
    state.scatterLayout = [];
    state.selectedDesign = null;
    state.activeRunLabel = payload.run_label || '';
    state.galleryAvailable = Boolean(payload.gallery_path);
    if (source === 'boltzgen') {
      const inferred = inferBoltzSpecName(payload) || requestedSpec || null;
      state.boltzActiveSpec = inferred;
      if (inferred && el.boltzSpec && el.boltzSpec.value.trim() !== inferred) {
        el.boltzSpec.value = inferred;
      }
    } else {
      state.boltzActiveSpec = null;
    }
    updateActiveRunDisplay();
    if (el.resultsRunLabel) {
      el.resultsRunLabel.value = state.activeRunLabel;
    }
    el.binderDetail.hidden = true;
    if (!contextMatches) {
      resetAnalysisPanel();
    }
    state.resultsSource = source;
    renderResults();
    if (!silent) {
      const label = state.activeRunLabel || 'latest';
      const suffix = source === 'boltzgen' ? 'BoltzGen metrics' : 'designs';
      setBadge(el.resultsMeta, `${payload.rows.length} ${suffix} · run ${label}`);
    }
    highlightRunChip(state.activeRunLabel);
    if (state.currentPdb) {
      fetchRunHistory(state.currentPdb);
      if (source === 'boltzgen') {
        fetchBoltzRunHistory(state.currentPdb);
      }
    }
    renderPresets();
    applyResultsEngineCapabilities();
    stopRankingsPolling();
    renderBoltzRunHistory(state.currentPdb);
  } catch (err) {
    if (!silent) {
      setBadge(el.resultsMeta, `Rankings failed: ${err.message || err}`, 'rgba(248, 113, 113, 0.2)');
    } else {
      console.warn('Rankings refresh failed', err);
    }
  } finally {
    if (!silent && el.refreshResultsBtn) {
      el.refreshResultsBtn.disabled = false;
    }
    applyResultsEngineCapabilities();
    state.rankingsFetching = false;
  }
}

function sortRankings() {
  const { key, dir } = state.tableSort;
  const rows = [...state.rankings];
  rows.sort((a, b) => {
    const va = a[key];
    const vb = b[key];
    const multiplier = dir === 'asc' ? 1 : -1;
    const numA = typeof va === 'number' ? va : parseFloat(va);
    const numB = typeof vb === 'number' ? vb : parseFloat(vb);
    if (!Number.isNaN(numA) && !Number.isNaN(numB)) {
      return (numA - numB) * multiplier;
    }
    return String(va || '').localeCompare(String(vb || '')) * multiplier;
  });
  state.rankings = rows;
  state.rankingsTotal = rows.length;
}

function renderResultsTable() {
  if (!el.resultsTable) return;
  const tbody = el.resultsTable.querySelector('tbody');
  if (!tbody) return;
  tbody.innerHTML = '';
  const limitValue = state.rankingsShowAll
    ? null
    : Math.max(0, state.rankingsDisplayLimit || RESULTS_TABLE_DISPLAY_LIMIT);
  const rowsToRender = limitValue ? state.rankings.slice(0, limitValue) : state.rankings;
  rowsToRender.forEach((row) => {
    const tr = document.createElement('tr');
    tr.dataset.design = row.design_name;
    const epitopeLabel = getEpitopeLabel(row);
    tr.innerHTML = `
      <td>${row.design_name}</td>
      <td>${row.iptm !== null && row.iptm !== undefined ? row.iptm.toFixed(3) : '—'}</td>
      <td>${row.rmsd_diego !== null && row.rmsd_diego !== undefined ? row.rmsd_diego.toFixed(3) : '—'}</td>
      <td>${row.ipsae_min !== null && row.ipsae_min !== undefined ? row.ipsae_min.toFixed(3) : '—'}</td>
      <td>${row.hotspot_min_distance !== null && row.hotspot_min_distance !== undefined ? row.hotspot_min_distance.toFixed(3) : '—'}</td>
      <td>${escapeHtml(epitopeLabel)}</td>
    `;
    tr.addEventListener('click', () => selectDesign(row.design_name));
    tbody.appendChild(tr);
  });
  updateResultsTruncateNotice();
  if (state.selectedDesign) {
    const selectedRow = Array.from(tbody.querySelectorAll('tr')).find(
      (tr) => tr.dataset.design === state.selectedDesign,
    );
    if (selectedRow) {
      selectedRow.classList.add('selected');
    } else if (!state.rankingsShowAll && el.binderDetail) {
      el.binderDetail.hidden = true;
    }
  }
}

function updateResultsTruncateNotice() {
  if (!el.resultsTruncate) return;
  const total = state.rankingsTotal || state.rankings.length || 0;
  const limit = state.rankingsDisplayLimit || RESULTS_TABLE_DISPLAY_LIMIT;
  if (!total || total <= limit) {
    el.resultsTruncate.hidden = true;
    if (el.resultsTruncateText) el.resultsTruncateText.textContent = '';
    if (el.resultsTruncateToggle) {
      el.resultsTruncateToggle.textContent = 'Show all rows';
      el.resultsTruncateToggle.setAttribute('aria-pressed', 'false');
    }
    return;
  }
  const showingAll = Boolean(state.rankingsShowAll);
  const showingCount = showingAll ? total : Math.min(limit, total);
  if (el.resultsTruncateText) {
    el.resultsTruncateText.textContent = showingAll
      ? `Showing all ${showingCount.toLocaleString()} rows.`
      : `Showing first ${showingCount.toLocaleString()} of ${total.toLocaleString()} rows.`;
  }
  if (el.resultsTruncateToggle) {
    el.resultsTruncateToggle.textContent = showingAll
      ? `Show first ${limit.toLocaleString()} rows`
      : 'Show all rows';
    el.resultsTruncateToggle.setAttribute('aria-pressed', showingAll ? 'true' : 'false');
  }
  el.resultsTruncate.hidden = false;
}

function computeScatterLayout(points) {
  if (!el.scatterCanvas) {
    state.scatterLayout = [];
    state.scatterStats = null;
    state.scatterColorMap = new Map();
    return;
  }
  const axes = getActiveScatterAxes();
  const xConfig = getScatterFieldConfig(axes.x);
  const yConfig = getScatterFieldConfig(axes.y);
  if (!xConfig || !yConfig) {
    state.scatterLayout = [];
    state.scatterStats = {
      xBins: [],
      yBins: [],
      total: 0,
      passingCount: 0,
      thresholds: { x: null, y: null },
      fields: { x: axes.x, y: axes.y },
      axis: {
        xMin: 0,
        xMax: 1,
        yMin: 0,
        yMax: 1,
        xRange: 1,
        yRange: 1,
        margin: { left: 72, right: 32, top: 48, bottom: 72 },
      },
    };
    state.scatterColorMap = new Map();
    return;
  }

  const source = Array.isArray(points) ? points : [];
  const filtered = source
    .map((p) => {
      const xValue = getScatterFieldValue(p, axes.x);
      const yValue = getScatterFieldValue(p, axes.y);
      return {
        ...p,
        xValue,
        yValue,
        epitopeLabel: getEpitopeLabel(p),
      };
    })
    .filter((p) => Number.isFinite(p.xValue) && Number.isFinite(p.yValue));

  const thresholdsRaw = state.scatterThresholds || {};
  const thresholds = {
    x: Number.isFinite(thresholdsRaw[axes.x]) ? thresholdsRaw[axes.x] : null,
    y: Number.isFinite(thresholdsRaw[axes.y]) ? thresholdsRaw[axes.y] : null,
  };

  const margin = { left: 72, right: 32, top: 48, bottom: 72 };

  if (filtered.length === 0) {
    state.scatterLayout = [];
    state.scatterStats = {
      xBins: [],
      yBins: [],
      total: 0,
      passingCount: 0,
      thresholds,
      fields: { x: axes.x, y: axes.y },
      axis: {
        xMin: 0,
        xMax: 1,
        yMin: 0,
        yMax: 1,
        xRange: 1,
        yRange: 1,
        margin,
      },
    };
    state.scatterColorMap = new Map();
    return;
  }

  const epitopeKeys = filtered.map((p) => (p.epitopeLabel ? p.epitopeLabel : '__default__'));
  const uniqueKeys = Array.from(new Set(epitopeKeys));
  uniqueKeys.sort((a, b) => {
    if (a === '__default__') return 1;
    if (b === '__default__') return -1;
    return a.localeCompare(b, undefined, { sensitivity: 'base' });
  });
  const colorMap = new Map();
  uniqueKeys.forEach((key, idx) => {
    const color = SCATTER_COLOR_PALETTE[idx % SCATTER_COLOR_PALETTE.length];
    colorMap.set(key, color);
  });
  state.scatterColorMap = colorMap;
  const width = el.scatterCanvas.width;
  const height = el.scatterCanvas.height;
  const plotWidth = Math.max(1, width - margin.left - margin.right);
  const plotHeight = Math.max(1, height - margin.top - margin.bottom);
  const xs = filtered.map((p) => p.xValue);
  const ys = filtered.map((p) => p.yValue);
  let xMin = Math.min(...xs);
  let xMax = Math.max(...xs);
  let yMin = Math.min(...ys);
  let yMax = Math.max(...ys);
  if (Number.isFinite(xConfig.defaultMin)) xMin = Math.min(xMin, xConfig.defaultMin);
  if (Number.isFinite(xConfig.defaultMax)) xMax = Math.max(xMax, xConfig.defaultMax);
  if (Number.isFinite(yConfig.defaultMin)) yMin = Math.min(yMin, yConfig.defaultMin);
  if (Number.isFinite(yConfig.defaultMax)) yMax = Math.max(yMax, yConfig.defaultMax);
  if (!Number.isFinite(xMin) || !Number.isFinite(xMax)) {
    xMin = 0;
    xMax = 1;
  }
  if (!Number.isFinite(yMin) || !Number.isFinite(yMax)) {
    yMin = 0;
    yMax = 1;
  }
  if (xMin === xMax) {
    const delta = Math.abs(xMin) * 0.1 || 1;
    xMin -= delta;
    xMax += delta;
  }
  if (yMin === yMax) {
    const delta = Math.abs(yMin) * 0.1 || 1;
    yMin -= delta;
    yMax += delta;
  }
  const xRange = xMax - xMin || 1;
  const yRange = yMax - yMin || 1;
  const xBinCount = 24;
  const yBinCount = 24;
  const xBins = Array.from({ length: xBinCount }, () => 0);
  const yBins = Array.from({ length: yBinCount }, () => 0);
  const activeThresholds = [];
  if (Number.isFinite(thresholds.x) && xConfig.threshold) {
    activeThresholds.push({ axis: 'x', value: thresholds.x, direction: xConfig.threshold.direction });
  }
  if (Number.isFinite(thresholds.y) && yConfig.threshold) {
    activeThresholds.push({ axis: 'y', value: thresholds.y, direction: yConfig.threshold.direction });
  }
  let passingCount = 0;

  state.scatterLayout = filtered.map((p) => {
    const normalizedX = (p.xValue - xMin) / xRange;
    const normalizedY = (p.yValue - yMin) / yRange;
    const x = margin.left + normalizedX * plotWidth;
    const y = height - margin.bottom - normalizedY * plotHeight;
    const colorKey = p.epitopeLabel ? p.epitopeLabel : '__default__';
    const color = colorMap.get(colorKey) || SCATTER_COLOR_PALETTE[0];
    if (Number.isFinite(normalizedX)) {
      const binX = Math.min(xBinCount - 1, Math.max(0, Math.floor(normalizedX * xBinCount)));
      xBins[binX] += 1;
    }
    if (Number.isFinite(normalizedY)) {
      const binY = Math.min(yBinCount - 1, Math.max(0, Math.floor(normalizedY * yBinCount)));
      yBins[binY] += 1;
    }
    const meetsThresholds = activeThresholds.length === 0 || activeThresholds.every((rule) => {
      const value = rule.axis === 'x' ? p.xValue : p.yValue;
      if (!Number.isFinite(value)) return false;
      if (rule.direction === 'max') return value <= rule.value;
      if (rule.direction === 'min') return value >= rule.value;
      return true;
    });
    if (meetsThresholds) {
      passingCount += 1;
    }
    return {
      ...p,
      x,
      y,
      color,
      xField: axes.x,
      yField: axes.y,
      xMin,
      xMax,
      yMin,
      yMax,
      xRange,
      yRange,
      margin,
    };
  });

  if (activeThresholds.length === 0) {
    passingCount = filtered.length;
  }

  state.scatterStats = {
    xBins,
    yBins,
    total: filtered.length,
    passingCount,
    thresholds,
    fields: { x: axes.x, y: axes.y },
    axis: {
      xMin,
      xMax,
      yMin,
      yMax,
      xRange,
      yRange,
      margin,
    },
  };
}

function formatLegendLabel(key) {
  if (!key || key === '__default__') return 'No epitope tag';
  return key;
}

function updateScatterLegend() {
  const container = el.scatterLegend;
  if (!container) return;
  const colorMap = state.scatterColorMap instanceof Map ? state.scatterColorMap : new Map();
  if (!colorMap || colorMap.size === 0) {
    container.innerHTML = '';
    container.hidden = true;
    return;
  }

  const counts = new Map();
  if (Array.isArray(state.scatterLayout)) {
    state.scatterLayout.forEach((point) => {
      const key = point && point.epitopeLabel ? point.epitopeLabel : '__default__';
      counts.set(key, (counts.get(key) || 0) + 1);
    });
  }

  container.innerHTML = '';
  const title = document.createElement('span');
  title.className = 'scatter-legend-title';
  title.textContent = 'Color by epitope:';
  container.appendChild(title);

  colorMap.forEach((color, key) => {
    const item = document.createElement('div');
    item.className = 'scatter-legend-item';

    const swatch = document.createElement('span');
    swatch.className = 'scatter-legend-swatch';
    swatch.style.backgroundColor = color;

    const label = document.createElement('span');
    label.className = 'scatter-legend-label';
    const count = counts.get(key) || 0;
    const text = formatLegendLabel(key);
    label.textContent = count > 0 ? `${text} (${count})` : text;

    item.appendChild(swatch);
    item.appendChild(label);
    container.appendChild(item);
  });

  container.hidden = false;
}

function renderScatter() {
  if (!el.scatterCanvas) return;
  const ctx = el.scatterCanvas.getContext('2d');
  ctx.clearRect(0, 0, el.scatterCanvas.width, el.scatterCanvas.height);
  renderScatterHistograms(null);
  updateScatterLegend();
  if (state.scatterLayout.length === 0) {
    updateScatterThresholdSummary();
    ctx.fillStyle = '#64748b';
    ctx.font = '13px Inter, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillText('No scatter data available', 20, 24);
    return;
  }

  const sample = state.scatterLayout[0];
  const {
    xMin, xMax, yMin, yMax, xRange, yRange, margin,
  } = sample;
  const width = el.scatterCanvas.width;
  const height = el.scatterCanvas.height;
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const stats = state.scatterStats || {};
  const thresholds = stats.thresholds || {};
  const fields = stats.fields || {};
  const axes = getActiveScatterAxes();
  const xField = fields.x || sample.xField || axes.x;
  const yField = fields.y || sample.yField || axes.y;
  const xConfig = getScatterFieldConfig(xField) || { label: xField };
  const yConfig = getScatterFieldConfig(yField) || { label: yField };
  const xLabel = xConfig.label || xField;
  const yLabel = yConfig.label || yField;
  const formatXTicks = (value) => {
    if (!Number.isFinite(value)) return '';
    const decimals = Number.isInteger(xConfig.tickDecimals) ? xConfig.tickDecimals : 2;
    return value.toFixed(decimals);
  };
  const formatYTicks = (value) => {
    if (!Number.isFinite(value)) return '';
    const decimals = Number.isInteger(yConfig.tickDecimals) ? yConfig.tickDecimals : 2;
    return value.toFixed(decimals);
  };

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = '#cbd5f5';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.rect(margin.left, margin.top, plotWidth, plotHeight);
  ctx.stroke();

  ctx.font = '13px Inter, sans-serif';
  ctx.fillStyle = '#0f172a';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(xLabel, margin.left + plotWidth / 2, height - margin.bottom + 36);
  ctx.save();
  ctx.translate(margin.left - 56, margin.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  ctx.strokeStyle = '#e2e8f0';
  ctx.fillStyle = '#0f172a';
  ctx.font = '12px Inter, sans-serif';
  const xTicks = 6;
  for (let i = 0; i < xTicks; i += 1) {
    const t = xTicks === 1 ? 0 : i / (xTicks - 1);
    const value = xMin + t * xRange;
    const x = margin.left + t * plotWidth;
    ctx.beginPath();
    ctx.moveTo(x, height - margin.bottom);
    ctx.lineTo(x, height - margin.bottom + 6);
    ctx.stroke();
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText(formatXTicks(value), x, height - margin.bottom + 10);
  }

  const yTicks = 6;
  for (let i = 0; i < yTicks; i += 1) {
    const t = yTicks === 1 ? 0 : i / (yTicks - 1);
    const value = yMin + t * yRange;
    const y = height - margin.bottom - t * plotHeight;
    ctx.beginPath();
    ctx.moveTo(margin.left - 6, y);
    ctx.lineTo(margin.left, y);
    ctx.stroke();
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(formatYTicks(value), margin.left - 10, y);
  }

  const hasXThreshold = Number.isFinite(thresholds.x) && xConfig.threshold;
  const hasYThreshold = Number.isFinite(thresholds.y) && yConfig.threshold;
  if (hasXThreshold && thresholds.x >= xMin && thresholds.x <= xMax) {
    const t = (thresholds.x - xMin) / xRange;
    const x = margin.left + t * plotWidth;
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = '#f97316';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(x, margin.top);
    ctx.lineTo(x, height - margin.bottom);
    ctx.stroke();
    ctx.restore();
    ctx.fillStyle = '#f97316';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    const directionSymbol = xConfig.threshold.direction === 'max' ? '≤' : '≥';
    ctx.fillText(
      `${xLabel} ${directionSymbol} ${formatScatterThreshold(xField, thresholds.x)}`,
      Math.min(x + 8, width - margin.right),
      margin.top + 18,
    );
  }

  if (hasYThreshold && thresholds.y >= yMin && thresholds.y <= yMax) {
    const t = (thresholds.y - yMin) / yRange;
    const y = height - margin.bottom - t * plotHeight;
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(margin.left, y);
    ctx.lineTo(width - margin.right, y);
    ctx.stroke();
    ctx.restore();
    ctx.fillStyle = '#0f172a';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    const directionSymbol = yConfig.threshold.direction === 'max' ? '≤' : '≥';
    ctx.fillText(
      `${yLabel} ${directionSymbol} ${formatScatterThreshold(yField, thresholds.y)}`,
      margin.left + 8,
      Math.max(margin.top + 8, y - 18),
    );
  }

  state.scatterLayout.forEach((point) => {
    ctx.beginPath();
    ctx.fillStyle = point.color || '#2563eb';
    ctx.globalAlpha = 0.82;
    ctx.arc(point.x, point.y, 4, 0, Math.PI * 2);
    ctx.fill();
    ctx.globalAlpha = 1;
  });

  renderScatterHistograms(sample);
  updateScatterThresholdSummary();
}

function renderScatterHistograms(sample) {
  const topCanvas = el.scatterHistX;
  const rightCanvas = el.scatterHistY;
  if (topCanvas) {
    const ctx = topCanvas.getContext('2d');
    ctx.clearRect(0, 0, topCanvas.width, topCanvas.height);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, topCanvas.width, topCanvas.height);
  }
  if (rightCanvas) {
    const ctx = rightCanvas.getContext('2d');
    ctx.clearRect(0, 0, rightCanvas.width, rightCanvas.height);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, rightCanvas.width, rightCanvas.height);
  }
  if (!sample || !state.scatterStats || state.scatterLayout.length === 0) return;

  const {
    margin, xRange, yRange, xMin, yMin,
  } = sample;
  const { thresholds = {}, xBins = [], yBins = [], fields = {} } = state.scatterStats;
  const axes = getActiveScatterAxes();
  const xField = fields.x || sample.xField || axes.x;
  const yField = fields.y || sample.yField || axes.y;
  const xConfig = getScatterFieldConfig(xField);
  const yConfig = getScatterFieldConfig(yField);
  const xThreshold = Number.isFinite(thresholds.x) && xConfig && xConfig.threshold ? thresholds.x : null;
  const yThreshold = Number.isFinite(thresholds.y) && yConfig && yConfig.threshold ? thresholds.y : null;

  if (topCanvas) {
    const ctx = topCanvas.getContext('2d');
    const marginTop = { left: margin.left, right: margin.right, top: 16, bottom: 28 };
    const plotWidth = Math.max(1, topCanvas.width - marginTop.left - marginTop.right);
    const plotHeight = Math.max(1, topCanvas.height - marginTop.top - marginTop.bottom);
    const maxBin = xBins.length ? Math.max(...xBins, 1) : 1;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, topCanvas.width, topCanvas.height);
    ctx.strokeStyle = '#cbd5f5';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(marginTop.left, topCanvas.height - marginTop.bottom);
    ctx.lineTo(topCanvas.width - marginTop.right, topCanvas.height - marginTop.bottom);
    ctx.stroke();
    const barWidth = plotWidth / (xBins.length || 1);
    ctx.fillStyle = '#3b82f6';
    xBins.forEach((count, idx) => {
      if (!Number.isFinite(count)) return;
      const heightRatio = maxBin === 0 ? 0 : count / maxBin;
      const barHeight = heightRatio * plotHeight;
      const x = marginTop.left + idx * barWidth;
      const y = topCanvas.height - marginTop.bottom - barHeight;
      ctx.fillRect(x + 0.5, y, Math.max(1, barWidth - 2), barHeight);
    });
    if (Number.isFinite(xThreshold) && xThreshold >= xMin && xThreshold <= xMin + xRange) {
      const t = (xThreshold - xMin) / xRange;
      const x = marginTop.left + t * plotWidth;
      ctx.save();
      ctx.setLineDash([6, 4]);
      ctx.strokeStyle = '#f97316';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(x, marginTop.top);
      ctx.lineTo(x, topCanvas.height - marginTop.bottom);
      ctx.stroke();
      ctx.restore();
    }
  }

  if (rightCanvas) {
    const ctx = rightCanvas.getContext('2d');
    const marginRight = { top: margin.top, bottom: margin.bottom, left: 16, right: 28 };
    const plotWidth = Math.max(1, rightCanvas.width - marginRight.left - marginRight.right);
    const plotHeight = Math.max(1, rightCanvas.height - marginRight.top - marginRight.bottom);
    const maxBin = yBins.length ? Math.max(...yBins, 1) : 1;
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, rightCanvas.width, rightCanvas.height);
    ctx.strokeStyle = '#cbd5f5';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(marginRight.left, marginRight.top);
    ctx.lineTo(marginRight.left, rightCanvas.height - marginRight.bottom);
    ctx.stroke();
    const barHeight = plotHeight / (yBins.length || 1);
    ctx.fillStyle = '#22c55e';
    yBins.forEach((count, idx) => {
      if (!Number.isFinite(count)) return;
      const widthRatio = maxBin === 0 ? 0 : count / maxBin;
      const barWidth = widthRatio * plotWidth;
      const y = marginRight.top + (yBins.length - 1 - idx) * barHeight;
      ctx.fillRect(marginRight.left, y + 0.5, barWidth, Math.max(1, barHeight - 2));
    });
    if (Number.isFinite(yThreshold) && yThreshold >= yMin && yThreshold <= yMin + yRange) {
      const t = (yThreshold - yMin) / yRange;
      const y = rightCanvas.height - marginRight.bottom - t * plotHeight;
      ctx.save();
      ctx.setLineDash([6, 4]);
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(marginRight.left, y);
      ctx.lineTo(rightCanvas.width - marginRight.right, y);
      ctx.stroke();
      ctx.restore();
    }
  }
}

function updateScatterThresholdSummary() {
  if (!el.scatterThresholdCount) return;
  const stats = state.scatterStats;
  if (!stats || !stats.total) {
    el.scatterThresholdCount.textContent = '';
    return;
  }
  const { thresholds = {}, passingCount, total, fields = {} } = stats;
  const axes = getActiveScatterAxes();
  const xField = fields.x || axes.x;
  const yField = fields.y || axes.y;
  const xConfig = getScatterFieldConfig(xField);
  const yConfig = getScatterFieldConfig(yField);
  const clauses = [];
  if (Number.isFinite(thresholds.x) && xConfig && xConfig.threshold) {
    const directionSymbol = xConfig.threshold.direction === 'max' ? '≤' : '≥';
    clauses.push(`${xConfig.label} ${directionSymbol} ${formatScatterThreshold(xField, thresholds.x)}`);
  }
  if (Number.isFinite(thresholds.y) && yConfig && yConfig.threshold) {
    const directionSymbol = yConfig.threshold.direction === 'max' ? '≤' : '≥';
    clauses.push(`${yConfig.label} ${directionSymbol} ${formatScatterThreshold(yField, thresholds.y)}`);
  }
  if (!clauses.length) {
    el.scatterThresholdCount.textContent = `${total} designs plotted.`;
    return;
  }
  const percent = total === 0 ? 0 : ((passingCount / total) * 100);
  const formatted = `${passingCount} / ${total} designs (${percent.toFixed(1)}%) meet ${clauses.join(' and ')}`;
  el.scatterThresholdCount.textContent = formatted;
}

function handleScatterThresholdChange(axis, rawValue) {
  if (!axis) return;
  const axes = getActiveScatterAxes();
  const fieldKey = axis === 'y' ? axes.y : axes.x;
  const config = getScatterFieldConfig(fieldKey);
  if (!config || !config.threshold) {
    updateScatterControls();
    return;
  }
  let value = Number.isFinite(rawValue) ? rawValue : null;
  if (value !== null) {
    value = clampScatterThreshold(fieldKey, value);
  }
  state.scatterThresholds = {
    ...state.scatterThresholds,
    [fieldKey]: value,
  };
  updateScatterThresholdInputs();
  computeScatterLayout(state.scatterPoints || []);
  renderScatter();
}

function analysisSignature(pdbId, runLabel, sourcePath, engineId = 'rfantibody') {
  return [pdbId || '', runLabel || '', sourcePath || '', engineId || ''].join('::');
}

function resetResultsPanel(options = {}) {
  const preserveMeta = Boolean(options.preserveMeta);
  const preserveRunLabel = Boolean(options.preserveRunLabel);
  state.rankings = [];
  state.rankingsResponse = null;
  state.rankingsTotal = 0;
  state.rankingsShowAll = false;
  state.rankingsDisplayLimit = RESULTS_TABLE_DISPLAY_LIMIT;
  state.scatterLayout = [];
  state.scatterPoints = [];
  state.scatterStats = null;
  state.scatterColorMap = new Map();
  state.selectedDesign = null;
  state.galleryAvailable = false;
  state.boltzActiveSpec = null;
  state.librarySummary = null;
  if (!preserveRunLabel && el.resultsRunLabel) {
    el.resultsRunLabel.value = '';
  }
  if (!preserveMeta && el.resultsMeta) {
    setBadge(el.resultsMeta, null);
  }
  if (el.binderDetail) {
    el.binderDetail.hidden = true;
  }
  renderResults();
  applyResultsEngineCapabilities();
}

function resetAnalysisPanel(options = {}) {
  const disableButton = Boolean(options.disableButton);
  const preserveContext = Boolean(options.preserveContext);
  state.analysisResults = null;
  if (!preserveContext) {
    state.analysisContext = null;
  }
  if (el.analysisPlots) el.analysisPlots.innerHTML = '';
  if (el.analysisMatrix) el.analysisMatrix.innerHTML = '';
  if (el.analysisLogs) el.analysisLogs.innerHTML = '';
  if (el.analysisPanel) el.analysisPanel.hidden = true;
  if (el.analysisStatus) setBadge(el.analysisStatus, null);
  if (el.analysisRun) {
    el.analysisRun.textContent = 'Analyze rankings';
    if (disableButton) {
      el.analysisRun.disabled = true;
    } else {
      el.analysisRun.disabled = state.rankings.length === 0 || state.analysisRunning;
    }
  }
}

function applyResultsEngineCapabilities() {
  const engine = state.rankingsResponse?.engine_id || 'rfantibody';
  const isBoltz = engine === 'boltzgen';
  const hasRows = state.rankings.length > 0;
  if (el.analysisRun) {
    el.analysisRun.textContent = 'Analyze rankings';
    el.analysisRun.disabled = isBoltz || !hasRows || state.analysisRunning;
  }
  if (el.exportButton) {
    el.exportButton.disabled = isBoltz || !hasRows;
  }
  if (el.exportOpen) {
    el.exportOpen.disabled = isBoltz || !hasRows;
  }
  if (el.libraryRun) {
    el.libraryRun.disabled = isBoltz || !hasRows;
  }
  if (el.pymolTop) {
    const boltzReady = isBoltz ? Boolean(hasRows && state.rankingsResponse?.source_path) : hasRows;
    el.pymolTop.disabled = !boltzReady;
    if (isBoltz) {
      el.pymolTop.textContent = 'PyMOL BoltzGen bundle';
    } else {
      el.pymolTop.textContent = state.galleryAvailable && hasRows ? 'Launch PyMOL Gallery' : 'PyMOL top 96';
    }
  }
  if (el.pymolMovie) {
    el.pymolMovie.disabled = isBoltz || !state.galleryAvailable;
  }
  if (el.syncResultsBtn) {
    el.syncResultsBtn.hidden = isBoltz;
  }
  if (el.syncBoltzResultsBtn) {
    el.syncBoltzResultsBtn.hidden = !isBoltz;
  }
  if (el.resultsMeta && isBoltz && state.rankingsResponse) {
    el.resultsMeta.hidden = false;
  }
}

function ensureAf3Results(featureDescription) {
  if (state.rankingsResponse?.engine_id === 'boltzgen') {
    showAlert(`${featureDescription} is only available for AF3 assessment runs.`);
    return false;
  }
  return true;
}

function similarityColor(value) {
  if (typeof value !== 'number' || Number.isNaN(value)) {
    return '#e2e8f0';
  }
  const t = Math.min(1, Math.max(0, value));
  const start = [37, 99, 235];
  const end = [220, 38, 38];
  const r = Math.round(start[0] + (end[0] - start[0]) * t);
  const g = Math.round(start[1] + (end[1] - start[1]) * t);
  const b = Math.round(start[2] + (end[2] - start[2]) * t);
  return `rgb(${r}, ${g}, ${b})`;
}

function renderSimilarityMatrix(data) {
  if (!el.analysisMatrix) return false;
  el.analysisMatrix.innerHTML = '';

  if (!data || !Array.isArray(data.designs) || !Array.isArray(data.matrix) || !data.designs.length) {
    const empty = document.createElement('div');
    empty.className = 'matrix-legend';
    empty.textContent = 'Sequence similarity matrix unavailable for this run.';
    el.analysisMatrix.appendChild(empty);
    return false;
  }

  const designs = data.designs;
  const matrix = data.matrix;
  const sequences = Array.isArray(data.sequences) ? data.sequences : [];
  const n = Math.min(designs.length, matrix.length);
  if (!n) {
    const empty = document.createElement('div');
    empty.className = 'matrix-legend';
    empty.textContent = 'Sequence similarity matrix unavailable for this run.';
    el.analysisMatrix.appendChild(empty);
    return false;
  }

  const cellSize = n <= 20 ? 24 : n <= 40 ? 16 : n <= 80 ? 12 : 8;
  const marginLeft = Math.max(120, Math.min(200, cellSize * 5));
  const marginTop = Math.max(110, Math.min(200, cellSize * 5));
  const marginBottom = 60;
  const marginRight = 40;
  const width = marginLeft + marginRight + cellSize * n;
  const height = marginTop + marginBottom + cellSize * n;

  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');

  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, width, height);

  ctx.save();
  ctx.translate(marginLeft, marginTop);
  for (let i = 0; i < n; i += 1) {
    const row = Array.isArray(matrix[i]) ? matrix[i] : [];
    for (let j = 0; j < n; j += 1) {
      const val = typeof row[j] === 'number' ? row[j] : Number.NaN;
      ctx.fillStyle = similarityColor(val);
      ctx.fillRect(j * cellSize, i * cellSize, cellSize, cellSize);
      if (!Number.isNaN(val) && cellSize >= 18) {
        const text = (val * 100).toFixed(0);
        ctx.fillStyle = val >= 0.65 ? '#f8fafc' : '#0f172a';
        ctx.font = `600 ${Math.max(10, Math.min(14, cellSize - 6))}px Inter, sans-serif`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(text, j * cellSize + cellSize / 2, i * cellSize + cellSize / 2);
      }
    }
  }
  ctx.restore();

  ctx.fillStyle = '#0f172a';
  ctx.font = `${Math.max(10, Math.min(12, cellSize - 2))}px Inter, sans-serif`;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (let i = 0; i < n; i += 1) {
    const label = designs[i] || `Design ${i + 1}`;
    const trimmed = label.length > 22 ? `${label.slice(0, 19)}…` : label;
    ctx.fillText(trimmed, marginLeft - 8, marginTop + i * cellSize + cellSize / 2);
  }

  ctx.save();
  ctx.translate(marginLeft, marginTop);
  for (let j = 0; j < n; j += 1) {
    const label = designs[j] || `Design ${j + 1}`;
    const trimmed = label.length > 22 ? `${label.slice(0, 19)}…` : label;
    ctx.save();
    ctx.translate(j * cellSize + cellSize / 2, -6);
    ctx.rotate(-Math.PI / 3);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(trimmed, 0, 0);
    ctx.restore();
  }
  ctx.restore();

  ctx.font = '12px Inter, sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('Binder index', marginLeft + (cellSize * n) / 2, height - marginBottom / 2);
  ctx.save();
  ctx.translate(24, marginTop + (cellSize * n) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('Binder index', 0, 0);
  ctx.restore();

  el.analysisMatrix.appendChild(canvas);

  const legend = document.createElement('div');
  legend.className = 'matrix-legend';
  const metricLabel = data.metric ? data.metric.replace(/_/g, ' ') : 'sequence match ratio';
  legend.textContent = `Sequence similarity across top ${n} binders (${metricLabel}). Values range from 0 (dissimilar) to 1 (identical).`;
  el.analysisMatrix.appendChild(legend);

  if (sequences.length) {
    const details = document.createElement('details');
    details.className = 'matrix-legend';
    const summary = document.createElement('summary');
    summary.textContent = 'Show binder sequences';
    details.appendChild(summary);
    const list = document.createElement('ol');
    list.className = 'analysis-sequence-list';
    sequences.slice(0, n).forEach((seq, idx) => {
      const item = document.createElement('li');
      item.innerHTML = `<strong>${designs[idx] || `Design ${idx + 1}`}</strong>: ${seq}`;
      list.appendChild(item);
    });
    details.appendChild(list);
    el.analysisMatrix.appendChild(details);
  }

  return true;
}

function renderAnalysis() {
  if (!el.analysisPanel) return;
  if (!state.analysisResults) {
    el.analysisPanel.hidden = true;
    if (el.analysisPlots) el.analysisPlots.innerHTML = '';
    if (el.analysisMatrix) el.analysisMatrix.innerHTML = '';
    if (el.analysisLogs) el.analysisLogs.innerHTML = '';
    return;
  }

  const { plots = [], similarity = null, logs = [] } = state.analysisResults;

  if (el.analysisPlots) {
    el.analysisPlots.innerHTML = '';
    if (Array.isArray(plots) && plots.length) {
      plots.forEach((plot) => {
        if (!plot || !plot.image_data) return;
        const figure = document.createElement('figure');
        figure.className = 'analysis-plot-card';
        const img = document.createElement('img');
        img.src = `data:image/png;base64,${plot.image_data}`;
        img.alt = plot.title || plot.name || 'Ranking diagnostic plot';
        const caption = document.createElement('figcaption');
        caption.textContent = plot.title || plot.name || 'Plot';
        figure.appendChild(img);
        figure.appendChild(caption);
        el.analysisPlots.appendChild(figure);
      });
    } else {
      const empty = document.createElement('div');
      empty.className = 'matrix-legend';
      empty.textContent = 'No plots generated yet. Run “Analyze rankings” to build diagnostics.';
      el.analysisPlots.appendChild(empty);
    }
  }

  let matrixVisible = false;
  if (similarity) {
    matrixVisible = renderSimilarityMatrix(similarity);
  } else if (el.analysisMatrix) {
    el.analysisMatrix.innerHTML = '';
    const empty = document.createElement('div');
    empty.className = 'matrix-legend';
    empty.textContent = 'Sequence similarity metrics were not generated (binder sequences missing).';
    el.analysisMatrix.appendChild(empty);
  }

  if (el.analysisLogs) {
    el.analysisLogs.innerHTML = '';
    if (logs && logs.length) {
      const heading = document.createElement('div');
      heading.textContent = 'Recent analysis messages:';
      el.analysisLogs.appendChild(heading);
      const list = document.createElement('ul');
      logs.forEach((line) => {
        const li = document.createElement('li');
        li.textContent = line;
        list.appendChild(li);
      });
      el.analysisLogs.appendChild(list);
    }
  }

  const hasPlots = Array.isArray(plots) && plots.length > 0;
  el.analysisPanel.hidden = !(hasPlots || matrixVisible);
}

async function runRankingsAnalysis() {
  if (!state.currentPdb) {
    showAlert('Initialize and load rankings before running the analysis.');
    return;
  }
  if (!ensureAf3Results('Rankings analysis')) {
    return;
  }
  if (state.analysisRunning) return;
  if (!state.rankings || state.rankings.length === 0) {
    showAlert('Load rankings first.');
    return;
  }

  state.analysisRunning = true;
  if (el.analysisRun) {
    el.analysisRun.disabled = true;
    el.analysisRun.textContent = 'Analyzing…';
  }
  if (el.analysisPanel) {
    el.analysisPanel.hidden = false;
  }
  if (el.analysisStatus) {
    setBadge(el.analysisStatus, 'Generating diagnostics…');
  }
  if (el.analysisPlots) {
    const loading = document.createElement('div');
    loading.className = 'matrix-legend';
    loading.textContent = 'Running plot_rankings.py…';
    el.analysisPlots.innerHTML = '';
    el.analysisPlots.appendChild(loading);
  }
  if (el.analysisMatrix) {
    el.analysisMatrix.innerHTML = '';
  }
  if (el.analysisLogs) {
    el.analysisLogs.innerHTML = '';
  }

  const runLabel = el.resultsRunLabel?.value.trim() || state.activeRunLabel || '';
  const payload = runLabel ? { run_label: runLabel } : {};

  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/rankings/analysis`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const data = await res.json();
    state.analysisResults = data;
    state.analysisContext = {
      pdbId: state.currentPdb,
      runLabel: state.activeRunLabel || '',
      sourcePath: state.rankingsResponse?.source_path || '',
      signature: analysisSignature(
        state.currentPdb,
        state.activeRunLabel || '',
        state.rankingsResponse?.source_path || '',
        state.rankingsResponse?.engine_id || 'rfantibody',
      ),
    };
    renderAnalysis();
    if (el.analysisStatus) {
      setBadge(el.analysisStatus, `Updated ${new Date().toLocaleString()}`, 'rgba(134, 239, 172, 0.25)');
    }
  } catch (err) {
    const message = err && err.message ? err.message : String(err);
    if (el.analysisStatus) {
      setBadge(el.analysisStatus, `Analysis failed: ${message}`, 'rgba(248, 113, 113, 0.25)');
    }
    showAlert(message);
  } finally {
    state.analysisRunning = false;
    if (el.analysisRun) {
      el.analysisRun.textContent = 'Analyze rankings';
      el.analysisRun.disabled = state.rankings.length === 0;
    }
  }
}

function renderResults() {
  sortRankings();
  renderResultsTable();
  const scatterPoints = resolveScatterPoints();
  state.scatterPoints = scatterPoints;
  updateScatterControls();
  computeScatterLayout(scatterPoints);
  renderScatter();
  const hasRows = state.rankings.length > 0;
  if (el.resultsTableWrapper) el.resultsTableWrapper.hidden = !hasRows;
  if (el.plotContainer) el.plotContainer.hidden = state.scatterLayout.length === 0;
  if (el.exportButton) el.exportButton.disabled = !hasRows;
  if (el.exportOpen) el.exportOpen.disabled = !hasRows;
  if (el.libraryRun) el.libraryRun.disabled = !hasRows;
  if (el.pymolTop) {
    const engine = state.rankingsResponse?.engine_id || 'rfantibody';
    const boltzReady = engine === 'boltzgen' ? Boolean(hasRows && state.rankingsResponse?.source_path) : hasRows;
    el.pymolTop.disabled = !boltzReady;
    if (engine === 'boltzgen') {
      el.pymolTop.textContent = 'PyMOL BoltzGen bundle';
    } else {
      el.pymolTop.textContent = state.galleryAvailable && hasRows ? 'Launch PyMOL Gallery' : 'PyMOL top 96';
    }
  }
  if (el.pymolMovie) {
    const canRenderMovie = state.galleryAvailable && hasRows;
    el.pymolMovie.disabled = !canRenderMovie;
    el.pymolMovie.title = canRenderMovie
      ? 'Render headless PyMOL gallery movie (this may take a minute)'
      : 'Load a gallery-enabled run to render a movie';
  }
  if (el.refreshResultsBtn) el.refreshResultsBtn.disabled = false;
  if (el.analysisRun) {
    el.analysisRun.disabled = !hasRows || state.analysisRunning;
  }
  if (!hasRows) {
    resetAnalysisPanel({ disableButton: true });
    state.librarySummary = null;
  } else {
    renderAnalysis();
  }
  renderLibrarySummary(state.librarySummary);
}

async function refreshScatterPlot() {
  if (!el.scatterRefresh) return;
  if (!state.currentPdb) {
    showAlert('Load results before refreshing the scatter plot.');
    return;
  }
  if (state.rankingsFetching) {
    showAlert('Results are still loading. Try again in a moment.');
    return;
  }
  const previousLabel = el.scatterRefresh.textContent || 'Refresh scatter';
  el.scatterRefresh.disabled = true;
  el.scatterRefresh.textContent = 'Refreshing…';
  try {
    const initialPoints = resolveScatterPoints();
    if (!state.rankingsResponse || initialPoints.length === 0) {
      throw new Error('Load rankings before refreshing the scatter plot.');
    }
    const scatterPoints = resolveScatterPoints();
    if (!scatterPoints.length) {
      throw new Error('No scatter metrics are available for this run.');
    }
    state.scatterPoints = scatterPoints;
    computeScatterLayout(scatterPoints);
    renderScatter();
    if (el.plotContainer) el.plotContainer.hidden = state.scatterLayout.length === 0;
  } catch (err) {
    const message = err && err.message ? err.message : String(err);
    showAlert(message);
  } finally {
    el.scatterRefresh.disabled = false;
    el.scatterRefresh.textContent = previousLabel;
  }
}

function renderLibrarySummary(summary) {
  if (!el.librarySummary) return;
  if (!summary) {
    el.librarySummary.innerHTML = '<em>Load rankings and generate the assembly plan to preview fragment layout.</em>';
    if (el.libraryDownloads) {
      el.libraryDownloads.hidden = true;
      el.libraryDownloads.innerHTML = '';
    }
    if (el.libraryAlignment) {
      el.libraryAlignment.hidden = true;
      if (el.libraryAlignmentTable) el.libraryAlignmentTable.innerHTML = '';
    }
    if (el.libraryAssembly) el.libraryAssembly.innerHTML = '';
    if (el.libraryPreview) el.libraryPreview.innerHTML = '';
    state.libraryAlignmentMode = 'aa';
    state.libraryJobId = null;
    return;
  }

  if (summary.job_id) {
    state.libraryJobId = summary.job_id;
  }

  const framework1 = summary.framework1 || {};
  const framework2 = summary.framework2 || {};
  const cdr = summary.cdr || {};
  const bsaI = summary.bsaI || {};
  const flanks = summary.flanks || {};

  const lines = [];
  lines.push(`<div><strong>Framework 1:</strong> ${framework1.length_aa || '—'} aa (${framework1.length_nt || '—'} nt)</div>`);
  lines.push(`<div><strong>Framework 2:</strong> ${framework2.length_aa || '—'} aa (${framework2.length_nt || '—'} nt)</div>`);
  const processed = summary.design_count || 0;
  const targetTop = summary.top_n || processed;
  const alignedAa = cdr.aligned_length_aa || cdr.length_aa || '—';
  const alignedNt = cdr.aligned_length_nt || cdr.length_nt || '—';
  lines.push(`<div><strong>CDR window:</strong> ${cdr.length_aa || '—'} aa (${cdr.length_nt || '—'} nt) · aligned ${alignedAa} aa / ${alignedNt} nt · ${processed} designs (top ${targetTop})</div>`);
  lines.push(`<div><strong>BsaI:</strong> ${bsaI.forward || 'GGTCTC'} / ${bsaI.reverse || 'GAGACC'} · flanks ${flanks.upstream || ''} … ${flanks.downstream || ''}</div>`);
  el.librarySummary.innerHTML = lines.join('');

  renderLibraryDownloads(summary.downloads);
  renderLibraryAlignment(summary.alignment);

  if (el.libraryAssembly) {
    const segments = Array.isArray(summary.example?.segments)
      ? summary.example.segments
      : [];
    if (segments.length) {
      const html = segments
        .map((seg) => {
          const typeClass = seg.type ? seg.type.toLowerCase() : 'framework';
          const label = seg.label ? `<strong>${escapeHtml(seg.label)}:</strong> ` : '';
          return `<span class="gg-segment ${typeClass}">${label}${escapeHtml(seg.sequence || '')}</span>`;
        })
        .join('');
      const finalInsert = summary.example?.final_insert || '';
      const finalLabel = summary.example?.design_name ? ` (${escapeHtml(summary.example.design_name)})` : '';
      const finalHtml = finalInsert
        ? `<div class="gg-final">Final insert${finalLabel}: ${escapeHtml(finalInsert)}</div>`
        : '';
      el.libraryAssembly.innerHTML = html + finalHtml;
    } else {
      el.libraryAssembly.innerHTML = '';
    }
  }

  if (el.libraryPreview) {
    const preview = Array.isArray(summary.preview) ? summary.preview : [];
    if (preview.length) {
      const rows = preview
        .map((row) => `
          <tr>
            <td>${escapeHtml(row.design_name || '')}</td>
            <td>${escapeHtml(row.cdr_aa || '')}</td>
            <td>${(row.cdr_dna || '').length}</td>
            <td>${(row.final_insert || '').length}</td>
            <td class="monospace">${escapeHtml(row.cdr_fragment || '')}</td>
          </tr>`)
        .join('');
      el.libraryPreview.innerHTML = `
        <table>
          <thead>
            <tr>
              <th>Design</th>
              <th>CDR (aa)</th>
              <th>CDR nt</th>
              <th>Final insert nt</th>
              <th>BsaI-flanked CDR fragment</th>
            </tr>
          </thead>
          <tbody>${rows}</tbody>
        </table>`;
    } else {
      el.libraryPreview.innerHTML = '<em>No preview available.</em>';
    }
  }
}

function renderLibraryDownloads(downloads) {
  if (!el.libraryDownloads) return;
  const entries = [];
  if (downloads && typeof downloads === 'object') {
    const order = [
      { key: 'csv', fallback: 'Download CSV' },
      { key: 'aa_fasta', fallback: 'Download AA FASTA' },
      { key: 'dna_fasta', fallback: 'Download DNA FASTA' },
    ];
    for (const item of order) {
      const info = downloads[item.key];
      if (info && info.path) {
        entries.push({
          key: item.key,
          label: info.label || item.fallback,
        });
      }
    }
  }

  if (!entries.length) {
    el.libraryDownloads.hidden = true;
    el.libraryDownloads.innerHTML = '';
    return;
  }

  const buttons = entries
    .map(
      (entry) =>
        `<button type="button" class="download-link" data-download="${entry.key}">${escapeHtml(entry.label)}</button>`,
    )
    .join('');
  el.libraryDownloads.innerHTML = `<div class="download-buttons">${buttons}</div>`;
  el.libraryDownloads.hidden = false;
}

function regionClass(region) {
  if (!region) return 'region-cdr';
  return `region-${String(region).toLowerCase().replace(/[^a-z0-9_-]+/g, '') || 'cdr'}`;
}

function buildAlignmentTable(view) {
  const columns = Array.isArray(view?.columns) ? view.columns : [];
  const rows = Array.isArray(view?.rows) ? view.rows : [];
  if (!columns.length || !rows.length) {
    return '<em>No alignment available.</em>';
  }

  const headerCells = columns
    .map((col) => {
      const refChar = col.reference ? escapeHtml(col.reference) : '&nbsp;';
      return `<th class="gg-col ${regionClass(col.region)}">${refChar}</th>`;
    })
    .join('');
  const positionCells = columns
    .map((col) => {
      const text = col.position !== null && col.position !== undefined ? String(col.position) : '';
      return `<th class="gg-col ${regionClass(col.region)}">${escapeHtml(text)}</th>`;
    })
    .join('');

  const bodyRows = rows
    .map((row) => {
      const seq = typeof row.sequence === 'string' ? row.sequence : '';
      const chars = seq.split('');
      const mutations = Array.isArray(row.mutations) ? row.mutations : [];
      const cells = columns
        .map((col, idx) => {
          const ch = idx < chars.length ? chars[idx] : '-';
          const mutated = mutations[idx] === true;
          const classes = ['gg-cell', regionClass(col.region)];
          if (mutated) classes.push('mutated');
          if (ch === '-') classes.push('gap');
          const displayChar = ch === '-' ? '&middot;' : escapeHtml(ch);
          return `<td class="${classes.join(' ')}" data-region="${escapeHtml(col.region || '')}">${displayChar}</td>`;
        })
        .join('');
      const label = row.is_reference ? `${row.design_name || 'Reference'} (ref)` : row.design_name || '';
      const rowClass = row.is_reference ? 'reference-row' : '';
      return `<tr class="${rowClass}"><th scope="row">${escapeHtml(label)}</th>${cells}</tr>`;
    })
    .join('');

  return `
    <div class="alignment-scroll">
      <table class="gg-alignment-table">
        <thead>
          <tr>
            <th scope="col">Design</th>${headerCells}
          </tr>
          <tr class="position-row">
            <th scope="col">#</th>${positionCells}
          </tr>
        </thead>
        <tbody>${bodyRows}</tbody>
      </table>
    </div>`;
}

function renderLibraryAlignment(alignment) {
  if (!el.libraryAlignment || !el.libraryAlignmentTable) return;
  if (!alignment || typeof alignment !== 'object') {
    el.libraryAlignment.hidden = true;
    el.libraryAlignmentTable.innerHTML = '';
    return;
  }

  const modes = [];
  if (alignment.aa?.columns?.length && alignment.aa?.rows?.length) modes.push('aa');
  if (alignment.dna?.columns?.length && alignment.dna?.rows?.length) modes.push('dna');
  if (!modes.length) {
    el.libraryAlignment.hidden = true;
    el.libraryAlignmentTable.innerHTML = '<em>No alignment available.</em>';
    return;
  }

  if (!modes.includes(state.libraryAlignmentMode)) {
    state.libraryAlignmentMode = modes[0];
  }

  if (el.libraryAlignmentToggle) {
    el.libraryAlignmentToggle.querySelectorAll('button[data-mode]').forEach((btn) => {
      const mode = btn.dataset.mode;
      btn.disabled = !modes.includes(mode);
      btn.classList.toggle('active', mode === state.libraryAlignmentMode);
    });
  }

  if (el.libraryAlignmentLegend) {
    const legendItems = Array.isArray(alignment.legends?.[state.libraryAlignmentMode])
      ? alignment.legends[state.libraryAlignmentMode]
      : [];
    if (legendItems.length) {
      const legendHtml = legendItems
        .map(
          (item) =>
            `<span class="legend-item"><span class="legend-swatch ${regionClass(item.type)}"></span>${escapeHtml(
              item.label || item.type || '',
            )}</span>`,
        )
        .join('');
      el.libraryAlignmentLegend.innerHTML = legendHtml;
      el.libraryAlignmentLegend.hidden = false;
    } else {
      el.libraryAlignmentLegend.innerHTML = '';
      el.libraryAlignmentLegend.hidden = true;
    }
  }

  const view = alignment[state.libraryAlignmentMode];
  el.libraryAlignmentTable.innerHTML = buildAlignmentTable(view);
  el.libraryAlignment.hidden = false;
}

function handleLibraryAlignmentToggle(event) {
  const btn = event.target.closest('button[data-mode]');
  if (!btn) return;
  const mode = btn.dataset.mode;
  if (!mode || btn.disabled || state.libraryAlignmentMode === mode) return;
  state.libraryAlignmentMode = mode;
  if (state.librarySummary?.alignment) {
    renderLibraryAlignment(state.librarySummary.alignment);
  }
}

function handleLibraryDownloadClick(event) {
  const btn = event.target.closest('button[data-download]');
  if (!btn) return;
  const kind = btn.dataset.download;
  if (!kind) return;
  downloadGoldenGateFile(kind);
}

async function downloadGoldenGateFile(kind) {
  if (!state.libraryJobId) {
    showAlert('Golden Gate plan is not available for download.');
    return;
  }
  try {
    const resp = await fetch(`/api/golden-gate/${state.libraryJobId}/download/${kind}`);
    if (!resp.ok) {
      throw new Error(`Download failed (${resp.status})`);
    }
    const blob = await resp.blob();
    const disposition = resp.headers.get('content-disposition') || '';
    let filename = state.librarySummary?.downloads?.[kind]?.filename || '';
    const match = disposition.match(/filename="?([^";]+)"?/i);
    if (match && match[1]) {
      filename = match[1];
    }
    if (!filename) {
      if (kind === 'csv') filename = 'golden_gate.csv';
      else if (kind === 'aa_fasta') filename = 'golden_gate_aa.fasta';
      else if (kind === 'dna_fasta') filename = 'golden_gate_dna.fasta';
      else filename = 'download.dat';
    }
    const url = window.URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = filename;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    window.URL.revokeObjectURL(url);
  } catch (err) {
    console.error(err);
    showAlert(err.message || 'Failed to download file.');
  }
}

function selectDesign(designName) {
  state.selectedDesign = designName;
  document.querySelectorAll('#results-table tbody tr').forEach((row) => {
    row.classList.toggle('selected', row.dataset.design === designName);
  });
  const row = state.rankings.find((r) => r.design_name === designName);
  if (!row) {
    el.binderDetail.hidden = true;
    return;
  }
  const metadata = row.metadata && typeof row.metadata === 'object' ? row.metadata : {};
  const info = [
    `<strong>${row.design_name}</strong>`,
    `ipTM: ${row.iptm !== null && row.iptm !== undefined ? row.iptm.toFixed(3) : 'NA'}`,
    `Binder RMSD: ${row.rmsd_diego !== null && row.rmsd_diego !== undefined ? row.rmsd_diego.toFixed(3) : 'NA'}`,
    `ipSAE_min: ${row.ipsae_min !== null && row.ipsae_min !== undefined ? row.ipsae_min.toFixed(3) : 'NA'}`,
    (() => {
      if (row.hotspot_min_distance !== null && row.hotspot_min_distance !== undefined && Number.isFinite(row.hotspot_min_distance)) {
        return `Hotspot min distance: ${row.hotspot_min_distance.toFixed(3)} Å`;
      }
      return '';
    })(),
    metadata.arm ? `Arm: ${metadata.arm}` : '',
    (() => {
      const epitopeLabel = getEpitopeLabel(row);
      return epitopeLabel ? `Epitope: ${epitopeLabel}` : '';
    })(),
    metadata.pymol_script_path ? `PyMOL script: ${metadata.pymol_script_path}` : '',
  ].filter(Boolean);
  el.binderDetail.innerHTML = info.join('<br>');
  el.binderDetail.hidden = false;
}

function handleTableHeaderClicks() {
  document.querySelectorAll('#results-table thead th').forEach((th) => {
    th.addEventListener('click', () => {
      const key = th.dataset.key;
      if (!key) return;
      if (state.tableSort.key === key) {
        state.tableSort.dir = state.tableSort.dir === 'asc' ? 'desc' : 'asc';
      } else {
        state.tableSort = { key, dir: th.dataset.type === 'number' ? 'desc' : 'asc' };
      }
      renderResults();
    });
  });
}

function stopJobPolling() {
  if (state.jobPoller) {
    clearInterval(state.jobPoller);
    state.jobPoller = null;
  }
}

async function fetchJob(jobId) {
  const res = await fetch(`/api/jobs/${jobId}`);
  if (!res.ok) throw new Error('Unable to retrieve job status');
  return res.json();
}

function startJobPolling(jobId, context, opts = {}) {
  state.currentJobId = jobId;
  state.jobContext = context;
  state.selectedJobId = jobId;
  renderJobList();
  stopJobPolling();
  const poll = async () => {
    try {
      const job = await fetchJob(jobId);
      updateJobUI(job);
    } catch (err) {
      appendLog(`Polling error: ${err.message || err}`);
      stopJobPolling();
      if (context === 'target') el.targetSubmit.disabled = false;
      if (context === 'design') el.designSubmit.disabled = false;
      if (context === 'export') {
        el.exportButton.disabled = false;
        if (el.exportOpen) el.exportOpen.disabled = false;
      }
      if (context === 'library') {
        if (el.libraryRun) el.libraryRun.disabled = false;
        if (el.libraryStatus) setBadge(el.libraryStatus, 'Polling failed', 'rgba(248, 113, 113, 0.25)');
      }
      if (context === 'sync' && el.syncResultsBtn) el.syncResultsBtn.disabled = false;
    }
  };
  poll();
  state.jobPoller = setInterval(poll, 2000);
}

function updateJobUI(job) {
  if (!job) return;
  resetJobLog(job.logs.join('\n'));

  if (job.details && typeof job.details === 'object' && job.details.run_label && el.resultsRunLabel) {
    el.resultsRunLabel.value = job.details.run_label;
    highlightRunChip(job.details.run_label);
  }

  const ctx = state.jobContext;
  if (ctx === 'target') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.targetBadge, job.status === 'running' ? 'Running…' : 'Queued…');
    } else if (job.status === 'success') {
      setBadge(el.targetBadge, 'Finished', 'rgba(134, 239, 172, 0.25)');
      el.targetSubmit.disabled = false;
      el.pymolHotspots.disabled = false;
      stopJobPolling();
      loadAlignment(state.currentPdb);
      el.designSubmit.disabled = false;
      setBadge(el.designStatus, 'Ready for submission');
      el.refreshResultsBtn.disabled = false;
      if (el.assessSubmit) el.assessSubmit.disabled = false;
      refreshRunLabel(true);
      if (state.currentPdb) {
        fetchRunHistory(state.currentPdb);
        fetchTargetStatus(state.currentPdb, { silent: true });
      }
      loadPresets();
      scheduleRankingsPolling();
    } else {
      setBadge(el.targetBadge, 'Failed', 'rgba(248, 113, 113, 0.25)');
      showAlert(job.message || 'Target workflow failed.');
      el.targetSubmit.disabled = false;
      stopJobPolling();
    }
  } else if (ctx === 'design') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.designStatus, job.status === 'running' ? 'Running…' : 'Queued…');
      el.designSubmit.disabled = true;
    } else if (job.status === 'success') {
      setBadge(el.designStatus, 'Cluster jobs submitted', 'rgba(134, 239, 172, 0.25)');
      el.designSubmit.disabled = false;
      el.refreshResultsBtn.disabled = false;
      if (el.targetSubmit) el.targetSubmit.disabled = false;
      if (el.assessSubmit) el.assessSubmit.disabled = false;
      const autoRunLabel = job.details?.assessment_run_label || '';
      if (autoRunLabel) {
        state.activeRunLabel = autoRunLabel;
        if (el.resultsRunLabel) el.resultsRunLabel.value = autoRunLabel;
        updateActiveRunDisplay();
        highlightRunChip(autoRunLabel);
      }
      const shouldAutoSync = Boolean(job.details?.run_assess_requested) && autoRunLabel;
      if (shouldAutoSync) {
        syncResultsFromCluster({
          runLabel: autoRunLabel,
          disableButton: false,
          useInputRunLabel: false,
        }).catch((err) => {
          console.error('Automatic cluster sync failed', err);
        });
      }
      stopJobPolling();
      refreshRunLabel(true);
      if (state.currentPdb) {
        fetchRunHistory(state.currentPdb);
      }
      scheduleRankingsPolling();
    } else {
      setBadge(el.designStatus, 'Failed', 'rgba(248, 113, 113, 0.25)');
      showAlert(job.message || 'Design pipeline failed.');
      el.designSubmit.disabled = false;
      if (el.targetSubmit) el.targetSubmit.disabled = false;
      if (el.assessSubmit) el.assessSubmit.disabled = false;
      stopJobPolling();
    }
  } else if (ctx === 'export') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.resultsMeta, 'Exporting…');
    } else if (job.status === 'success') {
      const outDir = job.details?.out_dir ? job.details.out_dir : '';
      const suffix = outDir ? ` → ${outDir}` : '';
      setBadge(el.resultsMeta, `Export complete${suffix}`, 'rgba(134, 239, 172, 0.25)');
      el.exportButton.disabled = false;
      if (el.exportOpen) el.exportOpen.disabled = false;
      stopJobPolling();
      if (outDir) {
        showAlert(`Exports saved in ${outDir}`, false);
      } else {
        showAlert('Exports finished successfully.', false);
      }
    } else {
      showAlert(job.message || 'Export failed.');
      el.exportButton.disabled = false;
      if (el.exportOpen) el.exportOpen.disabled = false;
      stopJobPolling();
    }
  } else if (ctx === 'library') {
    if (job.status === 'running' || job.status === 'pending') {
      if (el.libraryStatus) setBadge(el.libraryStatus, job.status === 'running' ? 'Planning…' : 'Queued…');
      if (el.libraryRun) el.libraryRun.disabled = true;
    } else if (job.status === 'success') {
      const summary = job.details?.summary || null;
      state.librarySummary = summary;
      state.libraryJobId = job.job_id;
      state.libraryAlignmentMode = 'aa';
      renderLibrarySummary(summary);
      if (el.libraryStatus) setBadge(el.libraryStatus, 'Golden Gate plan ready', 'rgba(134, 239, 172, 0.25)');
      if (el.libraryRun) el.libraryRun.disabled = false;
      stopJobPolling();
      const downloadCount = summary?.downloads
        ? Object.values(summary.downloads).filter((info) => info && info.path).length
        : 0;
      if (downloadCount > 0) {
        showAlert(`Golden Gate plan ready — ${downloadCount} download${downloadCount === 1 ? '' : 's'} available.`, false);
      } else {
        showAlert('Golden Gate plan finished.', false);
      }
    } else {
      if (el.libraryStatus) setBadge(el.libraryStatus, 'Planning failed', 'rgba(248, 113, 113, 0.25)');
      if (el.libraryRun) el.libraryRun.disabled = false;
      showAlert(job.message || 'Golden Gate planning failed.');
      state.librarySummary = null;
      state.libraryJobId = null;
      renderLibrarySummary(null);
      stopJobPolling();
    }
  } else if (ctx === 'sync') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.resultsMeta, job.status === 'running' ? 'Syncing…' : 'Queued…');
    } else if (job.status === 'success') {
      const skipped = Boolean(job.details?.skipped);
      const skippedReason = (job.details?.skipped_reason || '').trim();
      const skippedCode = (job.details?.skipped_code || '').trim();
      const localPath = job.details?.local_path ? String(job.details.local_path) : '';
      const baseStatus = skipped
        ? skippedReason || 'Assessments already synced'
        : 'Assessments synced';
      const suffix = localPath ? ` → ${localPath}` : '';
      setBadge(el.resultsMeta, `${baseStatus}${suffix}`, 'rgba(134, 239, 172, 0.25)');
      if (el.syncResultsBtn) el.syncResultsBtn.disabled = false;
      stopJobPolling();
      let alertMsg = localPath ? `${baseStatus}. Saved to ${localPath}` : baseStatus;
      let scheduledRetry = false;
      if (skipped && skippedCode === 'remote-missing') {
        state.lastClusterSyncAt = 0;
        if (state.clusterRetryTimer) {
          clearTimeout(state.clusterRetryTimer);
        }
        const retryRunLabel = job.details?.run_label || state.activeRunLabel || '';
        if (retryRunLabel) {
          state.clusterRetryTimer = setTimeout(() => {
            state.clusterRetryTimer = null;
            syncResultsFromCluster({
              runLabel: retryRunLabel,
              disableButton: false,
              useInputRunLabel: false,
              silent: true,
              force: true,
            }).catch((err) => {
              console.error('Retry cluster sync failed', err);
            });
          }, 60000);
          scheduledRetry = true;
        }
      }
      if (scheduledRetry) {
        alertMsg = `${alertMsg} Retrying shortly…`;
      }
      showAlert(alertMsg, false);
      if (state.currentPdb) {
        fetchRunHistory(state.currentPdb);
      }
    } else {
      if (el.syncResultsBtn) el.syncResultsBtn.disabled = false;
      showAlert(job.message || 'Cluster sync failed.');
      stopJobPolling();
    }
  } else if (ctx === 'assess') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.resultsMeta, 'Submitting assessment…');
      if (el.assessSubmit) el.assessSubmit.disabled = true;
    } else if (job.status === 'success') {
      setBadge(el.resultsMeta, 'Assessment submitted to cluster', 'rgba(134, 239, 172, 0.25)');
      if (el.assessSubmit) el.assessSubmit.disabled = false;
      stopJobPolling();
      if (state.currentPdb) {
        fetchRunHistory(state.currentPdb);
        const syncLabel = job.details?.run_label
          || el.resultsRunLabel?.value.trim()
          || state.activeRunLabel;
        syncResultsFromCluster({ runLabel: syncLabel, silent: true });
      }
    } else {
      setBadge(el.resultsMeta, 'Assessment failed', 'rgba(248, 113, 113, 0.25)');
      showAlert(job.message || 'Assessment submission failed.');
      if (el.assessSubmit) el.assessSubmit.disabled = false;
      stopJobPolling();
    }
  } else {
    const extra = job.message ? ` — ${job.message}` : '';
    const summary = `${job.label} (${job.kind}) ${job.status}${extra}`;
    const isError = job.status !== 'success' && job.status !== 'pending' && job.status !== 'running';
    showAlert(summary, isError);
  }
  fetchJobList();
}

async function loadAlignment(pdbId) {
  if (!pdbId) return;
  setBadge(el.alignmentMeta, `Loading alignment for ${pdbId}…`);
  el.alignmentMeta.hidden = false;
  el.alignmentView.hidden = true;
  try {
    const res = await fetch(`/api/targets/${pdbId}/alignment`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    if (!payload.chain_results || payload.chain_results.length === 0) {
      el.alignmentView.innerHTML = '<div>No alignment results available.</div>';
    } else {
      const fragments = payload.chain_results.map((result) => {
        const header = `Chains ${result.chain_ids.join('+')} · identity ${(result.identity * 100).toFixed(1)}% · coverage ${(result.coverage * 100).toFixed(1)}% · mismatches ${result.mismatches}`;
        const leftLen = result.left_unaligned_length || 0;
        const rightLen = result.right_unaligned_length || 0;
        const leftPreviewRaw = result.left_unaligned_preview || '';
        const rightPreviewRaw = result.right_unaligned_preview || '';
        const leftPreview = leftLen > leftPreviewRaw.length && leftPreviewRaw ? `…${leftPreviewRaw}` : leftPreviewRaw;
        const rightPreview = rightLen > rightPreviewRaw.length && rightPreviewRaw ? `${rightPreviewRaw}…` : rightPreviewRaw;
        const gapBits = [];
        if (leftLen) {
          gapBits.push(`Left gap ${leftLen} aa${leftPreview ? ` (${leftPreview})` : ''}`);
        }
        if (rightLen) {
          gapBits.push(`Right gap ${rightLen} aa${rightPreview ? ` (${rightPreview})` : ''}`);
        }
        const gapLine = gapBits.length ? `<div class="alignment-gap">${gapBits.join(' · ')}</div>` : '';
        const vendorLine = result.aligned_vendor
          .split('')
          .map((v, idx) => {
            const t = result.aligned_target[idx];
            const cls = v === t ? 'match' : 'mismatch';
            return `<span class="${cls}">${v}</span>`;
          })
          .join('');
        const targetLine = result.aligned_target
          .split('')
          .map((t, idx) => {
            const v = result.aligned_vendor[idx];
            const cls = v === t ? 'match' : 'mismatch';
            return `<span class="${cls}">${t}</span>`;
          })
          .join('');
        return `
          <div style="margin-bottom: 18px;">
            <div style="margin-bottom: 8px; font-size: 0.9rem; opacity: 0.85;">${header}</div>
            ${gapLine}
            <div class="monospace">Vendor&nbsp;&nbsp;&nbsp;| ${vendorLine}</div>
            <div class="monospace">Target&nbsp;&nbsp;&nbsp;| ${targetLine}</div>
          </div>`;
      });
      el.alignmentView.innerHTML = fragments.join('\n');
    }
    setBadge(el.alignmentMeta, `Vendor length ${payload.vendor_sequence_length} · ${payload.chain_results.length} solutions`);
    el.alignmentView.hidden = false;
  } catch (err) {
    setBadge(el.alignmentMeta, `Alignment failed: ${err.message || err}`, 'rgba(248, 113, 113, 0.2)');
  }
}

async function launchHotspots() {
  if (!state.currentPdb) {
    showAlert('Initialize target before launching PyMOL.');
    return;
  }
  if (!ensureAf3Results('PyMOL hotspot export')) {
    return;
  }
  el.pymolHotspots.disabled = true;
  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/pymol/hotspots`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ launch: true, bundle_only: false }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    if (payload.launched) {
      showAlert('Launching PyMOL hotspots viewer…', false);
    } else {
      showAlert(`Hotspot bundle ready at ${payload.bundle_path || 'temp directory'}`, false);
    }
  } catch (err) {
    showAlert(err.message || String(err));
  } finally {
    el.pymolHotspots.disabled = false;
  }
}

async function launchTopBinders() {
  if (!state.currentPdb || !state.rankingsResponse) {
    showAlert('Load rankings before launching PyMOL.');
    return;
  }
  const engineId = state.rankingsResponse.engine_id || 'rfantibody';
  if (engineId !== 'boltzgen' && !ensureAf3Results('PyMOL launcher')) {
    return;
  }
  el.pymolTop.disabled = true;
  try {
    const topN = Number(el.exportTopN?.value || 0) || 96;
    const requestPayload = {
      run_label: state.rankingsResponse.run_label,
      top_n: topN,
      launch: true,
      bundle_only: false,
      engine_id: engineId,
    };
    if (engineId === 'boltzgen') {
      const specName = state.boltzActiveSpec || el.boltzSpec?.value.trim() || '';
      if (specName) {
        requestPayload.spec = specName;
      }
    }
    const res = await fetch(`/api/targets/${state.currentPdb}/pymol/top-binders`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(requestPayload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const responsePayload = await res.json();
    showAlert(`PyMOL aggregate script: ${responsePayload.bundle_path || 'generated'}`, false);
  } catch (err) {
    showAlert(err.message || String(err));
  } finally {
    el.pymolTop.disabled = false;
  }
}

async function renderPymolMovie() {
  if (!state.currentPdb || !state.rankingsResponse) {
    showAlert('Load rankings before rendering a movie.');
    return;
  }
  if (!ensureAf3Results('PyMOL movie rendering')) {
    return;
  }
  if (!state.galleryAvailable) {
    showAlert('Gallery movie requires a run with PyMOL gallery outputs.');
    return;
  }

  const button = el.pymolMovie;
  const originalLabel = button ? button.textContent : '';
  if (button) {
    button.disabled = true;
    button.textContent = 'Rendering movie…';
  }

  showAlert('Rendering PyMOL gallery movie… this may take a minute.', false);
  appendLog('Rendering headless PyMOL gallery movie…');

  const topN = Number(el.exportTopN?.value || 0) || 96;
  const payload = {
    run_label: state.rankingsResponse.run_label || null,
    top_n: topN,
    fps: 10,
    interval_sec: 2.0,
    rotation_deg_per_sec: 30,
    rotation_axis: 'y',
    desired_states: 48,
  };

  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/pymol/gallery-movie`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const body = await res.json();
    const moviePath = body.movie_path || 'cache directory';
    showAlert(`PyMOL movie saved to ${moviePath}`, false);
    if (body.frames_pattern) {
      appendLog(`PyMOL movie frames: ${body.frames_pattern}`);
    }
    if (body.script_path) {
      appendLog(`PyMOL movie script: ${body.script_path}`);
    }
    if (body.log_path) {
      appendLog(`PyMOL movie log: ${body.log_path}`);
    }
  } catch (err) {
    showAlert(err.message || String(err));
    appendLog(`PyMOL movie render failed: ${err.message || err}`);
  } finally {
    if (button) {
      button.disabled = false;
      button.textContent = originalLabel || 'Render PyMOL movie';
    }
  }
}

async function syncResultsFromCluster(options = {}) {
  const {
    runLabel: explicitRunLabel,
    silent = false,
    disableButton = true,
    useInputRunLabel = true,
    force = false,
    allAssessments = false,
  } = options;
  if (!state.currentPdb) {
    if (!silent) showAlert('Initialize target first.');
    return;
  }
  if (state.clusterRetryTimer) {
    clearTimeout(state.clusterRetryTimer);
    state.clusterRetryTimer = null;
  }
  let runLabelInput = '';
  if (explicitRunLabel !== undefined) {
    runLabelInput = (explicitRunLabel || '').trim();
  } else if (useInputRunLabel && el.resultsRunLabel) {
    runLabelInput = el.resultsRunLabel.value.trim();
  }
  if (!allAssessments && !runLabelInput && state.activeRunLabel) {
    runLabelInput = state.activeRunLabel.trim();
  }

  const now = Date.now();
  if (!force && state.lastClusterSyncAt) {
    const elapsed = now - state.lastClusterSyncAt;
    if (elapsed < SYNC_RESULTS_MIN_INTERVAL_MS) {
      const remainingMs = SYNC_RESULTS_MIN_INTERVAL_MS - elapsed;
      if (!silent) {
        const remainingSeconds = Math.ceil(remainingMs / 1000);
        showAlert(`Cluster sync is throttled. Try again in ${remainingSeconds}s.`);
      }
      if (state.pendingClusterSyncTimer) {
        clearTimeout(state.pendingClusterSyncTimer);
      }
      state.pendingClusterSync = { ...options, force: true };
      state.pendingClusterSyncTimer = setTimeout(() => {
        state.pendingClusterSyncTimer = null;
        const pending = state.pendingClusterSync;
        state.pendingClusterSync = null;
        if (pending) {
          syncResultsFromCluster(pending);
        }
      }, remainingMs);
      return null;
    }
  }
  const params = runLabelInput ? `?run_label=${encodeURIComponent(runLabelInput)}` : '';
  if (state.pendingClusterSyncTimer) {
    clearTimeout(state.pendingClusterSyncTimer);
    state.pendingClusterSyncTimer = null;
  }
  state.pendingClusterSync = null;
  state.lastClusterSyncAt = now;
  if (disableButton && el.syncResultsBtn) {
    el.syncResultsBtn.disabled = true;
  }
  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/sync${params}`, {
      method: 'POST',
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    if (!silent) {
      const desc = [payload.message || 'Sync queued.'];
      if (payload.run_label) desc.push(`run ${payload.run_label}`);
      showAlert(desc.join(' · '), false);
    }
    if (payload.job_id) {
      startJobPolling(payload.job_id, 'sync');
    } else if (disableButton && el.syncResultsBtn) {
      el.syncResultsBtn.disabled = false;
    }
    return payload;
  } catch (err) {
    if (!silent) {
      showAlert(err.message || String(err));
    }
    if (disableButton && el.syncResultsBtn) {
      el.syncResultsBtn.disabled = false;
    }
    throw err;
  }
}

async function syncBoltzResultsFromCluster(options = {}) {
  const { runLabel: explicitRunLabel, silent = false } = options;
  if (!state.currentPdb) {
    if (!silent) showAlert('Initialize target first.');
    return;
  }
  const source = getResultsSource();
  if (source !== 'boltzgen') {
    if (!silent) showAlert('Switch the results source to BoltzGen to sync these outputs.');
    return;
  }
  let runLabel = '';
  if (explicitRunLabel !== undefined) {
    runLabel = (explicitRunLabel || '').trim();
  } else if (el.resultsRunLabel) {
    runLabel = el.resultsRunLabel.value.trim();
  }
  const params = new URLSearchParams();
  if (runLabel) params.append('run_label', runLabel);
  if (el.syncBoltzResultsBtn) {
    el.syncBoltzResultsBtn.disabled = true;
  }
  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/boltzgen/sync${params.size ? `?${params}` : ''}`, {
      method: 'POST',
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    if (!silent) {
      showAlert(payload.message || 'BoltzGen results synced.', false);
    }
    fetchBoltzRunHistory(state.currentPdb);
    return payload;
  } catch (err) {
    if (!silent) {
      showAlert(err.message || String(err));
    }
    throw err;
  } finally {
    if (el.syncBoltzResultsBtn) {
      el.syncBoltzResultsBtn.disabled = false;
    }
  }
}

function registerScatterClick() {
  if (!el.scatterCanvas) return;
  el.scatterCanvas.addEventListener('click', (event) => {
    if (state.scatterLayout.length === 0) return;
    event.preventDefault();
  });
}

function disableFutureSections() {
  stopRankingsPolling();
  el.designSubmit.disabled = true;
  setBadge(el.designStatus, 'Awaiting target prep');
  if (el.refreshResultsBtn) el.refreshResultsBtn.disabled = true;
  el.exportButton.disabled = true;
  if (el.exportOpen) el.exportOpen.disabled = true;
  el.pymolTop.disabled = true;
  el.pymolHotspots.disabled = true;
  if (el.pymolMovie) el.pymolMovie.disabled = true;
  if (el.assessSubmit) el.assessSubmit.disabled = true;
  if (el.libraryRun) el.libraryRun.disabled = true;
  state.librarySummary = null;
  renderLibrarySummary(null);
  state.targetStatus = null;
  state.targetDetails = null;
  renderTargetInsights({});
  renderCurrentSelection();
}

function initEventHandlers() {
  if (el.debugToggle) {
    setDebugMode(el.debugToggle.checked);
    el.debugToggle.addEventListener('change', (event) => {
      setDebugMode(event.target.checked);
    });
  }
  if (el.designEngineSelect) {
    el.designEngineSelect.addEventListener('change', (event) => {
      const value = (event.target.value || '').trim();
      state.selectedDesignEngine = value || 'rfantibody';
      updateDesignEngineHint();
      applyDesignEngineFieldConfig(getDesignEngineMeta(state.selectedDesignEngine));
    });
  }
  if (el.presetRefresh) {
    el.presetRefresh.addEventListener('click', () => loadPresets());
  }
  el.targetForm.addEventListener('submit', queueTargetInit);
  el.refreshAlignmentBtn.addEventListener('click', () => {
    if (state.currentPdb) loadAlignment(state.currentPdb);
  });
  el.designSubmit.addEventListener('click', queueDesignRun);
  if (el.exportOpen) el.exportOpen.addEventListener('click', () => {
    if (el.exportOpen.disabled) return;
    openExportModal();
  });
  if (el.exportClose) el.exportClose.addEventListener('click', () => closeExportModal());
  if (el.exportBackdrop) el.exportBackdrop.addEventListener('click', () => closeExportModal());
  el.exportButton.addEventListener('click', runExport);
  if (el.catalogOpen) {
    el.catalogOpen.addEventListener('click', () => openCatalogModal());
  }
  if (el.catalogClose) {
    el.catalogClose.addEventListener('click', () => closeCatalogModal());
  }
  if (el.catalogBackdrop) {
    el.catalogBackdrop.addEventListener('click', () => closeCatalogModal());
  }
  if (el.catalogRefresh) {
    el.catalogRefresh.addEventListener('click', () => loadCatalogFiles());
  }
  if (el.catalogLimit) {
    el.catalogLimit.addEventListener('change', () => {
      if (state.catalogSelectedFile) {
        selectCatalogFile(state.catalogSelectedFile);
      }
    });
  }
  if (el.catalogDownload) {
    el.catalogDownload.addEventListener('click', () => downloadCatalogFile());
  }
  if (el.libraryRun) el.libraryRun.addEventListener('click', runGoldenGatePlanner);
  if (el.counterSelectionMethod) {
    el.counterSelectionMethod.addEventListener('change', (event) => {
      state.counterSelectionMethod = event.target.value;
      updateCounterSelectionPlan();
    });
  }
  if (el.refreshResultsBtn) {
    el.refreshResultsBtn.addEventListener('click', () => fetchRankings({ silent: false }));
  }
  if (el.analysisRun) {
    el.analysisRun.addEventListener('click', () => runRankingsAnalysis());
  }
  el.pymolHotspots.addEventListener('click', launchHotspots);
  el.pymolTop.addEventListener('click', launchTopBinders);
  if (el.pymolMovie) {
    el.pymolMovie.addEventListener('click', renderPymolMovie);
  }
  if (el.syncResultsBtn) {
    el.syncResultsBtn.addEventListener('click', () =>
      syncResultsFromCluster({ useInputRunLabel: false, allAssessments: true })
    );
  }
  if (el.assessSubmit) {
    el.assessSubmit.addEventListener('click', queueAssessmentRun);
  }
  if (el.libraryDownloads) {
    el.libraryDownloads.addEventListener('click', handleLibraryDownloadClick);
  }
  if (el.libraryAlignmentToggle) {
    el.libraryAlignmentToggle.addEventListener('click', handleLibraryAlignmentToggle);
  }
  if (el.resultsRunLabel) {
    el.resultsRunLabel.addEventListener('input', (event) => {
      state.activeRunLabel = event.target.value.trim();
      updateActiveRunDisplay();
      highlightRunChip(event.target.value || '');
    });
  }
  if (el.resultsSource) {
    setResultsSource(el.resultsSource.value || state.resultsSource);
    el.resultsSource.addEventListener('change', (event) => {
      setResultsSource(event.target.value || 'af3');
      renderBoltzRunHistory(state.currentPdb);
      stopRankingsPolling();
    });
  } else {
    setResultsSource(state.resultsSource);
  }
  if (el.syncBoltzResultsBtn) {
    el.syncBoltzResultsBtn.addEventListener('click', () => syncBoltzResultsFromCluster({}));
  }
  // resultsLimit is read when reload is requested.
  if (el.resultsTruncateToggle) {
    el.resultsTruncateToggle.addEventListener('click', () => {
      const total = state.rankingsTotal || state.rankings.length;
      const limit = state.rankingsDisplayLimit || RESULTS_TABLE_DISPLAY_LIMIT;
      if (!total || total <= limit) return;
      state.rankingsShowAll = !state.rankingsShowAll;
      renderResultsTable();
    });
  }
  if (el.pdbInput) {
    el.pdbInput.addEventListener('blur', () => {
      const pdb = el.pdbInput.value.trim();
      if (pdb.length === 4) {
        setCurrentPdb(pdb);
      }
    });
    el.pdbInput.addEventListener('input', () => {
      const pdb = el.pdbInput.value.trim();
      if (pdb.length === 4) {
        syncActivePreset();
      }
    });
  }
  if (el.antigenInput) {
    el.antigenInput.addEventListener('blur', () => syncActivePreset());
    el.antigenInput.addEventListener('input', () => syncActivePreset());
  }
  if (el.refreshJobsBtn) {
    el.refreshJobsBtn.addEventListener('click', () => {
      fetchJobList();
    });
  }
  if (el.scatterXAxis) {
    el.scatterXAxis.addEventListener('change', (event) => {
      handleScatterAxisChange('x', event.target.value);
    });
  }
  if (el.scatterYAxis) {
    el.scatterYAxis.addEventListener('change', (event) => {
      handleScatterAxisChange('y', event.target.value);
    });
  }
  if (el.scatterThresholdX) {
    const handler = (event) => {
      const text = event.target.value != null ? String(event.target.value).trim() : '';
      if (text === '') {
        handleScatterThresholdChange('x', Number.NaN);
        return;
      }
      const value = Number.parseFloat(text);
      if (!Number.isNaN(value)) {
        handleScatterThresholdChange('x', value);
      }
    };
    el.scatterThresholdX.addEventListener('change', handler);
    el.scatterThresholdX.addEventListener('input', handler);
  }
  if (el.scatterThresholdY) {
    const handler = (event) => {
      const text = event.target.value != null ? String(event.target.value).trim() : '';
      if (text === '') {
        handleScatterThresholdChange('y', Number.NaN);
        return;
      }
      const value = Number.parseFloat(text);
      if (!Number.isNaN(value)) {
        handleScatterThresholdChange('y', value);
      }
    };
    el.scatterThresholdY.addEventListener('change', handler);
    el.scatterThresholdY.addEventListener('input', handler);
  }
  if (el.scatterRefresh) {
    el.scatterRefresh.addEventListener('click', () => {
      refreshScatterPlot();
    });
  }
  registerScatterClick();
  handleTableHeaderClicks();
  disableFutureSections();
}

function init() {
  initEventHandlers();
  renderDesignEngineOptions();
  loadDesignEngines();
  updateScatterControls();
  updateActiveRunDisplay();
  renderPresets();
  loadPresets();
  refreshRunLabel(true);
  updateClusterStatus();
  setInterval(updateClusterStatus, 15000);
  fetchJobList();
  if (!state.jobRefreshTimer) {
    state.jobRefreshTimer = setInterval(fetchJobList, 15000);
  }
}

document.addEventListener('DOMContentLoaded', init);
