const state = {
  currentPdb: null,
  jobPoller: null,
  jobContext: null,
  rankings: [],
  rankingsResponse: null,
  tableSort: { key: 'iptm', dir: 'desc' },
  scatterLayout: [],
  scatterPoints: [],
  scatterThresholds: { rmsd: 3.7, iptm: 0.8 },
  scatterStats: null,
  selectedDesign: null,
  jobs: [],
  selectedJobId: null,
  runHistory: {},
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
};

const SYNC_RESULTS_MIN_INTERVAL_MS = 600000;

function escapeHtml(value) {
  return String(value ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
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
  designTotal: document.querySelector('#design-total'),
  designNumSeq: document.querySelector('#design-num-seq'),
  designTemp: document.querySelector('#design-temp'),
  designRunLabel: document.querySelector('#design-run-label'),
  designBinderChain: document.querySelector('#design-binder-chain'),
  designAf3Seed: document.querySelector('#design-af3-seed'),
  refreshResultsBtn: document.querySelector('#refresh-results'),
  resultsMeta: document.querySelector('#results-meta'),
  binderDetail: document.querySelector('#binder-detail'),
  resultsTableWrapper: document.querySelector('#results-table-wrapper'),
  resultsTable: document.querySelector('#results-table'),
  scatterCanvas: document.querySelector('#scatter-plot'),
  scatterHistX: document.querySelector('#scatter-hist-x'),
  scatterHistY: document.querySelector('#scatter-hist-y'),
  scatterThresholdRmsd: document.querySelector('#scatter-threshold-rmsd'),
  scatterThresholdIptm: document.querySelector('#scatter-threshold-iptm'),
  scatterThresholdCount: document.querySelector('#scatter-threshold-count'),
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
  resultsLimit: document.querySelector('#results-limit'),
  pymolTop: document.querySelector('#pymol-top'),
  pymolMovie: document.querySelector('#pymol-movie'),
  syncResultsBtn: document.querySelector('#sync-results'),
  assessSubmit: document.querySelector('#assess-submit'),
  jobList: document.querySelector('#job-list'),
  refreshJobsBtn: document.querySelector('#refresh-jobs'),
  runHistory: document.querySelector('#run-history'),
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
};

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
      const title = document.createElement('span');
      title.className = 'insight-title';
      const chainId = chain.id || chain.chain_id || chain.chain || '—';
      title.textContent = `Chain ${chainId}`;
      li.appendChild(title);
      const descBits = [];
      if (chain.description) descBits.push(chain.description);
      if (chain.length) descBits.push(`${chain.length} aa`);
      if (chain.role) descBits.push(chain.role);
      const detail = document.createElement('div');
      detail.textContent = descBits.length ? descBits.join(' · ') : 'No additional description available.';
      li.appendChild(detail);
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
  fetchRankings({ silent: true }).catch(() => {});
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
  if (!state.currentPdb) return;
  state.rankingsPollTimer = setInterval(() => {
    if (!state.currentPdb || state.rankingsFetching) return;
    fetchRankings({ silent: true }).catch(() => {});
  }, 8000);
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
  if (options.updateInput !== false && el.pdbInput) {
    el.pdbInput.value = upper;
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
  fetchRunHistory(upper);
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
      setClusterStatus('error', `Cluster: not connected - ${msg}. Run ssh -o ControlPath=/Users/inagakit/.ssh/cm-initbinder-hpc3.rcic.uci.edu -o ControlMaster=auto -o ControlPersist=yes -MNf hpc3.rcic.uci.edu`);
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
      btn.addEventListener('click', async () => {
        state.activeRunLabel = run.run_label || '';
        updateActiveRunDisplay();
        if (el.resultsRunLabel) {
          el.resultsRunLabel.value = run.run_label || '';
        }
        highlightRunChip(run.run_label);
        try {
          if (run.available_local) {
            await fetchRankings();
          } else if (run.available_remote) {
            await syncResultsFromCluster({ runLabel: run.run_label, disableButton: false });
            await fetchRankings();
          } else {
            showAlert('No rankings available yet for this run.');
          }
          scheduleRankingsPolling();
        } catch (err) {
          console.warn('[run-history] action failed', err);
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
    total_designs: Number(el.designTotal.value) || 90,
    num_sequences: Number(el.designNumSeq.value) || 1,
    temperature: Number(el.designTemp.value) || 0.1,
    binder_chain_id: el.designBinderChain?.value.trim() || null,
    af3_seed: (() => {
      const raw = Number(el.designAf3Seed?.value ?? 1);
      const fallback = 1;
      if (!Number.isFinite(raw)) return fallback;
      return Math.max(0, Math.floor(raw));
    })(),
    run_label: el.designRunLabel.value.trim() || null,
    run_assess: true,
  };

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
  const limitVal = Number(el.resultsLimit?.value) || null;
  const params = new URLSearchParams();
  if (runLabel) params.append('run_label', runLabel);
  if (limitVal) params.append('limit', String(limitVal));

  if (!silent && el.refreshResultsBtn) {
    el.refreshResultsBtn.disabled = true;
  }
  if (!silent) {
    setBadge(el.resultsMeta, `Loading rankings for ${state.currentPdb}…`);
    el.resultsMeta.hidden = false;
  }

  const requestedSignature = analysisSignature(
    state.currentPdb,
    runLabel || '',
    state.rankingsResponse?.source_path || '',
  );
  const shouldPreserveAnalysis =
    state.analysisContext && state.analysisContext.signature === requestedSignature;

  if (!shouldPreserveAnalysis) {
    resetAnalysisPanel({ disableButton: true });
  } else if (el.analysisRun) {
    el.analysisRun.disabled = true;
  }

  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/rankings?${params.toString()}`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    const newSignature = analysisSignature(
      state.currentPdb,
      payload.run_label || '',
      payload.source_path || '',
    );
    const contextMatches =
      state.analysisContext && state.analysisContext.signature === newSignature;
    const previousSummary = state.librarySummary;
    const newSource = payload.source_path ? String(payload.source_path) : '';
    const keepSummary =
      Boolean(previousSummary && previousSummary.rankings_source && newSource) &&
      String(previousSummary.rankings_source) === newSource;
    state.rankings = payload.rows || [];
    state.rankingsResponse = payload;
    if (!keepSummary) {
      state.librarySummary = null;
    }
    state.tableSort = { key: 'iptm', dir: 'desc' };
    state.scatterLayout = [];
    state.selectedDesign = null;
    state.activeRunLabel = payload.run_label || '';
    state.galleryAvailable = Boolean(payload.gallery_path);
    updateActiveRunDisplay();
    if (el.resultsRunLabel) {
      el.resultsRunLabel.value = state.activeRunLabel;
    }
    el.binderDetail.hidden = true;
    if (!contextMatches) {
      resetAnalysisPanel();
    }
    renderResults();
    if (!silent) {
      const label = state.activeRunLabel || 'latest';
      setBadge(el.resultsMeta, `${payload.rows.length} designs · run ${label}`);
    }
    highlightRunChip(state.activeRunLabel);
    if (state.currentPdb) {
      fetchRunHistory(state.currentPdb);
    }
    renderPresets();
    scheduleRankingsPolling();
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
    if (el.analysisRun && !state.analysisRunning) {
      el.analysisRun.disabled = state.rankings.length === 0;
    }
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
}

function renderResultsTable() {
  if (!el.resultsTable) return;
  const tbody = el.resultsTable.querySelector('tbody');
  if (!tbody) return;
  tbody.innerHTML = '';
  state.rankings.forEach((row) => {
    const tr = document.createElement('tr');
    tr.dataset.design = row.design_name;
    tr.innerHTML = `
      <td>${row.design_name}</td>
      <td>${row.iptm !== null && row.iptm !== undefined ? row.iptm.toFixed(3) : '—'}</td>
      <td>${row.rmsd_diego !== null && row.rmsd_diego !== undefined ? row.rmsd_diego.toFixed(3) : '—'}</td>
      <td>${row.tm_score !== null && row.tm_score !== undefined ? row.tm_score.toFixed(3) : '—'}</td>
      <td>${row.metadata.epitope || row.metadata.arm || ''}</td>
    `;
    tr.addEventListener('click', () => selectDesign(row.design_name));
    tbody.appendChild(tr);
  });
}

function computeScatterLayout(points) {
  if (!el.scatterCanvas) {
    state.scatterLayout = [];
    state.scatterStats = null;
    return;
  }
  const filtered = points
    .filter((p) => typeof p.iptm === 'number' && typeof p.rmsd_diego === 'number')
    .map((p) => ({ ...p }));
  if (filtered.length === 0) {
    state.scatterLayout = [];
    state.scatterStats = null;
    return;
  }
  const margin = { left: 72, right: 32, top: 48, bottom: 72 };
  const width = el.scatterCanvas.width;
  const height = el.scatterCanvas.height;
  const xs = filtered.map((p) => p.rmsd_diego);
  const ys = filtered.map((p) => p.iptm);
  const xMin = Math.min(0, ...xs);
  const xMax = Math.max(5, ...xs);
  const yMin = Math.min(0, ...ys);
  const yMax = Math.max(1, ...ys);

  const xRange = xMax - xMin || 1;
  const yRange = yMax - yMin || 1;
  const thresholds = {
    rmsd: Number.isFinite(state.scatterThresholds?.rmsd) ? state.scatterThresholds.rmsd : 3.7,
    iptm: Number.isFinite(state.scatterThresholds?.iptm) ? state.scatterThresholds.iptm : 0.8,
  };
  const xBinCount = 24;
  const yBinCount = 24;
  const xBins = Array.from({ length: xBinCount }, () => 0);
  const yBins = Array.from({ length: yBinCount }, () => 0);
  let passingCount = 0;

  state.scatterLayout = filtered.map((p) => {
    const normalizedX = (p.rmsd_diego - xMin) / xRange;
    const normalizedY = (p.iptm - yMin) / yRange;
    const x = margin.left + normalizedX * (width - margin.left - margin.right);
    const y = height - margin.bottom - normalizedY * (height - margin.top - margin.bottom);
    if (Number.isFinite(thresholds.rmsd) && Number.isFinite(thresholds.iptm)) {
      if (p.rmsd_diego <= thresholds.rmsd && p.iptm >= thresholds.iptm) {
        passingCount += 1;
      }
    }
    if (Number.isFinite(normalizedX)) {
      const binX = Math.min(xBinCount - 1, Math.max(0, Math.floor(normalizedX * xBinCount)));
      xBins[binX] += 1;
    }
    if (Number.isFinite(normalizedY)) {
      const binY = Math.min(yBinCount - 1, Math.max(0, Math.floor(normalizedY * yBinCount)));
      yBins[binY] += 1;
    }
    return {
      ...p,
      x,
      y,
      xMin,
      xMax,
      yMin,
      yMax,
      xRange,
      yRange,
      margin,
    };
  });

  state.scatterStats = {
    xBins,
    yBins,
    total: filtered.length,
    passingCount,
    thresholds,
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

function renderScatter() {
  if (!el.scatterCanvas) return;
  const ctx = el.scatterCanvas.getContext('2d');
  ctx.clearRect(0, 0, el.scatterCanvas.width, el.scatterCanvas.height);
  renderScatterHistograms(null);
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
  ctx.fillText('Binder RMSD (Å)', margin.left + plotWidth / 2, height - margin.bottom + 36);
  ctx.save();
  ctx.translate(margin.left - 56, margin.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('ipTM', 0, 0);
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
    ctx.fillText(value.toFixed(1), x, height - margin.bottom + 10);
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
    ctx.fillText(value.toFixed(2), margin.left - 10, y);
  }

  const hasRmsdThreshold = Number.isFinite(thresholds.rmsd);
  const hasIptmThreshold = Number.isFinite(thresholds.iptm);
  if (hasRmsdThreshold && thresholds.rmsd >= xMin && thresholds.rmsd <= xMax) {
    const t = (thresholds.rmsd - xMin) / xRange;
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
    ctx.fillText(`RMSD ≤ ${thresholds.rmsd.toFixed(2)} Å`, Math.min(x + 8, width - margin.right), margin.top + 18);
  }

  if (hasIptmThreshold && thresholds.iptm >= yMin && thresholds.iptm <= yMax) {
    const t = (thresholds.iptm - yMin) / yRange;
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
    ctx.fillText(`ipTM ≥ ${thresholds.iptm.toFixed(2)}`, margin.left + 8, Math.max(margin.top + 8, y - 18));
  }

  state.scatterLayout.forEach((point) => {
    ctx.beginPath();
    ctx.fillStyle = '#2563eb';
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
  const { thresholds = {}, xBins = [], yBins = [] } = state.scatterStats;

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
    if (Number.isFinite(thresholds.rmsd) && thresholds.rmsd >= xMin && thresholds.rmsd <= xMin + xRange) {
      const t = (thresholds.rmsd - xMin) / xRange;
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
    if (Number.isFinite(thresholds.iptm) && thresholds.iptm >= yMin && thresholds.iptm <= yMin + yRange) {
      const t = (thresholds.iptm - yMin) / yRange;
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
  const { thresholds, passingCount, total } = stats;
  if (!thresholds || !Number.isFinite(thresholds.rmsd) || !Number.isFinite(thresholds.iptm)) {
    el.scatterThresholdCount.textContent = '';
    return;
  }
  const percent = total === 0 ? 0 : ((passingCount / total) * 100);
  const formatted = `${passingCount} / ${total} designs (${percent.toFixed(1)}%) meet RMSD ≤ ${thresholds.rmsd.toFixed(2)} Å and ipTM ≥ ${thresholds.iptm.toFixed(2)}`;
  el.scatterThresholdCount.textContent = formatted;
}

function handleScatterThresholdChange(key, rawValue) {
  if (!Number.isFinite(rawValue)) return;
  const thresholds = {
    rmsd: Number.isFinite(state.scatterThresholds?.rmsd) ? state.scatterThresholds.rmsd : 3.7,
    iptm: Number.isFinite(state.scatterThresholds?.iptm) ? state.scatterThresholds.iptm : 0.8,
  };
  let value = rawValue;
  if (key === 'rmsd') {
    value = Math.max(0, value);
  }
  if (key === 'iptm') {
    value = Math.min(1, Math.max(0, value));
  }
  if (!Number.isFinite(value)) return;
  const updated = { ...thresholds, [key]: value };
  state.scatterThresholds = updated;
  if (key === 'rmsd' && el.scatterThresholdRmsd) {
    el.scatterThresholdRmsd.value = value;
  }
  if (key === 'iptm' && el.scatterThresholdIptm) {
    el.scatterThresholdIptm.value = value;
  }
  computeScatterLayout(state.scatterPoints || []);
  renderScatter();
}

function analysisSignature(pdbId, runLabel, sourcePath) {
  return [pdbId || '', runLabel || '', sourcePath || ''].join('::');
}

function resetResultsPanel(options = {}) {
  const preserveMeta = Boolean(options.preserveMeta);
  const preserveRunLabel = Boolean(options.preserveRunLabel);
  state.rankings = [];
  state.rankingsResponse = null;
  state.scatterLayout = [];
  state.scatterPoints = [];
  state.scatterStats = null;
  state.selectedDesign = null;
  state.galleryAvailable = false;
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
  const scatterPoints = state.rankingsResponse?.scatter || [];
  state.scatterPoints = scatterPoints;
  const thresholds = state.scatterThresholds || {};
  if (el.scatterThresholdRmsd && Number.isFinite(thresholds.rmsd)) {
    el.scatterThresholdRmsd.value = thresholds.rmsd;
  }
  if (el.scatterThresholdIptm && Number.isFinite(thresholds.iptm)) {
    el.scatterThresholdIptm.value = thresholds.iptm;
  }
  computeScatterLayout(scatterPoints);
  renderScatter();
  const hasRows = state.rankings.length > 0;
  if (el.resultsTableWrapper) el.resultsTableWrapper.hidden = !hasRows;
  if (el.plotContainer) el.plotContainer.hidden = state.scatterLayout.length === 0;
  if (el.exportButton) el.exportButton.disabled = !hasRows;
  if (el.exportOpen) el.exportOpen.disabled = !hasRows;
  if (el.libraryRun) el.libraryRun.disabled = !hasRows;
  if (el.pymolTop) {
    el.pymolTop.disabled = !hasRows;
    el.pymolTop.textContent = state.galleryAvailable && hasRows ? 'Launch PyMOL Gallery' : 'PyMOL top 96';
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
  const info = [
    `<strong>${row.design_name}</strong>`,
    `ipTM: ${row.iptm !== null && row.iptm !== undefined ? row.iptm.toFixed(3) : 'NA'}`,
    `Binder RMSD: ${row.rmsd_diego !== null && row.rmsd_diego !== undefined ? row.rmsd_diego.toFixed(3) : 'NA'}`,
    row.metadata.arm ? `Arm: ${row.metadata.arm}` : '',
    row.metadata.epitope ? `Epitope: ${row.metadata.epitope}` : '',
    row.metadata.pymol_script_path ? `PyMOL script: ${row.metadata.pymol_script_path}` : '',
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
        if (state.activeRunLabel) {
          fetchRankings({ silent: true });
        }
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
  el.pymolTop.disabled = true;
  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/pymol/top-binders`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        run_label: state.rankingsResponse.run_label,
        top_n: Number(el.exportTopN.value) || 96,
        launch: true,
        bundle_only: false,
      }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    showAlert(`PyMOL aggregate script: ${payload.bundle_path || 'generated'}`, false);
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

function registerScatterClick() {
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
  if (el.libraryRun) el.libraryRun.addEventListener('click', runGoldenGatePlanner);
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
  if (el.resultsLimit) {
    el.resultsLimit.addEventListener('change', () => {
      fetchRankings({ silent: false });
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
  if (el.scatterThresholdRmsd) {
    const handler = (event) => {
      const value = Number.parseFloat(event.target.value);
      if (Number.isFinite(value)) {
        handleScatterThresholdChange('rmsd', value);
      }
    };
    el.scatterThresholdRmsd.addEventListener('change', handler);
    el.scatterThresholdRmsd.addEventListener('input', handler);
  }
  if (el.scatterThresholdIptm) {
    const handler = (event) => {
      const value = Number.parseFloat(event.target.value);
      if (Number.isFinite(value)) {
        handleScatterThresholdChange('iptm', value);
      }
    };
    el.scatterThresholdIptm.addEventListener('change', handler);
    el.scatterThresholdIptm.addEventListener('input', handler);
  }
  registerScatterClick();
  handleTableHeaderClicks();
  disableFutureSections();
}

function init() {
  initEventHandlers();
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
