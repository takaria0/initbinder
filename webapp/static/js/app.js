const state = {
  currentPdb: null,
  jobPoller: null,
  jobContext: null,
  rankings: [],
  rankingsResponse: null,
  tableSort: { key: 'iptm', dir: 'desc' },
  scatterLayout: [],
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
};

const el = {
  targetForm: document.querySelector('#target-form'),
  clusterStatus: document.querySelector('#cluster-status'),
  debugToggle: document.querySelector('#debug-toggle'),
  pdbInput: document.querySelector('#pdb-id'),
  antigenInput: document.querySelector('#antigen-url'),
  targetName: document.querySelector('#target-name'),
  targetEpitopes: document.querySelector('#target-epitopes'),
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
  syncResultsBtn: document.querySelector('#sync-results'),
  assessSubmit: document.querySelector('#assess-submit'),
  jobList: document.querySelector('#job-list'),
  refreshJobsBtn: document.querySelector('#refresh-jobs'),
  runHistory: document.querySelector('#run-history'),
  runLabelOptions: document.querySelector('#run-label-options'),
  activeRunName: document.querySelector('#active-run-name'),
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
  state.currentPdb = upper;
  if (options.updateInput !== false && el.pdbInput) {
    el.pdbInput.value = upper;
  }
  state.activeRunLabel = '';
  updateActiveRunDisplay();
  syncActivePreset();
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
      setClusterStatus('error', `Cluster: not connected - ${msg}. Run ssh hpc3.rcic.uci.edu -MNf`);
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

async function fetchRankings(options = {}) {
  const { silent = false } = options;
  if (!state.currentPdb) {
    if (!silent) showAlert('Initialize target first.');
    return;
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

  try {
    const res = await fetch(`/api/targets/${state.currentPdb}/rankings?${params.toString()}`);
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    state.rankings = payload.rows || [];
    state.rankingsResponse = payload;
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
  const tbody = el.resultsTable.querySelector('tbody');
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
  const filtered = points
    .filter((p) => typeof p.iptm === 'number' && typeof p.rmsd_diego === 'number')
    .map((p) => ({ ...p }));
  if (filtered.length === 0) {
    state.scatterLayout = [];
    return;
  }
  const margin = { left: 48, right: 12, top: 12, bottom: 40 };
  const width = el.scatterCanvas.width;
  const height = el.scatterCanvas.height;
  const xs = filtered.map((p) => p.iptm);
  const ys = filtered.map((p) => p.rmsd_diego);
  const xMin = Math.min(0, ...xs);
  const xMax = Math.max(1, ...xs);
  const yMin = Math.min(0, ...ys);
  const yMax = Math.max(...ys, 5);

  state.scatterLayout = filtered.map((p) => {
    const x = margin.left + ((p.iptm - xMin) / (xMax - xMin || 1)) * (width - margin.left - margin.right);
    const y = height - margin.bottom - ((p.rmsd_diego - yMin) / (yMax - yMin || 1)) * (height - margin.top - margin.bottom);
    return { ...p, x, y, xMin, xMax, yMin, yMax, margin };
  });
}

function renderScatter() {
  const ctx = el.scatterCanvas.getContext('2d');
  ctx.clearRect(0, 0, el.scatterCanvas.width, el.scatterCanvas.height);
  if (state.scatterLayout.length === 0) {
    ctx.fillStyle = '#64748b';
    ctx.fillText('No scatter data available', 20, 40);
    return;
  }

  const sample = state.scatterLayout[0];
  const { xMin, xMax, yMin, yMax, margin } = sample;
  const width = el.scatterCanvas.width;
  const height = el.scatterCanvas.height;

  ctx.strokeStyle = '#cbd5f5';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(margin.left, margin.top);
  ctx.lineTo(margin.left, height - margin.bottom);
  ctx.lineTo(width - margin.right, height - margin.bottom);
  ctx.stroke();

  ctx.fillStyle = '#0f172a';
  ctx.font = '12px sans-serif';
  ctx.fillText('ipTM', width / 2, height - 8);
  ctx.save();
  ctx.translate(16, height / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('RMSD (Å)', 0, 0);
  ctx.restore();

  state.scatterLayout.forEach((point) => {
    const selected = point.design_name === state.selectedDesign;
    ctx.beginPath();
    ctx.fillStyle = selected ? '#1e3a8a' : '#2563eb';
    ctx.globalAlpha = selected ? 0.95 : 0.8;
    ctx.arc(point.x, point.y, selected ? 6 : 4, 0, Math.PI * 2);
    ctx.fill();
    ctx.globalAlpha = 1;
  });
}

function renderResults() {
  sortRankings();
  renderResultsTable();
  const scatterPoints = state.rankingsResponse?.scatter || [];
  computeScatterLayout(scatterPoints);
  renderScatter();
  const hasRows = state.rankings.length > 0;
  if (el.resultsTableWrapper) el.resultsTableWrapper.hidden = !hasRows;
  if (el.plotContainer) el.plotContainer.hidden = state.scatterLayout.length === 0;
  if (el.exportButton) el.exportButton.disabled = !hasRows;
  if (el.exportOpen) el.exportOpen.disabled = !hasRows;
  if (el.pymolTop) {
    el.pymolTop.disabled = !hasRows;
    el.pymolTop.textContent = state.galleryAvailable && hasRows ? 'Launch PyMOL Gallery' : 'PyMOL top 96';
  }
  if (el.refreshResultsBtn) el.refreshResultsBtn.disabled = false;
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
    `RMSD Diego: ${row.rmsd_diego !== null && row.rmsd_diego !== undefined ? row.rmsd_diego.toFixed(3) : 'NA'}`,
    row.metadata.arm ? `Arm: ${row.metadata.arm}` : '',
    row.metadata.epitope ? `Epitope: ${row.metadata.epitope}` : '',
    row.metadata.pymol_script_path ? `PyMOL script: ${row.metadata.pymol_script_path}` : '',
  ].filter(Boolean);
  el.binderDetail.innerHTML = info.join('<br>');
  el.binderDetail.hidden = false;
  renderScatter();
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
       if (el.assessSubmit) el.assessSubmit.disabled = false;
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
  } else if (ctx === 'sync') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.resultsMeta, job.status === 'running' ? 'Syncing…' : 'Queued…');
    } else if (job.status === 'success') {
      const skipped = Boolean(job.details?.skipped);
      const localPath = job.details?.local_path ? String(job.details.local_path) : '';
      const statusText = skipped ? 'Assessments already synced' : 'Assessments synced';
      const suffix = localPath ? ` → ${localPath}` : '';
      setBadge(el.resultsMeta, `${statusText}${suffix}`, 'rgba(134, 239, 172, 0.25)');
      if (el.syncResultsBtn) el.syncResultsBtn.disabled = false;
      stopJobPolling();
      const alertMsg = localPath ? `${statusText}. Saved to ${localPath}` : statusText;
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

async function syncResultsFromCluster(options = {}) {
  const {
    runLabel: explicitRunLabel,
    silent = false,
    disableButton = true,
    useInputRunLabel = true,
  } = options;
  if (!state.currentPdb) {
    if (!silent) showAlert('Initialize target first.');
    return;
  }
  let runLabelInput = '';
  if (explicitRunLabel !== undefined) {
    runLabelInput = (explicitRunLabel || '').trim();
  } else if (useInputRunLabel && el.resultsRunLabel) {
    runLabelInput = el.resultsRunLabel.value.trim();
  }
  if (!runLabelInput && state.activeRunLabel) {
    runLabelInput = state.activeRunLabel.trim();
  }
  const params = runLabelInput ? `?run_label=${encodeURIComponent(runLabelInput)}` : '';
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
    const rect = el.scatterCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const hit = state.scatterLayout.find((point) => {
      const dx = point.x - x;
      const dy = point.y - y;
      return Math.sqrt(dx * dx + dy * dy) <= 6;
    });
    if (hit) selectDesign(hit.design_name);
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
  if (el.assessSubmit) el.assessSubmit.disabled = true;
  state.targetStatus = null;
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
  if (el.refreshResultsBtn) {
    el.refreshResultsBtn.addEventListener('click', () => fetchRankings({ silent: false }));
  }
  el.pymolHotspots.addEventListener('click', launchHotspots);
  el.pymolTop.addEventListener('click', launchTopBinders);
  if (el.syncResultsBtn) {
    el.syncResultsBtn.addEventListener('click', () => syncResultsFromCluster({ useInputRunLabel: false }));
  }
  if (el.assessSubmit) {
    el.assessSubmit.addEventListener('click', queueAssessmentRun);
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
