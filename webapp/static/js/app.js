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
};

const el = {
  targetForm: document.querySelector('#target-form'),
  clusterStatus: document.querySelector('#cluster-status'),
  pdbInput: document.querySelector('#pdb-id'),
  antigenInput: document.querySelector('#antigen-url'),
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
  refreshResultsBtn: document.querySelector('#refresh-results'),
  resultsMeta: document.querySelector('#results-meta'),
  binderDetail: document.querySelector('#binder-detail'),
  resultsTableWrapper: document.querySelector('#results-table-wrapper'),
  resultsTable: document.querySelector('#results-table'),
  scatterCanvas: document.querySelector('#scatter-plot'),
  plotContainer: document.querySelector('#plot-container'),
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
  jobList: document.querySelector('#job-list'),
  refreshJobsBtn: document.querySelector('#refresh-jobs'),
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
  } catch (err) {
    console.error('Failed to fetch jobs', err);
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
    const created = new Date(job.created_at * 1000).toLocaleString();
    meta.textContent = `${job.kind} - ${job.status} - ${created}`;
    li.appendChild(label);
    li.appendChild(meta);
    li.addEventListener('click', () => {
      state.selectedJobId = job.job_id;
      renderJobList();
      startJobPolling(job.job_id, null, { manual: true });
    });
    el.jobList.appendChild(li);
  });
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
  state.currentPdb = pdbId;
  setBadge(el.targetBadge, 'Queued…');
  el.targetSubmit.disabled = true;
  el.pymolHotspots.disabled = true;
  el.jobAlert.hidden = true;
  resetJobLog('Submitting init-target job…');

  const payload = {
    pdb_id: pdbId,
    antigen_url: el.antigenInput.value.trim() || null,
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
    binder_chain_id: el.designBinderChain.value.trim() || 'H',
    run_label: el.designRunLabel.value.trim() || null,
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
  } catch (err) {
    showAlert(err.message || String(err));
    el.exportButton.disabled = false;
  }
}

async function fetchRankings() {
  if (!state.currentPdb) {
    showAlert('Initialize target first.');
    return;
  }
  const runLabel = el.resultsRunLabel.value.trim();
  const limitVal = Number(el.resultsLimit.value) || null;
  const params = new URLSearchParams();
  if (runLabel) params.append('run_label', runLabel);
  if (limitVal) params.append('limit', String(limitVal));

  el.refreshResultsBtn.disabled = true;
  setBadge(el.resultsMeta, `Loading rankings for ${state.currentPdb}…`);
  el.resultsMeta.hidden = false;

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
    el.binderDetail.hidden = true;
    renderResults();
    setBadge(
      el.resultsMeta,
      `${payload.rows.length} designs · run ${payload.run_label || 'latest'}`,
    );
  } catch (err) {
    setBadge(el.resultsMeta, `Rankings failed: ${err.message || err}`, 'rgba(248, 113, 113, 0.2)');
  } finally {
    el.refreshResultsBtn.disabled = false;
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
  el.resultsTableWrapper.hidden = state.rankings.length === 0;
  el.plotContainer.hidden = state.scatterLayout.length === 0;
  el.exportButton.disabled = state.rankings.length === 0;
  el.pymolTop.disabled = state.rankings.length === 0;
  el.refreshResultsBtn.disabled = false;
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
      if (context === 'export') el.exportButton.disabled = false;
    }
  };
  poll();
  state.jobPoller = setInterval(poll, 2000);
}

function updateJobUI(job) {
  if (!job) return;
  resetJobLog(job.logs.join('\n'));

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
      refreshRunLabel(true);
    } else {
      setBadge(el.targetBadge, 'Failed', 'rgba(248, 113, 113, 0.25)');
      showAlert(job.message || 'Target workflow failed.');
      el.targetSubmit.disabled = false;
      stopJobPolling();
    }
  } else if (ctx === 'design') {
    if (job.status === 'running' || job.status === 'pending') {
      setBadge(el.designStatus, job.status === 'running' ? 'Running…' : 'Queued…');
    } else if (job.status === 'success') {
      setBadge(el.designStatus, 'Cluster jobs submitted', 'rgba(134, 239, 172, 0.25)');
      el.designSubmit.disabled = false;
      el.refreshResultsBtn.disabled = false;
      stopJobPolling();
      refreshRunLabel(true);
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
      const outDir = job.details?.out_dir ? ` → ${job.details.out_dir}` : '';
      setBadge(el.resultsMeta, `Export complete${outDir}`, 'rgba(134, 239, 172, 0.25)');
      el.exportButton.disabled = false;
      stopJobPolling();
    } else {
      showAlert(job.message || 'Export failed.');
      el.exportButton.disabled = false;
      stopJobPolling();
    }
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
    showAlert(`Hotspot bundle ready at ${payload.bundle_path || 'temp directory'}`, false);
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

async function syncResultsFromCluster() {
  if (!state.currentPdb) {
    showAlert('Initialize target first.');
    return;
  }
  el.syncResultsBtn.disabled = true;
  try {
    const runLabel = el.resultsRunLabel.value.trim();
    const params = runLabel ? `?run_label=${encodeURIComponent(runLabel)}` : '';
    const res = await fetch(`/api/targets/${state.currentPdb}/sync${params}`, {
      method: 'POST',
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Failed with ${res.status}`);
    }
    const payload = await res.json();
    showAlert(payload.message || 'Synced results.', false);
  } catch (err) {
    showAlert(err.message || String(err));
  } finally {
    el.syncResultsBtn.disabled = false;
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
  el.designSubmit.disabled = true;
  setBadge(el.designStatus, 'Awaiting target prep');
  el.refreshResultsBtn.disabled = true;
  el.exportButton.disabled = true;
  el.pymolTop.disabled = true;
  el.pymolHotspots.disabled = true;
}

function initEventHandlers() {
  el.targetForm.addEventListener('submit', queueTargetInit);
  el.refreshAlignmentBtn.addEventListener('click', () => {
    if (state.currentPdb) loadAlignment(state.currentPdb);
  });
  el.designSubmit.addEventListener('click', queueDesignRun);
  el.exportButton.addEventListener('click', runExport);
  el.refreshResultsBtn.addEventListener('click', fetchRankings);
  el.pymolHotspots.addEventListener('click', launchHotspots);
  el.pymolTop.addEventListener('click', launchTopBinders);
  el.syncResultsBtn.addEventListener('click', syncResultsFromCluster);
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
  refreshRunLabel(true);
  updateClusterStatus();
  setInterval(updateClusterStatus, 15000);
  fetchJobList();
}

document.addEventListener('DOMContentLoaded', init);
