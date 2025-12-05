const state = {
  bulkPreviewRows: [],
  jobPoller: null,
  currentJobId: null,
};

const el = {
  bulkStatus: document.querySelector('#bulk-status'),
  bulkCsvInput: document.querySelector('#bulk-csv-input'),
  bulkEpitopes: document.querySelector('#bulk-epitopes'),
  bulkEpitopePrompt: document.querySelector('#bulk-epitope-prompt'),
  bulkLaunchPymol: document.querySelector('#bulk-launch-pymol'),
  bulkExportInsights: document.querySelector('#bulk-export-insights'),
  bulkExportDesigns: document.querySelector('#bulk-export-designs'),
  bulkSubmitDesigns: document.querySelector('#bulk-submit-designs'),
  bulkForceInit: document.querySelector('#bulk-force-init'),
  designEngine: document.querySelector('#design-engine'),
  designTotal: document.querySelector('#design-total'),
  designNumSeq: document.querySelector('#design-num-seq'),
  designTemp: document.querySelector('#design-temp'),
  designBinderChain: document.querySelector('#design-binder-chain'),
  designAf3Seed: document.querySelector('#design-af3-seed'),
  designEnableRfdiffCrop: document.querySelector('#design-enable-rfdiff-crop'),
  designBoltzBinding: document.querySelector('#design-boltz-binding'),
  bulkRunPrefix: document.querySelector('#bulk-run-prefix'),
  bulkThrottle: document.querySelector('#bulk-throttle'),
  bulkPreviewBtn: document.querySelector('#bulk-preview'),
  bulkRunBtn: document.querySelector('#bulk-run'),
  bulkImportQueue: document.querySelector('#bulk-import-queue'),
  bulkImportInput: document.querySelector('#bulk-design-import'),
  bulkPreviewPanel: document.querySelector('#bulk-preview-panel'),
  bulkPreviewTable: document.querySelector('#bulk-preview-table tbody'),
  bulkPreviewSummary: document.querySelector('#bulk-preview-summary'),
  bulkPreviewRefresh: document.querySelector('#bulk-preview-refresh'),
  jobAlert: document.querySelector('#job-alert'),
  jobLog: document.querySelector('#job-log'),
  jobMeta: document.querySelector('#job-meta'),
  snapshotSection: document.querySelector('#snapshot-section'),
  snapshotGrid: document.querySelector('#snapshot-grid'),
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

function renderSnapshots(meta = []) {
  if (!el.snapshotGrid || !el.snapshotSection) return;
  el.snapshotGrid.innerHTML = '';
  const list = Array.isArray(meta) ? meta : [];
  if (!list.length) {
    el.snapshotSection.hidden = true;
    return;
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
    el.snapshotGrid.appendChild(card);
  });
  el.snapshotSection.hidden = false;
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
      row.antigen_url || '',
      row.resolved_pdb_id || row.pdb_id || '—',
      row.preset_id ? 'Matched preset' : '—',
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
    renderBulkPreview(state.bulkPreviewRows, body.message || '');
    setBadge(el.bulkStatus, body.message || 'Preview ready', 'rgba(134, 239, 172, 0.25)');
    if (!silent) showAlert(body.message || 'Preview ready.', false);
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
    setBadge(el.bulkStatus, 'Preview failed', 'rgba(248, 113, 113, 0.25)');
  }
}

function collectDesignSettings() {
  const settings = {
    model_engine: (el.designEngine?.value || 'rfantibody').trim() || 'rfantibody',
    total_designs: Number(el.designTotal?.value || 0) || 90,
    num_sequences: Number(el.designNumSeq?.value || 0) || 1,
    temperature: Number(el.designTemp?.value || 0.1),
    binder_chain_id: (el.designBinderChain?.value || '').trim() || null,
    af3_seed: Number(el.designAf3Seed?.value || 1) || 1,
    run_assess: true,
    rfdiff_crop_radius: null,
    run_label_prefix: (el.bulkRunPrefix?.value || '').trim() || null,
    boltz_binding: (el.designBoltzBinding?.value || '').trim() || null,
  };
  if (!Number.isFinite(settings.temperature)) settings.temperature = 0.1;
  if (el.designEnableRfdiffCrop && el.designEnableRfdiffCrop.checked) {
    settings.rfdiff_crop_radius = 14;
  }
  return settings;
}

async function startBulkRun() {
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
    const raw = Number(el.bulkThrottle?.value || 0);
    if (!Number.isFinite(raw) || raw < 0) return 0;
    return raw;
  })();
  const payload = {
    csv_text: csvText,
    num_epitopes: numEpitopes,
    decide_scope_prompt: (el.bulkEpitopePrompt?.value || '').trim() || null,
    launch_pymol: el.bulkLaunchPymol ? el.bulkLaunchPymol.checked : true,
    export_insights: el.bulkExportInsights ? el.bulkExportInsights.checked : true,
    export_designs: el.bulkExportDesigns ? el.bulkExportDesigns.checked : true,
    submit_designs: el.bulkSubmitDesigns ? el.bulkSubmitDesigns.checked : false,
    force_init: el.bulkForceInit ? el.bulkForceInit.checked : false,
    prepare_targets: true,
    design_settings: collectDesignSettings(),
    throttle_seconds: throttle,
  };
  if (el.bulkRunBtn) el.bulkRunBtn.disabled = true;
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
    } else if (el.bulkRunBtn) {
      el.bulkRunBtn.disabled = false;
    }
  } catch (err) {
    showAlert(err.message || String(err));
    if (el.bulkRunBtn) el.bulkRunBtn.disabled = false;
    setBadge(el.bulkStatus, 'Bulk submission failed', 'rgba(248, 113, 113, 0.25)');
  }
}

async function importBulkDesigns() {
  const csvText = el.bulkImportInput?.value || '';
  if (!csvText.trim()) {
    showAlert('Paste a design configuration CSV to import.');
    return;
  }
  const throttle = (() => {
    const raw = Number(el.bulkThrottle?.value || 0);
    if (!Number.isFinite(raw) || raw < 0) return 0;
    return raw;
  })();
  if (el.bulkImportQueue) el.bulkImportQueue.disabled = true;
  setBadge(el.bulkStatus, 'Queuing design submissions…');
  resetJobLog('Waiting for commands...');
  try {
    const res = await fetch('/api/bulk/designs/import', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ csv_text: csvText, throttle_seconds: throttle }),
    });
    if (!res.ok) {
      const detail = await res.json().catch(() => ({}));
      throw new Error(detail.detail || `Import failed with ${res.status}`);
    }
    const body = await res.json();
    showAlert(body.message || 'Bulk design submissions queued.', false);
    if (body.job_id) {
      startJobPolling(body.job_id);
    } else if (el.bulkImportQueue) {
      el.bulkImportQueue.disabled = false;
    }
  } catch (err) {
    showAlert(err.message || String(err));
    if (el.bulkImportQueue) el.bulkImportQueue.disabled = false;
    setBadge(el.bulkStatus, 'Bulk design import failed', 'rgba(248, 113, 113, 0.25)');
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
  resetJobLog(job.logs.join('\n'));
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
    if (Array.isArray(job.details?.snapshots) && job.details.snapshots.length) {
      job.details.snapshots.forEach((snap) => {
        const name = (snap || '').trim();
        if (name) addLink('Snapshot', name);
      });
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
    if (el.bulkRunBtn) el.bulkRunBtn.disabled = true;
    if (el.bulkImportQueue) el.bulkImportQueue.disabled = true;
  } else if (job.status === 'success') {
    setBadge(el.bulkStatus, 'Bulk complete', 'rgba(134, 239, 172, 0.25)');
    showAlert(job.message || 'Bulk job finished.', false);
    if (el.bulkRunBtn) el.bulkRunBtn.disabled = false;
    if (el.bulkImportQueue) el.bulkImportQueue.disabled = false;
    if (job.details?.design_config_filename) {
      const url = `/api/bulk/file?name=${encodeURIComponent(job.details.design_config_filename)}`;
      window.open(url, '_blank');
    }
    if (Array.isArray(job.details?.snapshots) && job.details.snapshots.length) {
      loadSnapshotMetadata(job.details.snapshots.filter(Boolean));
    } else {
      renderSnapshots([]);
    }
    stopJobPolling();
  } else {
    setBadge(el.bulkStatus, 'Bulk failed', 'rgba(248, 113, 113, 0.25)');
    showAlert(job.message || 'Bulk job failed.');
    if (el.bulkRunBtn) el.bulkRunBtn.disabled = false;
    if (el.bulkImportQueue) el.bulkImportQueue.disabled = false;
    stopJobPolling();
  }
}

function startJobPolling(jobId) {
  state.currentJobId = jobId;
  stopJobPolling();
  const poll = async () => {
    try {
      const job = await fetchJob(jobId);
      updateJobUI(job);
    } catch (err) {
      appendLog(`Polling error: ${err.message || err}`);
      stopJobPolling();
      if (el.bulkRunBtn) el.bulkRunBtn.disabled = false;
      if (el.bulkImportQueue) el.bulkImportQueue.disabled = false;
    }
  };
  poll();
  state.jobPoller = setInterval(poll, 2000);
}

function toggleBoltzBinding(show) {
  if (!el.designBoltzBinding) return;
  const row = el.designBoltzBinding.closest('.form-row');
  if (row) row.hidden = !show;
}

function init() {
  if (el.bulkPreviewBtn) el.bulkPreviewBtn.addEventListener('click', () => previewBulkCsv({ silent: false }));
  if (el.bulkPreviewRefresh) el.bulkPreviewRefresh.addEventListener('click', () => previewBulkCsv({ silent: true }));
  if (el.bulkRunBtn) el.bulkRunBtn.addEventListener('click', startBulkRun);
  if (el.bulkImportQueue) el.bulkImportQueue.addEventListener('click', importBulkDesigns);
  if (el.designEngine) {
    el.designEngine.addEventListener('change', () => {
      toggleBoltzBinding((el.designEngine.value || '').trim() === 'boltzgen');
    });
    toggleBoltzBinding((el.designEngine.value || '').trim() === 'boltzgen');
  }
}

document.addEventListener('DOMContentLoaded', init);
