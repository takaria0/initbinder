const state = {
  bulkPreviewRows: [],
  jobPoller: null,
  currentJobId: null,
  snapshotNames: [],
};

const el = {
  bulkStatus: document.querySelector('#bulk-status'),
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

function renderSnapshots(meta = []) {
  if (!el.snapshotGrid || !el.snapshotSection) return;
  el.snapshotGrid.innerHTML = '';
  const list = Array.isArray(meta) ? meta : [];
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
      row.accession || '',
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
    model_engine: (el.designEngine?.value || 'boltzgen').trim() || 'boltzgen',
    total_designs: Number(el.designTotal?.value || 0) || 90,
    run_assess: true,
    run_label_prefix: (el.bulkRunPrefix?.value || '').trim() || null,
  };
  return settings;
}

async function startBulkRun(options = {}) {
  const {
    submitDesignsOverride = null,
    triggerEl = null,
    llmDelaySeconds = null,
    decideScopeAttempts = null,
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
    setActionButtonsDisabled(true);
  } else if (job.status === 'success') {
    setBadge(el.bulkStatus, 'Bulk complete', 'rgba(134, 239, 172, 0.25)');
    showAlert(job.message || 'Bulk job finished.', false);
    setActionButtonsDisabled(false);
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
    setActionButtonsDisabled(false);
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
    llmDelaySeconds: 70,
    decideScopeAttempts: 3,
  });
  scrollToSnapshots();
}

function init() {
  if (el.bulkPreviewBtn) el.bulkPreviewBtn.addEventListener('click', () => previewBulkCsv({ silent: false }));
  if (el.bulkPreviewRefresh) el.bulkPreviewRefresh.addEventListener('click', () => previewBulkCsv({ silent: true }));
  if (el.bulkVisualizeEpitopes) el.bulkVisualizeEpitopes.addEventListener('click', handleVisualizeEpitopes);
  if (el.bulkRunBtn) el.bulkRunBtn.addEventListener('click', () => startBulkRun());
  if (el.snapshotDownloadPdf) el.snapshotDownloadPdf.addEventListener('click', downloadSnapshotPdf);
}

document.addEventListener('DOMContentLoaded', init);
