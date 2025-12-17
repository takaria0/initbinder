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
  expandedTargets: new Set(),
  diversityCsv: null,
  diversityHtml: null,
  diversityMessage: null,
  diversityFiles: [],
  diversityOutputDir: null,
  binderRows: [],
  binderTotal: 0,
  binderPage: 1,
  binderPageSize: 100,
  binderCsvName: null,
  binderMessage: null,
  binderPdbIds: [],
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
  epitopePlotSection: document.querySelector('#epitope-plot-section'),
  epitopePlotGrid: document.querySelector('#epitope-plot-grid'),
  epitopeReportDownload: document.querySelector('#epitope-report-download'),
  epitopeMetricsSection: document.querySelector('#epitope-metrics-section'),
  diversitySection: document.querySelector('#diversity-section'),
  diversityGrid: document.querySelector('#diversity-grid'),
  diversityRefresh: document.querySelector('#diversity-refresh'),
  diversityDownloadCsv: document.querySelector('#diversity-download-csv'),
  diversityDownloadHtml: document.querySelector('#diversity-download-html'),
  boltzPanel: document.querySelector('#boltz-config-panel'),
  boltzTable: document.querySelector('#boltz-config-table tbody'),
  boltzSummary: document.querySelector('#boltz-config-summary'),
  boltzDesignCount: document.querySelector('#boltz-design-count'),
  boltzTimeHours: document.querySelector('#boltz-time-hours'),
  boltzRefresh: document.querySelector('#boltz-config-refresh'),
  boltzShowRunAll: document.querySelector('#boltz-show-run-all'),
  binderPanel: document.querySelector('#boltz-binders-panel'),
  binderTable: document.querySelector('#boltz-binders-table tbody'),
  binderSummary: document.querySelector('#boltz-binders-summary'),
  binderCsvNote: document.querySelector('#boltz-binders-csv'),
  binderRefresh: document.querySelector('#boltz-binders-refresh'),
  binderDownload: document.querySelector('#boltz-binders-download'),
  binderPagination: document.querySelector('#boltz-binders-pagination'),
  binderPageLabel: document.querySelector('#boltz-binders-page-label'),
  boltzConfigModal: document.querySelector('#boltz-config-modal'),
  boltzConfigTitle: document.querySelector('#boltz-config-title'),
  boltzConfigBody: document.querySelector('#boltz-config-body'),
  boltzConfigClose: document.querySelector('#boltz-config-close'),
  boltzLogModal: document.querySelector('#boltz-log-modal'),
  boltzLogTitle: document.querySelector('#boltz-log-title'),
  boltzLogBody: document.querySelector('#boltz-log-body'),
  boltzLogClose: document.querySelector('#boltz-log-close'),
  boltzRunModal: document.querySelector('#boltz-run-modal'),
  boltzRunTitle: document.querySelector('#boltz-run-title'),
  boltzRunBody: document.querySelector('#boltz-run-body'),
  boltzRunClose: document.querySelector('#boltz-run-close'),
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

function buildLocalResultBadge(done = false) {
  const pill = document.createElement('span');
  pill.className = 'badge';
  pill.textContent = done ? 'DONE' : 'No results';
  pill.style.background = done ? 'rgba(134, 239, 172, 0.25)' : 'rgba(148, 163, 184, 0.25)';
  pill.style.color = done ? '#14532d' : '#475569';
  return pill;
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

function bindingTextFromConfig(cfg = {}) {
  return cfg.binding_label || cfg.include_label || '—';
}

function epitopeLabel(cfg = {}, fallbackIndex = null) {
  if (cfg.epitope_name) return cfg.epitope_name;
  if (cfg.epitope_id) return cfg.epitope_id;
  if (fallbackIndex !== null) return `Epitope ${fallbackIndex}`;
  return 'Epitope';
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
  return Boolean(state.boltzLocalRuns?.[key]);
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
  const list = Array.isArray(meta) ? meta : [];
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

function renderEpitopePlots(items = []) {
  if (!el.epitopePlotGrid || !el.epitopePlotSection) return;
  el.epitopePlotGrid.innerHTML = '';
  const list = Array.isArray(items) ? items.filter(Boolean) : [];
  if (!list.length) {
    state.epitopePlots = [];
    el.epitopePlotSection.hidden = true;
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
      'No BoltzGen diversity plots yet. Download metrics to targets/<PDB>/designs/boltzgen/*/*/final_ranked_designs/all_designs_metrics.csv, then click Refresh.';
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

async function refreshDiversity({ silent = false } = {}) {
  try {
    const res = await fetch('/api/bulk/boltzgen/diversity');
    if (!res.ok) throw new Error('Unable to load diversity report');
    const data = await res.json();
    state.diversityCsv = data.csv_name || null;
    state.diversityHtml = data.html_name || null;
    state.diversityMessage = data.message || null;
    state.diversityFiles = Array.isArray(data.metrics_files) ? data.metrics_files : [];
    state.diversityOutputDir = data.output_dir || null;
    state.diversityPlots = Array.isArray(data.plots) ? data.plots : [];
    const pdbIds = new Set();
    (state.diversityPlots || []).forEach((plot) => {
      const val = (plot.pdb_id || '').trim().toUpperCase();
      if (val) pdbIds.add(val);
    });
    (state.diversityFiles || []).forEach((file) => {
      const val = (file.pdb_id || '').trim().toUpperCase();
      if (val) pdbIds.add(val);
    });
    state.binderPdbIds = Array.from(pdbIds);
    if (state.diversityCsv) {
      state.binderCsvName = state.binderCsvName || state.diversityCsv;
    }
    renderDiversityPlots(state.diversityPlots);
    if (!silent && data.message) {
      showAlert(data.message, false);
    }
    await loadBinderTable({ page: 1, silent: true });
  } catch (err) {
    if (!silent) showAlert(err.message || String(err));
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

  const noteRow = () => {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = 7;
    td.style.color = '#64748b';
    td.textContent = state.binderTotal ? 'No binders on this page.' : (state.binderMessage || 'No BoltzGen binders found.');
    tr.appendChild(td);
    return tr;
  };

  if (!rows.length) {
    tbody.appendChild(noteRow());
  } else {
    rows.forEach((row, idx) => {
      const tr = document.createElement('tr');
      const values = [
        row.pdb_id || '—',
        row.epitope || '—',
        row.rank ?? '—',
        row.iptm !== null && row.iptm !== undefined ? Number(row.iptm).toFixed(3) : '—',
        row.rmsd !== null && row.rmsd !== undefined ? Number(row.rmsd).toFixed(3) : '—',
      ];
      values.forEach((val) => {
        const td = document.createElement('td');
        td.textContent = val;
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
      if (!row.design_path) pymolBtn.disabled = true;
      pymolTd.appendChild(pymolBtn);
      tr.appendChild(pymolTd);

      tbody.appendChild(tr);
    });
  }

  if (el.binderSummary) {
    const msg = state.binderMessage || '';
    const summaryText = msg || `Showing ${rows.length} of ${state.binderTotal || rows.length} binders`;
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

async function launchBinderPymol(index) {
  const row = Array.isArray(state.binderRows) ? state.binderRows[index] : null;
  if (!row) {
    showAlert('No binder selected.');
    return;
  }
  if (!row.design_path) {
    showAlert('Design path unavailable for this binder.');
    return;
  }
  try {
    const payload = {
      pdb_id: row.pdb_id,
      design_path: row.design_path,
      epitope_label: row.epitope,
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
    el.boltzPanel.hidden = true;
    if (el.boltzSummary) el.boltzSummary.hidden = true;
    return;
  }
  let configCount = 0;
  targets.forEach((target, idx) => {
    const configs = Array.isArray(target.configs) ? target.configs : [];
    configCount += configs.length;
    const key = (target.pdb_id || '').toUpperCase();
    const isExpanded = state.expandedTargets instanceof Set ? state.expandedTargets.has(key) : false;
    const tr = document.createElement('tr');
    tr.className = `target-row ${isExpanded ? 'expanded' : 'collapsed'}`;
    tr.dataset.pdbId = key;
    const cells = [
      idx + 1,
      target.preset_name || '—',
      target.pdb_id || '',
      'All epitopes',
      '—',
    ];
    cells.forEach((val, cellIdx) => {
      const td = document.createElement('td');
      if (cellIdx === 0) {
        td.innerHTML = `<span class="toggle-icon">▸</span>${val}`;
      } else {
        td.textContent = val;
      }
      tr.appendChild(td);
    });
    const statusCell = document.createElement('td');
    statusCell.className = 'boltz-status-cell';
    statusCell.appendChild(buildLocalResultBadge(hasLocalBoltzResults(target.pdb_id)));
    tr.appendChild(statusCell);

    const cmdCell = document.createElement('td');
    cmdCell.className = 'boltz-cmd-cell';
    const runBtn = document.createElement('button');
    runBtn.type = 'button';
    runBtn.textContent = 'Run all epitopes';
    runBtn.dataset.action = 'run-target';
    runBtn.dataset.pdbId = target.pdb_id || '';
    cmdCell.appendChild(runBtn);

    const cfgBtn = document.createElement('button');
    cfgBtn.type = 'button';
    cfgBtn.textContent = 'Show config';
    cfgBtn.className = 'ghost';
    cfgBtn.dataset.action = 'show-config-target';
    cfgBtn.dataset.pdbId = target.pdb_id || '';
    cmdCell.appendChild(cfgBtn);

    const logBtn = document.createElement('button');
    logBtn.type = 'button';
    logBtn.textContent = 'Show log';
    logBtn.className = 'ghost';
    logBtn.dataset.action = 'show-log';
    logBtn.dataset.pdbId = target.pdb_id || '';
    if (target.target_job_id) {
      logBtn.dataset.jobId = target.target_job_id;
      logBtn.dataset.logTitle = `BoltzGen config run · ${target.pdb_id}`;
    } else {
      logBtn.disabled = true;
    }
    cmdCell.appendChild(logBtn);
    const manualBtn = document.createElement('button');
    manualBtn.type = 'button';
    manualBtn.textContent = 'Show run command';
    manualBtn.className = 'ghost';
    manualBtn.dataset.action = 'show-run';
    manualBtn.dataset.pdbId = target.pdb_id || '';
    manualBtn.dataset.scope = 'target';
    cmdCell.appendChild(manualBtn);

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
        epitopeLabel(cfg, cfgIdx + 1),
        bindingTextFromConfig(cfg),
      ];
      epCells.forEach((val, cellIdx) => {
        const td = document.createElement('td');
        td.textContent = val || '';
        if (cellIdx === 3) td.style.paddingLeft = '14px';
        epRow.appendChild(td);
      });
      const epStatus = document.createElement('td');
      epStatus.className = 'boltz-status-cell';
      epStatus.appendChild(buildLocalResultBadge(hasLocalBoltzResults(target.pdb_id)));
      epRow.appendChild(epStatus);

      const epCmd = document.createElement('td');
      epCmd.className = 'boltz-cmd-cell';
      const epRun = document.createElement('button');
      epRun.type = 'button';
      epRun.textContent = 'Run';
      epRun.dataset.action = 'run-epitope';
      epRun.dataset.pdbId = target.pdb_id || '';
      epRun.dataset.configPath = cfg.config_path || '';
      epCmd.appendChild(epRun);

      const epCfgBtn = document.createElement('button');
      epCfgBtn.type = 'button';
      epCfgBtn.textContent = 'Show config';
      epCfgBtn.className = 'ghost';
      epCfgBtn.dataset.action = 'show-config-epitope';
      epCfgBtn.dataset.pdbId = target.pdb_id || '';
      epCfgBtn.dataset.configPath = cfg.config_path || '';
      epCmd.appendChild(epCfgBtn);

      const epLog = document.createElement('button');
      epLog.type = 'button';
      epLog.textContent = 'Show log';
      epLog.className = 'ghost';
      epLog.dataset.action = 'show-log';
      epLog.dataset.pdbId = target.pdb_id || '';
      if (cfg.job_id) {
        epLog.dataset.jobId = cfg.job_id;
        epLog.dataset.logTitle = `${epitopeLabel(cfg, cfgIdx + 1)} · ${target.pdb_id}`;
      } else {
        epLog.disabled = true;
      }
      epCmd.appendChild(epLog);
      const epManual = document.createElement('button');
      epManual.type = 'button';
      epManual.textContent = 'Show run command';
      epManual.className = 'ghost';
      epManual.dataset.action = 'show-run';
      epManual.dataset.pdbId = target.pdb_id || '';
      epManual.dataset.configPath = cfg.config_path || '';
      epManual.dataset.epitopeName = epitopeLabel(cfg, cfgIdx + 1);
      epCmd.appendChild(epManual);

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
  if (!ids.length) {
    state.boltzLocalRuns = results;
    return results;
  }
  await Promise.all(ids.map(async (id) => {
    try {
      const res = await fetch(`/api/targets/${id}/boltzgen/runs`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const body = await res.json();
      const runs = Array.isArray(body.runs) ? body.runs : [];
      const hasLocal = runs.some((run) => {
        if (Array.isArray(run.specs) && run.specs.some((spec) => spec && spec.has_metrics)) return true;
        return Boolean(run.local_path);
      });
      results[id] = hasLocal;
    } catch (err) {
      results[id] = false;
    }
  }));
  state.boltzLocalRuns = results;
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
    el.boltzPanel.hidden = true;
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
    if (el.boltzPanel) el.boltzPanel.hidden = true;
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
  const toolsRoot = workspaceRoot || '<remote_root>';
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

function buildRunCommandText(pdbId, specPaths = [], epitopeName = null) {
  const upper = (pdbId || '').toUpperCase();
  const {
    remoteRoot, targetRoot, toolsRoot, sshTarget, condaActivate, boltz, localTargetsRoot,
  } = clusterDefaults();
  const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${upper}`.replace(/\/+/g, '/');
  const localTarget = localTargetsRoot ? `${localTargetsRoot}/${upper}` : `targets/${upper}`;
  const runLabel = runLabelFor(upper);
  const designCount = getBoltzDesignCount() || 100;
  const timeHours = getBoltzTimeHours() || boltz.time_hours || 12;
  const partition = boltz.partition || '<partition>';
  const gpus = boltz.gpus || '<gpus>';
  const cpus = boltz.cpus || '<cpus>';
  const memVal = boltz.mem_gb ?? '<mem>';
  const memLabel = Number.isFinite(memVal) ? `${memVal}G` : `${memVal}`;
  const outputRoot = boltz.output_root
    ? `${boltz.output_root}`.replace(/\/$/, '')
    : `${remoteTarget}/designs/boltzgen`;
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
  const pipelinePath = `${remoteRoot || '<remote_root>'}/tools/boltzgen/pipeline.py`;
  return [
    `# Manual BoltzGen submission${epitopeName ? ` (${epitopeName})` : ''}`,
    '',
    '# 1) Sync target configs + tools to cluster',
    `rsync -az ${localTarget} ${sshTarget}:${targetRoot}`,
    `rsync -az tools ${sshTarget}:${toolsRoot || '<remote_root>'}/tools`,
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
    `rsync -az ${sshTarget}:${outputRoot}/ ${localTarget}/designs/boltzgen/`,
  ].filter(Boolean).join('\n');
}

function buildAllRunCommandsText(targets = []) {
  const {
    remoteRoot, targetRoot, toolsRoot, sshTarget, condaActivate, boltz, localTargetsRoot,
  } = clusterDefaults();
  const envRoot = remoteRoot || '<remote_root>';
  const envTargetRoot = targetRoot || `${envRoot}/targets`;
  const pipelinePath = `${remoteRoot || '<remote_root>'}/tools/boltzgen/pipeline.py`;
  const condaLine = condaActivate ? condaActivate : '# conda activate bg';

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
  lines.push(`rsync -az tools ${sshTarget}:${toolsRoot || '<remote_root>'}/tools`);

  lines.push('', '# 2) Verify boltzgen configs exist on cluster');
  targets.forEach((t) => {
    const pdb = (t?.pdb_id || '').toUpperCase();
    if (!pdb) return;
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const specs = Array.isArray(t?.configs) ? t.configs.map((cfg) => cfg.config_path).filter(Boolean) : [];
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
    const pdb = (t?.pdb_id || '').toUpperCase();
    if (!pdb) return;
    const runLabel = runLabelFor(pdb);
    const remoteTarget = `${targetRoot || '<remote_root>/targets'}/${pdb}`.replace(/\/+/g, '/');
    const outputRoot = boltz.output_root
      ? `${boltz.output_root}`.replace(/\/$/, '')
      : `${remoteTarget}/designs/boltzgen`;
    const specs = Array.isArray(t?.configs) ? t.configs.map((cfg) => cfg.config_path).filter(Boolean) : [];
    const remoteSpecs = (specs.length ? specs : ['configs/*/boltzgen_config.yaml'])
      .map((spec) => `${remoteTarget}/${spec}`.replace(/\/+/g, '/'));
    const specLines = remoteSpecs.map((spec) => `  --spec ${spec} \\`);
    const cacheLine = cacheDir ? `  --cache_dir ${cacheDir} \\` : null;
    const extraLine = extraArgs.length ? `  --extra_run_args ${extraArgs.join(' ')} \\` : null;
    const localTarget = localTargetsRoot ? `${localTargetsRoot}/${pdb}` : `targets/${pdb}`;

    lines.push(
      '',
      `# Target ${pdb}`,
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
    pullLines.push(`rsync -az ${sshTarget}:${outputRoot}/ ${localTarget}/designs/boltzgen/`);
  });

  if (monitorLines.length) {
    lines.push('', '# Monitor', ...monitorLines);
  }
  if (pullLines.length) {
    lines.push('', '# Pull results back', ...pullLines);
  }

  return lines.filter(Boolean).join('\n');
}

async function showRunCommand(pdbId, configPath = null, epitopeName = null) {
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
    specs = [configPath];
  } else if (target && Array.isArray(target.configs)) {
    specs = target.configs.map((cfg) => cfg.config_path).filter(Boolean);
  }
  const text = buildRunCommandText(pdbId, specs, epitopeName);
  if (el.boltzRunBody) el.boltzRunBody.textContent = text;
}

async function showRunCommandAll() {
  if (el.boltzRunBody) el.boltzRunBody.textContent = 'Preparing cluster commands...';
  if (el.boltzRunTitle) el.boltzRunTitle.textContent = 'Manual run · all targets';
  toggleModal(el.boltzRunModal, true);
  if (!state.clusterStatus) {
    await loadClusterStatus();
  }
  const targets = Array.isArray(state.boltzConfigs) ? state.boltzConfigs : [];
  if (!targets.length) {
    if (el.boltzRunBody) el.boltzRunBody.textContent = '# No BoltzGen targets detected.';
    return;
  }
  const text = buildAllRunCommandsText(targets);
  if (el.boltzRunBody) el.boltzRunBody.textContent = text;
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
    loadBinderTable({ page: current - 1, silent: true });
  } else if (dir === 'next' && current < totalPages) {
    loadBinderTable({ page: current + 1, silent: true });
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
    showBoltzConfig(pdbId, null);
  } else if (action === 'show-config-epitope') {
    showBoltzConfig(pdbId, configPath);
  } else if (action === 'show-log') {
    showBoltzLog(btn.dataset.jobId, btn.dataset.logTitle || 'Job log');
  } else if (action === 'show-run') {
    showRunCommand(pdbId, configPath, btn.dataset.epitopeName || null);
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
    const rawVal = el.bulkThrottle?.value;
    const raw = Number(rawVal);
    if (!Number.isFinite(raw)) return 3;
    if (raw < 0) return 0;
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
    llmDelaySeconds: 70,
    decideScopeAttempts: 3,
  });
  scrollToSnapshots();
}

function init() {
  loadClusterStatus();
  if (el.bulkPreviewBtn) el.bulkPreviewBtn.addEventListener('click', () => previewBulkCsv({ silent: false }));
  if (el.bulkPreviewRefresh) el.bulkPreviewRefresh.addEventListener('click', () => previewBulkCsv({ silent: true }));
  if (el.bulkVisualizeEpitopes) el.bulkVisualizeEpitopes.addEventListener('click', handleVisualizeEpitopes);
  if (el.bulkRunBtn) el.bulkRunBtn.addEventListener('click', () => startBulkRun());
  if (el.snapshotDownloadPdf) el.snapshotDownloadPdf.addEventListener('click', downloadSnapshotPdf);
  if (el.boltzRefresh) el.boltzRefresh.addEventListener('click', () => loadBoltzConfigs({ silent: true }));
  if (el.diversityRefresh) el.diversityRefresh.addEventListener('click', () => refreshDiversity());
  if (el.diversityDownloadCsv) el.diversityDownloadCsv.addEventListener('click', () => downloadDiversityFile('csv'));
  if (el.diversityDownloadHtml) el.diversityDownloadHtml.addEventListener('click', () => downloadDiversityFile('html'));
  if (el.boltzShowRunAll) el.boltzShowRunAll.addEventListener('click', showRunCommandAll);
  if (el.boltzTable) el.boltzTable.addEventListener('click', handleBoltzTableClick);
  if (el.binderTable) el.binderTable.addEventListener('click', handleBinderTableClick);
  if (el.binderPagination) el.binderPagination.addEventListener('click', handleBinderPagination);
  if (el.binderRefresh) el.binderRefresh.addEventListener('click', () => loadBinderTable({ page: state.binderPage || 1, force: true }));
  if (el.binderDownload) el.binderDownload.addEventListener('click', downloadBinderCsv);
  if (el.boltzConfigClose) el.boltzConfigClose.addEventListener('click', () => toggleModal(el.boltzConfigModal, false));
  if (el.boltzLogClose) el.boltzLogClose.addEventListener('click', () => toggleModal(el.boltzLogModal, false));
  if (el.boltzRunClose) el.boltzRunClose.addEventListener('click', () => toggleModal(el.boltzRunModal, false));
  if (el.epitopeReportDownload) el.epitopeReportDownload.addEventListener('click', downloadEpitopeReportHtml);
  document.addEventListener('click', (evt) => {
    const closeTarget = evt.target?.dataset?.close;
    if (closeTarget === 'boltz-config') toggleModal(el.boltzConfigModal, false);
    if (closeTarget === 'boltz-log') toggleModal(el.boltzLogModal, false);
    if (closeTarget === 'boltz-run') toggleModal(el.boltzRunModal, false);
  });

  refreshDiversity({ silent: true });
}

document.addEventListener('DOMContentLoaded', init);
