(function () {
  const catalogList = document.querySelector('#catalog-list');
  if (!catalogList) {
    return;
  }

  const state = {
    files: [],
    selectedFile: null,
    jobId: null,
    pollingTimer: null,
    previewText: '',
    sidebarWidth: null,
  };

  const catalogViewer = document.querySelector('#catalog-viewer');
  const selectedNameEl = document.querySelector('#catalog-selected-name');
  const selectedInfoEl = document.querySelector('#catalog-selected-info');
  const previewLimitInput = document.querySelector('#preview-limit');
  const refreshBtn = document.querySelector('#catalog-refresh');
  const downloadBtn = document.querySelector('#catalog-download');
  const copyBtn = document.querySelector('#catalog-copy');
  const avoidSelect = document.querySelector('#tg-avoid');
  const layoutEl = document.querySelector('#tg-catalog-layout');
  const resizerEl = document.querySelector('#catalog-resizer');

  const form = document.querySelector('#target-generation-form');
  const instructionInput = document.querySelector('#tg-instruction');
  const maxTargetsInput = document.querySelector('#tg-max-targets');
  const speciesInput = document.querySelector('#tg-species');
  const preferTagsInput = document.querySelector('#tg-prefer-tags');
  const prefixInput = document.querySelector('#tg-prefix');
  const extraInput = document.querySelector('#tg-extra');
  const noBrowserInput = document.querySelector('#tg-no-browser');
  const submitBtn = document.querySelector('#tg-submit');
  const statusBadge = document.querySelector('#tg-status');

  const jobIdEl = document.querySelector('#tg-job-id');
  const jobStatusEl = document.querySelector('#tg-job-status');
  const jobMessageEl = document.querySelector('#tg-job-message');
  const jobLogEl = document.querySelector('#tg-log');
  const refreshStatusBtn = document.querySelector('#tg-refresh');
  const clearLogBtn = document.querySelector('#tg-clear-log');
  const SIDEBAR_WIDTH_KEY = 'tgCatalogSidebarWidth';
  const SIDEBAR_MIN = 240;
  const SIDEBAR_MAX = 640;

  function applySidebarWidth(width, persist = false) {
    if (!layoutEl) {
      return;
    }
    const clamped = Math.max(SIDEBAR_MIN, Math.min(SIDEBAR_MAX, width));
    layoutEl.style.setProperty('--catalog-sidebar-width', `${clamped}px`);
    state.sidebarWidth = clamped;
    if (persist && window.localStorage) {
      try {
        window.localStorage.setItem(SIDEBAR_WIDTH_KEY, String(Math.round(clamped)));
      } catch (err) {
        console.warn('Unable to persist catalog width', err);
      }
    }
  }

  if (layoutEl && window.localStorage) {
    try {
      const savedWidth = Number.parseInt(window.localStorage.getItem(SIDEBAR_WIDTH_KEY) || '', 10);
      if (Number.isFinite(savedWidth)) {
        applySidebarWidth(savedWidth);
      }
    } catch (err) {
      console.warn('Unable to read saved catalog width', err);
    }
  }

  if (layoutEl && resizerEl) {
    let resizing = false;
    let startX = 0;
    let startWidth = 0;

    const stopResize = () => {
      if (!resizing) {
        return;
      }
      resizing = false;
      resizerEl.classList.remove('dragging');
      document.body.style.userSelect = '';
      document.body.style.cursor = '';
      document.removeEventListener('pointermove', handlePointerMove);
      document.removeEventListener('pointerup', stopResize);
      if (state.sidebarWidth) {
        applySidebarWidth(state.sidebarWidth, true);
      }
    };

    const handlePointerMove = (event) => {
      if (!resizing) {
        return;
      }
      const delta = event.clientX - startX;
      applySidebarWidth(startWidth + delta);
    };

    resizerEl.addEventListener('pointerdown', (event) => {
      if (window.innerWidth <= 960) {
        return;
      }
      const sidebar = layoutEl.querySelector('.catalog-column');
      if (!sidebar) {
        return;
      }
      event.preventDefault();
      resizing = true;
      startX = event.clientX;
      startWidth = sidebar.getBoundingClientRect().width;
      resizerEl.classList.add('dragging');
      document.body.style.userSelect = 'none';
      document.body.style.cursor = 'col-resize';
      document.addEventListener('pointermove', handlePointerMove);
      document.addEventListener('pointerup', stopResize);
    });

    resizerEl.addEventListener('keydown', (event) => {
      if (!['ArrowLeft', 'ArrowRight'].includes(event.key)) {
        return;
      }
      event.preventDefault();
      const sidebar = layoutEl.querySelector('.catalog-column');
      const currentWidth = state.sidebarWidth || (sidebar ? sidebar.getBoundingClientRect().width : SIDEBAR_MIN);
      const delta = event.key === 'ArrowLeft' ? -20 : 20;
      applySidebarWidth(currentWidth + delta, true);
    });
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

  function setStatusBadge(text, visible) {
    if (!statusBadge) {
      return;
    }
    if (visible) {
      statusBadge.hidden = false;
      statusBadge.textContent = text;
    } else {
      statusBadge.hidden = true;
    }
  }

  function clearViewer(resetPreview = true) {
    catalogViewer.innerHTML = '';
    if (resetPreview) {
      state.previewText = '';
      if (copyBtn) copyBtn.disabled = true;
    }
  }

  function renderTable(headers, rows, truncated) {
    const lines = [];
    if (headers && headers.length) {
      lines.push(headers.join('\t'));
    }
    rows.forEach((row) => {
      lines.push(row.join('\t'));
    });
    state.previewText = lines.join('\n');

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
    rows.forEach((row) => {
      const tr = document.createElement('tr');
      row.forEach((cell) => {
        const td = document.createElement('td');
        td.textContent = cell;
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);

    if (truncated) {
      const caption = document.createElement('caption');
      caption.textContent = 'Preview truncated. Increase the row limit to view more records.';
      table.appendChild(caption);
    }

    clearViewer(false);
    catalogViewer.appendChild(table);
    if (copyBtn) {
      copyBtn.disabled = !state.previewText;
    }
  }

  function renderCatalogList() {
    catalogList.innerHTML = '';
    if (!state.files.length) {
      const empty = document.createElement('li');
      empty.className = 'catalog-empty';
      empty.textContent = 'No catalog files found yet.';
      catalogList.appendChild(empty);
      downloadBtn.disabled = true;
      return;
    }

    state.files.forEach((file) => {
      const item = document.createElement('li');
      item.className = 'catalog-item';
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = file.name;
      button.dataset.filename = file.name;
      button.className = 'catalog-button';
      const meta = document.createElement('span');
      meta.className = 'catalog-item-meta';
      meta.textContent = `${formatBytes(file.size_bytes)} • ${formatTimestamp(file.modified_at)}`;
      if (state.selectedFile === file.name) {
        button.classList.add('active');
      }
      button.addEventListener('click', () => {
        selectFile(file.name);
      });
      item.appendChild(button);
      item.appendChild(meta);
      catalogList.appendChild(item);
    });
  }

  function populateAvoidOptions() {
    if (!avoidSelect) {
      return;
    }
    const currentSelection = new Set(Array.from(avoidSelect.selectedOptions).map((opt) => opt.value));
    avoidSelect.innerHTML = '';
    const tsvFiles = state.files.filter(
      (file) => file && typeof file.name === 'string' && file.name.toLowerCase().endsWith('.tsv'),
    );
    tsvFiles.forEach((file) => {
      const option = document.createElement('option');
      option.value = file.name;
      option.textContent = file.name;
      if (currentSelection.has(file.name)) {
        option.selected = true;
      }
      avoidSelect.appendChild(option);
    });
  }

  async function loadCatalog() {
    try {
      refreshBtn.disabled = true;
      const response = await fetch('/api/target-generation/catalog');
      if (!response.ok) {
        throw new Error(`Failed to load catalog (${response.status})`);
      }
      const payload = await response.json();
      state.files = Array.isArray(payload.files) ? payload.files : [];
      renderCatalogList();
      populateAvoidOptions();
      if (state.selectedFile && !state.files.some((file) => file.name === state.selectedFile)) {
        state.selectedFile = null;
        clearViewer();
        selectedNameEl.textContent = 'Select a file to preview';
        selectedInfoEl.hidden = true;
        downloadBtn.disabled = true;
      }
    } catch (err) {
      console.error(err);
    } finally {
      refreshBtn.disabled = false;
    }
  }

  async function selectFile(name) {
    const limit = Number.parseInt(previewLimitInput.value, 10);
    selectedNameEl.textContent = name;
    selectedInfoEl.hidden = true;
    downloadBtn.disabled = true;
    state.selectedFile = name;
    clearViewer();
    const spinner = document.createElement('div');
    spinner.className = 'catalog-spinner';
    spinner.textContent = 'Loading preview…';
    catalogViewer.appendChild(spinner);

    try {
      const response = await fetch(`/api/target-generation/catalog/${encodeURIComponent(name)}?limit=${Number.isFinite(limit) ? limit : 200}`);
      if (!response.ok) {
        throw new Error(`Preview failed (${response.status})`);
      }
      const payload = await response.json();
      renderTable(payload.headers || [], payload.rows || [], Boolean(payload.truncated));
      const extra = [];
      if (Number.isFinite(payload.total_rows)) {
        extra.push(`${payload.total_rows} rows`);
      }
      if (Number.isFinite(payload.displayed_rows)) {
        extra.push(`${payload.displayed_rows} shown`);
      }
      if (payload.filtered_by_biotin) {
        extra.push('biotinylated only');
      }
      if (payload.deduped_by_gene) {
        extra.push('unique genes');
      }
      if (extra.length) {
        selectedInfoEl.textContent = extra.join(' • ');
        selectedInfoEl.hidden = false;
      } else {
        selectedInfoEl.hidden = true;
      }
      if (downloadBtn) downloadBtn.disabled = false;
      if (copyBtn) {
        copyBtn.disabled = !state.previewText;
      }
    } catch (err) {
      console.error(err);
      clearViewer();
      const errorBox = document.createElement('div');
      errorBox.className = 'catalog-error';
      errorBox.textContent = err.message || 'Failed to load preview.';
      catalogViewer.appendChild(errorBox);
      if (downloadBtn) downloadBtn.disabled = true;
      if (copyBtn) copyBtn.disabled = true;
    }
  }

  async function submitRun(event) {
    event.preventDefault();
    if (!instructionInput.value.trim()) {
      instructionInput.focus();
      return;
    }

    const payload = {
      instruction: instructionInput.value.trim(),
      no_browser_popup: Boolean(noBrowserInput?.checked),
    };

    const maxTargets = Number.parseInt(maxTargetsInput.value, 10);
    if (Number.isFinite(maxTargets)) {
      payload.max_targets = maxTargets;
    }
    const species = speciesInput.value.trim();
    if (species) {
      payload.species = species;
    }
    const preferTags = preferTagsInput.value.trim();
    if (preferTags) {
      payload.prefer_tags = preferTags;
    }
    const prefix = prefixInput.value.trim();
    if (prefix) {
      payload.out_prefix = prefix;
    }
    const extra = extraInput.value.trim();
    if (extra) {
      payload.extra_args = extra;
    }
    if (avoidSelect && avoidSelect.selectedOptions.length) {
      payload.avoid_existing = Array.from(avoidSelect.selectedOptions).map((opt) => opt.value);
    }

    try {
      submitBtn.disabled = true;
      setStatusBadge('Submitting…', true);
      const response = await fetch('/api/target-generation/run', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.detail || `Request failed (${response.status})`);
      }
      const data = await response.json();
      state.jobId = data.job_id;
      jobIdEl.textContent = data.job_id;
      jobStatusEl.textContent = 'pending';
      jobMessageEl.textContent = data.message || 'Job queued';
      jobLogEl.textContent = '';
      refreshStatusBtn.disabled = false;
      startPolling();
    } catch (err) {
      console.error(err);
      jobMessageEl.textContent = err.message || 'Failed to submit job';
    } finally {
      submitBtn.disabled = false;
      setStatusBadge('Submitted', false);
    }
  }

  async function fetchJobStatus() {
    if (!state.jobId) {
      return;
    }
    try {
      const response = await fetch(`/api/jobs/${state.jobId}`);
      if (!response.ok) {
        throw new Error(`Job status failed (${response.status})`);
      }
      const data = await response.json();
      jobStatusEl.textContent = data.status;
      jobMessageEl.textContent = data.message || '—';
      if (Array.isArray(data.logs)) {
        jobLogEl.textContent = data.logs.join('\n');
        jobLogEl.scrollTop = jobLogEl.scrollHeight;
      }
      if (['success', 'failed', 'canceled'].includes(data.status)) {
        stopPolling();
      }
    } catch (err) {
      console.error(err);
    }
  }

  function startPolling() {
    stopPolling();
    fetchJobStatus();
    state.pollingTimer = window.setInterval(fetchJobStatus, 4000);
  }

  function stopPolling() {
    if (state.pollingTimer) {
      window.clearInterval(state.pollingTimer);
      state.pollingTimer = null;
    }
  }

  function handleDownload() {
    if (!state.selectedFile) {
      return;
    }
    window.open(`/api/target-generation/catalog/${encodeURIComponent(state.selectedFile)}/file`, '_blank');
  }

  async function handleCopy() {
    if (!state.previewText) {
      return;
    }
    try {
      await navigator.clipboard.writeText(state.previewText);
      selectedInfoEl.textContent = 'Preview copied to clipboard';
      selectedInfoEl.hidden = false;
    } catch (err) {
      console.error(err);
    }
  }

  refreshBtn.addEventListener('click', () => {
    loadCatalog();
  });

  downloadBtn.addEventListener('click', handleDownload);
  if (copyBtn) {
    copyBtn.addEventListener('click', handleCopy);
  }
  previewLimitInput.addEventListener('change', () => {
    if (state.selectedFile) {
      selectFile(state.selectedFile);
    }
  });

  form.addEventListener('submit', submitRun);
  refreshStatusBtn.addEventListener('click', fetchJobStatus);
  clearLogBtn.addEventListener('click', () => {
    jobLogEl.textContent = '';
  });

  window.addEventListener('beforeunload', () => {
    stopPolling();
  });

  loadCatalog();
})();
