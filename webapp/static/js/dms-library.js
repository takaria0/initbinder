const state = {
  resultId: null,
  downloadUrl: null,
};

const form = document.querySelector('#dms-form');
const statusBadge = document.querySelector('#dms-status');
const resultsCard = document.querySelector('#dms-results-card');
const summaryText = document.querySelector('#dms-summary');
const metaList = document.querySelector('#dms-meta');
const previewBody = document.querySelector('#dms-preview-table tbody');
const downloadLink = document.querySelector('#dms-download');
const pymolButton = document.querySelector('#dms-launch-pymol');
const logCard = document.querySelector('#dms-log-card');
const logArea = document.querySelector('#dms-log');
const mutatedContainer = document.querySelector('#dms-mutated-residues');

if (downloadLink) {
  downloadLink.classList.add('disabled');
  downloadLink.href = '#';
}

if (pymolButton) {
  pymolButton.disabled = true;
}

function setStatus(message, tone = 'idle') {
  if (!statusBadge) return;
  if (!message) {
    statusBadge.hidden = true;
    statusBadge.textContent = '';
    statusBadge.classList.remove('warn', 'error', 'success');
    return;
  }
  statusBadge.hidden = false;
  statusBadge.textContent = message;
  statusBadge.classList.remove('warn', 'error', 'success');
  if (tone === 'warn') statusBadge.classList.add('warn');
  if (tone === 'error') statusBadge.classList.add('error');
  if (tone === 'success') statusBadge.classList.add('success');
}

function clearLog() {
  if (logArea) {
    logArea.textContent = '';
  }
  if (logCard) {
    logCard.hidden = true;
  }
}

function appendLog(message) {
  if (!logArea || !logCard) return;
  const timestamp = new Date().toLocaleTimeString();
  logArea.textContent += `[${timestamp}] ${message}\n`;
  logCard.hidden = false;
  logArea.scrollTop = logArea.scrollHeight;
}

function createMetaRow(label, value) {
  const dt = document.createElement('dt');
  dt.textContent = label;
  const dd = document.createElement('dd');
  dd.textContent = value;
  metaList.append(dt, dd);
}

function renderMeta(data) {
  if (!metaList) return;
  metaList.innerHTML = '';
  createMetaRow('PDB ID', data.pdb_id || '—');
  createMetaRow('Chain', data.chain_id);
  createMetaRow('PDB path', data.pdb_path);
  createMetaRow('Total variants', data.total_variants.toLocaleString());
  createMetaRow('Residues mutated', data.candidate_residue_count.toLocaleString());
  createMetaRow('Surface residues ≥ threshold', data.surface_residue_count.toLocaleString());
  createMetaRow('Sequence length', data.sequence_length.toLocaleString());
  createMetaRow('Mutation menu', data.mutation_kind);
  createMetaRow('RSA threshold', data.rsa_threshold.toFixed(2));
  createMetaRow('Surface only', data.target_surface_only ? 'Yes' : 'No');
  const expressedState = data.restrict_to_expressed_region
    ? data.expressed_region_applied
      ? 'Applied'
      : 'Requested (not applied)'
    : 'Disabled';
  createMetaRow('Vendor expressed filter', expressedState);
  if (data.restrict_to_expressed_region) {
    createMetaRow('Vendor expressed range', data.expressed_region_vendor_range || '—');
    if (typeof data.expressed_region_sequence_length === 'number') {
      createMetaRow('Vendor expressed length', data.expressed_region_sequence_length.toLocaleString());
    }
    createMetaRow(
      'Residues overlapping vendor range',
      data.expressed_region_matched_residues.toLocaleString(),
    );
  }
}

function renderPreview(rows) {
  if (!previewBody) return;
  previewBody.innerHTML = '';
  rows.forEach((row) => {
    const tr = document.createElement('tr');
    const residueLabel = `${row.pdb_resnum}${row.icode || ''}`;
    tr.innerHTML = `
      <td>${residueLabel}</td>
      <td>${row.wt}</td>
      <td>${row.mut}</td>
      <td>${row.category}</td>
      <td>${row.rsa.toFixed(3)}</td>
      <td>${row.barcode_18nt || ''}</td>
    `;
    previewBody.appendChild(tr);
  });
}

function renderMutatedResidues(residues, surfaceCount) {
  if (!mutatedContainer) return;
  mutatedContainer.innerHTML = '';
  const header = document.createElement('p');
  header.className = 'mutated-summary-header';
  header.textContent = `Mutated residues (${residues.length} of ${surfaceCount} surface residues meeting RSA threshold)`;
  mutatedContainer.appendChild(header);

  if (!residues.length) {
    const empty = document.createElement('p');
    empty.className = 'help-text';
    empty.textContent = 'No residues satisfied the filter criteria. Lower the RSA threshold or disable the surface filter.';
    mutatedContainer.appendChild(empty);
    return;
  }

  const list = document.createElement('ul');
  list.className = 'chip-row';
  residues.forEach((res) => {
    const li = document.createElement('li');
    li.className = 'chip';
    const label = `${res.uid} (${res.wt}, RSA ${res.rsa.toFixed(2)})`;
    const categories = res.categories && res.categories.length ? res.categories.join(', ') : 'SSM';
    li.innerHTML = `<strong>${label}</strong><span>${categories}</span>`;
    list.appendChild(li);
  });
  mutatedContainer.appendChild(list);
}

async function handleSubmit(event) {
  event.preventDefault();
  clearLog();
  setStatus('Designing…', 'warn');
  pymolButton.disabled = true;
  downloadLink.classList.add('disabled');
  downloadLink.href = '#';
  resultsCard.hidden = true;

  const formData = new FormData(form);
  const payload = {
    pdb_id: (formData.get('pdb_id') || '').toString().trim().toUpperCase() || null,
    pdb_path: (formData.get('pdb_path') || '').toString().trim() || null,
    chain_id: (formData.get('chain_id') || '').toString().trim().toUpperCase(),
    mutation_kind: (formData.get('mutation_kind') || 'SSM').toString(),
    target_surface_only: document.querySelector('#dms-surface-only')?.checked ?? true,
    restrict_to_expressed_region: document.querySelector('#dms-vendor-range')?.checked ?? false,
    rsa_threshold: parseFloat(formData.get('rsa_threshold')) || 0.25,
    include_glycan_toggles: document.querySelector('#dms-glycan')?.checked ?? true,
    add_conservative_swaps: document.querySelector('#dms-conservative')?.checked ?? true,
    add_controls: document.querySelector('#dms-controls')?.checked ?? true,
    preview_limit: parseInt(formData.get('preview_limit'), 10) || 200,
  };

  const barcodeRaw = formData.get('add_barcodes');
  if (barcodeRaw) {
    const barcodeCount = parseInt(barcodeRaw, 10);
    if (Number.isFinite(barcodeCount) && barcodeCount > 0) {
      payload.add_barcodes = barcodeCount;
    }
  }

  try {
    const response = await fetch('/api/dms-library/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      const detail = errorData.detail || response.statusText;
      throw new Error(detail);
    }

    const data = await response.json();
    state.resultId = data.result_id;
    state.downloadUrl = data.download_url;

    summaryText.textContent = data.message;
    renderMeta(data);
    renderPreview(data.preview);
    renderMutatedResidues(data.mutated_residues, data.surface_residue_count);
    if (Array.isArray(data.expressed_region_notes) && data.expressed_region_notes.length) {
      data.expressed_region_notes.forEach((note) => appendLog(note));
    }

    downloadLink.href = data.download_url;
    downloadLink.classList.remove('disabled');
    pymolButton.disabled = false;
    resultsCard.hidden = false;

    setStatus('Ready', 'success');
  } catch (error) {
    console.error('Failed to generate DMS library', error);
    appendLog(error.message || String(error));
    setStatus('Failed', 'error');
  }
}

async function handleLaunchPyMol() {
  if (!state.resultId) return;
  setStatus('Preparing PyMOL…', 'warn');
  pymolButton.disabled = true;
  try {
    const response = await fetch(`/api/dms-library/${encodeURIComponent(state.resultId)}/pymol`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ launch: true, bundle_only: false }),
    });
    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      const detail = errorData.detail || response.statusText;
      throw new Error(detail);
    }
    const data = await response.json();
    appendLog(`PyMOL bundle: ${data.session_path}`);
    setStatus(data.launched ? 'PyMOL launched' : 'Bundle ready', 'success');
  } catch (error) {
    console.error('PyMOL launch failed', error);
    appendLog(error.message || String(error));
    setStatus('PyMOL failed', 'error');
  } finally {
    pymolButton.disabled = false;
  }
}

if (form) {
  form.addEventListener('submit', handleSubmit);
}

if (pymolButton) {
  pymolButton.addEventListener('click', handleLaunchPyMol);
}
