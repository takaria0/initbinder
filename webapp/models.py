"""Pydantic models shared across FastAPI endpoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field, validator


class TargetPreset(BaseModel):
    id: str
    name: str
    pdb_id: str
    antigen_url: Optional[str] = None
    num_epitopes: Optional[int] = Field(None, ge=1, le=32)
    created_at: float
    updated_at: float
    last_used: float


class TargetPresetRequest(BaseModel):
    name: Optional[str] = Field(None, max_length=120)
    pdb_id: str = Field(..., pattern=r"^[0-9A-Za-z]{4}$")
    antigen_url: Optional[str] = None
    num_epitopes: Optional[int] = Field(None, ge=1, le=32)
    preset_id: Optional[str] = None


class TargetPresetResponse(BaseModel):
    preset: TargetPreset


class TargetPresetListResponse(BaseModel):
    presets: List[TargetPreset]


class TargetInitRequest(BaseModel):
    pdb_id: str = Field(..., pattern=r"^[0-9A-Za-z]{4}$", description="4-character PDB accession")
    antigen_url: Optional[str] = Field(None, description="Vendor product URL")
    target_accession: str = Field(..., max_length=40, description="UniProt/RefSeq accession for target")
    target_vendor_range: Optional[str] = Field(
        None,
        max_length=64,
        description="Optional vendor expressed range (e.g., 1-167)",
    )
    preset_name: Optional[str] = Field(None, max_length=120, description="Friendly name for saved target preset")
    num_epitopes: Optional[int] = Field(None, ge=1, le=32, description="Desired number of epitopes to surface")
    decide_scope_prompt: Optional[str] = Field(
        None,
        max_length=2000,
        description="Optional natural language guidance for decide-scope epitope selection",
    )
    force_refresh: bool = Field(False, description="Redo init even if target folder exists")
    run_decide_scope: bool = Field(True, description="Also run decide-scope after initialization")
    run_prep: bool = Field(True, description="Also run prep-target after decide-scope")


class TargetInitResponse(BaseModel):
    job_id: str
    message: str


class JobStatusResponse(BaseModel):
    job_id: str
    kind: str
    label: str
    status: Literal["pending", "running", "success", "failed", "canceled"]
    message: Optional[str] = None
    progress: Optional[float] = None
    logs: List[str] = Field(default_factory=list)
    details: Dict[str, object] = Field(default_factory=dict)


class JobSummary(BaseModel):
    job_id: str
    kind: str
    label: str
    status: Literal["pending", "running", "success", "failed", "canceled"]
    created_at: float
    started_at: Optional[float] = None
    finished_at: Optional[float] = None
    pdb_id: Optional[str] = None
    run_label: Optional[str] = None


class AssessmentRunSummary(BaseModel):
    run_label: str
    updated_at: float
    rankings_path: Optional[str] = None
    total_rows: Optional[int] = None
    available_local: bool = False
    available_remote: bool = False
    local_path: Optional[str] = None
    remote_path: Optional[str] = None
    origin: Literal["local", "remote", "both"] = "local"


class AssessmentRunRequest(BaseModel):
    pdb_id: str
    run_label: str
    binder_chain_id: str = Field("H", min_length=1, max_length=1)
    include_keyword: Optional[str] = None


class AssessmentRunResponse(BaseModel):
    job_id: str
    message: str


class BoltzGenSyncResponse(BaseModel):
    message: str
    run_label: Optional[str] = None


class TargetCatalogFile(BaseModel):
    name: str
    size_bytes: int
    modified_at: float


class TargetCatalogListResponse(BaseModel):
    directory: str
    files: List[TargetCatalogFile] = Field(default_factory=list)


class TargetCatalogPreviewResponse(BaseModel):
    name: str
    headers: List[str] = Field(default_factory=list)
    rows: List[List[str]] = Field(default_factory=list)
    total_rows: int = 0
    displayed_rows: int = 0
    truncated: bool = False
    filtered_by_biotin: bool = False
    deduped_by_gene: bool = False


class TargetGenerationRequest(BaseModel):
    instruction: str = Field(..., min_length=3, max_length=4000)
    max_targets: Optional[int] = Field(None, ge=1, le=2000)
    species: Optional[str] = Field(None, max_length=120)
    prefer_tags: Optional[str] = Field(None, max_length=240)
    out_prefix: Optional[str] = Field(None, max_length=240)
    avoid_existing: List[str] = Field(default_factory=list)
    extra_args: Optional[str] = Field(None, max_length=500)
    no_browser_popup: bool = True

    @validator("avoid_existing", each_item=True)
    def _validate_avoid(cls, value: str) -> str:
        name = Path(value).name
        if name != value:
            raise ValueError("Invalid TSV file name")
        if not name.endswith(".tsv"):
            raise ValueError("Avoid list entries must end with .tsv")
        return name


class TargetGenerationResponse(BaseModel):
    job_id: str
    message: str


class BulkCsvRow(BaseModel):
    raw_index: int
    preset_name: str
    antigen_url: Optional[str] = None
    protein_name: Optional[str] = None
    expression_host: Optional[str] = Field(None, max_length=200)
    selection: Optional[str] = Field(None, max_length=64)
    biotinylated: Optional[bool] = None
    tags: Optional[str] = Field(None, max_length=500)
    pdb_id: Optional[str] = Field(None, max_length=32)
    accession: Optional[str] = Field(None, max_length=200)
    vendor_range: Optional[str] = Field(None, max_length=64)
    resolved_pdb_id: Optional[str] = None
    preset_id: Optional[str] = None
    warnings: List[str] = Field(default_factory=list)


class BulkDesignSettings(BaseModel):
    model_engine: Literal["rfantibody", "boltzgen"] = "boltzgen"
    total_designs: int = Field(90, ge=1, le=50000)
    num_sequences: int = Field(1, ge=1, le=32)
    temperature: float = Field(0.1, ge=0.0, le=1.0)
    binder_chain_id: Optional[str] = Field(None, min_length=1, max_length=1)
    af3_seed: int = Field(1, ge=0)
    run_assess: bool = True
    rfdiff_crop_radius: Optional[float] = Field(None, ge=0.0)
    run_label_prefix: Optional[str] = Field(None, max_length=80)
    boltz_binding: Optional[str] = Field(
        None,
        max_length=500,
        description="Optional BoltzGen binding_types override, e.g. A:142,151,229;B:10..30",
    )
    boltzgen_crop_radius: Optional[float] = Field(
        None,
        ge=0.0,
        description="Optional BoltzGen crop radius (Å) around hotspots for res_index.",
    )
    boltz_time_hours: Optional[int] = Field(
        None,
        ge=1,
        le=240,
        description="Optional per-job walltime (hours) for BoltzGen submissions",
    )


class BulkPreviewRequest(BaseModel):
    csv_text: str = Field(..., min_length=3, max_length=1000000)
    num_epitopes: Optional[int] = Field(None, ge=1, le=32)
    decide_scope_prompt: Optional[str] = Field(None, max_length=2000)


class BulkPreviewResponse(BaseModel):
    rows: List[BulkCsvRow] = Field(default_factory=list)
    total_rows: int
    resolved: int
    unresolved: int
    message: str


class BulkLlmMessage(BaseModel):
    role: Literal["user", "assistant"]
    content: str = Field(..., min_length=1, max_length=8000)


class BulkLlmCandidate(BaseModel):
    target_name: Optional[str] = Field(None, max_length=240)
    protein_name: Optional[str] = Field(None, max_length=240)
    gene: Optional[str] = Field(None, max_length=80)
    uniprot: Optional[str] = Field(None, max_length=80)
    pdb_id: Optional[str] = Field(None, max_length=32)
    antigen_catalog: Optional[str] = Field(None, max_length=120)
    accession: Optional[str] = Field(None, max_length=120)
    rationale: Optional[str] = Field(None, max_length=500)


class BulkCatalogMatch(BaseModel):
    candidate: BulkLlmCandidate
    row: BulkCsvRow
    match_type: str = Field(..., max_length=64)
    matched_field: Optional[str] = Field(None, max_length=64)
    matched_value: Optional[str] = Field(None, max_length=240)
    confidence: float = Field(0.0, ge=0.0, le=1.0)


class BulkUnmatchedSuggestion(BaseModel):
    candidate: BulkLlmCandidate
    reason: str = Field(..., max_length=240)
    nearest: List[BulkCatalogMatch] = Field(default_factory=list)


class BulkLlmTargetSuggestRequest(BaseModel):
    catalog_name: str = Field(..., min_length=3, max_length=255)
    prompt: str = Field(..., min_length=3, max_length=4000)
    history: List[BulkLlmMessage] = Field(default_factory=list)
    max_candidates: int = Field(24, ge=1, le=200)


class BulkLlmTargetSuggestResponse(BaseModel):
    catalog_name: str
    assistant_message: str
    candidates: List[BulkLlmCandidate] = Field(default_factory=list)
    matched_rows: List[BulkCsvRow] = Field(default_factory=list)
    matches: List[BulkCatalogMatch] = Field(default_factory=list)
    unmatched: List[BulkUnmatchedSuggestion] = Field(default_factory=list)
    message: str


class BulkLlmUnmatchedDiscoverRequest(BaseModel):
    catalog_name: str = Field(..., min_length=3, max_length=255)
    unmatched_key: str = Field(..., min_length=3, max_length=160)
    candidate: BulkLlmCandidate
    history: List[BulkLlmMessage] = Field(default_factory=list)
    vendor_scope: Literal["both", "sino", "acro"] = "both"
    planning_mode: Literal["balanced", "aggressive", "fast"] = "balanced"
    max_targets: int = Field(3, ge=1, le=10)
    launch_browser: bool = True


class BulkLlmUnmatchedDiscoverResponse(BaseModel):
    job_id: str
    unmatched_key: str
    message: str


class BulkLlmUnmatchedDiscoverStatusResponse(BaseModel):
    job_id: str
    status: Literal["pending", "running", "success", "failed", "canceled"]
    phase: Optional[
        Literal[
            "planning",
            "catalog_rematch",
            "vendor_search",
            "parsing",
            "appending",
            "rematching",
            "success",
            "failed",
        ]
    ] = None
    message: Optional[str] = None
    unmatched_key: Optional[str] = None
    catalog_name: Optional[str] = None
    resolved_species: Optional[str] = None
    planned_queries: List[str] = Field(default_factory=list)
    vendors_consulted: List[str] = Field(default_factory=list)
    llm_plan_summary: Optional[str] = None
    attempts: List[Dict[str, Any]] = Field(default_factory=list)
    matched_row: Optional[BulkCsvRow] = None
    match: Optional[BulkCatalogMatch] = None
    failure_reason: Optional[str] = None


class BulkRunRequest(BaseModel):
    csv_text: str = Field(..., min_length=3, max_length=1000000)
    num_epitopes: Optional[int] = Field(None, ge=1, le=32)
    decide_scope_prompt: Optional[str] = Field(None, max_length=2000)
    launch_pymol: bool = True
    render_pymol_snapshots: Optional[bool] = None
    export_insights: bool = True
    export_designs: bool = True
    submit_designs: bool = False
    prepare_targets: bool = True
    force_init: bool = False
    design_settings: BulkDesignSettings = Field(default_factory=BulkDesignSettings)
    throttle_seconds: float = Field(3.0, ge=0.0, le=300.0)
    limit: Optional[int] = Field(None, ge=1, le=2000)
    llm_delay_seconds: float = Field(0.0, ge=0.0, le=600.0, description="Optional delay before each LLM decide-scope call")
    decide_scope_attempts: int = Field(1, ge=1, le=5, description="Number of decide-scope retries")


class BulkDesignImportRequest(BaseModel):
    csv_text: str = Field(..., min_length=3, max_length=200000)
    throttle_seconds: float = Field(3.0, ge=0.0, le=300.0)


class BulkRunResponse(BaseModel):
    job_id: str
    message: str


class BoltzgenDiversityPlot(BaseModel):
    pdb_id: str
    png_name: Optional[str] = None
    png_path: Optional[str] = None
    svg_name: Optional[str] = None
    svg_path: Optional[str] = None
    epitope_colors: Dict[str, str] = Field(default_factory=dict)


class BoltzgenMetricsFile(BaseModel):
    pdb_id: str
    epitope_name: str
    datetime: str
    path: str


class BoltzgenDiversityResponse(BaseModel):
    csv_name: Optional[str] = None
    output_dir: Optional[str] = None
    html_name: Optional[str] = None
    plots: List[BoltzgenDiversityPlot] = Field(default_factory=list)
    metrics_files: List[BoltzgenMetricsFile] = Field(default_factory=list)
    message: Optional[str] = None
    binder_rows: List[BoltzgenBinderRow] = Field(default_factory=list)
    binder_total: int = 0
    binder_page: int = 1
    binder_page_size: int = 0
    binder_message: Optional[str] = None
    binder_counts: Dict[str, int] = Field(default_factory=dict)


class AntigenDiversityPlot(BaseModel):
    name: str
    svg_name: str
    svg_path: str


class AntigenDiversityRequest(BaseModel):
    pdb_ids: List[str] = Field(default_factory=list)


class AntigenDiversityResponse(BaseModel):
    output_dir: Optional[str] = None
    plots: List[AntigenDiversityPlot] = Field(default_factory=list)
    message: Optional[str] = None


class EpitopeDiversityPlot(BaseModel):
    title: str
    png_name: Optional[str] = None
    png_path: Optional[str] = None
    svg_name: Optional[str] = None
    svg_path: Optional[str] = None


class EpitopeDiversityRequest(BaseModel):
    selections: List[str] = Field(
        default_factory=list,
        description="Entries like 5WT9:epitope_1, 3J8F:4",
    )


class EpitopeDiversityResponse(BaseModel):
    output_dir: Optional[str] = None
    csv_name: Optional[str] = None
    hotspot_csv_name: Optional[str] = None
    plots: List[EpitopeDiversityPlot] = Field(default_factory=list)
    message: Optional[str] = None


class BoltzgenBinderRow(BaseModel):
    pdb_id: str
    epitope: Optional[str] = None
    epitope_id: Optional[str] = None
    engine: Optional[str] = None
    rank: int
    iptm: Optional[float] = None
    rmsd: Optional[float] = None
    hotspot_dist: Optional[float] = None
    ipsae_min: Optional[float] = None
    binder_seq: Optional[str] = None
    design_path: Optional[str] = None
    metrics_path: Optional[str] = None
    run_label: Optional[str] = None
    config_path: Optional[str] = None
    binding_label: Optional[str] = None
    include_label: Optional[str] = None
    target_path: Optional[str] = None


class BoltzgenBinderResponse(BaseModel):
    rows: List[BoltzgenBinderRow] = Field(default_factory=list)
    total_rows: int = 0
    page: int = 1
    page_size: int = 100
    csv_name: Optional[str] = None
    message: Optional[str] = None


class BoltzgenBinderExportRequest(BaseModel):
    selections: List[str] = Field(
        default_factory=list,
        description="PDB:epitope entries (e.g. 5WT9:epitope_1).",
    )
    per_group: int = Field(48, ge=1, description="Binders to export per antigen:epitope group.")
    include_summary: bool = Field(True, description="Whether to emit a summary CSV.")
    upstream_flank: str = Field(
        "GGAG",
        min_length=1,
        description="Deprecated compatibility field; ignored by binder export adapter generation.",
    )
    downstream_flank: str = Field(
        "CGCT",
        min_length=1,
        description="Deprecated compatibility field; ignored by binder export adapter generation.",
    )


class BoltzgenBinderExportPlot(BaseModel):
    engine: str
    png_name: Optional[str] = None
    svg_name: Optional[str] = None
    map_csv_name: Optional[str] = None
    point_count: int = 0
    skipped_missing_metrics: int = 0


class BoltzgenBinderExportResponse(BaseModel):
    csv_name: Optional[str] = None
    summary_csv_name: Optional[str] = None
    selection_count: int = 0
    invalid: List[str] = Field(default_factory=list)
    plot_exports: List[BoltzgenBinderExportPlot] = Field(default_factory=list)
    message: Optional[str] = None


class PipelineRefreshRequest(BaseModel):
    pdb_id: str
    force: bool = True
    expected_epitopes: Optional[int] = 3
    decide_scope_attempts: int = Field(3, ge=1, le=5)
    decide_scope_prompt: Optional[str] = Field(None, max_length=2000)
    antigen_url: Optional[str] = None
    target_accession: Optional[str] = Field(None, max_length=40)
    target_vendor_range: Optional[str] = Field(None, max_length=64)
    design_count: Optional[int] = Field(
        None,
        ge=1,
        le=50000,
        description="Designs per epitope to store in BoltzGen configs after prep",
    )


class PipelineRefreshResponse(BaseModel):
    job_id: str
    message: str


class BulkCommandBoltzgenDefaults(BaseModel):
    partition: Optional[str] = None
    account: Optional[str] = None
    gpus: Optional[str] = None
    cpus: Optional[int] = None
    mem_gb: Optional[int] = None
    time_hours: Optional[int] = None
    cache_dir: Optional[str] = None
    output_root: Optional[str] = None
    conda_activate: Optional[str] = None
    extra_args: List[str] = Field(default_factory=list)


class BulkCommandDefaultsResponse(BaseModel):
    ssh_target: Optional[str] = None
    remote_root: Optional[str] = None
    target_root: Optional[str] = None
    local_root: Optional[str] = None
    conda_activate: Optional[str] = None
    boltzgen: BulkCommandBoltzgenDefaults = Field(default_factory=BulkCommandBoltzgenDefaults)


class BulkUiClusterConfig(BaseModel):
    mock: bool = False
    ssh_config_alias: Optional[str] = None
    remote_root: Optional[str] = None
    target_root: Optional[str] = None
    conda_activate: Optional[str] = None
    pymol_path: Optional[str] = None
    pymol_conda_env: Optional[str] = None


class BulkUiBoltzgenConfig(BaseModel):
    partition: Optional[str] = None
    account: Optional[str] = None
    gpus: Optional[str] = None
    cpus: Optional[int] = Field(None, ge=1, le=256)
    mem_gb: Optional[int] = Field(None, ge=1, le=2048)
    time_hours: Optional[int] = Field(None, ge=1, le=240)
    default_num_designs: Optional[int] = Field(None, ge=1, le=50000)
    nanobody_scaffolds: List[str] = Field(default_factory=list)


class BulkUiInputConfig(BaseModel):
    default_input_path: Optional[str] = None
    auto_load_default_input: bool = False


class BulkUiLlmConfig(BaseModel):
    openai_api_key: Optional[str] = None
    openai_model: Optional[str] = None


class BulkUiConfigResponse(BaseModel):
    local_config_path: str
    cluster: BulkUiClusterConfig = Field(default_factory=BulkUiClusterConfig)
    boltzgen: BulkUiBoltzgenConfig = Field(default_factory=BulkUiBoltzgenConfig)
    input: BulkUiInputConfig = Field(default_factory=BulkUiInputConfig)
    llm: BulkUiLlmConfig = Field(default_factory=BulkUiLlmConfig)


class BulkUiConfigUpdateRequest(BaseModel):
    cluster: BulkUiClusterConfig = Field(default_factory=BulkUiClusterConfig)
    boltzgen: BulkUiBoltzgenConfig = Field(default_factory=BulkUiBoltzgenConfig)
    input: BulkUiInputConfig = Field(default_factory=BulkUiInputConfig)
    llm: BulkUiLlmConfig = Field(default_factory=BulkUiLlmConfig)


class BulkDefaultInputResponse(BaseModel):
    path: str
    size_bytes: int
    text: str


class BulkGuiReadmeResponse(BaseModel):
    path: str
    text: str


class BoltzgenEpitopeConfig(BaseModel):
    epitope_id: Optional[str] = None
    epitope_name: Optional[str] = None
    config_path: str
    binding_label: Optional[str] = None
    include_label: Optional[str] = None
    hotspot_count: Optional[int] = None
    hotspot_surface_ok: Optional[bool] = None
    hotspot_surface_exposed_count: Optional[int] = None
    hotspot_surface_total: Optional[int] = None
    hotspot_surface_missing: Optional[int] = None
    hotspot_surface_cutoff: Optional[float] = None
    job_id: Optional[str] = None
    job_status: Optional[str] = None
    run_label: Optional[str] = None
    submitted_at: Optional[float] = None
    parent_job_id: Optional[str] = None


class BoltzgenTargetConfig(BaseModel):
    pdb_id: str
    preset_name: Optional[str] = None
    configs: List[BoltzgenEpitopeConfig] = Field(default_factory=list)
    target_job_id: Optional[str] = None
    target_job_status: Optional[str] = None
    antigen_url: Optional[str] = None
    has_prep: Optional[bool] = None
    antigen_expressed_range: Optional[str] = None
    antigen_expressed_length: Optional[int] = None
    allowed_epitope_range: Optional[str] = None
    allowed_epitope_length: Optional[int] = None
    epitope_count: Optional[int] = None


class BoltzgenConfigListResponse(BaseModel):
    targets: List[BoltzgenTargetConfig] = Field(default_factory=list)


class BoltzgenConfigContent(BaseModel):
    pdb_id: str
    config_path: str
    epitope_name: Optional[str] = None
    yaml_text: str


class RfaPipelineScript(BaseModel):
    path: str
    name: str
    stage: Optional[str] = None
    epitope: Optional[str] = None
    variant: Optional[str] = None


class RfaPipelineTargetScripts(BaseModel):
    pdb_id: str
    scripts: List[RfaPipelineScript] = Field(default_factory=list)
    launcher_path: Optional[str] = None
    launcher_name: Optional[str] = None


class RfaPipelineConfigListResponse(BaseModel):
    targets: List[RfaPipelineTargetScripts] = Field(default_factory=list)


class RfaPipelineScriptContent(BaseModel):
    pdb_id: str
    script_path: str
    script_name: str
    script_text: str


class TargetYamlContent(BaseModel):
    pdb_id: str
    path: str
    yaml_text: str


class BoltzgenConfigRunRequest(BaseModel):
    pdb_id: str
    design_count: int = Field(90, ge=1, le=50000, description="Designs to generate per epitope")
    config_path: Optional[str] = Field(
        None,
        description="Optional relative path to a boltzgen_config.yaml to restrict run to one epitope",
    )
    run_label_prefix: Optional[str] = Field(None, max_length=80)
    throttle_seconds: float = Field(0.0, ge=0.0, le=120.0)
    time_hours: Optional[int] = Field(
        None,
        ge=1,
        le=240,
        description="Optional walltime (hours) for BoltzGen config runs",
    )


class BoltzgenConfigRunResponse(BaseModel):
    job_id: str
    message: str


class BoltzgenConfigRegenerateRequest(BaseModel):
    pdb_ids: List[str] = Field(default_factory=list, description="PDB IDs to regenerate configs for")
    design_count: int = Field(100, ge=1, le=50000, description="Designs per epitope to store in configs")
    boltzgen_crop_radius: Optional[float] = Field(
        None,
        ge=0.0,
        description="Optional BoltzGen target crop radius (Å) around hotspots for res_index.",
    )


class BoltzgenConfigRegenerateResult(BaseModel):
    pdb_id: str
    status: Literal["ok", "skipped", "error"] = "ok"
    configs_written: int = 0
    message: Optional[str] = None


class BoltzgenConfigRegenerateResponse(BaseModel):
    results: List[BoltzgenConfigRegenerateResult] = Field(default_factory=list)


class BoltzgenBinderPymolRequest(BaseModel):
    pdb_id: str
    design_path: str
    epitope_label: Optional[str] = None
    binding_label: Optional[str] = None
    include_label: Optional[str] = None
    target_path: Optional[str] = None
    config_path: Optional[str] = None


class BoltzgenBinderPymolResponse(BaseModel):
    session_dir: str
    script_path: str
    launched: bool


class AlignmentResponse(BaseModel):
    pdb_id: str
    antigen_url: Optional[str]
    vendor_range: Optional[List[int]] = None
    vendor_range_label: Optional[str] = None
    vendor_sequence_length: int
    chain_results: List[Dict[str, object]]


class DesignRunRequest(BaseModel):
    pdb_id: str
    model_engine: Literal["rfantibody", "boltzgen"] = Field(
        "rfantibody",
        description="Binder design engine to use (rfantibody default, boltzgen optional)",
    )
    total_designs: int = Field(90, ge=1, le=50000)
    num_sequences: int = Field(1, ge=1, le=32)
    temperature: float = Field(0.1, ge=0.0, le=1.0)
    binder_chain_id: Optional[str] = Field(
        None,
        min_length=1,
        max_length=1,
        description="Optional binder chain override; defaults to H when omitted",
    )
    af3_seed: int = Field(1, ge=0)
    run_label: Optional[str] = Field(None)
    run_assess: bool = Field(True, description="Auto-schedule assess-rfa-all after AF3")
    submit: bool = Field(False, description="Submit jobs to scheduler immediately")
    rfdiff_crop_radius: Optional[float] = Field(
        None,
        ge=0.0,
        description=(
            "Optional crop radius (Å) for RFdiffusion; when omitted the full prepared target is used."
        ),
    )
    boltzgen_crop_radius: Optional[float] = Field(
        None,
        ge=0.0,
        description="Optional crop radius (Å) around hotspots for BoltzGen res_index.",
    )
    boltz_binding: Optional[str] = Field(
        None,
        max_length=500,
        description="Optional chain:residue list for BoltzGen binding_types (e.g., A:142,151,229;B:10..30)",
    )
    boltz_time_hours: Optional[int] = Field(
        None,
        ge=1,
        le=240,
        description="Optional per-job walltime (hours) for BoltzGen runs",
    )


class DesignRunResponse(BaseModel):
    job_id: str
    message: str


class DesignEngineFieldInfo(BaseModel):
    field_id: str
    label: str
    description: Optional[str] = None
    visible: bool = True
    debug_only: bool = False


class DesignEngineInfo(BaseModel):
    engine_id: str
    label: str
    description: str
    is_default: bool = False
    fields: List[DesignEngineFieldInfo] = Field(default_factory=list)


class DesignEngineListResponse(BaseModel):
    engines: List[DesignEngineInfo]


class RankingRow(BaseModel):
    index: int
    design_name: str
    iptm: Optional[float]
    rmsd_diego: Optional[float]
    tm_score: Optional[float]
    ipsae_min: Optional[float]
    hotspot_min_distance: Optional[float]
    metadata: Dict[str, object] = Field(default_factory=dict)


class ScatterPoint(BaseModel):
    design_name: str
    iptm: Optional[float]
    rmsd_diego: Optional[float]
    ipsae_min: Optional[float]
    metadata: Dict[str, object] = Field(default_factory=dict)


class RankingResponse(BaseModel):
    pdb_id: str
    run_label: Optional[str]
    rows: List[RankingRow]
    scatter: List[ScatterPoint] = Field(default_factory=list)
    source_path: Optional[str] = None
    gallery_path: Optional[str] = None
    engine_id: Optional[str] = None


class RankingAnalysisRequest(BaseModel):
    run_label: Optional[str] = None


class RankingPlot(BaseModel):
    name: str
    title: str
    image_data: str


class SequenceSimilarityMatrix(BaseModel):
    designs: List[str]
    sequences: List[str]
    matrix: List[List[Optional[float]]]
    metric: str = "sequence_match_ratio"


class RankingAnalysisResponse(BaseModel):
    plots: List[RankingPlot] = Field(default_factory=list)
    similarity: Optional[SequenceSimilarityMatrix] = None
    logs: List[str] = Field(default_factory=list)


class BoltzGenSpecSummary(BaseModel):
    name: str
    has_metrics: bool = False
    metrics_path: Optional[str] = None
    design_count: Optional[int] = None


class BoltzGenRunSummary(BaseModel):
    run_label: str
    updated_at: float
    specs: List[BoltzGenSpecSummary] = Field(default_factory=list)
    local_path: Optional[str] = None


class BoltzGenRunListResponse(BaseModel):
    runs: List[BoltzGenRunSummary] = Field(default_factory=list)


class DMSLibraryDesignRow(BaseModel):
    chain: str
    uid: str
    pdb_resnum: int
    icode: Optional[str] = ""
    wt: str
    mut: str
    category: str
    rsa: float
    barcode_18nt: Optional[str] = None


class DMSResidueSummaryModel(BaseModel):
    uid: str
    pdb_resnum: int
    icode: Optional[str] = ""
    wt: str
    rsa: float
    categories: List[str] = Field(default_factory=list)


class AntigenDMSRequest(BaseModel):
    pdb_id: Optional[str] = Field(
        None,
        pattern=r"^[0-9A-Za-z]{4}$",
        description="Optional PDB identifier; required when pdb_path is omitted",
    )
    pdb_path: Optional[str] = Field(
        None,
        description="Optional override path to a prepared PDB file",
    )
    chain_id: str = Field(..., min_length=1, max_length=2)
    target_surface_only: bool = True
    restrict_to_expressed_region: bool = Field(
        False,
        description="Limit mutations to residues overlapping the vendor expressed/soluble construct",
    )
    rsa_threshold: float = Field(0.25, ge=0.0, le=1.0)
    mutation_kind: str = Field("SSM", description="Mutation menu to apply (SSM, alanine, charge)")
    include_glycan_toggles: bool = True
    add_conservative_swaps: bool = True
    add_controls: bool = True
    add_barcodes: Optional[int] = Field(None, ge=0, le=20000)
    preview_limit: int = Field(200, ge=1, le=1000)

    @validator("mutation_kind")
    def _normalize_menu(cls, value: str) -> str:
        menu = value.strip()
        if not menu:
            raise ValueError("mutation_kind cannot be empty")
        upper = menu.upper()
        mapping = {
            "SSM": "SSM",
            "ALANINE": "alanine",
            "CHARGE": "charge",
            "CHARGEFLIP": "charge",
        }
        if upper not in mapping:
            raise ValueError(f"Unsupported mutation_kind: {value}")
        return mapping[upper]

    @validator("pdb_path")
    def _trim_path(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        trimmed = value.strip()
        return trimmed or None

    @validator("chain_id")
    def _upper_chain(cls, value: str) -> str:
        return value.strip().upper()

    @validator("pdb_id")
    def _upper_pdb(cls, value: Optional[str]) -> Optional[str]:
        return value.upper() if value else value

    @validator("pdb_id", always=True)
    def _require_identifier(cls, value: Optional[str], values: Dict[str, object]) -> Optional[str]:
        if not value and not values.get("pdb_path"):
            raise ValueError("Either pdb_id or pdb_path must be provided")
        return value


class AntigenDMSResponse(BaseModel):
    result_id: str
    pdb_id: Optional[str]
    pdb_path: str
    chain_id: str
    total_variants: int
    preview: List[DMSLibraryDesignRow] = Field(default_factory=list)
    preview_count: int
    truncated: bool
    mutated_residues: List[DMSResidueSummaryModel] = Field(default_factory=list)
    surface_residue_count: int
    candidate_residue_count: int
    sequence_length: int
    target_surface_only: bool
    rsa_threshold: float
    mutation_kind: str
    restrict_to_expressed_region: bool
    expressed_region_applied: bool
    expressed_region_vendor_range: Optional[str]
    expressed_region_sequence_length: Optional[int]
    expressed_region_matched_residues: int
    expressed_region_notes: List[str] = Field(default_factory=list)
    download_url: str
    created_at: float
    message: str


class ExportRequest(BaseModel):
    pdb_id: str
    rankings_path: str
    top_n: int = Field(48, ge=1)
    codon_host: str = Field("yeast")
    gc_target: Optional[float] = Field(0.45, ge=0.0, le=1.0)
    gc_window: Optional[int] = Field(100, ge=1, le=1000)
    prefix_raw: Optional[str] = None
    suffix_raw: Optional[str] = None
    use_dnachisel: bool = Field(True)


class ExportResponse(BaseModel):
    job_id: str
    message: str


class GoldenGateRequest(BaseModel):
    pdb_id: str
    rankings_path: str
    top_n: int = Field(48, ge=1)
    codon_host: str = Field("yeast")
    sequence_column: Optional[str] = Field(
        None,
        description="Optional column override for binder amino acid sequence",
    )
    upstream_flank: str = Field("GGAG", min_length=2, description="5' flanking bases")
    downstream_flank: str = Field("CGCT", min_length=2, description="3' flanking bases")


class GoldenGateResponse(BaseModel):
    job_id: str
    message: str


class PyMolHotspotRequest(BaseModel):
    launch: bool = Field(True, description="Invoke PyMOL after bundle generation")
    bundle_only: bool = Field(False, description="Skip launching PyMOL; just return bundle path")
    epitope_name: Optional[str] = Field(None, description="Optional epitope label to highlight")


class PyMolHotspotResponse(BaseModel):
    bundle_path: Optional[str]
    launched: bool
    message: str


class PyMolDMSRequest(BaseModel):
    launch: bool = True
    bundle_only: bool = False


class PyMolDMSResponse(BaseModel):
    session_path: str
    script_path: str
    launched: bool
    message: str


class PyMolTopBindersRequest(BaseModel):
    run_label: Optional[str] = None
    top_n: int = Field(96, ge=1, le=500)
    launch: bool = Field(True)
    bundle_only: bool = Field(False)
    engine_id: Literal["rfantibody", "boltzgen"] = "rfantibody"
    spec: Optional[str] = None


class PyMolTopBindersResponse(BaseModel):
    bundle_path: Optional[str]
    launched: bool
    message: str


class PyMolGalleryMovieRequest(BaseModel):
    run_label: Optional[str] = None
    top_n: int = Field(96, ge=1, le=500)
    fps: int = Field(10, ge=1, le=120)
    interval_sec: float = Field(2.0, ge=0.1, le=60.0)
    rotation_deg_per_sec: float = Field(30.0, ge=0.0, le=360.0)
    rotation_axis: str = Field("y", pattern="^[xyzXYZ]$")
    desired_states: int = Field(48, ge=1, le=400)


class PyMolGalleryMovieResponse(BaseModel):
    bundle_path: Optional[str]
    movie_path: Optional[str]
    frames_prefix: Optional[str]
    frames_pattern: Optional[str]
    script_path: Optional[str]
    log_path: Optional[str]
    message: str
