"""Pydantic models shared across FastAPI endpoints."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Literal, Optional

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


class AssessmentSyncResponse(BaseModel):
    job_id: str
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


class AlignmentResponse(BaseModel):
    pdb_id: str
    antigen_url: Optional[str]
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
