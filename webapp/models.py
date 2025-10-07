"""Pydantic models shared across FastAPI endpoints."""

from __future__ import annotations

from typing import Dict, List, Literal, Optional

from pydantic import BaseModel, Field


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


class AlignmentResponse(BaseModel):
    pdb_id: str
    antigen_url: Optional[str]
    vendor_sequence_length: int
    chain_results: List[Dict[str, object]]


class DesignRunRequest(BaseModel):
    pdb_id: str
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


class DesignRunResponse(BaseModel):
    job_id: str
    message: str


class RankingRow(BaseModel):
    index: int
    design_name: str
    iptm: Optional[float]
    rmsd_diego: Optional[float]
    tm_score: Optional[float]
    ipsae_min: Optional[float]
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
