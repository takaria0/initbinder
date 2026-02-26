from pathlib import Path

from fastapi.testclient import TestClient

from webapp.bulk import (
    _append_generated_row_to_catalog,
    _build_discovery_instruction,
    _select_best_generated_row,
    discover_unmatched_bulk_target,
)
from webapp.main import app
from webapp.models import BulkLlmCandidate, BulkLlmUnmatchedDiscoverRequest


def test_build_discovery_instruction_uses_candidate_fields():
    candidate = BulkLlmCandidate(
        target_name="PD-L1",
        gene="CD274",
        uniprot="Q9NZQ7",
        antigen_catalog="10084-H08H-B",
    )
    instruction = _build_discovery_instruction(candidate)
    assert "CD274" in instruction
    assert "Q9NZQ7" in instruction
    assert "10084-H08H-B" in instruction


def test_select_best_generated_row_prefers_exact_identity():
    candidate = BulkLlmCandidate(gene="TNFRSF9", uniprot="Q07011")
    rows = [
        {
            "gene": "RANDOM1",
            "uniprot": "P00001",
            "chosen_pdb": "1AAA",
            "antigen_catalog": "X-1",
        },
        {
            "gene": "TNFRSF9",
            "uniprot": "Q07011",
            "chosen_pdb": "8GYE",
            "antigen_catalog": "10041-H08H-B",
        },
    ]
    chosen = _select_best_generated_row(candidate, rows)
    assert chosen["gene"] == "TNFRSF9"
    assert chosen["uniprot"] == "Q07011"


def test_append_generated_row_to_catalog_dedupes_and_increments_rank(tmp_path: Path):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "rank\tselection\tuniprot\tgene\tprotein_name\tchosen_pdb\tvendor_accession\tantigen_catalog\tantigen_url\n"
        "1\tbiotin\tQ07011\tTNFRSF9\tTNFRSF9 receptor\t8GYE\tQ07011\t10041-H08H-B\thttps://example.org/a\n",
        encoding="utf-8",
    )

    duplicate_row = {
        "selection": "any",
        "uniprot": "Q07011",
        "gene": "TNFRSF9",
        "protein_name": "TNFRSF9 receptor",
        "chosen_pdb": "8GYE",
        "vendor_accession": "Q07011",
        "antigen_catalog": "10041-H08H-B",
        "antigen_url": "https://example.org/a",
    }
    assert _append_generated_row_to_catalog(catalog, duplicate_row) is False

    new_row = {
        "selection": "any",
        "uniprot": "Q9NZQ7",
        "gene": "CD274",
        "protein_name": "PD-L1",
        "chosen_pdb": "5O45",
        "vendor_accession": "Q9NZQ7",
        "antigen_catalog": "10084-H08H-B",
        "antigen_url": "https://example.org/b",
    }
    assert _append_generated_row_to_catalog(catalog, new_row) is True

    lines = catalog.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3
    assert lines[-1].startswith("2\t")
    assert "CD274" in lines[-1]


def test_discover_unmatched_bulk_target_rematches_and_appends(tmp_path: Path, monkeypatch):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "rank\tselection\tuniprot\tgene\tprotein_name\tchosen_pdb\tvendor_accession\tantigen_catalog\tantigen_url\n"
        "1\tbiotin\tQ07011\tTNFRSF9\tTNFRSF9 receptor\t8GYE\tQ07011\t10041-H08H-B\thttps://example.org/a\n",
        encoding="utf-8",
    )

    monkeypatch.setattr("webapp.bulk._catalog_dir", lambda: tmp_path)

    def fake_run_target_generation_discovery(
        *,
        instruction,
        out_prefix,
        max_targets,
        launch_browser=True,
        log_hook=None,
    ):
        generated = tmp_path / f"{out_prefix}_all.tsv"
        generated.write_text(
            "rank\tselection\tuniprot\tgene\tprotein_name\tchosen_pdb\tvendor_accession\tantigen_catalog\tantigen_url\n"
            "1\tany\tQ9NZQ7\tCD274\tPD-L1\t5O45\tQ9NZQ7\t10084-H08H-B\thttps://example.org/b\n",
            encoding="utf-8",
        )

    monkeypatch.setattr("webapp.bulk._run_target_generation_discovery", fake_run_target_generation_discovery)

    request = BulkLlmUnmatchedDiscoverRequest(
        catalog_name="catalog.tsv",
        unmatched_key="unm_cd274",
        candidate=BulkLlmCandidate(gene="CD274", uniprot="Q9NZQ7", target_name="PD-L1"),
        max_targets=3,
    )
    result = discover_unmatched_bulk_target(request, job_id="job_test_1")
    matched_row = result["matched_row"]

    assert matched_row.resolved_pdb_id == "5O45"
    assert result["catalog_appended"] is True

    lines = catalog.read_text(encoding="utf-8").splitlines()
    assert any("CD274" in line and "5O45" in line for line in lines)


def test_unmatched_discover_endpoint_requires_existing_catalog():
    client = TestClient(app)
    payload = {
        "catalog_name": "missing_catalog.tsv",
        "unmatched_key": "unm_missing",
        "candidate": {"target_name": "PD-L1", "gene": "CD274"},
        "max_targets": 3,
    }
    response = client.post("/api/bulk/llm-targets/unmatched/discover", json=payload)
    assert response.status_code == 404
