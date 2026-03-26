import csv
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from fastapi.testclient import TestClient

from webapp import bulk as bulk_mod
from webapp import main as main_mod
from webapp.models import (
    BoltzgenBinderExportRequest,
    BoltzgenBinderRow,
    BoltzgenDiversityResponse,
)


def _patch_bulk_config(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    targets_dir = tmp_path / "targets"
    targets_dir.mkdir(parents=True, exist_ok=True)
    cfg = SimpleNamespace(
        paths=SimpleNamespace(
            targets_dir=targets_dir,
            workspace_root=tmp_path,
            project_root=tmp_path,
            cache_dir=tmp_path / "cache",
            antigen_category_map=None,
        ),
        log_dir=tmp_path / "logs",
        cluster=SimpleNamespace(
            rfantibody=SimpleNamespace(binder_chain_id="H"),
        ),
    )
    monkeypatch.setattr(bulk_mod, "load_config", lambda: cfg)
    return targets_dir


def _write_target_with_metrics(targets_dir: Path, pdb_id: str, *, iptm: float, rmsd: float) -> None:
    target_dir = targets_dir / pdb_id.upper()
    target_dir.mkdir(parents=True, exist_ok=True)
    (target_dir / "target.yaml").write_text(
        yaml.safe_dump(
            {
                "target_name": f"Target {pdb_id.upper()}",
                "sequences": {"pdb": {"A": "ACDEFG"}, "pdb_residue_numbers": {"A": [1, 2, 3, 4, 5, 6]}},
                "epitopes": [{"name": "epitope_1", "residues": ["A:1-3"], "hotspots": ["A:2"]}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    metrics_dir = target_dir / "designs" / "boltzgen" / "epitope_1" / "20260101_0000" / "final_ranked_designs"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = metrics_dir / "all_designs_metrics.csv"
    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["file_name", "designed_chain_sequence", "filter_rmsd", "design_to_target_iptm", "rank"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "file_name": "rank_001_model.cif",
                "designed_chain_sequence": "ACDEFG",
                "filter_rmsd": rmsd,
                "design_to_target_iptm": iptm,
                "rank": 1,
            }
        )


def _patch_diversity_lightweight(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bulk_mod, "build_epitope_diversity_plots", lambda **kwargs: ([], None, 0))
    monkeypatch.setattr(bulk_mod, "_collect_rfa_scatter_metrics", lambda *args, **kwargs: ({}, {}))
    monkeypatch.setattr(bulk_mod, "_collect_rfa_design_rows", lambda *args, **kwargs: [])


def test_build_diversity_report_scopes_to_requested_pdb_ids(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_bulk_config(monkeypatch, tmp_path)
    _patch_diversity_lightweight(monkeypatch)
    _write_target_with_metrics(targets_dir, "1AAA", iptm=0.82, rmsd=1.2)
    _write_target_with_metrics(targets_dir, "2BBB", iptm=0.74, rmsd=1.8)
    _write_target_with_metrics(targets_dir, "3CCC", iptm=0.65, rmsd=2.4)

    response = bulk_mod.build_boltzgen_diversity_report(
        include_binders=True,
        binder_page=1,
        binder_page_size=50,
        pdb_ids=["1AAA", "2BBB"],
        epitope_min_designs=1,
    )

    assert response.csv_name is not None
    metrics_pdbs = {item["pdb_id"] for item in (response.metrics_files or []) if item.get("pdb_id")}
    binder_pdbs = {row.pdb_id for row in (response.binder_rows or [])}
    assert metrics_pdbs <= {"1AAA", "2BBB"}
    assert binder_pdbs <= {"1AAA", "2BBB"}
    assert "3CCC" not in metrics_pdbs
    assert "3CCC" not in binder_pdbs


def test_diversity_cache_respects_scope(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    targets_dir = _patch_bulk_config(monkeypatch, tmp_path)
    _patch_diversity_lightweight(monkeypatch)
    _write_target_with_metrics(targets_dir, "1AAA", iptm=0.82, rmsd=1.2)
    _write_target_with_metrics(targets_dir, "2BBB", iptm=0.74, rmsd=1.8)

    original = bulk_mod._response_from_diversity_cache
    cache_hits = {"count": 0}

    def _wrapped_cache(*args, **kwargs):
        cache_hits["count"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(bulk_mod, "_response_from_diversity_cache", _wrapped_cache)

    first = bulk_mod.build_boltzgen_diversity_report(pdb_ids=["1AAA"], epitope_min_designs=1)
    assert first.csv_name is not None

    second = bulk_mod.build_boltzgen_diversity_report(pdb_ids=["1AAA"], epitope_min_designs=1)
    assert second.csv_name is not None
    assert cache_hits["count"] == 1

    third = bulk_mod.build_boltzgen_diversity_report(pdb_ids=["2BBB"], epitope_min_designs=1)
    assert third.csv_name is not None
    assert cache_hits["count"] == 1
    third_pdbs = {item["pdb_id"] for item in (third.metrics_files or []) if item.get("pdb_id")}
    assert third_pdbs <= {"2BBB"}


def test_list_binders_requires_cached_csv_without_rebuild(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    _patch_bulk_config(monkeypatch, tmp_path)
    monkeypatch.setattr(
        bulk_mod,
        "build_boltzgen_diversity_report",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("diversity build should not run")),
    )

    response = bulk_mod.list_boltzgen_binders(["1AAA"], page=1, page_size=20)

    assert response.rows == []
    assert response.total_rows == 0
    assert response.message is not None
    assert "Refresh" in response.message


def test_export_requires_cached_csv_without_rebuild(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    _patch_bulk_config(monkeypatch, tmp_path)
    monkeypatch.setattr(
        bulk_mod,
        "build_boltzgen_diversity_report",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("diversity build should not run")),
    )

    response = bulk_mod.export_selected_binders(
        BoltzgenBinderExportRequest(
            selections=["1AAA:epitope_1"],
            per_group=1,
            include_summary=False,
        )
    )

    assert response.csv_name is None
    assert response.message is not None
    assert "Refresh" in response.message


def test_diversity_refresh_endpoint_passes_scope_pdb_ids(monkeypatch: pytest.MonkeyPatch):
    captured: dict = {}

    def _fake_build(**kwargs):
        captured["pdb_ids"] = kwargs.get("pdb_ids")
        return BoltzgenDiversityResponse(
            csv_name="all_design_metrics_20260101_000000.csv",
            output_dir="/tmp",
            message="ok",
            binder_rows=[
                BoltzgenBinderRow(
                    pdb_id="1AAA",
                    epitope="epitope_1",
                    epitope_id="epitope_1",
                    rank=1,
                    iptm=0.8,
                    rmsd=1.2,
                    engine="boltzgen",
                )
            ],
            binder_total=1,
            binder_page=1,
            binder_page_size=100,
            binder_message="Showing 1 of 1 binders",
            binder_counts={"1AAA": 1},
        )

    monkeypatch.setattr(main_mod, "build_boltzgen_diversity_report", _fake_build)
    client = TestClient(main_mod.app)
    response = client.post(
        "/api/bulk/boltzgen/diversity/refresh",
        params={"pdb_ids": "1AAA,2BBB", "page": 1, "page_size": 100},
    )
    assert response.status_code == 200
    payload = response.json()
    assert captured.get("pdb_ids") == ["1AAA", "2BBB"]
    assert payload.get("binder_total") == 1
    assert payload.get("binder_rows", [])[0].get("pdb_id") == "1AAA"
