import pathlib
import sys
import types
import pytest

sys.path.append(str(pathlib.Path(__file__).resolve().parents[1]))

# Provide lightweight stubs so decide_scope imports cleanly in the test env.
sys.modules.setdefault("requests", types.SimpleNamespace())
sys.modules.setdefault("google", types.SimpleNamespace())
sys.modules.setdefault("google.generativeai", types.SimpleNamespace())
sys.modules.setdefault(
    "yaml",
    types.SimpleNamespace(
        safe_load=lambda *args, **kwargs: {}, safe_dump=lambda *args, **kwargs: ""
    ),
)
sys.modules.setdefault(
    "jsonschema", types.SimpleNamespace(validate=lambda *args, **kwargs: None)
)
sys.modules.setdefault(
    "utils",
    types.SimpleNamespace(
        _ensure_dir=lambda *args, **kwargs: None,
        ROOT=pathlib.Path("."),
        TARGETS_ROOT_LOCAL=pathlib.Path("."),
        SCHEMA={},
        RCSB_ENTRY="",
        RCSB_ASSEM="",
        RCSB_PDB="",
        UNIPROT_IDMAPPING_RUN_API="",
        UNIPROT_IDMAPPING_STATUS_API="",
        UNIPROT_API="",
    ),
)
sys.modules.setdefault(
    "env",
    types.SimpleNamespace(GOOGLE_API_KEY="", MODEL="", USE_LLM=False, GROQ_API_KEY=""),
)

from decide_scope import _ensure_epitopes_within_target_chains


def test_ensure_epitopes_trims_invalid_chains_and_residues(capsys):
    cfg = {
        "chains": ["A", "B", "X"],
        "epitopes": [
            {"name": "keep", "residues": ["A:1-3"]},
            {"name": "trim", "residues": ["B:2-6", "B:10-11"]},
        ],
    }

    allowed = {"A", "B"}
    valid_residues = {"A": {1, 2, 3}, "B": {2, 3, 4, 10}}

    _ensure_epitopes_within_target_chains(cfg, allowed, valid_residue_numbers=valid_residues)

    assert cfg["chains"] == ["A", "B"]

    trimmed_epitope = next(ep for ep in cfg["epitopes"] if ep["name"] == "trim")
    assert trimmed_epitope["residues"] == ["B:2-4", "B:10"]

    kept_epitope = next(ep for ep in cfg["epitopes"] if ep["name"] == "keep")
    assert kept_epitope["residues"] == ["A:1-3"]

    out = capsys.readouterr().out
    assert "Removing chain(s) [X]" in out
