from webapp.bulk import (
    _CatalogRecord,
    _build_catalog_lookup,
    _coerce_candidate,
    _match_candidates_to_catalog,
    preview_bulk_targets,
)
from webapp.models import BulkCsvRow, BulkLlmCandidate, BulkPreviewRequest


def _record(raw_index, *, preset_name, protein_name, pdb_id, accession, gene=None, uniprot=None, catalog=None):
    row = BulkCsvRow(
        raw_index=raw_index,
        preset_name=preset_name,
        protein_name=protein_name,
        pdb_id=pdb_id,
        resolved_pdb_id=pdb_id,
        accession=accession,
        antigen_url=None,
        vendor_range=None,
        warnings=[],
    )
    values = {}
    if gene:
        values["gene"] = gene
    if uniprot:
        values["uniprot"] = uniprot
    if catalog:
        values["antigen_catalog"] = catalog
    return _CatalogRecord(row=row, values=values, name_key=protein_name.lower())


def test_exact_match_precedence_over_fuzzy():
    records = [
        _record(
            1,
            preset_name="PD-L1",
            protein_name="Programmed cell death 1 ligand 1",
            pdb_id="5O45",
            accession="Q9NZQ7",
            gene="CD274",
            uniprot="Q9NZQ7",
            catalog="10084-H08H-B",
        ),
        _record(
            2,
            preset_name="CD137",
            protein_name="Tumor necrosis factor receptor superfamily member 9",
            pdb_id="8GYE",
            accession="Q07011",
            gene="TNFRSF9",
            uniprot="Q07011",
            catalog="10041-H08H-B",
        ),
    ]
    lookup = _build_catalog_lookup(records)
    candidate = BulkLlmCandidate(
        target_name="TNFRSF9 receptor target",
        protein_name="Tumor necrosis factor receptor 9",
        gene="TNFRSF9",
    )

    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, [candidate])
    assert not unmatched
    assert len(matches) == 1
    assert matches[0].match_type == "exact_gene"
    assert matches[0].row.resolved_pdb_id == "8GYE"
    assert len(matched_rows) == 1


def test_fuzzy_name_match_when_exact_fields_missing():
    records = [
        _record(
            1,
            preset_name="CD137",
            protein_name="Tumor necrosis factor receptor superfamily member 9",
            pdb_id="8GYE",
            accession="Q07011",
            gene="TNFRSF9",
        ),
    ]
    lookup = _build_catalog_lookup(records)
    candidate = BulkLlmCandidate(target_name="Tumor necrosis factor receptor member 9")
    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, [candidate])
    assert not unmatched
    assert len(matches) == 1
    assert matches[0].match_type == "fuzzy_name"
    assert matched_rows[0].resolved_pdb_id == "8GYE"


def test_unmatched_candidate_returns_reason():
    records = [
        _record(
            1,
            preset_name="CD137",
            protein_name="Tumor necrosis factor receptor superfamily member 9",
            pdb_id="8GYE",
            accession="Q07011",
            gene="TNFRSF9",
        ),
    ]
    lookup = _build_catalog_lookup(records)
    candidate = BulkLlmCandidate(target_name="Completely unrelated enzyme family target")
    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, [candidate])
    assert not matches
    assert not matched_rows
    assert len(unmatched) == 1
    assert "No catalog match" in unmatched[0].reason


def test_candidate_coercion_normalizes_pdb_and_alias_fields():
    candidate = _coerce_candidate(
        {
            "name": "PD-L1",
            "chosen_pdb": " 5o45 ",
            "uniprot_id": "Q9NZQ7",
            "catalog": "10084-H08H-B",
            "vendor_accession": "Q9NZQ7",
        }
    )
    assert candidate is not None
    assert candidate.target_name == "PD-L1"
    assert candidate.pdb_id == "5O45"
    assert candidate.uniprot == "Q9NZQ7"
    assert candidate.antigen_catalog == "10084-H08H-B"


def test_semantic_alias_merges_equivalent_target_into_matched_group():
    records = [
        _record(
            1,
            preset_name="PD-L1",
            protein_name="Programmed cell death 1 ligand 1",
            pdb_id="5O45",
            accession="Q9NZQ7",
            gene="CD274",
            uniprot="Q9NZQ7",
            catalog="10084-H08H-B",
        ),
    ]
    lookup = _build_catalog_lookup(records)
    candidates = [
        BulkLlmCandidate(gene="CD274"),
        BulkLlmCandidate(target_name="PD-L1"),
    ]
    matches, unmatched, matched_rows = _match_candidates_to_catalog(lookup, candidates)
    assert len(matches) == 2
    assert any(match.match_type == "semantic_alias" for match in matches)
    assert len(unmatched) == 0
    assert len(matched_rows) == 1
    assert matched_rows[0].resolved_pdb_id == "5O45"


def test_preview_parses_selection_biotinylated_and_tags_columns():
    csv_text = (
        "rank\tselection\tbiotinylated\ttags\tgene\tprotein_name\tchosen_pdb\tvendor_accession\tantigen_url\n"
        "1\tbiotin\ttrue\tAviTag;His\tCD274\tPD-L1\t5O45\tQ9NZQ7\thttps://example.org/pd-l1\n"
    )
    preview = preview_bulk_targets(BulkPreviewRequest(csv_text=csv_text))
    assert len(preview.rows) == 1
    row = preview.rows[0]
    assert row.selection == "biotin"
    assert row.biotinylated is True
    assert row.tags == "AviTag;His"
