# DMS Soluble Region Plan

We need to let the antigen DMS builder focus on soluble or extracellular regions that vendors express for membrane proteins.

## Context
- Vendor metadata in `targets/<PDB>/target.yaml` keeps `sequences.accession.expressed_range` (start-end) and optionally `expressed_aa` for the soluble construct.
- The DMS generator already loads PDB sequences but ignores vendor context, so it mutates the full chain including transmembrane segments.
- We can align the vendor construct back to the prepared PDB chain to figure out which residues correspond to the soluble fragment.

## TODOs
1. [x] **Identify metadata loader** – Add a backend helper that reads the expressed range / sequence from `target.yaml` and aligns it to a PDB chain to recover the matching residue UIDs.
2. [x] **Wire backend option** – Extend `DMSLibraryOptions` / generator / API models so a flag like `restrict_to_expressed_region` filters candidates to those UIDs, with informative messages when the range is missing.
3. [x] **Expose in UI + docs** – Surface the toggle in the Antigen DMS form, include it in the request payload, and show the restriction in the returned summary metadata.
4. [x] **Test coverage** – Add unit tests for the vendor-region mapping helper and overall DMS generation to prove the filter works and degrades gracefully when metadata is absent.
