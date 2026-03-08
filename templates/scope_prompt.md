You are a protein design expert. Your task is to analyze the provided PDB metadata and propose a design scope for creating a binder against this target.

Based on the metadata below, identify potential epitopes. For each epitope, provide a list of residue ranges (e.g., "A:10-20") and a brief rationale for your choice.
Please ensure the proposed epitopes are distinct regions on the target protein that are likely to be accessible and relevant for binding.

The final output must be a single YAML block.

--- PDB METADATA ---
{meta}

{uniprot_context}

{chain_constraints}
--- END METADATA ---

Propose a scope in the following YAML format:
```yaml
id: "PDB_ID_HERE"
assembly_id: "1"
target_name: "A descriptive name for the target protein"
chains: ["A", "B"] # List the chains that make up the antigen target (validated against the commercial antigen)
target_chains: ["A", "B"] # Must exactly match the validated antigen-supported chains
epitopes:
  - name: "Epitope 1 Site"
    residues: ["A:50-65", "A:80-88"] # Residue ranges for the epitope
    rationale: "Rationale for choosing this epitope, based on function, literature, or structural features."
  - name: "Epitope 2 Site"
    residues: ["B:25-40"]
    rationale: "Another rationale."
  - name: "Epitope 3 Site"
    residues: ["B:60-70"]
    rationale: "Another rationale."
