# 4 · Evaluation against TF knockouts

> Benchmark the three annotation strategies (network-based / common-target / common-function) against DEGs from published TF-knockout experiments. **Manuscript:** Fig 2.

**Pipeline position:** **⬅ prev** [`3_Functional_annotation`](../3_Functional_annotation) · [⌂ overview](../README.md) · **next ➡** [`6_Prioritization_and_conditions`](../6_Prioritization_and_conditions)

| Script | What it does | Inputs | Output |
|---|---|---|---|
| `Fig_MethodsComparison.R` | compare strategy predictions vs knockout DEGs; assemble the per-TF annotation summary | per-strategy enrichment tables from step 3 (`CommonTarg_…`, `CommonFunction_…`, `NetworkBased_…` PWY/GO enrichment), knockout-DEG database | Fig 2 panels; `Summary.Total.Annotation.txt` (→ step 6, Table S8 ⭐) |
| `DEGdb_vs_random/DEGdb_vs_random.R` | null control: knockout-DEG recovery vs random gene sets | knockout-DEG database | Fig 2 support |
