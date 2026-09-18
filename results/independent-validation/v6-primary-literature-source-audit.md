# ACE v6 Primary-Literature Source Audit

Date: 2026-08-29

## Panel

- Panel file: `docs/publication-engine/papers/01-ace-rediscovery/known-positive-v6-independent-primary-panel.csv`
- Peptides: `52`
- Primary-reference groups: `18`
- Source-material groups: `18`
- Length <= 5 residues: `39`
- Length > 5 residues: `13`
- C-terminal proline: `4`
- Non-C-terminal proline: `48`

## Inclusion Rule

The v6 panel was curated only after the v6 candidate recipes were frozen under `docs/publication-engine/papers/01-ace-rediscovery/frozen-recipes/2026-08-29-v6-freeze-001/`. Included rows required a specific peptide sequence with reported ACE-inhibitory activity from a primary source, preferably synthesized peptide IC50 assays or direct peptide activity assays.

Pure docking-only hits, database-only predictions, and exact-sequence overlaps with previous training or holdout panels were excluded.

## Independence Audit

The final v6 panel passed exact-sequence exclusion against every prior ACE seed, expanded, post-freeze H12, v3, v4, v5, and AHTPDB-derived training panel.

- Audit JSON: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/v6-exclusion-audit.json`
- Overlap sequences: `0`
- Overlap rows: `0`

During curation, `VSLPEW` and `GVSLPEW` were removed because they overlapped AHTPDB-derived training panels. Additional candidate rows from chicken, whey, monkfish, and soybean sources were checked and excluded where exact overlaps with v5 or AHTPDB-derived panels were found.

## Source Mix

| Source material | Peptides | Reference |
| --- | ---: | --- |
| Takifugu flavidus hydrolysate | 8 | https://doi.org/10.3390/md21100522 |
| broccoli protein hydrolysate | 5 | https://doi.org/10.1021/acs.jafc.9b01137 |
| porcine liver and placenta hydrolysate | 4 | https://doi.org/10.3390/molecules30030754 |
| Tenebrio molitor protein hydrolysate | 4 | https://doi.org/10.26599/FSHW.2025.9250609 |
| maize germ protein hydrolysate | 4 | https://doi.org/10.1016/j.lwt.2023.115254 |
| alpha-lactalbumin hydrolysate | 3 | https://doi.org/10.1016/j.lwt.2021.112984 |
| hazelnut protein hydrolysate | 3 | https://www.sciencedirect.com/science/article/abs/pii/S1756464616303942 |
| hemp seed protein hydrolysate | 3 | https://doi.org/10.1021/acs.jafc.7b04522 |
| monkfish swim bladder collagen hydrolysate | 3 | https://doi.org/10.3389/fnut.2022.957778 |
| tuna meat hydrolysate | 3 | https://doi.org/10.1021/acs.jafc.6c03183 |
| Cangkuk fermented beef | 2 | https://doi.org/10.5713/ab.23.0433 |
| Salmo salar protein hydrolysate | 2 | https://doi.org/10.1002/jsfa.8908 |
| pearl oyster meat hydrolysate | 2 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6723713/ |
| peony seed protein hydrolysate | 2 | https://www.sciencedirect.com/science/article/pii/S1756464622002213 |
| Torreya grandis protein hydrolysate | 1 | https://doi.org/10.3390/nu15102374 |
| oyster hydrolysate | 1 | https://doi.org/10.3389/fnut.2022.981163 |
| soybean protein hydrolysate | 1 | https://doi.org/10.3390/foods11172667 |
| tree bean seed protein hydrolysate | 1 | https://doi.org/10.1007/s00217-025-05009-0 |

## Curation Caveats

The v6 panel is exact-sequence independent, not literature-domain independent: it remains intentionally focused on food-derived and natural-source ACE-inhibitory peptides. The panel is also still short-peptide heavy, with `39/52` sequences at five residues or fewer. This makes the v6 readout appropriate for repeated blinded rediscovery and prioritization claims, not for a general ACE activity predictor claim across all peptide classes.
