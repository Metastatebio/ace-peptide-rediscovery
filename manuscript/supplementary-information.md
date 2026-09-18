# Supplementary Information

Manuscript: ProteinLM-Augmented Blinded Rediscovery of Food-Derived ACE-Inhibitory Peptides

Authors: Oğuzcan Ünver

Affiliation: Metastate Bio Inc, Wilmington, Delaware, USA

## Supplementary Methods

### Positive Panels

The public `v0.1.0` release stores curated positive-panel CSV files under `results/core-results/`. Each row records sequence, target, activity class, source material, generation context, evidence type, reference URL, and notes. The claim-bearing independent holdouts are:

- v4: `results/core-results/known-positive-v4-independent-primary-panel.csv`
- v5: `results/core-results/known-positive-v5-failure-mode-primary-panel.csv`
- v6: `results/core-results/known-positive-v6-independent-primary-panel.csv`

### Exclusion Audits

Exact-sequence overlap audits were run with `scripts/ace_audit_panel_overlap.py`. v6 was audited against seed, expanded, post-freeze H12, v3, v4, v5, and all AHTPDB-derived training panels. The final v6 audit reported zero overlapping sequences (`results/independent-validation/v6-exclusion-audit.json`).

### Motif-Neighborhood Audit

The v6 motif-neighborhood audit was run with `scripts/ace_motif_neighborhood_audit.py`. It reports nearest prior positive sequence, edit distance, normalized identity, one-edit and two-edit neighbor counts, and shared C-terminal dipeptide/tripeptide motifs for every v6 positive (`results/independent-validation/v6-motif-neighborhood-audit.csv`).

### Frozen Recipes

Frozen v4 recipes are stored under `results/frozen-recipes/2026-08-29-v4-freeze-001/`. Frozen v6 recipes are stored under `results/frozen-recipes/2026-08-29-v6-freeze-001/`. Recipe JSON files specify gate components, rank components, gate percentages, required score inputs, and fusion weights.

### External Comparator

The PeptideRanker-style comparator was run with `scripts/ace_run_peptideranker_style_baseline.py`. It is a transparent local composition and length baseline inspired by Mooney et al. 2012, not the official PeptideRanker neural-network server.

### ProteinLM Fusion

ProteinLM scenarios were run with `scripts/ace_protein_lm_fusion_scenario.py` using ESM2 checkpoints. Parent-grouped repeated cross-validation was run with `scripts/ace_plm_cross_validated_fusion_audit.py`. The full-universe frozen 0.65 ACE / 0.35 length-calibrated PLM rank was generated with `scripts/ace_plm_full_universe_rank.py`. Cross-checkpoint masked-variant consensus was generated with `scripts/ace_plm_variant_consensus.py`.

## Supplementary Tables

Supplementary Table S1. Structured assay metadata for v3, v4, v5, and v6 panels: `results/tables/ace-positive-panel-assay-metadata.csv`.

Supplementary Table S2. Frozen-validation summary: `results/tables/ace-frozen-validation-summary.csv`.

Supplementary Table S3. External PeptideRanker-style baseline summary: `results/tables/ace-external-baseline-summary.csv`.

Supplementary Table S4. Panel composition summary: `results/tables/ace-panel-composition-summary.csv`.

Supplementary Table S5. Release verifier and reproducibility metadata: `scripts/verify_release.py`, `results/`, and `data/labels-post-unblinding/`.

Supplementary Table S6. v6 motif-neighborhood audit: `results/independent-validation/v6-motif-neighborhood-audit.csv`.

Supplementary Table S7. ProteinLM parent-grouped cross-validation summary: `results/proteinlm/cv-fusion-audit/ace_plm_cv_fusion_summary.csv`.

Supplementary Table S8. ProteinLM full-universe frozen fusion summary: `results/proteinlm/full-universe-rank/ace_plm_lengthz_full_universe_summary.csv`.

Supplementary Table S9. ProteinLM full-universe digestome discovery shortlist: `results/proteinlm/full-universe-rank/ace_plm_lengthz_digestome_discovery_shortlist.csv`.

Supplementary Table S10. ProteinLM full-universe computational discovery shortlist including synthetic background: `results/proteinlm/full-universe-rank/ace_plm_lengthz_discovery_shortlist.csv`.

Supplementary Table S11. Cross-checkpoint masked-variant consensus: `results/proteinlm/variant-consensus/ace_plm_consensus_variant_candidates.csv`.

## Supplementary Figures

Supplementary Figure S1. Frozen holdout recovery: `results/figures/frozen-validation-recovery.svg`.

Supplementary Figure S2. Matched-decoy AUROC: `results/figures/frozen-validation-matched-auroc.svg`.

Supplementary Figure S3. Panel length composition: `results/figures/panel-composition-length.svg`.

## Current Archive Status

The public reproducibility release is tagged `v0.1.0` at https://github.com/Metastatebio/ace-peptide-rediscovery/releases/tag/v0.1.0. The automated verification script checks the corrected headline metrics against released labels and results. An archival DOI may be added later but is not required for submission.
