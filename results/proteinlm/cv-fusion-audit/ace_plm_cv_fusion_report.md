# ACE ProteinLM Cross-Validated Fusion Audit

Date: 2026-09-18

## Inputs

| Scenario | Scored panel |
| --- | --- |
| esm2_8m | /tmp/ace-qc.jXpMGi/no-zero-plm/cv/all/esm2_8m.csv |
| esm2_35m | /tmp/ace-qc.jXpMGi/no-zero-plm/cv/all/esm2_35m.csv |
| esm2_150m | /tmp/ace-qc.jXpMGi/no-zero-plm/cv/all/esm2_150m.csv |

## Parent-Grouped Cross-Validation

Parent grouping keeps each literature positive and its composition-shuffled or length-matched decoys in the same fold. Fusion weights are selected on each training fold, then evaluated on held-out parent groups.

| Scenario | PLM feature | CV AUROC | ACE-only AUROC | Delta | Win rate | Median ACE weight | Sign-test p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| esm2_8m | length_z | 0.7835 | 0.7643 | 0.0192 | 0.800 | 0.65 | 1.12e-09 |
| esm2_8m | raw | 0.7813 | 0.7643 | 0.0170 | 0.800 | 0.70 | 1.12e-09 |
| esm2_150m | raw | 0.7702 | 0.7643 | 0.0059 | 0.800 | 0.80 | 1.12e-09 |
| esm2_150m | length_z | 0.7694 | 0.7643 | 0.0051 | 0.800 | 0.75 | 1.12e-09 |
| esm2_35m | raw | 0.7648 | 0.7643 | 0.0005 | 0.400 | 0.85 | 0.0569 |
| esm2_35m | length_z | 0.7647 | 0.7643 | 0.0004 | 0.400 | 0.75 | 0.0569 |

Verdict: `strengthens`.

Best cross-validated scenario: `esm2_8m` with `length_z` PLM feature. Mean held-out AUROC `0.7835` versus ACE-only `0.7643`, delta `0.0192`.

## Full-Panel Exploratory Readout

The full-panel table is exploratory because it selects a weight after seeing all labels. It is useful for prioritizing the next preregistered run, not as a stand-alone validation claim.

| Scenario | PLM feature | ACE wt | PLM wt | Top 1% | AUROC | AP | Median positive rank |
| --- | --- | --- | --- | --- | --- | --- | --- |
| esm2_150m | length_z | 0.70 | 0.30 | 4 | 0.7761 | 0.0399 | 734.5 |
| esm2_35m | length_z | 0.75 | 0.25 | 4 | 0.7715 | 0.0398 | 828.0 |
| esm2_8m | length_z | 0.65 | 0.35 | 4 | 0.7850 | 0.0430 | 758.0 |

## Files

- Summary: `/tmp/ace-qc.jXpMGi/no-zero-plm/cv/all-results/ace_plm_cv_fusion_summary.csv`
- Fold-level rows: `/tmp/ace-qc.jXpMGi/no-zero-plm/cv/all-results/ace_plm_cv_fusion_folds.csv`
- Full-panel exploratory metrics: `/tmp/ace-qc.jXpMGi/no-zero-plm/cv/all-results/ace_plm_full_panel_exploratory_metrics.csv`
- Manifest: `/tmp/ace-qc.jXpMGi/no-zero-plm/cv/all-results/ace_plm_cv_fusion_manifest.json`
