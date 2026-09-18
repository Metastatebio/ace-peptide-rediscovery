# ACE Frozen Rerank Recipe Application

Date: 2026-09-18

## Frozen Recipe

- Recipe: `ace_v6_candidate_from_v5_balanced_v1`
- Recipe file: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/frozen-recipes/2026-08-29-v6-freeze-001/balanced-v6-candidate.recipe.json`
- Gate component: `fusion:v5dev_mode_dirichlet_1353`
- Rank component: `lenrank:ml:ahtpdb-plus-expanded-plus-h12:logreg_l2_c0_05`
- Gate percent: `1.5`
- Include length-stratified components: `True`
- Candidate rows: `58192`
- Label counts: `{'digestome_background': 3980, 'digestome_known_positive': 3, 'known_positive': 49, 'matched_decoy': 4160, 'random_background': 50000}`

## Inputs

- Model score directories: `/tmp/ace-qc.jXpMGi/model-scores/ahtpdb-plus-expanded-plus-h12`
- Baseline sweep directory: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/score-mode-sweep`
- Components available after frozen expansion: `15`

## Validation Metrics

| Top 1% | Top 10% | Top 100 | AUROC matched | AUROC all | Mean positive rank |
| --- | --- | --- | --- | --- | --- |
| 6 | 18 | 0 | 0.7647 | 0.7926 | 12086.3 |

## Interpretation

This run applies a locked recipe and does not search gate components, rank components, gate percentages, or fusion weights on the current panel. Any validation labels are used only after ranking to compute readout metrics.

## Outputs

- Ranked candidates: `/tmp/ace-qc.jXpMGi/no-zero-frozen-rerank/frozen-rerank-ranked-candidates.csv`
- Summary: `/tmp/ace-qc.jXpMGi/no-zero-frozen-rerank/frozen-rerank-summary.csv`
- Positive ranks: `/tmp/ace-qc.jXpMGi/no-zero-frozen-rerank/frozen-rerank-positive-ranks.csv`
