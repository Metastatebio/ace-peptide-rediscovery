# ACE PLM Variant Consensus

Date: 2026-08-30

## Inputs

| Model | Variant table |
| --- | --- |
| esm2_8m | /home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-30-ace-protein-lm-fusion-001/ace_plm_masked_variant_candidates.csv |
| esm2_35m | /home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-30-ace-protein-lm-fusion-002-esm35m/ace_plm_masked_variant_candidates.csv |
| esm2_150m | /home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-30-ace-protein-lm-fusion-003-esm150m/ace_plm_masked_variant_candidates.csv |

## Consensus Candidates

Candidates are ranked first by checkpoint agreement, then by mean rank within the model-specific variant tables. These are not validated inhibitors; they are the strongest current computational candidates for a wet-lab or receptor-conditioned follow-up batch.

| Rank | Sequence | Models | Model IDs | Mean model rank | Parents | Nearest known positive | Distance | Developability |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | MSSP | 3 | esm2_150m;esm2_35m;esm2_8m | 2 | IYSP | AGSP | 2 | screen_pass |
| 2 | MLSP | 3 | esm2_150m;esm2_35m;esm2_8m | 6 | IYSP | LSP | 1 | screen_pass |
| 3 | MFLLQ | 3 | esm2_150m;esm2_35m;esm2_8m | 7 | MFGPQ | MFGPQ | 2 | screen_pass |
| 4 | MYSL | 3 | esm2_150m;esm2_35m;esm2_8m | 9.333 | IYSP | IYLL | 2 | screen_pass |
| 5 | MFLPL | 3 | esm2_150m;esm2_35m;esm2_8m | 11.667 | MFGPQ | MFGPQ | 2 | screen_pass |
| 6 | MPLLQ | 3 | esm2_150m;esm2_35m;esm2_8m | 13 | YPGLQ | YPGLQ | 2 | screen_pass |
| 7 | MLLPQ | 3 | esm2_150m;esm2_35m;esm2_8m | 15.667 | PGLPQ | MFGPQ | 2 | screen_pass |
| 8 | MFGLL | 3 | esm2_150m;esm2_35m;esm2_8m | 20 | MFGPQ | MFGPQ | 2 | screen_pass |
| 9 | MYSS | 3 | esm2_150m;esm2_35m;esm2_8m | 20.667 | IYSP | AGSS | 2 | screen_pass |
| 10 | MFFLQ | 3 | esm2_150m;esm2_35m;esm2_8m | 21 | MFGPQ | MDFLI | 2 | screen_pass |
| 11 | MGLPL | 3 | esm2_150m;esm2_35m;esm2_8m | 21.667 | PGLPQ | PGLPQ | 2 | screen_pass |
| 12 | MYSI | 3 | esm2_150m;esm2_35m;esm2_8m | 22 | IYSP | IYSP | 2 | screen_pass |
| 13 | MPGTAA | 3 | esm2_150m;esm2_35m;esm2_8m | 22.667 | SPGTAF | SPGTAF | 2 | screen_pass |
| 14 | PGPPG | 3 | esm2_150m;esm2_35m;esm2_8m | 22.667 | PGLPQ | PGLPQ | 2 | screen_pass |
| 15 | MRWPR | 3 | esm2_150m;esm2_35m;esm2_8m | 23.667 | LHWPR | LHWPR | 2 | screen_pass |
| 16 | MPGLL | 3 | esm2_150m;esm2_35m;esm2_8m | 24.333 | YPGLQ | YPGLQ | 2 | screen_pass |
| 17 | MPSLQ | 3 | esm2_150m;esm2_35m;esm2_8m | 27.333 | YPGLQ | YPGLQ | 2 | screen_pass |
| 18 | MFFPL | 3 | esm2_150m;esm2_35m;esm2_8m | 29.667 | MFGPQ | MFGPQ | 2 | screen_pass |
| 19 | MGLLQ | 3 | esm2_150m;esm2_35m;esm2_8m | 30 | PGLPQ | PGLPQ | 2 | screen_pass |
| 20 | PPLPL | 3 | esm2_150m;esm2_35m;esm2_8m | 33.333 | PGLPQ | PGLPQ | 2 | screen_pass |

## Output

- Consensus CSV: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-30-ace-plm-variant-consensus-001/ace_plm_consensus_variant_candidates.csv`
- Manifest: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-30-ace-plm-variant-consensus-001/ace_plm_variant_consensus_manifest.json`
