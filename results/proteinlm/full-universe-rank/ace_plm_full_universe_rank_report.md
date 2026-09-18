# ACE PLM Full-Universe Rank

Date: 2026-09-18

## Frozen Fusion

- Fusion: `ace_0.65_plm_length_z_0.35`
- Total candidates: `58192`
- Positive labels: `52`
- Input v6 top 1% positives: `6`
- Input v6 top 100 positives: `0`
- Input v6 best positive rank: `129`
- Input v6 matched-decoy AUROC: `0.7647`
- Input v6 all-negative AUROC: `0.7926`

## Full-Universe Result

| Metric | Value |
| --- | --- |
| Top 1% positives | 6 |
| Top 1% enrichment | 11.54x |
| Top 1% hypergeometric p | 0.000013 |
| Top 100 positives | 1 |
| Best positive rank | 45 |
| Matched-decoy AUROC | 0.7817 |
| All-negative AUROC | 0.8129 |
| Median positive rank | 8360.0 |

## Top Positive Recoveries

| Rank | Sequence | Label | v6 rank | v6 score | PLM z | Fused |
| --- | --- | --- | --- | --- | --- | --- |
| 45 | MFGPQ | known_positive | 195 | 0.9445 | 1.568 | 0.7578 |
| 133 | SPGTAF | known_positive | 227 | 0.9343 | 0.787 | 0.7381 |
| 162 | IYSP | known_positive | 149 | 0.9627 | 0.277 | 0.7329 |
| 215 | YPGLQ | known_positive | 129 | 0.9676 | -0.062 | 0.7265 |
| 301 | LHWPR | known_positive | 288 | 0.9021 | 0.115 | 0.7162 |
| 440 | IPPAYTK | known_positive | 417 | 0.8210 | -0.126 | 0.6932 |
| 1403 | LVLPGELAK | known_positive | 1507 | -1.0120 | 1.741 | 0.3377 |
| 1513 | LLRP | known_positive | 2639 | -1.0325 | 1.807 | 0.3347 |
| 1558 | LVLPGE | known_positive | 1370 | -1.0094 | 1.520 | 0.3333 |
| 1633 | VLGA | known_positive | 3586 | -1.0497 | 1.839 | 0.3317 |
| 1736 | PGLPQ | known_positive | 1073 | -1.0038 | 1.295 | 0.3295 |
| 1972 | LPGP | digestome_known_positive | 1646 | -1.0145 | 1.200 | 0.3251 |

## Computational Digestome Discovery Shortlist

This shortlist excludes known positives, matched decoys, and synthetic random peptides. It is the most relevant computational prioritization list for a food-peptide validation paper.

| Rank | Sequence | Origin | v6 rank | PLM z | Fused |
| --- | --- | --- | --- | --- | --- |
| 1 | AGPAGP | collagen:P02453:thermolysin_like | 89 | 4.191 | 0.8242 |
| 3 | LHLPLP | milk:P02666:thermolysin_like | 1 | 3.299 | 0.8084 |
| 4 | HLPLPL | milk:P02666:chymotrypsin | 102 | 3.495 | 0.8079 |
| 6 | APGPAG | collagen:P02453:thermolysin_like | 205 | 3.520 | 0.8006 |
| 8 | AGPPGFPGA | collagen:P02453:thermolysin_like | 167 | 3.063 | 0.7937 |
| 9 | AGPPGFPG | collagen:P02453:thermolysin_like | 150 | 2.901 | 0.7914 |
| 10 | AGAPGFPG | collagen:P02453:thermolysin_like | 148 | 2.726 | 0.7876 |
| 11 | MGPPGL | collagen:P02453:thermolysin_like | 323 | 3.475 | 0.7871 |
| 15 | AGLPG | collagen:P02453:thermolysin_like | 59 | 2.239 | 0.7823 |
| 16 | MFPPQSV | milk:P02666:thermolysin_like | 146 | 2.419 | 0.7808 |
| 20 | PGTAGL | collagen:P02453:pepsin_like | 139 | 1.938 | 0.7704 |
| 21 | APGAPG | collagen:P02453:thermolysin_like | 502 | 3.775 | 0.7690 |
| 22 | APGLQG | collagen:P02453:thermolysin_like | 118 | 1.792 | 0.7684 |
| 27 | MFPPQS | milk:P02666:thermolysin_like | 301 | 2.380 | 0.7652 |
| 32 | AAGP | collagen:P02453:thermolysin_like | 275 | 2.097 | 0.7625 |

## Computational Discovery Shortlist

This shortlist excludes known positives and matched decoys. It is a computational prioritization list, not a binding or potency claim.

| Rank | Sequence | Label | Origin | v6 rank | PLM z | Fused |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | AGPAGP | digestome_background | collagen:P02453:thermolysin_like | 89 | 4.191 | 0.8242 |
| 2 | MLPGR | random_background | random_background | 31 | 3.560 | 0.8130 |
| 3 | LHLPLP | digestome_background | milk:P02666:thermolysin_like | 1 | 3.299 | 0.8084 |
| 4 | HLPLPL | digestome_background | milk:P02666:chymotrypsin | 102 | 3.495 | 0.8079 |
| 5 | MRVSPWPL | random_background | random_background | 57 | 3.191 | 0.8036 |
| 6 | APGPAG | digestome_background | collagen:P02453:thermolysin_like | 205 | 3.520 | 0.8006 |
| 8 | AGPPGFPGA | digestome_background | collagen:P02453:thermolysin_like | 167 | 3.063 | 0.7937 |
| 9 | AGPPGFPG | digestome_background | collagen:P02453:thermolysin_like | 150 | 2.901 | 0.7914 |
| 10 | AGAPGFPG | digestome_background | collagen:P02453:thermolysin_like | 148 | 2.726 | 0.7876 |
| 11 | MGPPGL | digestome_background | collagen:P02453:thermolysin_like | 323 | 3.475 | 0.7871 |
| 12 | MGITPYP | random_background | random_background | 80 | 2.450 | 0.7860 |
| 13 | MKGPPL | random_background | random_background | 130 | 2.537 | 0.7843 |
| 14 | MKP | random_background | random_background | 45 | 2.230 | 0.7829 |
| 15 | AGLPG | digestome_background | collagen:P02453:thermolysin_like | 59 | 2.239 | 0.7823 |
| 16 | MFPPQSV | digestome_background | milk:P02666:thermolysin_like | 146 | 2.419 | 0.7808 |

## Files

- Ranked full universe: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_lengthz_full_universe_ranked_candidates.csv`
- Positive ranks: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_lengthz_positive_ranks.csv`
- Summary: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_lengthz_full_universe_summary.csv`
- Digestome discovery shortlist: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_lengthz_digestome_discovery_shortlist.csv`
- Discovery shortlist: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_lengthz_discovery_shortlist.csv`
- Manifest: `/tmp/ace-qc.jXpMGi/no-zero-plm/results/ace_plm_full_universe_rank_manifest.json`
