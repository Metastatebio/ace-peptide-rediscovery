# ACE Rediscovery Stratified Report

Date: 2026-08-29

This report stratifies the unblinded known positives because the expanded literature panel contains heterogeneous evidence: short tripeptides, longer hydrolysate-derived peptides, reported antihypertensive peptides, and peptides with IC50 values reported in different units.

Total candidates: `16305`  
Known positives: `4`  
Top 1 percent cutoff rank: `164`

| Subset | N | Top 1% | Top 1% frac | Enrichment | Hypergeom p | Top 100 | Mean rank | Median rank | Best rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_literature_positives | 4 | 4 | 1.0 | 99.420732 | 9.87e-09 | 3 | 31.25 | 8.5 | 4 |
| primary_four_anchors | 4 | 4 | 1.0 | 99.420732 | 9.87e-09 | 3 | 31.25 | 8.5 | 4 |
| reported_antihypertensive | 2 | 2 | 1.0 | 99.420732 | 0.000101 | 2 | 7.5 | 7.5 | 4 |
| reported_ic50_uM_le_30 | 0 | 0 |  |  |  | 0 |  |  |  |
| reported_ic50_uM_le_10 | 0 | 0 |  |  |  | 0 |  |  |  |
| short_len_le_5 | 4 | 4 | 1.0 | 99.420732 | 9.87e-09 | 3 | 31.25 | 8.5 | 4 |
| short_len_le_5_ic50_uM_le_30 | 0 | 0 |  |  |  | 0 |  |  |  |
| naturally_released_digestome_positive | 2 | 2 | 1.0 | 99.420732 | 0.000101 | 2 | 7.5 | 7.5 | 4 |

## Interpretation

The cleanest near-term claim should focus on recovery and enrichment of short ACE-like literature positives and the four primary anchors. The broader 81-sequence panel is valuable as a stress test, but it includes long and weak positives that the current sequence-only score does not yet model well.
