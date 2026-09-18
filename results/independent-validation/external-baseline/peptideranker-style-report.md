# ACE PeptideRanker-Style Baseline

Date: 2026-08-30

## Baseline

- Baseline ID: `peptideranker_style_mooney2012_composition_v1`
- Family: general bioactive-peptide comparator inspired by Mooney et al. 2012 PeptideRanker.
- Literature reference: https://doi.org/10.1371/journal.pone.0045012
- Evaluation directory: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001`
- Candidate rows: `58192`
- Known positives: `52`

This is not the official PeptideRanker neural-network server. It is a transparent local comparator using peptide length and amino-acid composition only. It does not use ACE-specific motifs, ACE labels, target structure, or holdout labels during ranking.

## Results

| Top 1% | Top 5% | Top 10% | Top-1 p | AUROC matched | AUROC all |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | 3 | 0.407 | 0.5585 | 0.4597 |

## Outputs

- Ranked candidates: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/external-baselines/peptideranker-style/peptideranker-style-ranked-candidates.csv`
- Positive ranks: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/external-baselines/peptideranker-style/peptideranker-style-positive-ranks.csv`
- Summary: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/external-baselines/peptideranker-style/peptideranker-style-summary.csv`
