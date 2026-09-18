# ACE v6 Motif-Neighborhood Audit

Date: 2026-08-30

## Purpose

Check whether the funding-grade ACE v6 balanced rediscovery result is just exact reuse or near-neighbor reuse from prior positive panels.

## Inputs

- Candidate panel: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/known-positive-v6-independent-primary-panel.csv`
- Frozen positive ranks: `/home/openclaw/.openclaw/workspace/peptide-lab-github/docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/frozen-rerank/balanced-v6-candidate/frozen-rerank-positive-ranks.csv`
- Prior positive panels audited: `11`
- Prior unique positive sequences: `1356`

## Summary

- v6 positives audited: `52`
- Exact overlaps with prior panels: `0`
- Top-1% positives: `6`
- Top-1% positives with a one-edit prior neighbor: `1`
- Top-1% positives with a high-similarity prior neighbor, identity >= 0.8: `0`
- Top-1% positives classified only by shared C-terminal motif: `5`
- Top-1% positives distant from prior positive panels: `0`
- Median nearest-prior identity: `0.6`

## Neighborhood Classes

| Class | All positives | Top 1% positives |
| --- | --- | --- |
| distant_from_prior_positive_panels | 3 | 0 |
| shared_cterm2_motif | 29 | 5 |
| shared_cterm3_motif | 3 | 0 |
| single_edit_neighbor | 17 | 1 |

## Interpretation

The exact-sequence leakage result remains clean. The right claim boundary is narrower: this is a frozen ACE rediscovery and prioritization result under exact-sequence exclusion, not a proof of de novo molecular understanding. Sequence-neighborhood structure is now explicit and must travel with the manuscript.

For a funding-grade upgrade, the next holdout should be source-family and motif-neighborhood pre-registered before scoring, with success reported separately for close-neighbor and distant positives.
