# ACE v6 Frozen Independent Validation Report

Date: 2026-08-29

## Purpose

v6 tests the two candidate recipes frozen after the v5 failure-mode development cycle. The v6 panel was curated after recipe freeze, exact-sequence-audited before scoring, and then ranked once with the stored recipe components, gate percentages, and fusion weights unchanged.

## Blinded Universe

- Universe: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/blinded-candidate-universe.csv`
- Candidates: `58,192`
- Hidden positives: `52`
- Matched decoys: `4,160`
- Random background: `50,000`
- Digestome background: `3,980`
- Digestome-known positives: `3`
- Seed: `2026082906`
- Exact-sequence overlap audit: `0` overlaps against seed, expanded, H12, v3, v4, v5, and all AHTPDB-derived training panels.

## Frozen Validation Results

| Frozen recipe | Purpose | Top 1% | Top 5% | Top 10% | Best rank | AUROC matched | Interpretation |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `ace_v6_candidate_from_v5_broad_shortlist_v1` | broad top-10 recovery | 6/52 | 13/52 | 20/52 | 93 | 0.6868 | positive early and broad enrichment, weaker hard-decoy separation |
| `ace_v6_candidate_from_v5_balanced_v1` | early recovery with stronger hard-decoy discrimination | 6/52 | 16/52 | 18/52 | 121 | 0.7642 | positive early enrichment and stronger matched-decoy AUROC |

For both recipes, top-1-percent recovery was `6/52`, equal to `11.54x` over random expectation with hypergeometric p = `1.34e-05`. The broad recipe recovered `20/52` positives in the top 10 percent, `3.85x` over random expectation with p = `5.10e-08`. The balanced recipe recovered `16/52` positives in the top 5 percent, `6.15x` over random expectation with p = `2.73e-09`, and retained matched-decoy AUROC above `0.75`.

## Baseline Score Modes

| Mode | Top 1% | Top 10% | Best rank | AUROC matched | Role |
| --- | ---: | ---: | ---: | ---: | --- |
| `generic_zero_shot` | 5/52 | 16/52 | 74 | 0.5788 | target-agnostic control |
| `oral_stability_zero_shot` | 4/52 | 17/52 | 144 | 0.5972 | developability control |
| `no_position_priors` | 4/52 | 16/52 | 82 | 0.6138 | position-prior control |
| `ace_no_exact_proline_motif` | 0/52 | 12/52 | 1835 | 0.6128 | ACE target-informed, no direct terminal-proline motif |
| `ace_pharmacophore_zero_label` | 2/52 | 11/52 | 380 | 0.6265 | ACE target-informed, label-free |
| `literature_sar` | 2/52 | 11/52 | 419 | 0.6279 | ACE SAR upper-bound heuristic |

## Top Recovered Positives

The broad recipe's first twelve recovered positives were `LHWPR` rank `93`, `IPPAYTK` rank `130`, `YPGLQ` rank `302`, `MFGPQ` rank `318`, `SPGTAF` rank `411`, `IYSP` rank `546`, `SLPEW` rank `789`, `LLRP` rank `860`, `PGLPQ` rank `997`, `IGPR` rank `1256`, `NLGPR` rank `2113`, and `LPGP` rank `2152`.

The balanced recipe's first twelve recovered positives were `YPGLQ` rank `121`, `IYSP` rank `150`, `MFGPQ` rank `206`, `SPGTAF` rank `225`, `LHWPR` rank `292`, `IPPAYTK` rank `381`, `PGLPQ` rank `1010`, `SLPEW` rank `1136`, `LVLPGE` rank `1306`, `IGPR` rank `1309`, `LVLPGELAK` rank `1394`, and `LPGP` rank `1708`.

## Interpretation

v6 converts the v5-derived candidate recipes into a positive independent validation result. The cleanest claim is the balanced recipe: it prospectively recovers new exact-sequence-excluded primary-literature ACE peptides at the top of a blinded universe while maintaining matched-decoy AUROC above `0.75`.

The result still needs conservative language. v6 does not erase v5's failure-mode finding, and the v6 panel is short-peptide heavy. The manuscript should be positioned as leakage-controlled blinded rediscovery and prioritization with repeated frozen-holdout evidence, not as novel peptide discovery or a universal ACE activity predictor.

## Outputs

- Source audit: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/v6-primary-literature-source-audit.md`
- Baseline sweep: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/score-mode-sweep/score-mode-sweep-report.md`
- Broad frozen recipe report: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/frozen-rerank/broad-v6-candidate/frozen-rerank-report.md`
- Balanced frozen recipe report: `docs/publication-engine/papers/01-ace-rediscovery/experiments/2026-08-29-v6-independent-validation-001/frozen-rerank/balanced-v6-candidate/frozen-rerank-report.md`
