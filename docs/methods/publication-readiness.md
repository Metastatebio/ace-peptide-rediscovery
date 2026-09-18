# ACE Rediscovery Publication Readiness

Date: 2026-08-29

## Current Experimental State

Ten computational experiment families have been run: pilot rediscovery, expanded-panel rediscovery, score-mode leakage controls, ACE docking probes, a 12-row post-freeze primary-literature holdout, a 51-row v3 primary-literature development stress set, v3 model/fusion iteration, a 31-row v4 frozen independent validation, a 44-row v5 failure-mode validation with post-unblinding v6 candidate freeze, and a 52-row v6 frozen independent validation.

### Pilot v1

- Blinded universe: `16,305` candidates.
- Known positives: `4` primary anchors: `VPP`, `IPP`, `LKPNM`, `LKP`.
- Top 1 percent cutoff: rank `164`.
- Top 1 percent recovery: `4/4`.
- AUROC versus matched decoys: `0.9930`.
- Primary endpoint: passed.

### Expanded v2

- Blinded universe: `63,023` candidates.
- Known positives: `81` literature-reported ACE-inhibitory peptides.
- Top 1 percent cutoff: rank `631`.
- Top 1 percent recovery: `11/81`.
- Top 1 percent enrichment over random expectation: `13.56x`.
- Top 1 percent hypergeometric p-value: `5.98e-10`.
- AUROC versus all negatives: `0.7877`.
- AUROC versus matched decoys: `0.7231`.
- Strict composite endpoint: not passed because matched-decoy AUROC did not reach the configured `0.75` gate.

### Score-Mode Leakage Controls

The expanded universe was rerun with fixed scoring modes that separate target-agnostic scoring from ACE-target-informed scoring.

| Mode | Claim tier | Broad top 1 percent recovery | Broad enrichment | Matched-decoy AUROC | Key result |
| --- | --- | --- | --- | --- | --- |
| `generic_zero_shot` | target-agnostic zero-shot | `7/81` | `8.63x` | `0.5776` | Enriches broad positives but weak hard-decoy separation. |
| `oral_stability_zero_shot` | target-agnostic developability zero-shot | `8/81` | `9.86x` | `0.6575` | Recovers `VPP`, `IPP`, and `LKP` in the top 1 percent without ACE labels or ACE target priors. |
| `ace_pharmacophore_zero_label` | ACE-target-informed, label-free | `11/81` | `13.56x` | `0.7202` | Recovers all four primary anchors in the top 1 percent, including `VPP`, `IPP`, `LKP`, and `LKPNM`. |
| `ace_no_exact_proline_motif` | ACE-target-informed, no exact terminal-proline motif | `4/81` | `4.93x` | `0.6464` | Stress test: only `LKP` remains in the top 1 percent among the four primary anchors. |
| `literature_sar` | ACE SAR label-free but literature-informed | `11/81` | `13.56x` | `0.7231` | Useful upper-bound heuristic, not a clean discovery claim. |

### Molecular Docking Probes

Two Vina docking probes were run against human ACE structure `1O86` with the same 8-positive/16-hard-decoy panel.

| Probe | Ligand setup | Raw Vina AUROC | Pose-QC AUROC | Interpretation |
| --- | --- | --- | --- | --- |
| `2026-08-29-docking-probe-001` | PeptideBuilder output without explicit terminal `OXT` | `0.1719` | `0.3047` | Negative control; ligand chemistry was incomplete. |
| `2026-08-29-docking-probe-oxt-001` | Explicit C-terminal `OXT`; VPP PDBQT confirmed as `C(=O)O` | `0.1250` | `0.1250` | Corrected chemistry but still anti-discriminative against hard decoys. |

The current docking tier should not be used as supporting evidence for ACE inhibition. It is useful because it shows that PepLab is not overclaiming from molecular simulation. A paper-grade physics tier needs corrected protonation states, peptide conformer ensembles, receptor ensemble docking, zinc-aware constraints, and independent pose/free-energy controls before it becomes a main result.

### Post-Freeze Primary-Literature Holdout

After freezing the current scoring modes, a new 12-peptide holdout panel was curated from primary literature sources not included in the original 81-peptide panel: spinach Rubisco peptides, ovotransferrin `EWL`, Spirulina `VTY`/`LGVP`, and Sardina peptides.

- Fresh universe: `54,954` blinded candidates.
- Hidden positives: `12`.
- Top 1 percent cutoff: rank `550`.

| Mode | Claim tier | Holdout top 1 percent recovery | Holdout enrichment | Matched-decoy AUROC | Interpretation |
| --- | --- | --- | --- | --- | --- |
| `generic_zero_shot` | target-agnostic zero-shot | `0/12` | `0.00x` | `0.4769` | Failed on this holdout. |
| `oral_stability_zero_shot` | target-agnostic developability zero-shot | `1/12` | `8.33x` | `0.5322` | Weak partial recovery; `LGVP` rank `199`. |
| `no_position_priors` | negative control | `0/12` | `0.00x` | `0.5127` | Negative control. |
| `ace_no_exact_proline_motif` | ACE-target-informed, no exact terminal-proline motif | `2/12` | `16.65x` | `0.6784` | Partial holdout signal; `MRW` rank `76`, `VTY` rank `330`. |
| `ace_pharmacophore_zero_label` | ACE-target-informed, label-free | `0/12` | `0.00x` | `0.6685` | Failed top 1 percent recovery despite modest AUROC. |
| `literature_sar` | ACE SAR label-free but literature-informed | `0/12` | `0.00x` | `0.6706` | Failed top 1 percent recovery. |

This holdout does not support a broad general ACE peptide predictor claim. It does support a narrower claim that the frozen engine can recover some unseen primary-literature ACE peptides, with strongest signal in the no-exact-terminal-proline target-informed mode.

### Post-Freeze v3 Primary-Literature Panel

A larger 51-peptide holdout was curated from open primary-literature sources after the current scoring modes existed. It uses no sequences from the earlier 81-row expanded panel or the 12-row post-freeze holdout.

- Fresh universe: `58,114` blinded candidates.
- Hidden positives: `51`.
- Top 1 percent cutoff: rank `582`.

| Mode | Claim tier | Holdout top 1 percent recovery | Holdout enrichment | Matched-decoy AUROC | Interpretation |
| --- | --- | --- | --- | --- | --- |
| `generic_zero_shot` | target-agnostic zero-shot | `1/51` | `1.96x` | `0.5199` | Failed broad holdout gate. |
| `oral_stability_zero_shot` | target-agnostic developability zero-shot | `1/51` | `1.96x` | `0.6070` | Failed broad holdout gate. |
| `no_position_priors` | negative control | `1/51` | `1.96x` | `0.5597` | Negative control did not separate meaningfully. |
| `ace_no_exact_proline_motif` | ACE-target-informed, no exact terminal-proline motif | `2/51` | `3.92x` | `0.7467` | Moderate matched-decoy separation but too little top-rank recovery. |
| `ace_pharmacophore_zero_label` | ACE-target-informed, label-free | `1/51` | `1.96x` | `0.7530` | Moderate AUROC but failed top-rank recovery. |
| `literature_sar` | ACE SAR label-free but literature-informed | `1/51` | `1.96x` | `0.7514` | Upper-bound heuristic still fails top-rank recovery. |

The v3 result blocks a broad claim that the current sequence-only PepLab scorer generally rediscovers heterogeneous ACE-inhibitory peptides across primary literature. It supports a more useful development conclusion: the current engine has coarse ACE-target signal but insufficient top-rank resolution.

### Frozen v4 Independent Validation

After v3 was converted into a development stress set, two v3-selected recipes were frozen before v4 literature curation. A new exact-sequence-excluded v4 panel was curated from five primary-literature sources: Larimichthys crocea, sesame protein, Flammulina velutipes, camel casein, and Pelodiscus sinensis. The reproducible overlap audit found `0` exact sequence overlaps against seed, expanded, H12, v3, and all AHTPDB-derived training panels.

- Fresh universe: `56,494` blinded candidates.
- Hidden positives: `31`.
- Matched decoys: `2,480`.
- Random background: `50,000`.
- Digestome background: `3,983`.

| Frozen endpoint | Purpose | Top 1 percent recovery | Top 5 percent recovery | Top 10 percent recovery | Matched-decoy AUROC | Interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| `early-enrichment-v1` | aggressive early recovery | `1/31`; p `0.268` | `6/31`; p `0.00389` | `7/31`; p `0.0305` | `0.6586` | Failed early endpoint. |
| `broad-shortlist-v1` | broad shortlist validation | `2/31`; p `0.0384` | `9/31`; p `1.43e-05` | `16/31`; p `6.76e-09` | `0.7575` | Passed broad-shortlist endpoint. |

The v4 result supports a narrow but real claim that a prospectively frozen two-stage PepLab recipe enriched known active food-derived ACE-inhibitory peptides from an exact-sequence-excluded blinded search space. It does not support a claim that PepLab discovered novel peptides or generally predicts all ACE-inhibitory peptide classes.

### v5 Failure-Mode Validation

The unchanged v4 frozen recipes were applied to a new 44-peptide exact-sequence-excluded v5 panel. v5 was deliberately harder: 42 of 44 peptides do not end in proline, 18 of 44 are longer than five residues, and the panel spans 11 source-material groups.

- Fresh universe: `57,547` blinded candidates.
- Hidden positives: `44`.
- Matched decoys: `3,520`.
- Random background: `50,000`.
- Digestome background: `3,983`.
- Exact-sequence overlap audit: `0` overlaps against seed, expanded, H12, v3, v4, and all AHTPDB-derived training panels.

| Frozen endpoint | Purpose | Top 1 percent recovery | Top 5 percent recovery | Top 10 percent recovery | Matched-decoy AUROC | Interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| `early-enrichment-v1` | aggressive early recovery | `1/44`; p `0.358` | not primary | `16/44` | `0.6217` | Failed early endpoint. |
| `broad-shortlist-v1` | broad shortlist validation | `2/44`; p `0.0718` | `11/44`; p `7.95e-06` | `18/44`; p `7.70e-08` | `0.7023` | Mixed: broad enrichment persisted, but top-1 recovery and matched-decoy AUROC were weaker than v4. |

After unblinding, v5 was converted into development data. A v5 two-stage search selected two candidate v6 recipes:

- broad candidate: `8/44` top 1 percent, `20/44` top 10 percent, matched-decoy AUROC `0.6313`;
- balanced candidate: `8/44` top 1 percent, `17/44` top 10 percent, matched-decoy AUROC `0.7306`.

These are not validation results. They are frozen under `docs/publication-engine/papers/01-ace-rediscovery/frozen-recipes/2026-08-29-v6-freeze-001/` and require a new v6 holdout.

### Frozen v6 Independent Validation

After the v5 development cycle, two candidate v6 recipes were frozen before v6 literature curation. The new 52-peptide panel spans 18 primary-reference groups and passed exact-sequence exclusion against seed, expanded, H12, v3, v4, v5, and all AHTPDB-derived training panels.

- Fresh universe: `58,192` blinded candidates.
- Hidden positives: `52`.
- Matched decoys: `4,160`.
- Random background: `50,000`.
- Digestome background: `3,980`.
- Exact-sequence overlap audit: `0` overlaps.

| Frozen endpoint | Purpose | Top 1 percent recovery | Top 5 percent recovery | Top 10 percent recovery | Matched-decoy AUROC | Interpretation |
| --- | --- | --- | --- | --- | --- | --- |
| `ace_v6_candidate_from_v5_broad_shortlist_v1` | broad top-10 recovery | `6/52`; p `1.34e-05` | `13/52`; p `1.21e-06` | `20/52`; p `5.10e-08` | `0.6868` | Positive broad enrichment, but hard-decoy AUROC remains weaker. |
| `ace_v6_candidate_from_v5_balanced_v1` | early recovery with stronger hard-decoy discrimination | `6/52`; p `1.34e-05` | `16/52`; p `2.73e-09` | `18/52`; p `1.46e-06` | `0.7642` | Cleanest v6 validation readout: positive early enrichment plus matched-decoy AUROC above `0.75`. |

v6 supports a positive repeated-rediscovery claim after the v5 failure-mode cycle. The caveat remains important: v6 is exact-sequence independent but short-peptide heavy, so it supports blinded prioritization and leakage-controlled rediscovery rather than a general ACE activity predictor claim.

## Best Claims Today

Strongest repeated frozen-validation claim:

> PepLab prospectively froze ACE peptide prioritization recipes, tested them on exact-sequence-excluded primary-literature holdouts, observed broad top-10 enrichment on v4, diagnosed a harder-panel failure mode on v5, and then validated a post-v5 balanced frozen recipe on v6 with `6/52` positives in the top 1 percent of a `58,192`-candidate blinded universe, an `11.54x` enrichment over random expectation with hypergeometric p `1.34e-05` and matched-decoy AUROC `0.7642`.

Strongest single broad-shortlist claim:

> PepLab prospectively froze a two-stage ACE peptide prioritization recipe and, on an exact-sequence-excluded independent primary-literature v4 panel, recovered `16/31` known active food-derived ACE-inhibitory peptides in the top 10 percent of a `56,494`-candidate blinded universe, a `5.16x` enrichment over random expectation with hypergeometric p `6.76e-09` and matched-decoy AUROC `0.7575`.

Claim-safe, target-agnostic:

> In a blinded 63,023-candidate search space, PepLab's target-agnostic oral peptide prioritization mode recovered three known working ACE-inhibitory food-derived peptides, `VPP`, `IPP`, and `LKP`, in the top 1 percent without using ACE activity labels or ACE target priors; the broader 81-peptide ACE literature panel was enriched `9.86x` in the top 1 percent.

Stronger, but target-informed:

> In the same blinded search space, an ACE-target-informed but label-free pharmacophore mode recovered all four primary canonical ACE-inhibitory peptide anchors, `VPP`, `IPP`, `LKP`, and `LKPNM`, in the top 1 percent and enriched the broader ACE literature panel `13.56x`.

Do not merge these two claims. They answer different leakage questions.

Important caveat:

> The all-four-anchor ACE-target-informed result depends materially on terminal-proline/proline-position priors. When direct terminal-proline motif priors are removed, broad enrichment remains statistically positive but the four-anchor recovery falls to `1/4`.

Post-freeze caveat:

> On a newly curated 12-peptide primary-literature holdout, the broad claim weakens. The strongest holdout mode recovers `MRW` and `VTY` in the top 1 percent, while the original ACE pharmacophore mode fails top 1 percent recovery. This makes the current manuscript a rigorous rediscovery/enrichment paper, not a general predictor paper.

v3 caveat:

> On a larger 51-peptide primary-literature holdout, all frozen scoring modes failed the broad top-1-percent recovery gate. The ACE-target-informed modes reached matched-decoy AUROC around `0.75`, but top-rank recall remained too low for a general predictor claim.

v3 development iteration:

> After v3 was converted into a development stress set, two-stage reranking improved early enrichment to `9/51` known positives in the top 1 percent (`17.62x`, hypergeometric p `1.99e-09`) and produced a separate broad-shortlist candidate with `33/51` positives in the top 10 percent. These are model-selection results only; they require a freshly curated v4 holdout before use as claim-bearing evidence.

v4 frozen-validation caveat:

> The broad-shortlist endpoint validated on v4, but the early top-1-percent endpoint did not. Use broad-shortlist language and avoid top-rank discovery language.

v5 failure-mode caveat:

> A harder exact-sequence-excluded v5 panel retained statistically significant top-10 broad enrichment but failed top-1 recovery. v5 should be reported as a failure-mode stress test and development set, not as a repeated validation pass.

v6 repeated-validation caveat:

> A post-v5 frozen balanced recipe passed on v6 with statistically significant top-1 enrichment and matched-decoy AUROC above `0.75`, but the panel remains short-peptide heavy. The manuscript can claim repeated leakage-controlled rediscovery and prioritization, not general ACE activity prediction.

## Claim Boundary

Do not claim:

- Novel ACE-inhibitory peptide discovery.
- Antihypertensive efficacy.
- General ACE activity prediction across all peptide classes.
- Human vascular benefit.
- Molecular docking validation of ACE inhibition.

The current evidence supports rediscovery and prioritization only.

## Why v2 Is Useful Despite the Failed Composite Gate

The expanded panel deliberately includes weak, long, review-table, and mixed-unit positives. That makes it a stress test. The first sequence-only ACE score strongly recovers short ACE-like peptides but does not yet separate all literature positives from composition-matched decoys.

This is a good reviewer-facing limitation. It shows the system is not simply declaring victory on a hand-picked easy panel.

## Next Gates Before Preprint

1. Manually verify every expanded-panel row against the primary paper where possible.
2. Add structured fields for IC50 value, IC50 unit, assay substrate, species, and in vivo evidence.
3. Freeze a v3 endpoint before running: primary anchors, potent short peptides, broad panel, and matched-decoy discrimination as separate hypotheses.
4. Freeze the target-agnostic and ACE-target-informed modes before any further positive-panel expansion.
5. Convert the v3 panel into structured assay metadata fields rather than notes-only rows.
6. Add source-protein split controls so naturally generated candidates from the same parent protein cannot drive apparent rediscovery.
7. Rebuild the ACE structural tier with corrected protonation, conformer ensembles, and zinc-aware controls before using it as evidence.
8. Benchmark against at least one public bioactive-peptide predictor or simple published baseline.
9. Produce manuscript figures from versioned artifacts.
10. Archive code and outputs with a fixed commit hash.
11. Freeze v4/v5/v6 artifact hashes against a commit and refresh all manuscript figures and summary tables from versioned outputs.

## Submission Positioning

Near-term realistic target:

- Bioinformatics Advances Original Article built around ProteinLM-augmented blinded ACE rediscovery.
- Fallback journal fit: BMC Bioinformatics, then Scientific Reports.

Nature-family path:

- Expand this into a resource plus benchmark paper: broad ACE/vascular peptide atlas, leakage-controlled benchmark, and omics-linked graph.
- Or connect it to replicated human metabolomics/genomics through GEMS.
