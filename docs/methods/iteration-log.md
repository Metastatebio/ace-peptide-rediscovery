# ACE Rediscovery Iteration Log

Date: 2026-08-29

## Iteration 1: Four-Anchor Pilot

- Universe: `16,305` candidates.
- Hidden positives: `VPP`, `IPP`, `LKPNM`, `LKP`.
- Result: `4/4` in top 1 percent; matched-decoy AUROC `0.9930`.
- Use: feasibility only. This panel is too small for a publication-grade claim.

## Iteration 2: Expanded Literature Panel

- Universe: `63,023` candidates.
- Hidden positives: `81` ACE literature peptides.
- Result: `11/81` in top 1 percent; `13.56x` enrichment; hypergeometric p `5.98e-10`; matched-decoy AUROC `0.7231`.
- Use: honest expanded stress test. It supports enrichment, not broad ACE activity prediction.

## Iteration 3: Leakage-Controlled Score Modes

- Target-agnostic generic mode: `7/81` in top 1 percent; matched-decoy AUROC `0.5776`.
- Target-agnostic oral-stability mode: `8/81` in top 1 percent; `VPP`, `IPP`, and `LKP` recovered in top 1 percent; matched-decoy AUROC `0.6575`.
- ACE-target-informed label-free mode: `11/81` in top 1 percent; `VPP`, `IPP`, `LKP`, and `LKPNM` recovered in top 1 percent; matched-decoy AUROC `0.7202`.
- Use: the strongest current publication line. Keep target-agnostic and target-informed claims separate.

## Iteration 4: Raw Docking Probe

- Structure: human ACE `1O86`; zinc-centered Vina box.
- Panel: 8 positives and 16 hard matched decoys.
- Result: raw Vina AUROC `0.1719`; pose-QC AUROC `0.3047`.
- Protocol issue: PeptideBuilder ligand PDBs lacked explicit terminal `OXT`, causing incomplete terminal chemistry.
- Use: negative control only.

## Iteration 5: Corrected-Terminal Docking Probe

- Change: add explicit C-terminal `OXT` before PDBQT conversion.
- Chemistry check: VPP converted as terminal carboxylic acid `C(=O)O`.
- Result: raw Vina AUROC `0.1250`; pose-QC AUROC `0.1250`.
- Use: strong evidence that the current simple docking tier is not paper-grade. Do not use docking as supporting evidence until the physics protocol is rebuilt.

## Iteration 6: No Exact Terminal-Proline Motif Stress Test

- Change: add `ace_no_exact_proline_motif`, an ACE-target-informed mode that removes direct terminal-proline motif priors.
- Result: `4/81` broad positives in the top 1 percent; `4.93x` enrichment; matched-decoy AUROC `0.6464`.
- Primary anchors: only `LKP` remains in the top 1 percent; `LKPNM` rank `1108`, `IPP` rank `5604`, and `VPP` rank `5657`.
- Use: leakage-boundary test. The all-four-anchor target-informed claim is materially driven by proline-position priors and must be described as target-informed, not target-agnostic discovery.

## Iteration 7: Post-Freeze Primary-Literature Holdout

- Change: curate a new 12-peptide primary-literature holdout after freezing the current scoring code.
- Holdout peptides: `MRWRD`, `MRW`, `IAYKPAG`, `EWL`, `VTY`, `LGVP`, `FIGR`, `FILR`, `FQRL`, `FRAL`, `KFL`, `KLF`.
- Universe: `54,954` candidates with private labels withheld from scoring.
- Results:
  - `generic_zero_shot`: `0/12` top 1 percent; matched-decoy AUROC `0.4769`.
  - `oral_stability_zero_shot`: `1/12` top 1 percent; `LGVP` rank `199`; matched-decoy AUROC `0.5322`.
  - `ace_no_exact_proline_motif`: `2/12` top 1 percent; `MRW` rank `76`, `VTY` rank `330`; matched-decoy AUROC `0.6784`.
  - `ace_pharmacophore_zero_label`: `0/12` top 1 percent; matched-decoy AUROC `0.6685`.
- Use: the current engine shows partial held-out rediscovery, but not enough for a broad predictor claim.

## Iteration 8: v3 Primary-Literature Panel

- Change: expand the post-freeze holdout to `51` non-overlapping ACE-inhibitory peptides from open primary-literature sources without changing the existing scoring modes.
- Sources: Ulva intestinalis, walnut glutelin-1, green coffee, Ziziphus jujuba, Gigantidas platifrons, skipjack tuna, areca nut kernel globulin, and Pyropia pseudolinearis.
- Universe: `58,114` candidates with private labels withheld from scoring.
- Results:
  - `generic_zero_shot`: `1/51` top 1 percent; matched-decoy AUROC `0.5199`.
  - `oral_stability_zero_shot`: `1/51` top 1 percent; matched-decoy AUROC `0.6070`.
  - `ace_no_exact_proline_motif`: `2/51` top 1 percent; matched-decoy AUROC `0.7467`.
  - `ace_pharmacophore_zero_label`: `1/51` top 1 percent; matched-decoy AUROC `0.7530`.
  - `literature_sar`: `1/51` top 1 percent; matched-decoy AUROC `0.7514`.
- Use: the broad holdout gate failed. The current engine has coarse ACE-target signal but insufficient top-rank resolution for a broad predictor claim.

## Iteration 9: v3 Development Model Iteration

- Change: treat v3 as a development stress set and iterate supervised/fusion/two-stage models. Added AHTPDB IC50 training panels, potency thresholds, multi-model fusion, XGBoost, physchem descriptors, length-stratified components, and two-stage reranking.
- Best early-enrichment recipe: gate with within-length `expanded-v2-plus-h12:logreg_l2_c0_05` at `4.85%`, then rerank with `ace_no_exact_proline_motif`.
- Best early-enrichment result on v3 development set: `9/51` top 1 percent; `17.62x` enrichment; hypergeometric p `1.99e-09`; `25/51` top 10 percent; matched-decoy AUROC `0.7977`.
- Best broad-shortlist recipe: gate with physchem XGBoost at `10%`, then rerank with the min3/physchem fusion.
- Best broad-shortlist result on v3 development set: `8/51` top 1 percent; `33/51` top 10 percent; matched-decoy AUROC `0.8016`; all-negative AUROC `0.8766`.
- Use: model selection only. Because v3 was used for threshold and recipe selection, these are not independent validation results. Freeze the selected recipes before curating v4.

## Iteration 10: v4 Frozen Independent Validation

- Change: freeze the two selected v3-derived recipes before curating the v4 panel, then curate a new 31-peptide primary-literature validation set from Larimichthys crocea, sesame protein, Flammulina velutipes, camel casein, and Pelodiscus sinensis sources.
- Independence audit: `31/31` v4 sequences have zero exact-sequence overlap with seed, expanded, H12, v3, and all AHTPDB-derived training panels.
- Universe: `56,494` blinded candidates with `31` hidden positives, `2,480` matched decoys, `50,000` random background peptides, and `3,983` digestome background peptides.
- Frozen early-enrichment endpoint: `1/31` top 1 percent; `7/31` top 10 percent; matched-decoy AUROC `0.6586`. Use: failed early endpoint.
- Frozen broad-shortlist endpoint: `2/31` top 1 percent; `9/31` top 5 percent; `16/31` top 10 percent; top-10 enrichment `5.16x`; top-10 hypergeometric p `6.76e-09`; matched-decoy AUROC `0.7575`.
- Use: claim-bearing broad-shortlist evidence. The clean statement is that PepLab prospectively froze a two-stage ACE peptide prioritization recipe and enriched exact-sequence-new known active food-derived ACE peptides in a blinded universe. Do not claim novel discovery or general ACE prediction.

## Iteration 11: v5 Failure-Mode Validation And v6 Candidate Freeze

- Change: curate a new 44-peptide v5 primary-literature stress panel after v4, emphasizing longer, non-C-terminal-proline, source-diverse, mixed-potency ACE-inhibitory peptides.
- Independence audit: `44/44` v5 sequences have zero exact-sequence overlap with seed, expanded, H12, v3, v4, and all AHTPDB-derived training panels.
- Universe: `57,547` blinded candidates with `44` hidden positives, `3,520` matched decoys, `50,000` random background peptides, and `3,983` digestome background peptides.
- Frozen early-enrichment endpoint: `1/44` top 1 percent; `16/44` top 10 percent; matched-decoy AUROC `0.6217`. Use: failed early endpoint.
- Frozen broad-shortlist endpoint: `2/44` top 1 percent; `11/44` top 5 percent; `18/44` top 10 percent; top-10 hypergeometric p `7.70e-08`; matched-decoy AUROC `0.7023`. Use: mixed stress result, not a repeated-validation pass.
- v5 development iteration after unblinding: best broad two-stage candidate reached `8/44` top 1 percent and `20/44` top 10 percent, but matched-decoy AUROC was only `0.6313`.
- v5 balanced development candidate: `8/44` top 1 percent and `17/44` top 10 percent with matched-decoy AUROC `0.7306`.
- Use: v5 is now development data. The v6 candidate recipes are frozen under `docs/publication-engine/papers/01-ace-rediscovery/frozen-recipes/2026-08-29-v6-freeze-001/` and require a new exact-sequence-excluded holdout before any validation claim.

## Iteration 12: v6 Frozen Independent Validation

- Change: curate a new 52-peptide v6 primary-literature panel after the v6 recipe freeze, excluding exact sequences from seed, expanded, H12, v3, v4, v5, and all AHTPDB-derived training panels.
- Independence audit: `52/52` v6 sequences have zero exact-sequence overlap with all listed exclusion panels.
- Universe: `58,192` blinded candidates with `52` hidden positives, `4,160` matched decoys, `50,000` random background peptides, and `3,980` digestome background peptides.
- Frozen broad v6 candidate: `6/52` top 1 percent; `13/52` top 5 percent; `20/52` top 10 percent; top-1 hypergeometric p `1.34e-05`; top-10 hypergeometric p `5.10e-08`; matched-decoy AUROC `0.6868`.
- Frozen balanced v6 candidate: `6/52` top 1 percent; `16/52` top 5 percent; `18/52` top 10 percent; top-1 hypergeometric p `1.34e-05`; top-5 hypergeometric p `2.73e-09`; matched-decoy AUROC `0.7642`.
- Use: positive repeated frozen-holdout evidence for blinded ACE peptide rediscovery and prioritization. The balanced recipe is the cleaner validation readout because it preserves hard-decoy AUROC above `0.75`; the broad recipe gives slightly higher top-10 recovery but weaker hard-decoy separation.

## Iteration 13: ProteinLM-Augmented v6 Rescoring

- Change: add an ESM2 masked pseudo-log-likelihood layer on top of the frozen v6 balanced recipe. Scores were calibrated within peptide length and fused with the frozen v6 score.
- Local checkpoint scenarios:
  - ESM2-8M raw fusion: best matched-panel AUROC `0.7835` versus v6 `0.7642`.
  - ESM2-35M: did not materially strengthen the result.
  - ESM2-150M: modest support only, best matched-panel AUROC `0.7744`.
- Parent-grouped repeated cross-validation: ESM2-8M length-calibrated fusion averaged held-out AUROC `0.7832` versus ACE-only `0.7624`, delta `+0.0208`, with `100/100` held-out folds beating ACE-only.
- Full-universe frozen application: scored all `58,192` v6 candidates with ESM2-8M and applied the cross-validated `0.65` ACE / `0.35` length-calibrated PLM fusion without refitting.
- Full-universe result: retained `6/52` top-1-percent positives, improved matched-decoy AUROC from `0.7642` to `0.7813`, improved all-negative AUROC from `0.7919` to `0.8124`, moved best positive rank from `121` to `48`, and moved top-100 positives from `0` to `1`.
- Discovery outputs: digestome shortlist led by `AGPAGP`, `LHLPLP`, `HLPLPL`, `APGPAG`, and `AGPPGFPGA`; cross-checkpoint masked-variant consensus led by `MSSP`, `MLSP`, `MFLLQ`, `MYSL`, and `MFLPL`.
- Use: stronger claim-bearing computational evidence. The paper can now be framed as blinded rediscovery plus PLM-augmented prioritization. It still cannot claim receptor binding, potency, or therapeutic efficacy without wet-lab or receptor-conditioned structural validation.

## Current Readout

PepLab is closest to a defensible paper around leakage-controlled sequence rediscovery, benchmark design, and PLM-augmented peptide prioritization. The strongest current claim is now repeated frozen-holdout prioritization plus an orthogonal ESM2 layer: v4 validated broad top-10 enrichment, v5 exposed a harder-panel failure mode, v6 validated a post-v5 frozen balanced recipe with `6/52` exact-sequence-new positives in the top 1 percent, and ESM2-8M length-calibrated fusion improved v6 hard-decoy AUROC to `0.7813` on the full universe while preserving top-1-percent recovery. The paper should be positioned as rigorous blinded rediscovery and PLM-augmented computational prioritization with explicit failure-mode reporting, not as a wet-lab discovery or receptor-binding paper.

## Next Iteration

Do not tune on v6 if it remains claim-bearing. The next technical iteration should be reviewer-facing controls and external validation: source-protein split controls, official public-predictor baselines where executable access is available, preregistered v7 validation for the `0.65/0.35` PLM fusion, and wet-lab or receptor-conditioned structural evaluation of the digestome and masked-variant shortlists.
