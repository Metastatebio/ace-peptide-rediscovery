# Manuscript Skeleton: ACE Rediscovery

## Working Title

Blinded rediscovery of antihypertensive food-derived ACE-inhibitory peptides with an evidence-aware peptide AI workflow

## Abstract Structure

Background:
Food-derived ACE-inhibitory peptides are a mature literature test case for computational peptide prioritization.

Methods:
We constructed a blinded candidate universe from food-protein digestomes and matched decoys, froze a literature-positive set, ranked candidates with an ACE-specific PepLab module, and unblinded after ranking.

Results:
Pilot v1 recovered all four primary anchors in the top 1 percent of a 16,305-candidate blinded universe. Expanded v2 recovered 11 of 81 literature positives in the top 1 percent of a 63,023-candidate universe, a 13.56-fold enrichment over random expectation by top-1-percent recovery. After a v3 development cycle, a two-stage broad-shortlist recipe was frozen before v4 curation and recovered 16 of 31 exact-sequence-new primary-literature positives in the top 10 percent of a 56,494-candidate blinded universe, a 5.16-fold enrichment over random expectation with hypergeometric p = 6.76e-09 and matched-decoy AUROC = 0.7575. A harder v5 panel produced mixed broad enrichment and was treated as development data. A post-v5 frozen balanced recipe then recovered 6 of 52 exact-sequence-new v6 positives in the top 1 percent of a 58,192-candidate blinded universe, an 11.54-fold enrichment with hypergeometric p = 1.34e-05 and matched-decoy AUROC = 0.7642.

Conclusions:
PepLab recovered known working ACE peptide biology from blinded search spaces and provides a reproducible route from food proteins to vascular-function peptide hypotheses. The validated result is leakage-controlled rediscovery and prioritization, not novel peptide discovery or general ACE activity prediction.

## Introduction

- Food-derived bioactive peptides are difficult to prioritize because source protein, proteolysis, stability, target biology, and assay evidence are scattered.
- ACE inhibition is a strong first benchmark because multiple food-derived peptides have reproducible literature support.
- The study tests rediscovery, not new efficacy.

## Methods

- Literature-positive freeze.
- Digestome and decoy universe construction.
- Blinding procedure.
- ACE-specific feature extraction.
- Ranking and statistical endpoints.
- Unblinding and robustness analyses.
- Reproducibility and code release.

## Results

- Candidate universe composition.
- Known-positive rank recovery.
- Enrichment versus decoys.
- Feature contributions and ablations.
- Evidence graph for recovered peptides.

Current result language:

- Four-anchor pilot: `VPP`, `IPP`, `LKPNM`, and `LKP` all recovered in the top 1 percent.
- Expanded literature panel: 81 known positives tested; 11 recovered in the top 1 percent; top-1-percent hypergeometric p-value `5.98e-10`.
- Strongest subset result: the four primary anchors recovered 4/4 in the top 1 percent; short peptides recovered 11/48 in the top 1 percent.
- Limitation: the expanded run does not yet meet the stricter matched-decoy AUROC gate.
- Frozen v4 validation: 31 exact-sequence-new primary-literature positives tested; broad-shortlist recipe recovered 16/31 in the top 10 percent; top-10 hypergeometric p-value `6.76e-09`; matched-decoy AUROC `0.7575`.
- Limitation: the v4 early-enrichment endpoint recovered only 1/31 in the top 1 percent and failed.
- v5 failure-mode validation: 44 exact-sequence-new primary-literature positives tested; the unchanged v4 broad-shortlist recipe recovered 18/44 in the top 10 percent but only 2/44 in the top 1 percent, so v5 became development data rather than a repeated-validation pass.
- Frozen v6 validation: 52 exact-sequence-new primary-literature positives tested; the post-v5 frozen balanced recipe recovered 6/52 in the top 1 percent and 18/52 in the top 10 percent with matched-decoy AUROC `0.7642`.

## Discussion

- What PepLab demonstrated.
- Why rediscovery is a meaningful validation step.
- Why the result does not establish novel peptide activity.
- How this leads to wet-lab validation and omics-linked peptide systems.

## Claim Language

Use:

> PepLab rediscovered literature-validated ACE-inhibitory peptides from a blinded food-derived peptide search space.

For the frozen v4 result, use:

> PepLab prospectively froze a two-stage ACE peptide prioritization recipe and enriched exact-sequence-new known active food-derived ACE-inhibitory peptides in the top 10 percent of an independent blinded validation universe.

For the repeated v6 result, use:

> After a harder v5 stress panel exposed top-rank failure modes, PepLab froze a revised balanced recipe before v6 curation and recovered exact-sequence-new primary-literature ACE-inhibitory peptides in the top 1 percent of a new blinded universe while preserving matched-decoy AUROC above 0.75.

For the expanded v2 result, use:

> PepLab significantly enriched known ACE-inhibitory literature positives in the top 1 percent of a blinded 63,023-candidate sequence universe, with strongest recovery for short ACE-like peptides and the four primary anchors.

Avoid:

> PepLab discovered antihypertensive peptides.

Also avoid:

> PepLab generally predicts ACE inhibition across all reported ACE-inhibitory peptides.

## Target Journals

- Scientific Reports or Communications Biology for a focused methods validation.
- Nature Food or Nature Communications only if scaled into a broader digestome/resource/benchmark package.
