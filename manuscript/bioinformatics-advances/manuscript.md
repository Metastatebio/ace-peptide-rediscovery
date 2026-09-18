# ProteinLM-Augmented Blinded Rediscovery of Food-Derived ACE-Inhibitory Peptides

Authors: Oğuzcan Ünver

Affiliation: Metastate Bio Inc, Wilmington, Delaware, USA

Corresponding author: Oğuzcan Ünver; can@metastate.bio; ORCID: 0009-0007-2023-5084

Keywords: bioactive peptides; ACE inhibition; protein language models; blinded validation; peptide prioritization

## Abstract

Food-derived angiotensin-converting enzyme (ACE)-inhibitory peptides are a useful benchmark for computational peptide prioritization because many active sequences have been reported across heterogeneous source materials. We built leakage-controlled blinded candidate universes containing hidden literature positives, matched decoys, random peptides, and food-protein digestome background. Ranking recipes were frozen before independent holdout curation and evaluated after exact-sequence overlap audits. A v4 frozen broad-shortlist recipe recovered 16 of 31 positives in the top 10 percent of a 56,494-candidate universe. A harder v5 stress panel exposed reduced top-rank recovery. After a post-v5 freeze, the v6 balanced recipe recovered 6 of 52 positives in the top 1 percent of a 58,192-candidate universe, with 11.54-fold enrichment and matched-decoy AUROC 0.7647. ESM2 ProteinLM rescoring improved parent-grouped held-out matched-decoy AUROC from 0.7643 to 0.7835, and a frozen full-universe 0.65 ACE / 0.35 length-calibrated PLM fusion improved matched-decoy AUROC to 0.7817 while preserving 6 of 52 top-1-percent recovery. These results support computational prioritization and rediscovery, not wet-lab potency or receptor-binding claims.

## Introduction

Bioactive peptides derived from food proteins are an important class of functional food and nutraceutical hypotheses. Angiotensin-converting enzyme (ACE)-inhibitory peptides are especially useful for computational validation because decades of food-peptide research have produced many sequence-level reports across dairy, plant, marine, meat, insect, and fermentation sources. This makes ACE inhibition a practical benchmark for testing whether a peptide workflow can recover known biology from a realistic search space.

The central validation problem is leakage. A model that ranks known literature peptides highly may be exploiting exact-sequence reuse, closely related source tables, post-hoc threshold tuning, or target-specific motifs selected after seeing the test set. A useful benchmark must therefore separate model development from validation, blind labels during scoring, freeze recipes before holdout curation, and audit exact-sequence overlap against every panel used for training or endpoint selection.

Here we evaluate PepLab, a peptide analysis workflow for prioritizing food-derived peptide hypotheses. We use ACE-inhibitory peptides not to claim clinical efficacy, but to test rediscovery under controlled conditions. The contribution is threefold: a reproducible blinded candidate-universe construction procedure, a series of exact-sequence-excluded primary-literature holdouts, and frozen two-stage ranking recipes extended with ESM2 ProteinLM rescoring.

## Results

### Blinded Benchmark Construction

Each validation universe included hidden known positives, matched decoys, random peptide background, and food-protein digestome background. The public scoring file contained candidate ID, sequence, and length only. Private key files contained labels, origins, and provenances and were used only after ranking.

The claim-bearing holdouts were v4, v5, and v6. The v4 panel contained 31 exact-sequence-new primary-literature peptides from five source-material groups. The v5 panel contained 44 peptides from 11 source-material groups and intentionally stressed harder peptide classes. The v6 panel contained 52 peptides from 18 source-material groups and was curated only after the v6 recipes were frozen. All three panels passed exact-sequence overlap auditing.

### Frozen Validation And Failure-Mode Reporting

The v4 broad-shortlist recipe recovered 16 of 31 positives in the top 10 percent of a 56,494-candidate universe. This corresponded to 5.16-fold enrichment over random expectation, hypergeometric p = 6.76e-09, and matched-decoy AUROC = 0.7575.

The unchanged v4 broad recipe was then applied to a harder v5 panel. The recipe recovered 18 of 44 positives in the top 10 percent, p = 7.70e-08, but only 2 of 44 positives in the top 1 percent, with matched-decoy AUROC = 0.7023. v5 therefore demonstrated retained broad enrichment but weaker early-rank discrimination. After unblinding, v5 was converted into development data rather than presented as another success-only validation.

After v5 development, two candidate v6 recipes were frozen before v6 curation. The balanced v6 recipe recovered 6 of 52 positives in the top 1 percent, 16 of 52 in the top 5 percent, and 18 of 52 in the top 10 percent, with matched-decoy AUROC = 0.7647. The top-1-percent recovery was 11.54-fold above random expectation with hypergeometric p = 1.34e-05.

The v6 motif-neighborhood audit compared the 52 v6 positives to 1,356 unique prior positive sequences from 11 prior panels. Exact overlaps remained zero. Among the six top-1-percent v6 positives, one had a one-edit prior neighbor and five were classified by shared C-terminal motif; none was distant from prior positive panels. This narrows the claim to exact-sequence-excluded rediscovery and prioritization rather than de novo molecular-mechanism discovery.

### ProteinLM Rescoring

ESM2 rescoring added an independent sequence-language layer to the v6 benchmark. In parent-grouped repeated five-fold cross-validation, ESM2-8M with length calibration improved mean held-out AUROC from 0.7643 for ACE-only to 0.7835 for the fused score, with 80 of 100 held-out folds outperforming ACE-only. ESM2-35M did not materially improve the result, and ESM2-150M gave only a modest improvement.

Using the cross-validated weight without refitting, the full 58,192-candidate v6 universe was rescored with ESM2-8M and ranked by a frozen 0.65 ACE / 0.35 length-calibrated PLM fusion. The fused rank preserved 6 of 52 positives in the top 1 percent, improved matched-decoy AUROC from 0.7647 to 0.7817, improved all-negative AUROC from 0.7926 to 0.8129, and moved the best positive from rank 129 to rank 45.

The full-universe PLM fusion also prioritized a digestome-first computational shortlist led by AGPAGP, LHLPLP, HLPLPL, APGPAG, and AGPPGFPGA. A separate cross-checkpoint masked-variant consensus across ESM2-8M, ESM2-35M, and ESM2-150M prioritized MSSP, MLSP, MFLLQ, MYSL, and MFLPL. These are candidate hypotheses for wet-lab or receptor-conditioned evaluation, not validated inhibitors.

### Baseline Comparator

The PeptideRanker-style general bioactivity comparator recovered 0 of 31 v4 positives, 0 of 44 v5 positives, and 1 of 52 v6 positives in the top 1 percent. Its matched-decoy AUROC values were 0.6235, 0.4937, and 0.5585. These results show that the frozen PepLab recipes are not explained by generic bioactive-peptide composition alone.

## Discussion

This study demonstrates that peptide prioritization can be evaluated under a disciplined prospective protocol: freeze the ranking recipe, curate the holdout after freeze, hide labels during scoring, audit exact-sequence overlap, and unblind only after ranking. Under that protocol, PepLab repeatedly enriched independently reported ACE-inhibitory peptides from large blinded search spaces.

The strongest evidence is the combination of v4, v6, and the ESM2 extension. v4 validated broad-shortlist enrichment. v5 exposed a meaningful hard-panel failure mode. v6 showed that a post-v5 balanced recipe, frozen before holdout curation, could recover exact-sequence-new positives in the top 1 percent while preserving matched-decoy AUROC above 0.75. ESM2 length-calibrated rescoring then improved hard-decoy discrimination without changing the frozen top-1-percent endpoint.

The validated claim is computational prioritization and rediscovery. The study does not establish new ACE-inhibitory activity, antihypertensive efficacy, universal ACE peptide prediction, or molecular-mechanism discovery. Exact-sequence independence also does not eliminate every possible relationship among source proteins, peptide families, assay types, motif neighborhoods, or literature domains. The current work provides a reproducible foundation for prospective wet-lab and receptor-conditioned follow-up.

## Methods

### Positive-Panel Curation

Known-positive panels were curated from primary-literature reports of ACE-inhibitory peptide sequences. Rows record sequence, target, activity class, source material, generation context, evidence type, reference URL, and notes. v4, v5, and v6 panels were curated after the relevant recipe freeze and were exact-sequence-audited before scoring.

### Development Training Data And Features

Development-only ACE training positives were derived from the AHTPDB IC50 export (Kumar et al., 2015). The acquisition script retained canonical peptides of length 2--16 residues with a positive, parseable reported IC50, deduplicated sequences by the lowest parsed IC50, and excluded evaluation-panel sequences. The final v6 rank component was a length-stratified L2-regularized logistic-regression score trained on AHTPDB positives together with the pre-v6 expanded and H12 development panels and generated background negatives. Its feature families were frozen ACE-heuristic scores, amino-acid composition, N- and C-terminal categories, edge n-grams, internal dipeptide and tripeptide frequencies, and simple motif flags. AHTPDB data were used only for development; all reported v4--v6 positive panels were exact-sequence-excluded from every training panel.

### Candidate Universe Construction

Validation universes were built with `scripts/ace_build_blinded_universe.py`. Each universe contained hidden known positives, matched decoys, random peptide background, and food-protein digestome background. The public scoring table exposed only candidate identifiers, peptide sequences, and sequence lengths during ranking.

### Frozen Ranking Recipes

Frozen recipes specify a gate component, gate percentage, rank component, required model-score directories, required baseline modes, and fusion weights. Recipe application used `scripts/ace_apply_frozen_rerank_recipe.py` and did not search thresholds, weights, or model identities on validation panels.

### ProteinLM Scoring And Fusion

ProteinLM scenarios used ESM2 (Lin et al., 2023) and were run with `scripts/ace_protein_lm_fusion_scenario.py`. Masked pseudo-log-likelihood was computed by masking each residue position, scoring the observed amino acid under the model, and averaging log-likelihood across peptide positions. Raw pseudo-log-likelihood and length-calibrated z scores were evaluated. Fusion weights were selected inside parent-grouped cross-validation with `scripts/ace_plm_cross_validated_fusion_audit.py`, keeping each literature positive and its composition-shuffled or length-matched decoys in the same fold. The full-universe ESM2-8M run then used the cross-validated 0.65 ACE / 0.35 length-calibrated PLM weight without refitting.

### Statistical Analysis

Top-1-percent, top-5-percent, and top-10-percent recovery were computed after unblinding. Hypergeometric survival probabilities measured the chance of observing at least the recovered number of hidden positives within each top-ranked subset under random ranking. AUROC versus all negatives and AUROC versus matched decoys were computed from unblinded labels and scores.

### Use Of AI-Assisted Tools

Language-model tools were used to assist code implementation, artifact generation, and manuscript drafting. The author remains responsible for the study design, curation, code, analysis, interpretation, and final manuscript.

## Data Availability

Analysis code, frozen recipes, post-unblinding evaluation labels, non-third-party derived results, and a reproducibility manifest are available at https://github.com/Metastatebio/ace-peptide-rediscovery. AHTPDB-derived training inputs are not redistributed; they can be reconstructed from the cited upstream database using the included acquisition script and source manifest. Private keys were unavailable during scoring and are released only as post-unblinding evaluation artifacts.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Competing Interests

Oğuzcan Ünver is affiliated with Metastate Bio. This affiliation is disclosed as a potential competing interest; no other competing interests are declared.

## Author Contributions

Oğuzcan Ünver conceived the study, curated the benchmark panels, implemented the analyses, interpreted the results, and wrote the manuscript.

## Acknowledgements

Not applicable.

## References

1. Kumar R, Chaudhary K, Sharma M, Nagpal G, Chauhan JS, Singh S, Gautam A, Raghava GPS. AHTPDB: a comprehensive platform for analysis and presentation of antihypertensive peptides. Nucleic Acids Research. 2015;43:D956--D962. https://doi.org/10.1093/nar/gku1141
2. Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science. 2023;379:1123--1130. https://doi.org/10.1126/science.ade2574
3. Mooney C, Haslam NJ, Pollastri G, Shields DC. Towards the improved discovery and design of functional peptides: common features of diverse classes permit generalized prediction of bioactivity. PLOS ONE. 2012;7:e45012. https://doi.org/10.1371/journal.pone.0045012
4. Minkiewicz P, Iwaniak A, Darewicz M. BIOPEP-UWM database of bioactive peptides: current opportunities. International Journal of Molecular Sciences. 2019;20:5978. https://pmc.ncbi.nlm.nih.gov/articles/PMC6928608/
5. Yang et al. ACE-inhibitory peptides from Larimichthys crocea protein. Molecules. 2024. https://doi.org/10.3390/molecules29051134
6. Zhang et al. ACE-inhibitory peptides from Flammulina velutipes. Foods. 2025. https://doi.org/10.3390/foods14152619
7. Takifugu flavidus ACE-inhibitory peptide study. Marine Drugs. 2023. https://doi.org/10.3390/md21100522
8. Broccoli protein ACE-inhibitory peptide study. Journal of Agricultural and Food Chemistry. 2019. https://doi.org/10.1021/acs.jafc.9b01137
9. Cangkuk fermented beef ACE-inhibitory peptide study. Animal Bioscience. 2024. https://doi.org/10.5713/ab.23.0433
10. Porcine liver and placenta ACE-inhibitory peptide study. Molecules. 2025. https://doi.org/10.3390/molecules30030754
11. Tenebrio molitor protein ACE-inhibitory peptide study. Food Science and Human Wellness. 2026. https://doi.org/10.26599/FSHW.2025.9250609

## Figure Legends

Figure 1. Frozen ACE holdout recovery. Top-1-percent, top-5-percent, and top-10-percent recovery for the claim-bearing frozen validation runs.

Figure 2. Matched-decoy discrimination. AUROC versus matched decoys for v4 broad, v5 broad, v6 broad, v6 balanced, and v6 balanced plus ESM2-8M length-calibrated PLM fusion.

Figure 3. Primary-literature panel length mix. Counts of short peptides, defined as length <= 5 residues, and longer peptides across v3, v4, v5, and v6 panels.
