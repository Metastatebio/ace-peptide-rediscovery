# Bioinformatics Advances Cover Letter

Dear Editor,

We are submitting "ProteinLM-Augmented Blinded Rediscovery of Food-Derived ACE-Inhibitory Peptides" for consideration as an Original Article in Bioinformatics Advances.

The manuscript presents a computational biology benchmark for food-derived ACE-inhibitory peptide prioritization. It addresses a recurrent validation problem in peptide AI: apparent performance can be driven by exact-sequence reuse, post-hoc recipe selection, easy motif recovery, or weak negative controls. The study therefore uses blinded candidate universes, exact-sequence-excluded primary-literature holdouts, matched decoys, frozen ranking recipes, motif-neighborhood auditing, and explicit post-unblinding statistics.

The central results combine prospective frozen validation with ProteinLM rescoring. A v4 frozen broad-shortlist recipe recovered 16 of 31 positives in the top 10 percent of a 56,494-candidate universe, with hypergeometric p = 6.76e-09 and matched-decoy AUROC = 0.7575. After a harder v5 stress panel exposed top-rank failure modes, a post-v5 balanced recipe recovered 6 of 52 v6 positives in the top 1 percent of a 58,192-candidate universe, with p = 1.34e-05 and matched-decoy AUROC = 0.7647. ESM2 ProteinLM rescoring improved parent-grouped held-out AUROC from 0.7643 to 0.7835 and improved full-universe matched-decoy AUROC to 0.7817 under a frozen 0.65 ACE / 0.35 PLM fusion.

The manuscript fits Bioinformatics Advances because it contributes an auditable sequence-analysis benchmark, evaluates machine-learning and ProteinLM-supported ranking of biological sequence data, and reports failure modes alongside positive validations. The validated claim is computational prioritization and rediscovery, not wet-lab potency, receptor binding, or therapeutic efficacy.

Funding: none. Competing interests: Oğuzcan Ünver is affiliated with Metastate Bio; this affiliation is disclosed as a potential competing interest. Code and non-third-party derived study artifacts are available in the public v0.1.0 release at https://github.com/Metastatebio/ace-peptide-rediscovery/releases/tag/v0.1.0.

Sincerely,

Oğuzcan Ünver

Metastate Bio Inc, Wilmington, Delaware, USA

can@metastate.bio
