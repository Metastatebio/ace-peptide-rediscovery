# Bioinformatics Advances Submission Metadata

Date prepared: 2026-08-30

## Journal

Bioinformatics Advances

## Article Type

Original Article

## Title

ProteinLM-Augmented Blinded Rediscovery of Food-Derived ACE-Inhibitory Peptides

## Author

Oğuzcan Ünver

## Affiliation

Metastate Bio Inc, Wilmington, Delaware, USA

## Corresponding Author

Oğuzcan Ünver; can@metastate.bio; ORCID: 0009-0007-2023-5084

## Abstract

Food-derived angiotensin-converting enzyme (ACE)-inhibitory peptides are a useful benchmark for computational peptide prioritization because many active sequences have been reported across heterogeneous source materials. We built leakage-controlled blinded candidate universes containing hidden literature positives, matched decoys, random peptides, and food-protein digestome background. Ranking recipes were frozen before independent holdout curation and evaluated after exact-sequence overlap audits. A v4 frozen broad-shortlist recipe recovered 16 of 31 positives in the top 10 percent of a 56,494-candidate universe. A harder v5 stress panel exposed reduced top-rank recovery. After a post-v5 freeze, the v6 balanced recipe recovered 6 of 52 positives in the top 1 percent of a 58,192-candidate universe, with 11.54-fold enrichment and matched-decoy AUROC 0.7647. ESM2 ProteinLM rescoring improved parent-grouped held-out matched-decoy AUROC from 0.7643 to 0.7835, and a frozen full-universe 0.65 ACE / 0.35 length-calibrated PLM fusion improved matched-decoy AUROC to 0.7817 while preserving 6 of 52 top-1-percent recovery. These results support computational prioritization and rediscovery, not wet-lab potency or receptor-binding claims.

## Keywords

bioactive peptides; ACE inhibition; protein language models; blinded validation; peptide prioritization

## Declarations

- Funding: This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.
- Competing interests: Oğuzcan Ünver is affiliated with Metastate Bio. This affiliation is disclosed as a potential competing interest; no other competing interests are declared.
- Data availability: Analysis code, frozen recipes, post-unblinding evaluation labels, non-third-party derived results, and a reproducibility manifest are available in the public `v0.1.0` release at https://github.com/Metastatebio/ace-peptide-rediscovery/releases/tag/v0.1.0. AHTPDB-derived training inputs are reconstructed from the cited upstream source using the included acquisition script and source manifest.
- Ethics: Not applicable; this is a computational study based on literature-derived peptide evidence and generated analysis artifacts.
- AI-assisted tools: Language-model tools were used to assist code implementation, artifact generation, and manuscript drafting. The author remains responsible for study design, curation, code, analysis, interpretation, and final manuscript.

## Submission Blockers

- Author attestation of authorship, funding, competing-interest, ethics, and AI-use declarations.
- Assay-metadata QC attestation after review of the documented mixed-unit and fixed-concentration records.
- Public `v0.1.0` release verified; an archival DOI is optional rather than a submission requirement.
