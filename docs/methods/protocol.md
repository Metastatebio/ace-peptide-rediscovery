# Paper 1 Protocol: ACE Rediscovery

## Thesis

PepLab can recover known working food-derived ACE-inhibitory peptides from a blinded candidate universe. This is the cleanest non-DPP4 "real known peptide" proof because the endpoint is simple, the literature positives are strong, and the assay class is well established.

## Primary Claim

After successful completion:

> PepLab rediscovered literature-validated ACE-inhibitory food-derived peptides from a blinded food-digestome and matched-decoy search space, prioritizing them as vascular-function peptide systems.

## Known Positives

The frozen positive set for the first run is:

| Sequence | Target | Source | Evidence |
| --- | --- | --- | --- |
| `VPP` | ACE | Sour milk / casein | Antihypertensive peptide report in PubMed `7673515` |
| `IPP` | ACE | Sour milk / casein | Antihypertensive peptide report in PubMed `7673515` |
| `LKPNM` | ACE | Dried bonito muscle | Prodrug-type ACE-inhibitory peptide in PubMed `10604535` |
| `LKP` | ACE | Dried bonito muscle | ACE-generated active peptide in PubMed `10604535` |

The positive seed file is `known-positive-seed.csv`.

## Search Space

Build a blinded candidate universe before scoring:

- Food protein sources: milk caseins and whey proteins, bonito and other fish proteins, egg proteins, collagen, soy, and neutral background proteins.
- Processing routes: gastrointestinal enzymes, thermolysin, fermentation-inspired proteases, and matched random cleavage controls.
- Candidate length: 2 to 12 amino acids for the first pass.
- Decoys: shuffled sequence decoys, composition-matched decoys, source-protein matched negatives, and length-matched random peptides.
- Blinding: labels, PubMed URLs, and source-positive flags must not be available to the scoring module.

## ACE-Specific Scoring Module

Add a vascular/ACE scoring module rather than reusing unrelated scoring logic. The module should use interpretable features:

- Short peptide preference.
- C-terminal proline, aromatic, and hydrophobic residue patterns.
- Proline-rich motifs compatible with ACE resistance and binding.
- Basic/hydrophobic balance.
- Absence of reactive cysteine unless literature-supported.
- Protease-release plausibility from source proteins.
- Peptidase-stability and GI survival flags.
- Synthesis ease and safety/developability filters.
- Optional structural tier: ACE docking and binding-site contact consistency.

The scoring manifest must list every feature and prove labels were unavailable.

## Primary Endpoint

The primary endpoint is frozen before unblinding:

- Success if at least two of the four known positives rank in the top 1 percent of the blinded universe, or if all four are significantly enriched versus matched decoys under a pre-specified rank-enrichment test.

Secondary endpoints:

- Top-k recall at 10, 25, 50, and 100.
- AUROC and AUPRC versus decoys.
- Median rank of known positives versus matched decoys.
- Source/protease enrichment.
- Stability of recovery under ablations of individual scoring features.

## Controls

Required controls:

- Composition-matched decoys for every known positive.
- Source-protein matched decoys from the same proteins.
- Scrambled peptides with identical length and residue counts where possible.
- A no-literature-prior run that removes exact motif priors.
- A random-ranking baseline with fixed seeds.

## Figures

1. PepLab ACE rediscovery workflow.
2. Blinded universe composition and source/protease coverage.
3. Rank distribution of known positives versus decoys.
4. Feature-ablation robustness.
5. Evidence graph for `VPP`, `IPP`, `LKPNM`, and `LKP`.

## Failure Criteria

The paper should not be written as a positive result if:

- Known positives are not enriched versus decoys.
- Recovery depends on exact sequence lookup or leaked labels.
- The search space is so narrow that rediscovery is trivial.
- The ranking cannot be reproduced from a frozen commit and manifest.

## Next Implementation Tasks

1. Implement `scripts/ace_build_blinded_universe.py`.
2. Implement `backend/app/analysis/ace_features.py`.
3. Implement `scripts/ace_score_blinded_universe.py`.
4. Implement `scripts/ace_unblind_rediscovery_report.py`.
5. Implement `scripts/ace_run_ablation_experiment.py`.
6. Add tests for blinding, feature values, deterministic ranking, and enrichment metrics.
