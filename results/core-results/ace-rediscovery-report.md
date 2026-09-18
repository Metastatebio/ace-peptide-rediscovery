# ACE Rediscovery Pilot Report

Date: 2026-08-29

## Claim Boundary

This is a computational pilot. It tests whether a sequence-only ACE-specific PepLab score can recover known literature-positive ACE inhibitory peptides after blinding. It does not demonstrate new wet-lab activity or clinical efficacy.

## Primary Endpoint

Success criterion: all configured gates below must pass.

Configured gates:

- Minimum known positives in top 1 percent: `2`
- Minimum top 1 percent recovery fraction: `0.0`
- Minimum AUROC versus matched decoys: `0.0`

- Total candidates: `16305`
- Top 1 percent cutoff rank: `164`
- Known positives: `4`
- Known positives in top 1 percent: `4`
- Top 1 percent recovery fraction: `1.0000`
- Top 1 percent enrichment over random: `99.42x`
- Top 1 percent hypergeometric p-value: `9.87e-09`
- Primary endpoint passed: `yes`

## Recovery Metrics

- Top 10 recall: `2/4`
- Top 25 recall: `3/4`
- Top 50 recall: `3/4`
- Top 100 recall: `3/4`
- Mean positive rank: `31.25`
- Mean matched-decoy rank: `5897.92`
- AUROC versus all negatives: `0.9981`
- AUROC versus matched decoys: `0.9930`
- Average precision versus all negatives: `0.7596`

## Label Counts

| Label | Count |
| --- | --- |
| digestome_background | 3981 |
| digestome_known_positive | 2 |
| known_positive | 2 |
| matched_decoy | 320 |
| random_background | 12000 |

## Known Positive Ranks

| Sequence | Rank | ACE score | Length | Label | Origin |
| --- | --- | --- | --- | --- | --- |
| IPP | 4 | 1.0 | 3 | digestome_known_positive | milk:P02666:thermolysin_like |
| LKP | 6 | 1.0 | 3 | known_positive | blinded_known_positive_spike |
| VPP | 11 | 1.0 | 3 | digestome_known_positive | milk:P02666:thermolysin_like |
| LKPNM | 104 | 0.836647 | 5 | known_positive | blinded_known_positive_spike |

## Top 25 After Unblinding

| Rank | Sequence | ACE score | Length | Label |
| --- | --- | --- | --- | --- |
| 1 | AKP | 1.0 | 3 | digestome_background |
| 2 | IHP | 1.0 | 3 | digestome_background |
| 3 | IKP | 1.0 | 3 | digestome_background |
| 4 | IPP | 1.0 | 3 | digestome_known_positive |
| 5 | LHP | 1.0 | 3 | random_background |
| 6 | LKP | 1.0 | 3 | known_positive |
| 7 | LPP | 1.0 | 3 | digestome_background |
| 8 | LRP | 1.0 | 3 | matched_decoy |
| 9 | MHP | 1.0 | 3 | random_background |
| 10 | VKP | 1.0 | 3 | random_background |
| 11 | VPP | 1.0 | 3 | digestome_known_positive |
| 12 | VRP | 1.0 | 3 | matched_decoy |
| 13 | AKPA | 1.0 | 4 | digestome_background |
| 14 | IHPF | 1.0 | 4 | digestome_background |
| 15 | IKPL | 1.0 | 4 | digestome_background |
| 16 | LRPW | 1.0 | 4 | random_background |
| 17 | MKPW | 1.0 | 4 | digestome_background |
| 18 | VKPA | 1.0 | 4 | random_background |
| 19 | VKPW | 1.0 | 4 | random_background |
| 20 | VPPF | 1.0 | 4 | digestome_background |
| 21 | IAPP | 0.935424 | 4 | random_background |
| 22 | VLPP | 0.935424 | 4 | random_background |
| 23 | VVPP | 0.935424 | 4 | digestome_background |
| 24 | WLPP | 0.935424 | 4 | random_background |
| 25 | ARPKHP | 0.93291 | 6 | digestome_background |

## Interpretation

If the primary endpoint passed, the defensible claim is:

> PepLab rediscovered literature-validated ACE-inhibitory peptides from a blinded sequence search space in a computational pilot.

Do not claim novel bioactivity. The next gate is a larger universe with stricter source-protein split controls and a no-motif-prior ablation.
