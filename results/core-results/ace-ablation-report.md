# ACE Rediscovery Ablation Report

Date: 2026-08-29

This is a post-hoc sensitivity analysis over the same blinded ranking table. It asks whether known-positive recovery depends on specific ACE-position priors.

Total candidates: `16305`  
Top 1 percent cutoff rank: `164`

| Ablation | Top 1% | Top 10 | Top 25 | Top 100 | Mean positive rank | AUROC all | AUROC matched | AUPRC all |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| full | 4 | 2 | 3 | 3 | 31.25 | 0.9981 | 0.9930 | 0.7596 |
| no_prodrug_triad | 3 | 0 | 0 | 3 | 109.00 | 0.9933 | 0.9812 | 0.0646 |
| no_cterm_proline_bonus | 4 | 0 | 3 | 4 | 27.75 | 0.9985 | 0.9953 | 0.1129 |
| no_position_priors | 0 | 0 | 0 | 0 | 3983.75 | 0.7557 | 0.4965 | 0.0016 |
| composition_only | 1 | 0 | 0 | 1 | 343.50 | 0.9794 | 0.8797 | 0.0108 |
| random_hash_baseline | 0 | 0 | 0 | 0 | 9671.25 | 0.4069 | 0.3875 | 0.0003 |

Interpretation:

- `full` is the registered pilot score.
- `no_prodrug_triad` removes the N-terminal hydrophobic/basic/proline pattern that helps `LKPNM`.
- `no_cterm_proline_bonus` removes only the terminal-proline bonus.
- `no_position_priors` removes most sequence-position priors and leaves mainly global composition/developability features.
- `composition_only` is a weaker baseline using length, hydrophobicity, proline fraction, charge, and penalties.
- `random_hash_baseline` is a deterministic random-ranking control.
