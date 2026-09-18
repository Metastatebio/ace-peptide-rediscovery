# Evidence guide

This directory is the manuscript-facing evidence bundle.

- `core-results/` contains the original panel files, blinded universe, primary rediscovery and
  ablation reports, and their associated result tables.
- `independent-validation/` contains the v6 holdout universe, exact-overlap exclusion audit,
  motif-neighborhood audit, frozen-validation report, two frozen rerank summaries, and a
  PeptideRanker-style baseline summary.
- `proteinlm/` contains the frozen full-universe rank, parent-grouped cross-validation audit, and
  masked-variant consensus outputs.
- `frozen-recipes/` captures the v4 and v6 frozen ranking definitions.
- `figures/` and `tables/` contain manuscript figures and summary tables.

The frozen ProteinLM full-universe ranked candidate file is retained because it is a direct
manuscript-facing follow-up artifact. Redundant ranked-intermediate files for alternate score modes
and transient model-score matrices are intentionally not copied from the roughly 1.2 GB source
experiment tree. This is a compact inspection and submission bundle, not a guarantee that every
intermediate used during development is included.

Blinding keys are intentionally outside this directory in `../restricted/`; they are internal-only
and must never be included in a public export.

