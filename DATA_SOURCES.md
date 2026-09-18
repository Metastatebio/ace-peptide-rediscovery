# Data Sources and Redistribution Boundaries

## Included material

The public release includes Metastate's code, frozen recipes, manuscript-facing derived results, sequence-level curation, source citations, and labels released after unblinding. The label files contain only candidate identifiers, peptide sequences, generated class labels, origins, and provenance strings; they contain no participant-level or personal data.

## Upstream sources not redistributed

### AHTPDB training source

The training-panel builder retrieves AHTPDB's IC50 export from:

`http://crdd.osdd.net/raghava/ahtpdb/downloads/pepic50.txt`

AHTPDB describes itself as a manually curated database of experimentally validated antihypertensive peptides and requests citation of Kumar et al. (2015), *Nucleic Acids Research* 43:D956-D965. Raw AHTPDB exports and processed AHTPDB panel CSVs are excluded from this repository because no explicit redistribution licence was identified during release preparation. The included script, documented source URL, parser, upstream citation, and manifest allow a user with lawful source access to reconstruct the panel.

### Primary literature

The curation tables preserve peptide-level annotations and source URLs, not article full text. Each source remains subject to its publisher's terms.

### UniProt sequences and ESM2 weights

Source-protein accessions and peptide candidates are documented in the analysis artifacts. Users who fetch UniProt records or ESM2 weights must follow the relevant upstream terms. No model weights are included here.

## Release-specific QC change

The prior AHTPDB parser treated a reported `0 uM` as a valid value. In this release, non-positive reported IC50 values are rejected. A sensitivity rerun excluding 22 such records preserved the v6 top-1-percent endpoint and produced the corrected values reported in the manuscript and `results/`.
