# ProteinLM-Augmented Blinded Rediscovery of Food-Derived ACE-Inhibitory Peptides

This repository is the public reproducibility release for the accompanying computational manuscript. It contains the analysis code, frozen ranking recipes, post-unblinding v4/v5/v6 evaluation labels, manuscript-facing result tables, figures, and a release-verification script.

The paper's claim is deliberately limited to leakage-controlled computational prioritization and rediscovery of literature-reported ACE-inhibitory peptides. It does not report new inhibition experiments, receptor binding, antihypertensive efficacy, or validated therapeutic candidates.

## Release contents

- `scripts/`: ACE benchmark construction, scoring, frozen-reranking, ProteinLM, and audit code.
- `backend/app/analysis/`: minimal shared feature code required by the core scoring utilities.
- `results/`: manuscript-facing derived results, frozen recipes, figures, tables, and the corrected QC outputs.
- `data/labels-post-unblinding/`: labels released after the frozen analyses were complete, allowing reported v4/v5/v6 metrics to be checked.
- `manuscript/`: the Bioinformatics Advances submission source.
- `docs/ASSAY_METADATA_QC.md`: the assay-metadata scope and remaining human-review gate.

## Fast verification

No model download or third-party data is needed to check the headline reported results:

```bash
python3 scripts/verify_release.py
```

The verifier checks the corrected v6 balanced endpoint (6/52 top-1-percent recovery; matched-decoy AUROC 0.7647), the frozen ESM2 fusion endpoint (matched-decoy AUROC 0.7817), and the hashes of the released post-unblinding label files.

## Data boundaries and provenance

This release includes original curation and derived analysis outputs. It does not redistribute raw AHTPDB exports or model weights. The AHTPDB-derived training panel can be reconstructed from the documented upstream download using `scripts/ace_build_ahtpdb_training_panel.py`; the script now rejects zero or negative reported IC50 values. See [DATA_SOURCES.md](DATA_SOURCES.md).

The post-unblinding label files were withheld from the scoring inputs during the frozen analyses. Releasing them now supports independent metric verification; it does not retroactively change the blinded scoring protocol.

## Installation for rerunning code

The lightweight verification script uses the Python standard library only. The broader analysis requires the packages listed in `requirements-core.txt`. ProteinLM scenarios additionally use the pinned optional runtime in `requirements-proteinlm.txt` and require users to obtain model weights under their upstream terms.

## Licensing

Code is licensed under Apache-2.0. Original documentation, figures, curated annotations, and derived non-third-party result tables are made available under CC BY 4.0. Third-party sources and model weights retain their own terms; this release does not grant rights to them.

## Citation

See [CITATION.cff](CITATION.cff). The repository should be cited by its tagged release and archival DOI, once created.
