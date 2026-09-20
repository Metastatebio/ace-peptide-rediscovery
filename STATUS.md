# Status of v0.1.x

Date: 2026-09-20

## Archived exploratory analysis — not submission-ready

The v0.1.x repository and release tags are preserved as historical analysis
artifacts. A subsequent artifact-level review reproduced several issues that
materially narrow their interpretation:

- archived repeated-CV result blocks are numerically duplicated and cannot
  support the reported fold-win inference;
- the selected v6 PLM fusion was selected and described on the same v6 panel,
  so it is exploratory rather than an independent validation;
- pooled matched-control AUROC is dominated by length-matched random controls
  and does not demonstrate strong discrimination against composition-
  preserving controls;
- some v6 background candidates occur in earlier released positive panels, so
  they must not be interpreted as experimentally inactive peptides.

The original files remain available so the audit can be reproduced. The
release verifier confirms selected hashes and released headline calculations;
it is not a validation of the broader scientific claims above.

Future work will use a separately versioned recovery analysis with explicit
label taxonomy, source-level curation, control-stratified reporting, and a
new pre-specified holdout before making any prospective-validation claim.

