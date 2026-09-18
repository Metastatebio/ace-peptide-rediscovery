# ACE Assay-Metadata QC Record

Date of automated review: 2026-09-18

## Scope

The structured assay table covers all 178 peptides in the v3, v4, v5, and v6 literature-curated panels. It is a provenance and claim-boundary aid; the manuscript does not assert cross-study potency comparability.

## Automated checks completed

- 178/178 rows were assigned a panel, source URL, assay-readout class, and evidence-strength class.
- 138 rows contain an IC50 expression; 107 are directly comparable after conversion from uM or mM.
- 31 mass-concentration IC50 records (30 mg/mL and 1 ug/mL) are preserved in their reported units and are **not** converted to uM, because a molecular-weight and assay-context review would be required.
- 11 records report fixed-concentration inhibition and remain separately classified from IC50 measurements.
- 29 records report activity without a parsed IC50 and remain separately classified.
- The AHTPDB training-panel parser was corrected to reject zero or negative reported IC50 values. The 22 invalid records were excluded in a frozen-recipe sensitivity rerun.

## Impact assessment

The corrected v6 balanced recipe retained 6/52 top-1-percent recovery, 11.54-fold enrichment, and p = 1.34e-05. Matched-decoy AUROC was 0.7647. The corrected ESM2-8M full-universe fusion had matched-decoy AUROC 0.7817. The full result files and checksums are in the release bundle.

## Required human attestation

Before submission, the corresponding author or a qualified assay-data reviewer must review the 31 mass-concentration rows and the 11 fixed-concentration rows against their cited sources, confirm that the manuscript keeps them out of cross-assay potency claims, and record the reviewer name and date below.

- Reviewer: pending
- Date: pending
- Outcome: pending
- Notes / exceptions: pending
