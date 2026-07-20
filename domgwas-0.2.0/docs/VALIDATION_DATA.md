# Validation data and reproducibility

This repository includes the aggregate validation report, figures, audit
summary, source tests, and a fully synthetic runnable example. It intentionally
does not include the private or potentially restricted rat genotype and
phenotype files used for the real-data demonstration.

## Included

- `VALIDATION_REPORT.md`: complete statistical validation and ADDO comparison.
- `publication_result_audit.json`: machine-readable counts and invariants.
- `figures/`: figures referenced by the validation report.
- `../tests/`: direct GLS and core behavioral tests.
- `../examples/quickstart.py`: public synthetic end-to-end example.

## Not included

- `r11.2.1.bed/.bim/.fam` genotype data.
- `ALLTRAITSALLNORMALIZES_07apr26 .csv` phenotype data.
- Large per-variant Parquet outputs.
- The locally patched ADDO installation and third-party executables.

The real-data inputs must only be deposited if their owners and the applicable
animal-study approvals permit public redistribution. A manuscript should give
an accession, DOI, or controlled-access procedure for those inputs.

The validation report distinguishes actual ADDO executions from an in-Python
experiment that reproduces an ADDO-like covariance structure. It also records
the limits of the experimental low-rank approximation at large sample sizes.
