# domgwas


`domgwas` is a pure-Python, dominance-aware GWAS pipeline for quantitative
traits. It streams PLINK 1 binary genotypes, constructs additive
leave-one-chromosome-out (LOCO) relationship matrices, estimates a variance
component by profile REML, and tests additive and dominance terms using
complete-design generalized least squares (GLS).

Version 0.2.0 corrects the phenotype-only transformation used in version
0.1.0. The phenotype, intercept, covariates, additive dosage, and dominance
code are now transformed with the same chromosome-specific covariance
operator.

## Features

- PLINK 1 `.bed/.bim/.fam` input with block-wise genotype streaming.
- Additive and conditional dominance association tests.
- Additive LOCO mixed-model correction with profile REML.
- Numeric and automatically one-hot-encoded categorical covariates.
- Per-trait phenotype missingness and genotype-cell count filters.
- Signed `d/a` estimates and additive, partial-dominance, complete-dominance,
  overdominance, underdominance, and unstable classifications.
- Per-chromosome compressed Parquet output for large scans.
- No PLINK, GCTA, R, or project-specific compiled extension required for
  association testing.

## Installation

Clone or download this repository, enter its root directory, and run:

```bash
python -m pip install .
```

For development and tests:

```bash
python -m pip install -e ".[test,plot]"
python -m pytest -q
```

Python 3.9 or newer is required. Direct dependencies are NumPy, pandas,
SciPy, scikit-learn, and PyArrow.

## Quick start

Phenotype and covariate DataFrames must be indexed by the same IID values used
in the PLINK `.fam` file.

```python
import domgwas
import pandas as pd

data = pd.read_csv("phenotypes.tsv", sep="\t").set_index("iid")

results = domgwas.run_gwas(
    "genotypes/study",                  # study.bed/.bim/.fam
    data[["trait"]],
    traits=["trait"],
    covar=data[["sex", "batch", "PC1", "PC2", "PC3"]],
    model="add-dom",
    min_genotype_count=5,
    out="results/trait",
)
```

The complete runnable example at [`examples/quickstart.py`](examples/quickstart.py)
generates synthetic PLINK data, runs domgwas, and writes results without
requiring external data.

```bash
python examples/quickstart.py
```

For a scan too large to retain in memory, use `return_results=False` together
with `out`. Each chromosome is written as `{out}.chrN.parquet`.

## Statistical model

For chromosome `c`, domgwas estimates

```text
V_c = sigma2 [h2 K_add,-c + (1-h2) I],
```

where `K_add,-c` is the additive LOCO relationship matrix. If
`L_c = V_c^(-1/2)`, the joint marker model is fitted as

```text
L_c y = L_c C beta + a L_c A + d L_c D + error,
```

where `C` contains the intercept and covariates, `A` is dosage 0/1/2, and `D`
is the heterozygote indicator 0/1/0. OLS in this transformed space is
equivalent to GLS in the original space.

ADDO uses a different covariance model with whole-genome additive and
dominance relationship matrices. domgwas uses additive-only LOCO covariance.
The methods are therefore empirically comparable but not mathematically
equivalent.

## Principal outputs

Results include raw and standardized additive/dominance effect estimates,
standard errors, marginal and conditional p-values, `-log10(p)` scores,
AA/AB/BB counts, MAF, genotype-filter status, chromosome-specific heritability,
sample size, signed `d/a`, and inheritance labels. See
[`docs/EXPECTED_OUTPUTS.md`](docs/EXPECTED_OUTPUTS.md) for the output contract.

## Validation

The 0.2.0 implementation was checked against direct GLS to machine precision,
tested in LD-aware null and power simulations, compared with PLINK 1.9,
PLINK 2, GCTA, and ADDO, and benchmarked through one million variants. The
included validation report records methods, results, figures, and limitations:

- [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md)
- [`docs/publication_result_audit.json`](docs/publication_result_audit.json)
- [`docs/COVARIATE_SUPPORT.md`](docs/COVARIATE_SUPPORT.md)

Private or potentially restricted real genotype and phenotype files are not
included in this repository. The report retains aggregate results only.

## Scope and limitations

- Quantitative traits only; binary outcomes are rejected.
- Autosomal PLINK 1 variants; sex-aware X/Y/mitochondrial models are not
  implemented.
- Missing genotypes are mean-imputed within marker for association.
- Use `min_genotype_count` to avoid unstable sparse-cell dominance tests.
- Exact dense eigendecomposition is the validated path for up to 4,000
  samples. The truncated large-sample representation is experimental and its
  approximation error must be assessed before inference.
- Associations do not establish biological causality.

## Reproducibility and citation

Run the test suite with `python -m pytest -q`. The release history is in
[`CHANGELOG.md`](CHANGELOG.md), and citation metadata is in
[`CITATION.cff`](CITATION.cff). For a paper, archive the exact GitHub release
in Zenodo and add the resulting DOI to the citation metadata.

## License

Released under the [MIT License](LICENSE).
