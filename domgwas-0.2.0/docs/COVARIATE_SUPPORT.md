# Covariate support in domgwas 0.2

`domgwas.run_gwas(..., covar=...)` accepts a pandas `DataFrame` indexed by IID
or an ordered numeric array. A DataFrame is recommended because its rows are
aligned to genotype and phenotype IDs.

Numeric columns are used directly. Object, string, category, and Boolean
columns are one-hot encoded with one reference level. The design always
contains an intercept, and redundant columns are removed before model fitting.

Version 0.2 applies the same `V^(-1/2)` transformation to phenotype,
intercept, covariates, additive dosage, and dominance code. Direct-reference
validation summarized in `VALIDATION_REPORT.md` verifies
coefficients and p-values to machine precision. This supersedes the
phenotype-only adjustment in version 0.1.0.

Recommended fixed effects include sex, batch, cohort/site, prespecified genetic
PCs, age, and measured technical variables. Small categorical levels should be
merged or excluded before analysis. Missing covariates exclude a sample from
that trait; phenotype missingness is handled separately for each trait.

```bash
domgwas \
  --bfile data/cohort \
  --pheno data/phenotypes.tsv \
  --id-column iid \
  --traits trait1,trait2 \
  --covariates PC1,PC2,PC3,age \
  --categorical-covariates sex,batch,cohort \
  --min-genotype-count 5 \
  --out results/traits
```

The implementation is restricted to quantitative traits. Binary covariates
are supported; binary phenotypes are not.
