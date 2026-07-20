# Changelog

All notable changes to domgwas are documented here.

## 0.2.0 - 2026-07-20

- Fit complete-design GLS by applying the same covariance transform to the
  phenotype, intercept, covariates, additive dosage, and dominance code.
- Include the complete fixed-effect design in profile REML estimation.
- Add numeric and categorical covariates with redundant-column removal.
- Handle missing phenotypes separately by trait.
- Add AA/AB/BB counts and configurable minimum genotype-cell filtering.
- Add signed `d/a`, allele-direction-aware inheritance modes, and unstable
  classification when the additive effect is near zero.
- Add a separate relationship-panel input with `grm_bfile`.
- Add selected-sample PLINK decoding and streamed chromosome output.
- Reject binary phenotypes because logistic mixed models are not implemented.
- Add direct-GLS, covariate, missingness, decoding, filter, classification, and
  binary-outcome tests.

## 0.1.0

- Initial dominance-aware GWAS implementation.
- This release transformed the phenotype without transforming the complete
  marker design and should not be used for mixed-model inference.
