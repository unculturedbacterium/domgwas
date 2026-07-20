# Association output contract

`domgwas.run_gwas` returns one row per trait and variant. If `out` is supplied,
it also writes compressed `{out}.chrN.parquet` files. The installed `domgwas`
command writes the same chromosome files, an optional combined
`{out}.all.parquet`, and `{out}.metadata.json`.

Principal columns include:

- Variant identity: `chrom`, `snp`, `pos`, alleles, and `trait`.
- Genotype summary: `n_AA`, `n_AB`, `n_BB`, `maf`, and
  `genotype_filter_pass`.
- Marginal tests: additive and dominance effects, standard errors, p-values,
  and `-log10(p)` scores from separate one-term models.
- Joint tests: `beta_add_joint_raw`, `beta_dom_joint_raw`, corresponding
  standard errors, p-values, and `neglog_p_add_joint` /
  `neglog_p_dom_joint` from the conditional two-term model.
- Nested-model test: `p_avsad` and `neglog_p_avsad` for additive versus
  additive-plus-dominance.
- Interpretation: signed `degree_of_dominance`, absolute ratio,
  `dominance_class`, and `inheritance_mode`.
- Model metadata: chromosome-specific `h2` and `n_samples`.

Markers that fail `min_genotype_count` remain in the output with
`genotype_filter_pass=False`; non-estimable test statistics are missing. Count
eligible markers explicitly before multiple-testing correction.

Column names can be inspected programmatically:

```python
results = domgwas.run_gwas(...)
print(results.columns.tolist())
```
