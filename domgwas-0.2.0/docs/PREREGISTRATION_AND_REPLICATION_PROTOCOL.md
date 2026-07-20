# domgwas confirmatory analysis and replication protocol

Status: draft for investigator approval. This document is not a registered
timestamp and does not substitute for independent data.

## Primary analysis

- Population: animals with genotype, normalized phenotype, sex, cohort, and
  PC1-PC5 available.
- Primary phenotype:
  `u01_olivier_george_cocaine:regressedlr_pr_max`.
- Software: domgwas 0.2.0 release archived with a permanent DOI.
- Trait type: quantitative. No binary-outcome claims.
- Association panel: autosomal variants in `r11.2.1.bed/.bim/.fam`.
- Relationship panel: the deterministic 35,723-marker QC panel used by the
  publication validation.
- Fixed effects: intercept, sex, cohort, and PC1-PC5.
- Covariance: additive LOCO GRM; profile REML separately by chromosome.
- Marker model: joint additive dosage and heterozygote deviation.
- Genotype eligibility: at least five observed AA, AB, and BB animals.
- Missing genotypes: within-marker mean imputation after eligibility counts.
- Primary tests: conditional additive and conditional dominance coefficients
  from the joint model.
- Secondary tests: marginal additive, marginal dominance, and the nested
  additive-versus-additive-plus-dominance F test.
- Multiple testing: `0.05 / M`, where `M` is the number of eligible autosomal
  variants, reported separately for the two primary tests. A stricter joint
  family threshold of `0.05 / (2M)` will be included as sensitivity analysis.

## Quality-control reporting

Report initial and final sample counts; phenotype and covariate exclusions;
variant counts by chromosome; missingness; MAF; AA/AB/BB counts; estimated h2;
genomic lambda; QQ plots; and all software/environment versions. Investigate
lambda outside 0.95-1.05 before biological interpretation, without changing
the primary model solely to improve lambda.

## Replication cohort

An independent cohort must not overlap the discovery animals or their genotype
records. Apply the same phenotype definition when available. Otherwise define
the closest construct before examining genotype associations and document the
difference.

For each discovery locus, define the replication window and lead variant before
opening replication association results. Replication requires:

1. the same effect direction for the prespecified additive or dominance term;
2. nominal one-sided `p < 0.05 / R`, where `R` is the number of independent
   discovery loci advanced; and
3. a consistent local signal after allele harmonization and cohort covariates.

Report exact failures as well as successes. Meta-analysis may follow the
cohort-specific replication test but must not replace it.

## Exploratory analyses

Alternative dominance thresholds, sex-stratified effects, cohort interactions,
regional conditioning, alternative PC counts, and other cocaine phenotypes are
exploratory unless separately preregistered. Label all post hoc model changes.

## Items requiring investigator sign-off

- Final phenotype transformation and exclusion rules.
- Whether the two primary tests share a `0.05` family or each receive `0.05`.
- Discovery-locus clumping/window rule.
- Available independent cohort and harmonized phenotype.
- Biological follow-up criteria and responsible investigators.
