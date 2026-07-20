"""End-to-end validation of the domgwas core on simulated data."""
import time
import numpy as np
import pandas as pd

import domgwas
from domgwas.testutils import simulate
from domgwas import PlinkReader, compute_loco_grms, estimate_h2, whiten_y

np.set_printoptions(suppress=True)
PREFIX = "/tmp/domsim"

# ---------------------------------------------------------------- simulate
meta = simulate(PREFIX, n_fam=200, sibs=4, n_chrom=3, snps_per_chrom=1000,
                Va=0.25, Vd=0.35, Ve=0.40, seed=1)
pheno = meta["pheno"]
print(f"simulated n={meta['n']} m={meta['m']}  big over-dominant SNP = {meta['big_dom_name']}")

# ---------------------------------------------------------------- 1. I/O round-trip
reader = PlinkReader(PREFIX)
blk = reader.read_block(0, meta["m"])          # (n, m)
print(f"\n[1] I/O: reader shape {blk.shape}, "
      f"genotype values present = {sorted(pd.unique(blk[~np.isnan(blk)]))[:5]}")

# ---------------------------------------------------------------- run GWAS
t0 = time.time()
res = domgwas.run_gwas(PREFIX, pheno, model="add-dom", verbose=True)
elapsed = time.time() - t0
print(f"\nGWAS finished in {elapsed:.2f}s  -> {len(res)} (variant x trait) rows")

# ---------------------------------------------------------------- 2. correctness vs manual OLS
# Reproduce the whitened phenotype for chr 1 / y_mix using the library, then
# fit the joint additive+dominance model by explicit OLS on a few SNPs.
grms = compute_loco_grms(reader, chroms=["1"], verbose=False)
U, s = grms["1"]["U"], grms["1"]["s"]
y = pheno["y_mix"].values
hh = estimate_h2(y, U, s)
ystar = whiten_y(y, U, s, hh["h2"], hh["yvar"]); ystar -= ystar.mean()

idx1 = reader.chrom_variant_index("1")
G1 = reader.read_block(int(idx1.min()), len(idx1))
rng = np.random.default_rng(0)
check = rng.choice(len(idx1), 8, replace=False)

max_beta_diff = 0.0
max_p_diff = 0.0
from scipy.stats import t as tdist
sub = res[(res["chrom"] == "1") & (res["trait"] == "y_mix")]
for j in check:
    g = G1[:, j].astype(float)          # raw 0/1/2 (no missing in the sim)
    add = g
    het = (g == 1).astype(float)         # raw 0/1/0
    X = np.column_stack([np.ones_like(add), add, het])
    b, *_ = np.linalg.lstsq(X, ystar, rcond=None)
    resid = ystar - X @ b
    dof = len(ystar) - 3
    sigma2 = (resid @ resid) / dof
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    p_dom = 2 * tdist.sf(abs(b[2] / se[2]), df=dof)
    row = sub[sub["snp"] == f"rs{idx1[j]}"].iloc[0]
    max_beta_diff = max(max_beta_diff, abs(b[2] - row["beta_dom_joint_raw"]))
    max_p_diff = max(max_p_diff, abs((-np.log10(p_dom)) - row["neglog_p_dom_joint"]))
print(f"\n[2] correctness vs manual OLS (joint dominance term, 8 SNPs):")
print(f"    max |beta| diff = {max_beta_diff:.2e}   max |-log10 p| diff = {max_p_diff:.2e}")

# ---------------------------------------------------------------- 3. null calibration
MED = 0.45493642
null_snps = np.setdiff1d(np.arange(meta["m"]),
                         np.r_[meta["causal_add"], meta["causal_dom"], meta["big_dom_snp"]])
null_names = {f"rs{i}" for i in null_snps}
nn = res[res["snp"].isin(null_names)]
print("\n[3] null calibration (lambda_GC, should be ~1):")
for trait in ["y_add", "y_dom", "y_mix", "y_null"]:
    d = nn[nn["trait"] == trait]
    for col, lab in [("neglog_p_additive", "add"),
                     ("neglog_p_dom_joint", "dom"),
                     ("neglog_p_avsad", "avsad")]:
        chi = (10 ** (-d[col].values)) 
        from scipy.stats import chi2 as _c
        stat = _c.isf(chi, 1)
        lam = np.nanmedian(stat) / MED
        print(f"    {trait:7s} {lab:5s} lambda = {lam:.3f}", end="   ")
    print()

# ---------------------------------------------------------------- 4. power + classification
print("\n[4] power & classification for the injected over-dominant SNP "
      f"({meta['big_dom_name']}) on trait y_dom:")
dom_res = res[res["trait"] == "y_dom"].copy()
dom_res = dom_res.sort_values("neglog_p_dom_joint", ascending=False).reset_index(drop=True)
rank = dom_res.index[dom_res["snp"] == meta["big_dom_name"]][0] + 1
big = dom_res[dom_res["snp"] == meta["big_dom_name"]].iloc[0]
print(f"    rank by joint dominance p-value : {rank} of {len(dom_res)}")
print(f"    -log10 p (dom joint)            : {big['neglog_p_dom_joint']:.1f}")
print(f"    -log10 p (add vs add+dom F)     : {big['neglog_p_avsad']:.1f}")
print(f"    degree of dominance |d/a|       : {big['degree_of_dominance']:.2f}")
print(f"    inheritance class               : {big['dominance_class']}")
print("\n    class counts among genome-wide-significant dominance hits (y_dom):")
sig = dom_res[dom_res["neglog_p_dom_joint"] > 5]
print("    ", sig["dominance_class"].value_counts().to_dict())
