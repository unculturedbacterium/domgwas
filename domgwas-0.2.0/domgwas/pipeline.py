"""One-call dominance-aware GWAS driver.

    import domgwas
    res = domgwas.run_gwas("mydata", pheno_df, model="add-dom")

Steps: align samples -> build LOCO additive GRMs (streamed) -> per chromosome,
estimate h2 and whiten each trait -> stream the chromosome's variants through
the association engine -> concatenate (and optionally write parquet per chrom).
"""
import warnings

import numpy as np
import pandas as pd

from .io import PlinkReader
from .grm import compute_loco_grms
from .vc import estimate_h2, covariance_matrix
from .assoc import scan_chromosome_gls


def _load_pheno(pheno):
    if isinstance(pheno, pd.DataFrame):
        return pheno
    df = pd.read_csv(pheno, sep=r"\s+")
    cols = {c.lower(): c for c in df.columns}
    id_col = cols.get("iid") or cols.get("id") or df.columns[1]
    return df.set_index(id_col)


def _prepare_covariates(covar, ids):
    if covar is None:
        return pd.DataFrame(index=ids)
    if isinstance(covar, pd.DataFrame):
        covar = covar.copy()
        covar.index = covar.index.astype(str)
        frame = covar.loc[ids].copy()
        categorical = frame.select_dtypes(include=["object", "category", "string", "bool"]).columns
        numeric = frame.columns.difference(categorical, sort=False)
        parts = []
        if len(numeric):
            parts.append(frame[numeric].apply(pd.to_numeric, errors="coerce"))
        if len(categorical):
            parts.append(pd.get_dummies(frame[categorical], drop_first=True, dtype=float))
        return pd.concat(parts, axis=1) if parts else pd.DataFrame(index=ids)
    values = np.asarray(covar, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    if len(values) != len(ids):
        raise ValueError("covariate array must have one row per aligned phenotype sample")
    return pd.DataFrame(values, index=ids)


def _full_rank_design(covariates):
    columns = [np.ones(len(covariates), dtype=float)]
    rank = 1
    for column in np.asarray(covariates, dtype=float).T:
        candidate = np.column_stack(columns + [column])
        new_rank = np.linalg.matrix_rank(candidate)
        if new_rank > rank:
            columns.append(column)
            rank = new_rank
    return np.column_stack(columns)


def run_gwas(bfile, pheno, traits=None, covar=None, model="add-dom",
             chroms=None, block_size=8192, n_components=None,
             y_correction="ystar", out=None, da_thresholds=(0.25, 0.75, 1.25),
             min_genotype_count=0, stability_z=1.0,
             per_trait_missing=True, grm_bfile=None, return_results=True,
             verbose=True):
    """Run a dominance-aware (or additive-only) GWAS over a PLINK1 fileset.

    Parameters
    ----------
    bfile : str
        PLINK prefix (``.bed/.bim/.fam``).
    pheno : DataFrame or path
        Phenotypes indexed by IID (DataFrame) or a whitespace table with an
        IID column. Missing values are handled separately for each trait.
    model : {'add-dom', 'additive'}
    y_correction : {'ystar', None}
        ``ystar`` performs exact complete-design GLS. ``None`` performs OLS.
    out : str or None
        If given, write ``{out}.chr{c}.parquet`` per chromosome.
    """
    reader = PlinkReader(bfile)
    grm_reader = reader if grm_bfile is None else PlinkReader(grm_bfile)
    pheno = _load_pheno(pheno)
    if traits is None:
        traits = [c for c in pheno.columns if c != "FID"]

    if model not in {"add-dom", "additive"}:
        raise ValueError("model must be 'add-dom' or 'additive'")
    if y_correction not in {"ystar", None}:
        raise ValueError("only exact GLS ('ystar') and uncorrected OLS (None) are supported")

    pheno.index = pheno.index.astype(str)
    fam_iid = reader.fam["iid"].astype(str).values
    fam_pos = pd.Index(fam_iid)
    pset = set(pheno.index)
    common = [i for i in fam_iid if i in pset]
    sub = pheno.loc[common, traits].apply(pd.to_numeric, errors="coerce")
    binary = [trait for trait in traits if sub[trait].dropna().nunique() == 2]
    if binary:
        raise ValueError(
            f"binary phenotypes are not supported by the quantitative-trait model: {binary}"
        )
    covar_frame = _prepare_covariates(covar, sub.index)
    covar_complete = np.isfinite(covar_frame.to_numpy(dtype=float)).all(axis=1)
    trait_masks = {
        trait: (sub[trait].notna().to_numpy() & covar_complete)
        for trait in traits
    }
    if not per_trait_missing:
        joint = np.logical_and.reduce(list(trait_masks.values()))
        trait_masks = {trait: joint for trait in traits}

    groups = {}
    for trait, mask in trait_masks.items():
        groups.setdefault(mask.tobytes(), {"mask": mask, "traits": []})["traits"].append(trait)
    if not return_results and out is None:
        raise ValueError("out is required when return_results=False")
    if not return_results and len(groups) != 1:
        raise ValueError("streamed output currently requires traits with one shared sample mask")

    chrom_results = {}
    for group in groups.values():
        mask = group["mask"]
        group_traits = group["traits"]
        analysis_ids = sub.index[mask].tolist()
        if len(analysis_ids) < 5:
            warnings.warn(f"skipping {group_traits}: fewer than five complete samples")
            continue
        sample_index = fam_pos.get_indexer(analysis_ids)
        grm_fam_pos = pd.Index(grm_reader.fam["iid"].astype(str).values)
        grm_sample_index = grm_fam_pos.get_indexer(analysis_ids)
        if np.any(grm_sample_index < 0):
            missing = [analysis_ids[i] for i in np.where(grm_sample_index < 0)[0][:5]]
            raise ValueError(f"analysis samples missing from grm_bfile: {missing}")
        Y = sub.loc[analysis_ids, group_traits].to_numpy(dtype=np.float64)
        C = _full_rank_design(covar_frame.loc[analysis_ids].to_numpy(dtype=float))
        if verbose:
            print(f"analysis set: {len(analysis_ids)} samples x {len(group_traits)} traits; "
                  f"{reader.n_variants} variants; model={model}; fixed-effect rank={C.shape[1]}")
        grms = compute_loco_grms(
            grm_reader, sample_index=grm_sample_index, block_size=block_size,
            n_components=n_components, chroms=chroms, verbose=verbose
        )
        run_chroms = [str(c) for c in chroms] if chroms is not None else list(grms.keys())
        for c in run_chroms:
            U, s = grms[c]["U"], grms[c]["s"]
            whiteners = []
            h2s = []
            for j in range(Y.shape[1]):
                if y_correction is None:
                    whiteners.append(np.eye(len(Y), dtype=np.float64))
                    h2s.append(np.nan)
                else:
                    hh = estimate_h2(Y[:, j], U, s, covar=C)
                    h2s.append(hh["h2"])
                    whiteners.append(covariance_matrix(U, s, hh["h2"], hh["yvar"], power=-0.5))
            if verbose and y_correction is not None:
                print(f"chr {c}: h2 = {np.round(h2s, 3).tolist()}")
            df = scan_chromosome_gls(
                reader, c, Y, group_traits, sample_index, whiteners,
                [C] * len(group_traits), block_size=block_size, model=model,
                da_thresholds=da_thresholds,
                min_genotype_count=min_genotype_count,
                stability_z=stability_z,
            )
            df["h2"] = df["trait"].map(dict(zip(group_traits, h2s)))
            df["n_samples"] = len(analysis_ids)
            if not return_results:
                df.to_parquet(f"{out}.chr{c}.parquet", index=False, compression="zstd")
                continue
            chrom_results.setdefault(c, []).append(df)

    if not return_results:
        return pd.DataFrame()
    if not chrom_results:
        raise ValueError("no traits had enough complete samples")
    results = []
    for c, frames in chrom_results.items():
        df = pd.concat(frames, ignore_index=True)
        if out is not None:
            df.to_parquet(f"{out}.chr{c}.parquet", index=False, compression="zstd")
        results.append(df)
    return pd.concat(results, ignore_index=True)
