"""Additive genomic relatedness matrices with leave-one-chromosome-out (LOCO).

Two-pass streaming design to bound memory:

1. One pass over every variant accumulates the genome-wide cross-products
   ``Z Z^T`` and the pairwise observed-SNP counts.
2. For each chromosome to be tested, a second short pass over only that
   chromosome's variants gives its cross-products; the LOCO GRM is
   ``(ZZ_all - ZZ_chr) / (w_all - w_chr)``, which is eigendecomposed and the
   full matrix discarded.

Peak memory is a small number of n*n matrices, independent of variant count.
"""
import numpy as np
from scipy.sparse.linalg import eigsh
from sklearn.utils.extmath import randomized_svd

from ._utils import standardize_block


def _eig(G, n_components=None):
    """Symmetric eigendecomposition, truncated to the top ``n_components``."""
    n = G.shape[0]
    if n_components is None or n_components >= n:
        s, U = np.linalg.eigh(G)
        order = np.argsort(s)[::-1]
        return U[:, order], np.abs(s[order])
    U, s, _ = randomized_svd(G, n_components=n_components, random_state=0)
    return U, np.abs(s)


def compute_loco_grms(reader, sample_index=None, block_size=8192,
                      n_components=None, chroms=None, verbose=True):
    """Return ``{chrom: {'U':U, 's':s}}`` of LOCO GRM eigendecompositions.

    ``sample_index`` selects rows (analysis samples) from the fileset.
    ``n_components`` truncates each eigendecomposition (defaults to full when
    n <= 4000, else 2000).
    """
    n_all = reader.n_samples
    sidx = np.arange(n_all) if sample_index is None else np.asarray(sample_index)
    n = len(sidx)
    if n_components is None:
        n_components = n if n <= 4000 else 2000
    chrom_arr = reader.bim["chrom"].values
    run_chroms = ([str(c) for c in chroms] if chroms is not None
                  else list(dict.fromkeys(chrom_arr)))

    # Pass 1: genome-wide cross-products.
    zzt_all = np.zeros((n, n))
    w_all = np.zeros((n, n))
    for _, g in reader.iter_blocks(block_size, sample_index=sidx):
        z, _, mask = standardize_block(g)
        zzt_all += z @ z.T
        w_all += mask @ mask.T

    out = {}
    for c in run_chroms:
        idx_c = np.where(chrom_arr == c)[0]
        zzt_c = np.zeros((n, n))
        w_c = np.zeros((n, n))
        for _, g in reader.iter_blocks(block_size, variant_index=idx_c, sample_index=sidx):
            z, _, mask = standardize_block(g)
            zzt_c += z @ z.T
            w_c += mask @ mask.T
        loco = (zzt_all - zzt_c) / np.maximum(w_all - w_c, 1.0)
        U, s = _eig(loco, n_components=n_components)
        out[c] = {"U": U, "s": s}
        if verbose:
            print(f"  LOCO GRM chr {c}: {len(idx_c)} variants left out, "
                  f"rank {len(s)}")
    return out


def compute_grm(reader, sample_index=None, block_size=8192):
    """Genome-wide additive GRM (dense n*n), NaN-aware."""
    n_all = reader.n_samples
    sidx = np.arange(n_all) if sample_index is None else np.asarray(sample_index)
    n = len(sidx)
    zzt = np.zeros((n, n))
    w = np.zeros((n, n))
    for _, g in reader.iter_blocks(block_size, sample_index=sidx):
        z, _, mask = standardize_block(g)
        zzt += z @ z.T
        w += mask @ mask.T
    return zzt / np.maximum(w, 1.0)
