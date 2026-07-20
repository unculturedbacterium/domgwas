"""Variance-component estimation and covariance transformations.

The default association path uses the same ``V^(-1/2)`` transformation for the
phenotype, fixed-effect design, and tested genotypes.  This is generalized
least squares (GLS), unlike the phenotype-only transformation used by domgwas
0.1.0.  Truncated eigendecompositions use the residual covariance eigenvalue in
the omitted orthogonal complement.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def _apply_eigen_operator(x, U, s, h2, scale=1.0, power=-0.5, tol=1e-10):
    """Apply ``[scale * (h2 G + (1-h2) I)]**power`` to vectors/matrices."""
    x = np.asarray(x, dtype=np.float64)
    was_vector = x.ndim == 1
    if was_vector:
        x = x[:, None]
    U = np.asarray(U, dtype=np.float64)
    s = np.maximum(np.asarray(s, dtype=np.float64), 0.0)
    d_top = np.maximum(scale * (h2 * s + (1.0 - h2)), tol)
    d_res = max(scale * (1.0 - h2), tol)
    utx = U.T @ x
    residual = x - U @ utx
    result = U @ (utx * (d_top ** power)[:, None]) + residual * (d_res ** power)
    return result[:, 0] if was_vector else result


def covariance_matrix(U, s, h2, scale=1.0, power=-0.5, tol=1e-10):
    """Materialize a covariance operator, primarily for blockwise GLS scans."""
    identity = np.eye(np.asarray(U).shape[0], dtype=np.float64)
    return _apply_eigen_operator(identity, U, s, h2, scale, power, tol)


def estimate_h2(y, U, s, method="REML", tol=1e-10, covar=None):
    """Profile ML/REML estimate of additive heritability.

    ``covar`` is the complete fixed-effect design.  An intercept is used when
    it is omitted.  The returned ``yvar`` is the profiled covariance scale, so
    ``V = yvar * (h2 G + (1-h2) I)``.
    """
    y = np.asarray(y, dtype=np.float64).ravel()
    n = y.size
    C = np.ones((n, 1), dtype=np.float64) if covar is None else np.asarray(covar, dtype=np.float64)
    if C.ndim == 1:
        C = C[:, None]
    if len(C) != n or not np.isfinite(y).all() or not np.isfinite(C).all():
        raise ValueError("y and covar must be finite and have matching rows")
    rank_c = int(np.linalg.matrix_rank(C))
    if rank_c == 0 or rank_c >= n:
        raise ValueError("fixed-effect design must have rank between 1 and n-1")
    U = np.asarray(U, dtype=np.float64)
    s = np.maximum(np.asarray(s, dtype=np.float64), 0.0)
    k_resid = max(n - len(s), 0)

    def objective(h2):
        vinv_y = _apply_eigen_operator(y, U, s, h2, power=-1.0, tol=tol)
        vinv_c = _apply_eigen_operator(C, U, s, h2, power=-1.0, tol=tol)
        info = C.T @ vinv_c
        sign, logdet_info = np.linalg.slogdet(info)
        if sign <= 0:
            return np.inf
        beta = np.linalg.solve(info, C.T @ vinv_y)
        resid = y - C @ beta
        q = max(float(resid @ _apply_eigen_operator(resid, U, s, h2, power=-1.0, tol=tol)), tol)
        df = n - rank_c if method.upper() == "REML" else n
        sigma2 = max(q / df, tol)
        d_top = np.maximum(h2 * s + (1.0 - h2), tol)
        d_res = max(1.0 - h2, tol)
        logdet = float(np.log(d_top).sum() + k_resid * np.log(d_res))
        value = df * np.log(sigma2) + logdet
        if method.upper() == "REML":
            value += logdet_info
        return 0.5 * value

    res = minimize_scalar(objective, bounds=(1e-4, 1 - 1e-4), method="bounded")
    h2 = float(res.x)
    vinv_y = _apply_eigen_operator(y, U, s, h2, power=-1.0, tol=tol)
    vinv_c = _apply_eigen_operator(C, U, s, h2, power=-1.0, tol=tol)
    beta = np.linalg.solve(C.T @ vinv_c, C.T @ vinv_y)
    resid = y - C @ beta
    df = n - rank_c if method.upper() == "REML" else n
    sigma2 = max(float(resid @ _apply_eigen_operator(resid, U, s, h2, power=-1.0, tol=tol)) / df, tol)
    return {"h2": h2, "yvar": sigma2, "converged": bool(res.success), "objective": float(res.fun)}


def whiten_y(y, U, s, h2, yvar, tol=1e-10):
    """Return V^(-1/2) y for V = yvar * (h2 * G + (1 - h2) I).

    Correct even when ``U, s`` are truncated: the retained eigenspace is scaled
    by its eigenvalues and the orthogonal complement by yvar * (1 - h2).
    """
    return _apply_eigen_operator(y, U, s, h2, yvar, power=-0.5, tol=tol)


def whiten_matrix(x, U, s, h2, yvar, tol=1e-10):
    """Return ``V^(-1/2) x`` without centering any columns."""
    return _apply_eigen_operator(x, U, s, h2, yvar, power=-0.5, tol=tol)


def blup_resid(y, U, s, h2, yvar, tol=1e-10):
    """GRAMMAR residual: y minus the polygenic BLUP (alternative correction)."""
    y = np.asarray(y, dtype=np.float64).ravel()
    ybar = y.mean()
    yc = y - ybar
    s = np.maximum(np.asarray(s, dtype=np.float64), 0.0)
    w = (h2 * s) / (h2 * s + (1.0 - h2))
    Uty = U.T @ yc
    blup = U @ (w * Uty)
    return yc - blup
