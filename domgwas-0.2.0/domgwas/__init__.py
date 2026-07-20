"""domgwas: a pure-Python, pip-installable, dominance-aware GWAS pipeline.

Public API
----------
    run_gwas         one-call dominance-aware GWAS over a PLINK1 fileset
    PlinkReader      streaming PLINK1 reader
    compute_loco_grms, compute_grm
    estimate_h2, whiten_y, blup_resid
    scan_chromosome
    degree_of_dominance, classify_da
"""
from .io import PlinkReader, read_bim, read_fam
from .grm import compute_loco_grms, compute_grm
from .vc import estimate_h2, whiten_y, whiten_matrix, covariance_matrix, blup_resid
from .assoc import scan_chromosome
from .classify import degree_of_dominance, classify_da, classify_inheritance, DEFAULT_THRESHOLDS
from .pipeline import run_gwas

__version__ = "0.2.0"
__all__ = [
    "run_gwas", "PlinkReader", "read_bim", "read_fam",
    "compute_loco_grms", "compute_grm",
    "estimate_h2", "whiten_y", "whiten_matrix", "covariance_matrix", "blup_resid",
    "scan_chromosome", "degree_of_dominance", "classify_da", "classify_inheritance",
    "DEFAULT_THRESHOLDS", "__version__",
]
