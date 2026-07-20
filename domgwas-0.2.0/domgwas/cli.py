"""Command-line interface for quantitative-trait domgwas scans."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

from . import __version__, run_gwas


def _names(value: str | None) -> list[str]:
    return [] if not value else [item.strip() for item in value.split(",") if item.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run a quantitative-trait dominance-aware GWAS")
    parser.add_argument("--version", action="version", version=f"domgwas {__version__}")
    parser.add_argument("--bfile", required=True, help="PLINK 1 prefix without .bed/.bim/.fam")
    parser.add_argument("--grm-bfile", help="Optional separate PLINK 1 relationship-panel prefix")
    parser.add_argument("--pheno", required=True, type=Path, help="Tab-separated phenotype/covariate table")
    parser.add_argument("--id-column", default="iid", help="Sample ID column in --pheno")
    parser.add_argument("--traits", required=True, help="Comma-separated quantitative phenotype columns")
    parser.add_argument("--covariates", help="Comma-separated numeric covariate columns")
    parser.add_argument("--categorical-covariates", help="Comma-separated categorical covariate columns")
    parser.add_argument("--chromosomes", help="Comma-separated chromosome labels")
    parser.add_argument("--model", choices=["add-dom", "additive"], default="add-dom")
    parser.add_argument("--block-size", type=int, default=8192)
    parser.add_argument("--min-genotype-count", type=int, default=5)
    parser.add_argument("--complete-case-across-traits", action="store_true")
    parser.add_argument("--stream", action="store_true", help="Write chromosomes without retaining all results")
    parser.add_argument("--out", required=True, type=Path, help="Output prefix")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    frame = pd.read_csv(args.pheno, sep="\t", dtype={args.id_column: "string"})
    frame[args.id_column] = frame[args.id_column].str.strip()
    frame = frame.set_index(args.id_column)
    traits = _names(args.traits)
    numeric = _names(args.covariates)
    categorical = _names(args.categorical_covariates)
    requested = traits + numeric + categorical
    missing = sorted(set(requested) - set(frame.columns))
    if missing:
        raise ValueError(f"columns missing from phenotype table: {missing}")

    covar = None
    if numeric or categorical:
        if numeric:
            frame[numeric] = frame[numeric].apply(pd.to_numeric, errors="raise")
        covar = frame[numeric + categorical]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    result = run_gwas(
        args.bfile,
        frame[traits],
        traits=traits,
        covar=covar,
        model=args.model,
        chroms=_names(args.chromosomes) or None,
        block_size=args.block_size,
        min_genotype_count=args.min_genotype_count,
        per_trait_missing=not args.complete_case_across_traits,
        grm_bfile=args.grm_bfile,
        out=str(args.out),
        return_results=not args.stream,
        verbose=True,
    )
    elapsed = time.perf_counter() - started
    if not args.stream:
        result.to_parquet(Path(f"{args.out}.all.parquet"), index=False)
    metadata = {
        "domgwas_version": __version__,
        "bfile": str(Path(args.bfile).resolve()),
        "grm_bfile": str(Path(args.grm_bfile).resolve()) if args.grm_bfile else None,
        "phenotype": str(args.pheno.resolve()),
        "traits": traits,
        "numeric_covariates": numeric,
        "categorical_covariates": categorical,
        "model": args.model,
        "block_size": args.block_size,
        "min_genotype_count": args.min_genotype_count,
        "per_trait_missingness": not args.complete_case_across_traits,
        "streamed_output": args.stream,
        "trait_type": "quantitative_only",
        "runtime_s": elapsed,
        "n_output_rows_in_memory": len(result),
    }
    Path(f"{args.out}.metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
