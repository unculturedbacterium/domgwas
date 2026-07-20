"""Check that the source tree is suitable for a GitHub release."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import domgwas


ROOT = Path(__file__).resolve().parent
EXPECTED_VERSION = "0.2.0"
FORBIDDEN_DIRECTORIES = {"__pycache__", ".pytest_cache", "build", "dist", "domgwas.egg-info"}


def main() -> None:
    assert domgwas.__version__ == EXPECTED_VERSION, (
        f"package version is {domgwas.__version__}, expected {EXPECTED_VERSION}"
    )
    required = [
        "README.md",
        "LICENSE",
        "CHANGELOG.md",
        "CITATION.cff",
        "pyproject.toml",
        ".github/workflows/tests.yml",
        "examples/quickstart.py",
        "docs/VALIDATION_REPORT.md",
        "docs/publication_result_audit.json",
    ]
    missing_required = [path for path in required if not (ROOT / path).is_file()]
    assert not missing_required, f"missing required files: {missing_required}"

    forbidden = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_dir() and path.name in FORBIDDEN_DIRECTORIES
    ]
    assert not forbidden, f"generated directories must be removed: {forbidden}"

    oversized = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file() and path.stat().st_size > 50 * 1024 * 1024
    ]
    assert not oversized, f"files larger than 50 MiB: {oversized}"

    missing_links: list[str] = []
    markdown_files = list(ROOT.rglob("*.md"))
    for source in markdown_files:
        text = source.read_text(encoding="utf-8")
        for target in re.findall(r"!?(?:\[[^]]*\])\(([^)]+)\)", text):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            path_text = target.split("#", 1)[0]
            if path_text and not (source.parent / path_text).resolve().exists():
                missing_links.append(f"{source.relative_to(ROOT)} -> {target}")
    assert not missing_links, f"broken local Markdown links: {missing_links}"

    audit = json.loads((ROOT / "docs/publication_result_audit.json").read_text(encoding="utf-8"))
    assert audit["status"] == "PASS", "publication validation audit did not pass"
    print(
        f"PASS: domgwas {domgwas.__version__}; {len(markdown_files)} Markdown files; "
        f"{sum(1 for path in ROOT.rglob('*') if path.is_file())} repository files"
    )


if __name__ == "__main__":
    main()
