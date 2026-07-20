# Uploading domgwas 0.2.0 to GitHub

This directory is the repository root. Upload its contents so that `README.md`,
`pyproject.toml`, `domgwas/`, `tests/`, and `.github/` appear at the top level
of the GitHub repository. Do not upload the parent directory as a ZIP file if
you want GitHub to display and index the source.

## GitHub website

1. Create a new empty repository named `domgwas` under the
   `unculturedbacterium` account. Do not initialize it with a README, license,
   or `.gitignore`; those files are already present here.
2. Open the empty repository and select **uploading an existing file**.
3. Drag all contents of this directory into the upload page, including the
   `.github` directory.
4. Commit the upload to the `main` branch.
5. Open the Actions tab and confirm that the test workflow passes.
6. Create a release tagged `v0.2.0` and attach the tested wheel if desired.

## Command line

Run these commands from this directory after creating an empty GitHub
repository:

```bash
git init
git add .
git commit -m "Release domgwas 0.2.0"
git branch -M main
git remote add origin https://github.com/unculturedbacterium/domgwas.git
git push -u origin main
git tag -a v0.2.0 -m "domgwas 0.2.0"
git push origin v0.2.0
```

Before uploading, run:

```bash
python release_check.py
python -m pip install -e ".[test]"
python -m pytest -q
```

The real rat genotype and phenotype files are intentionally excluded. Confirm
data-owner permission and animal-study requirements before depositing those
inputs in any public repository.
