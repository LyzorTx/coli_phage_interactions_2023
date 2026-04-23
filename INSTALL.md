# Getting Started

This repository uses Conda-compatible environment manifests (via
[Miniconda](https://docs.anaconda.com/miniconda/) or `micromamba`) to manage environments.

Install Miniconda using the official instructions for your platform. If `conda activate` is not available yet,
run `conda init <your-shell>` once and restart your shell. Then create and activate the main development environment:

```bash
conda env create -f environment.yml
conda activate phage_env
```

If you prefer `micromamba`, the same manifest works there:

```bash
micromamba create -n phage_env -f environment.yml
micromamba activate phage_env
```

To enable running Jupyter notebooks using this environment, create a dedicated Jupyter kernel for it:

```bash
python -m ipykernel install --user --name phage_env --display-name phage_env
```

Now when you run notebooks in `jupyter lab`, make sure to select the `phage_env` kernel.

## Optional specialized environments

Most local development should stay in `phage_env`, but a few heavyweight or conflict-prone bioinformatics toolchains
now live in dedicated manifests at the repository root so developers can reproduce the same boundaries used by the CI
images:

- `environment.phylogroup-caller.yml`
- `environment.serotype-caller.yml`
- `environment.sequence-type-caller.yml`
- `environment.phage-annotation-tools.yml`

Create them only when you need those tools locally:

```bash
conda env create -f environment.phylogroup-caller.yml
conda env create -f environment.serotype-caller.yml
conda env create -f environment.sequence-type-caller.yml
conda env create -f environment.phage-annotation-tools.yml
```

These manifests are intentionally checked in next to `environment.yml` rather than buried under CI-only paths so local
developers can run the same split toolchains as GitHub Actions.

## Optional: automatic env activation with direnv

This repository includes a tracked `.envrc` that runs `conda activate phage_env` when you enter the repository
directory.

1. Install `direnv` and enable its shell integration.
2. For `zsh`, add the hook to `~/.zshrc` if needed:
   ```bash
   eval "$(direnv hook zsh)"
   ```
3. Reload your shell and allow the repo env file:
   ```bash
   direnv allow
   ```

## Git hooks

This repository uses `pre-commit` for automated checks. Install both hook types once per clone:

```bash
pre-commit install
pre-commit install --hook-type pre-push
```

This activates:

- **pre-commit stage:** ruff linting/formatting, gitignore enforcement, pymarkdown fixes.
- **pre-push stage:** verifies your branch is rebased on `origin/main` before allowing `git push`.

## Markdown linting

Commands for manual runs:

```bash
pymarkdown --config .pymarkdown.yaml fix -r .
pre-commit run pymarkdown --all-files
```

After auto-fixes, stage updated files before committing.
