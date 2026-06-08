# Environment Notes

This project keeps the runnable source separate from local data and generated
outputs. A new machine only needs the tracked files plus the input data placed
in the expected local folders.

## Conda Setup

```powershell
conda env create -f environment.yml
conda activate gplw-rc
python code\run_all_make_two_figures.py --help
```

## Pip Setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r code\requirements.txt
python code\run_all_make_two_figures.py --help
```

The pip path may require extra system packages for `geopandas` and `cartopy`.
Use the conda environment when reproducing maps on a clean workstation.

## Required Local Inputs

- Station summary table: pass with `--main-csv`, or place it in `Data/` using
  the default filename expected by `run_all_make_two_figures.py`.
- Climate-zone shapefile: pass with `--climate-shp`. Use `--no-climate` to skip
  this optional panel during quick checks.

## Expected Outputs

- Intermediate CSVs are written to `code/RC_figure_csvs_1minus_beta/`.
- Final figures are written to `code/RC_final_figures_1minus_beta/`.

Both output folders are ignored by Git.
