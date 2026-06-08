# GPLW RC Figure Workspace

This repository is organized as a reproducible workspace for the GPLW runoff
generation and confluence figure pipeline.

The tracked project core is the Python code under `code/`, the environment
definition, and documentation. Large raw data, generated figures, notebooks,
GIS layers, archives, and temporary files are intentionally kept out of Git by
the root `.gitignore`.

## Project Layout

- `code/`: RC figure export and plotting scripts.
- `code/run_all_make_two_figures.py`: one-command entry point for the main
  two-figure pipeline.
- `code/01_export_csvs_for_two_figures.py`: converts station data into
  figure-ready CSV tables.
- `code/02_plot_two_figures.py`: renders the final figures from exported CSVs.
- `Data/`: local input data folder, ignored by Git.
- `code/RC_figure_csvs_1minus_beta/`: generated intermediate CSVs, ignored.
- `code/RC_final_figures_1minus_beta/`: generated PNG/PDF/SVG figures, ignored.
- Other top-level folders are legacy analysis outputs or local data products;
  keep them local unless they are explicitly curated into the source workflow.

## Recreate The Environment

Conda is recommended because optional map and GIS dependencies are easier to
install from conda-forge.

```powershell
conda env create -f environment.yml
conda activate gplw-rc
```

If you prefer pip:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r code\requirements.txt
```

`geopandas` is optional for climate-zone matching and `cartopy` is optional for
map rendering. Install both when reproducing the full publication figures.

## Run The Pipeline

From the project root:

```powershell
cd code
python run_all_make_two_figures.py
```

Useful variants:

```powershell
python run_all_make_two_figures.py --no-climate
python run_all_make_two_figures.py --skip-export
python run_all_make_two_figures.py --main-csv ..\Data\filtered_station_summary_.xlsx
python run_all_make_two_figures.py --csv-dir .\RC_figure_csvs_1minus_beta --fig-dir .\RC_final_figures_1minus_beta
```

The default station table is expected under `Data/`. To run climate-zone
matching, provide the local Koppen climate shapefile with `--climate-shp` or
place it at the default path used by `code/run_all_make_two_figures.py`.

## Git Policy

Track source code, environment files, and documentation. Do not commit raw
data, generated figures, GIS files, notebooks, large archives, or temporary
sync files. If a result must be published, export it separately or add a small
curated artifact with a clear note in the commit message.
