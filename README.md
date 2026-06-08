# global-flood-scaling-organization

Code and data-processing workflows for analysing scale-dependent runoff
generation and flood organization using global streamflow observations and
ISIMIP experiments.

This repository is organized as a reproducible workspace for the GPLW runoff
generation and confluence figure pipeline. The tracked project core is the
Python code under `code/`, the environment definition, and documentation. Large
raw data, generated figures, notebooks, GIS layers, archives, and temporary
files are intentionally kept out of Git by the root `.gitignore`.

## Project Layout

- `code/`: RC figure export and plotting scripts.
- `code/run_all_make_two_figures.py`: one-command entry point for the main
  two-figure pipeline.
- `code/01_export_csvs_for_two_figures.py`: converts station data into
  figure-ready CSV tables.
- `code/02_plot_two_figures.py`: renders the final figures from exported CSVs.
- `code/compute_event_rc_from_flood_events.py`: computes event-level `R` and
  `C` from `flood_events__*.nc` files.
- `docs/`: environment notes and reproduction cases.
- `Data/`: local input data folder, ignored by Git.
- `code/RC_figure_csvs_1minus_beta/`: generated intermediate CSVs, ignored by
  Git.
- `code/RC_final_figures_1minus_beta/`: generated PNG/PDF/SVG figures, ignored
  by Git.
- Other top-level folders are legacy analysis outputs or local data products;
  keep them local unless they are explicitly curated for publication.

## Recreate The Environment

Create a new conda environment from the project root:

```powershell
conda env create -f environment.yml
conda activate gplw-rc
```

If you prefer pip, install the packages listed in `environment.yml` into a clean
Python environment. See `docs/ENVIRONMENT.md` for additional setup notes.

## Run The RC Figure Pipeline

From the project root:

```powershell
cd code
python run_all_make_two_figures.py
```

Optional variants:

```powershell
python run_all_make_two_figures.py --no-show
python run_all_make_two_figures.py --export-script 01_export_csvs_for_two_figures.py --plot-script 02_plot_two_figures.py
```

Default local paths used by the figure scripts:

- Data root: `../Data`
- CMIP6 station table: `../Data/csv_from_o3_processed/df_results_modified_split.csv`
- GRDC station table: `../Data/GRDC_processed/df_results_modified_split.csv`
- Output CSV folder: `RC_figure_csvs_1minus_beta`
- Output figure folder: `RC_final_figures_1minus_beta`

Outputs:

- `code/RC_figure_csvs_1minus_beta/figure1_data.csv`
- `code/RC_figure_csvs_1minus_beta/figure2_data.csv`
- `code/RC_final_figures_1minus_beta/figure1_comprehensive.png`
- `code/RC_final_figures_1minus_beta/figure2_comprehensive.{png,pdf,svg}`

## Reproduction Cases

For the two-figure RC reproduction workflow, see:

- `docs/CASE_RC_FIGURES.md`

For the upstream flood-event extraction and event-level `R`/`C` computation
workflow, see:

- `docs/CASE_FLOOD_RC_PIPELINE.md`
- `code/compute_event_rc_from_flood_events.py`

## Git Policy

The repository tracks code, documentation, and small environment files. It does
not track large input data, generated figures, notebooks, archives, temporary
files, or local helper shims such as `git.cmd`.
