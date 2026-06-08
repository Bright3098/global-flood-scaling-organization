# RC Figure Pipeline

This folder contains the code used to generate the two RC figures:

1. `Fig1_RC_trends_map_bar_climate.*`
2. `Fig2_area_threshold_moving_window.*`

The confluence metric is defined as `C = 1 - beta`, so beta trend directions are
converted before plotting C.

## Main Files

- `run_all_make_two_figures.py`: one-command entry point.
- `01_export_csvs_for_two_figures.py`: exports standardized station tables and figure-ready CSV files.
- `02_plot_two_figures.py`: renders the two final figures from exported CSV files.
- `Fig3_A4_final_clean_1minus_beta_*.py`: separate Figure 3 variants.
- `111.py` and `2.py`: legacy robustness-analysis scripts; they are not part of the two-figure pipeline.
- `requirements.txt`: Python dependencies for this folder.

## Run

From this folder:

```powershell
python run_all_make_two_figures.py
```

Useful variants:

```powershell
python run_all_make_two_figures.py --no-climate
python run_all_make_two_figures.py --skip-export
python run_all_make_two_figures.py --main-csv ..\Data\filtered_station_summary_.xlsx
python run_all_make_two_figures.py --csv-dir .\RC_figure_csvs_1minus_beta --fig-dir .\RC_final_figures_1minus_beta
```

Default local paths:

- main table: `..\Data\filtered_station_summary_.xlsx`
- climate shapefile: the local climate-zone shapefile folder under the project root
- intermediate CSVs: `.\RC_figure_csvs_1minus_beta`
- final figures: `.\RC_final_figures_1minus_beta`

## Outputs

The pipeline writes intermediate CSV files to `RC_figure_csvs_1minus_beta` and
PNG/PDF/SVG figure files to `RC_final_figures_1minus_beta`.

For GitHub, keep source files and documentation tracked. Large data files,
generated figures, notebooks checkpoints, and temporary files are excluded by
the root `.gitignore`.
