# Case: Reproduce RC Figures In A New Environment

## Goal

Use this repository on a clean workstation to regenerate the two main RC
figures:

- `Fig1_RC_trends_map_bar_climate.*`
- `Fig2_area_threshold_moving_window.*`

The case is designed for a GitHub checkout where source code is tracked, while
large data and generated outputs are supplied locally and ignored by Git.

## What This Case Proves

- The project environment can be rebuilt from `environment.yml`.
- The main station table can be converted into figure-ready CSV outputs.
- The plotting step can regenerate the publication figure files.
- Optional climate-zone matching can be enabled when the local Koppen shapefile
  is available, or skipped for a quick smoke test.

## Tracked Repository Inputs

These files are committed in Git and should exist after cloning:

- `environment.yml`
- `code/requirements.txt`
- `code/run_all_make_two_figures.py`
- `code/01_export_csvs_for_two_figures.py`
- `code/02_plot_two_figures.py`

## Local Data Inputs

Place local data under the project root. These files are intentionally ignored
by Git.

Required:

- `Data/filtered_station_summary_.xlsx`

Optional for Figure 1 climate-zone panel:

- `气候区划/气候区划.shp`
- `气候区划/气候区划.dbf`
- `气候区划/气候区划.shx`
- `气候区划/气候区划.prj`
- `气候区划/气候区划.cpg`, if available

## Setup

From the project root:

```powershell
conda env create -f environment.yml
conda activate gplw-rc
python code\run_all_make_two_figures.py --help
```

The `--help` command is the first smoke test. It should print the pipeline
options without requiring data files.

## Quick Smoke Test Without Climate Matching

Use this when the station table is available but the climate shapefile is not:

```powershell
cd code
python run_all_make_two_figures.py --no-climate
```

Expected outputs:

- `code/RC_figure_csvs_1minus_beta/00_standardized_station_RC.csv`
- `code/RC_figure_csvs_1minus_beta/11_method_metadata.json`
- `code/RC_final_figures_1minus_beta/Fig1_RC_trends_map_bar_climate.png`
- `code/RC_final_figures_1minus_beta/Fig2_area_threshold_moving_window.png`

In this mode, Figure 1 panel d is expected to show that climate-zone data is
missing or skipped.

## Full Figure Reproduction

Use this when both the station table and Koppen climate-zone shapefile are
available:

```powershell
cd code
python run_all_make_two_figures.py ^
  --main-csv ..\Data\filtered_station_summary_.xlsx ^
  --climate-shp ..\气候区划\气候区划.shp
```

Expected outputs:

- Intermediate CSV files under `code/RC_figure_csvs_1minus_beta/`
- PNG, PDF, and SVG figure files under `code/RC_final_figures_1minus_beta/`
- Climate-zone summary CSV:
  `code/RC_figure_csvs_1minus_beta/04_climate_zone_net_trend_summary.csv`

## Plot From Existing CSVs

After CSV export has succeeded once, regenerate only the figures:

```powershell
cd code
python run_all_make_two_figures.py --skip-export
```

## Pass Criteria

The case passes when:

- `python code\run_all_make_two_figures.py --help` works in the rebuilt
  environment.
- The export step creates the expected CSV folder.
- The plot step creates Figure 1 and Figure 2 files.
- Generated outputs remain untracked because `.gitignore` excludes them.

## Troubleshooting

- If `cartopy` fails to install with pip, use the conda environment.
- If climate matching fails, rerun with `--no-climate` to check the rest of the
  pipeline first.
- If the station table is stored elsewhere, pass its path with `--main-csv`.
- If map fonts differ on a new computer, matplotlib will fall back to an
  available font; this may slightly change text layout.
