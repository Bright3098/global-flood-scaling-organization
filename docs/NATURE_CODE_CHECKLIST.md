# Nature Research Code and Software Submission Checklist

This repository is prepared as a minimal code and case-data package for review.

## Required Content

### Source Code

The repository includes source code for four workflows:

- `code/01_flood_extraction/extract_flood_events.py`
- `code/02_rc_calculation/compute_event_rc_from_flood_events.py`
- `code/03_trend_calculation/compute_qpeak_trends.py`
- `code/04_otsu_like_method/otsu_like_area_threshold.py`

No compiled standalone software is required.

### Demo Dataset

Small real-data subsets are included under `examples/case_data/`:

- Daily streamflow CSV for flood extraction.
- A small `flood_events__*.nc` NetCDF file for `R` and `C` calculation.
- A flood-event summary CSV for trend calculation.
- A station area/trend CSV for the Otsu-like threshold method.

These files are real GPLW workspace subsets, not synthetic data. Full research
datasets are intentionally excluded from GitHub.

## README Requirements

### System Requirements

Tested operating system:

- Windows 10

Expected compatible systems:

- Windows, Linux, or macOS with Python and the listed scientific Python
  packages installed.

Software dependencies:

- Python 3.11 or newer
- numpy
- pandas
- scipy
- netCDF4
- openpyxl

Environment files:

- `environment.yml`
- `code/requirements.txt`

Tested version snapshot from the local machine:

- Python 3.13.5
- numpy 2.1.3
- pandas 2.2.3
- scipy 1.15.3
- netCDF4 1.7.2

No non-standard hardware is required. The code runs on CPU.

### Installation Guide

Recommended conda installation:

```powershell
conda env create -f environment.yml
conda activate gplw-core
```

Alternative pip installation:

```powershell
pip install -r code/requirements.txt
```

Typical install time on a normal desktop computer is a few minutes with conda or
pip, depending on network speed.

### Demo Instructions

Print the commands:

```powershell
python examples/example_four_methods.py
```

Run all four demos:

```powershell
python examples/example_four_methods.py --execute
```

Expected runtime for the included demo on a normal desktop computer is under one
minute.

Expected output folder:

```text
outputs/case_data_run/
```

Expected key files:

- `outputs/case_data_run/01_flood_events_csv/ALL_flood_events_summary.csv`
- `outputs/case_data_run/01_flood_events_csv/ALL_flood_events_timeseries.csv`
- `outputs/case_data_run/02_event_rc/scaling_RC__barcoo_blackall_case.nc`
- `outputs/case_data_run/03_qpeak_trends/qpeak_trend_summary.csv`
- `outputs/case_data_run/04_otsu_like_threshold/otsu_like_threshold_summary.csv`

### Instructions For Use On New Data

Flood event extraction:

```powershell
python code/01_flood_extraction/extract_flood_events.py --input-root Data/daily_streamflow_csv --out-dir outputs/flood_events_csv --workers 4
```

`R` and `C` calculation:

```powershell
python code/02_rc_calculation/compute_event_rc_from_flood_events.py --in-dir Data/flood_event_nc --out-dir outputs/event_rc
```

Flood-peak trend calculation:

```powershell
python code/03_trend_calculation/compute_qpeak_trends.py --events-csv outputs/flood_events_csv/ALL_flood_events_summary.csv --out-dir outputs/qpeak_trends
```

Otsu-like area threshold:

```powershell
python code/04_otsu_like_method/otsu_like_area_threshold.py --input Data/station_area_trend_table.csv --out-dir outputs/otsu_like_threshold --area-col area_km2 --trend-col trend_direction
```

## Additional Information

### License

The repository uses the MIT License. See `LICENSE`.

### Open Source Repository

Repository link:

https://github.com/Bright3098/global-flood-scaling-organization

### Code Functionality

High-level pseudocode:

1. Flood extraction
   - Read daily streamflow CSV files.
   - Detect date and discharge columns.
   - For each station-year, identify the annual peak.
   - Expand backward to the rising limb and forward to recession.
   - Save event summaries and event time series.

2. `R` and `C` calculation
   - Read `flood_events__*.nc` files.
   - For each event, calculate maximum accumulated water volume `W(T)` over
     multiple durations `T`.
   - Fit `log10(W) = log10(R) + C * log10(T_hours)`.
   - Save event-level `R`, `C`, fit `r2`, and point counts.

3. Flood-peak trend calculation
   - Read flood-event summaries.
   - Convert multiple events per station-year to annual maximum peak flow.
   - Normalize each station by its baseline mean.
   - Estimate Theil-Sen slopes and Kendall trend significance.
   - Save station-level and annual trend summaries.

4. Otsu-like area threshold
   - Read a station area/trend table.
   - Convert trend direction to a binary increase/decrease indicator.
   - Scan candidate thresholds along `log10(area)`.
   - Select the threshold maximizing between-group trend separation.
   - Save prepared input, full threshold scan, and best-threshold summary.
