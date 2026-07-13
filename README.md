# global-flood-scaling-organization

Minimal code repository for four hydrology workflows:

1. Flood event extraction from daily streamflow CSV files.
2. Event-level runoff generation `R` and confluence `C` calculation.
3. Flood-peak trend calculation.
4. Otsu-like drainage-area threshold detection.

Large data, generated figures, notebooks, GIS layers, and local workspace
outputs are intentionally excluded from GitHub.

## Repository Layout

- `code/01_flood_extraction/extract_flood_events.py`
  Extracts annual flood events around station-year peaks and writes event
  summary/time-series CSV files.
- `code/02_rc_calculation/compute_event_rc_from_flood_events.py`
  Computes event-level `R` and `C` from `flood_events__*.nc` NetCDF files.
- `code/03_trend_calculation/compute_qpeak_trends.py`
  Computes station and global trends for annual flood peak discharge.
- `code/04_otsu_like_method/otsu_like_area_threshold.py`
  Finds an Otsu-like area threshold that maximizes trend-direction separation.
- `examples/example_four_methods.py`
  Example command file showing how to run the four workflows.
- `examples/case_data/`
  Small real-data subsets for runnable examples.

## Environment

Create a clean environment:

```powershell
conda env create -f environment.yml
conda activate gplw-core
```

Or install the Python requirements:

```powershell
pip install -r code/requirements.txt
```

## Example

Print the example commands for the included real case data:

```powershell
python examples/example_four_methods.py
```

Run the included case data:

```powershell
python examples/example_four_methods.py --execute
```

The case data are small subsets copied or converted from the local GPLW
workspace. They are included only to demonstrate file formats and command usage;
they are not the full research datasets.

Expected demo outputs are written to `outputs/case_data_run/`:

- `01_flood_events_csv/ALL_flood_events_summary.csv`
- `01_flood_events_csv/ALL_flood_events_timeseries.csv`
- `02_event_rc/scaling_RC__barcoo_blackall_case.nc`
- `03_qpeak_trends/qpeak_trend_summary.csv`
- `04_otsu_like_threshold/otsu_like_threshold_summary.csv`

On a normal desktop or laptop, the included case-data demo should finish in
under one minute.

## Individual Commands

Flood extraction:

```powershell
python code/01_flood_extraction/extract_flood_events.py --input-root Data/daily_streamflow_csv --out-dir outputs/flood_events_csv --workers 4
```

R/C calculation:

```powershell
python code/02_rc_calculation/compute_event_rc_from_flood_events.py --in-dir Data/flood_event_nc --out-dir outputs/event_rc
```

Trend calculation:

```powershell
python code/03_trend_calculation/compute_qpeak_trends.py --events-csv outputs/flood_events_csv/ALL_flood_events_summary.csv --out-dir outputs/qpeak_trends
```

Otsu-like threshold:

```powershell
python code/04_otsu_like_method/otsu_like_area_threshold.py --input Data/station_area_trend_table.csv --out-dir outputs/otsu_like_threshold --area-col area_km2 --trend-col trend_direction
```

## Submission Checklist

For a Nature Research-style code/software submission checklist, see:

- `docs/NATURE_CODE_CHECKLIST.md`
