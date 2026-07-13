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

Print the example commands:

```powershell
python examples/example_four_methods.py
```

After editing the input paths in `examples/example_four_methods.py`, run:

```powershell
python examples/example_four_methods.py --execute
```

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
