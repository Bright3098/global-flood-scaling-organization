# Code

This folder keeps only the four core workflows intended for GitHub.

## 01 Flood Extraction

`01_flood_extraction/extract_flood_events.py`

Reads daily streamflow CSV files, extracts annual flood events around each
station-year peak, and writes merged summary/time-series CSV outputs.

## 02 R/C Calculation

`02_rc_calculation/compute_event_rc_from_flood_events.py`

Reads `flood_events__*.nc` files and computes event-level `R` and `C` from the
volume scaling relation:

```text
log10(W) = log10(R) + C * log10(T_hours)
```

## 03 Trend Calculation

`03_trend_calculation/compute_qpeak_trends.py`

Reads a flood-event summary table, converts multiple events per station-year to
annual maximum Qpeak, normalizes stations by baseline mean Qpeak, and estimates
Theil-Sen/Kendall trends.

## 04 Otsu-Like Method

`04_otsu_like_method/otsu_like_area_threshold.py`

Scans log10 drainage-area thresholds and selects the split that maximizes
between-group trend-direction separation.
