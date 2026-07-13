#!/usr/bin/env python3
"""Compute station and global trends for annual flood peak discharge.

Input is an event summary table with at least station, year, and peak-flow
columns. If a station has multiple events in one year, the annual maximum peak
is used. Each station is normalized by its own baseline mean before global
aggregation.

Outputs:
    station_year_qpeak_normalized.csv
    annual_qpeak_index.csv
    station_qpeak_trends.csv
    qpeak_trend_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


STATION_CANDIDATES = ("station_id", "station_name", "folder", "site_no", "gauge_id")
PEAK_CANDIDATES = ("peak_runoff", "peak_value", "peak_value(ML/day)", "annual_daily_qpeak", "qpeak")


def find_column(columns: list[str], requested: str | None, candidates: tuple[str, ...]) -> str:
    if requested:
        if requested not in columns:
            raise KeyError(f"Column not found: {requested}")
        return requested
    lower = {str(c).lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    raise KeyError(f"Cannot find any of {candidates}. Existing columns: {columns}")


def station_id_from_series(series: pd.Series) -> pd.Series:
    station = (
        series.astype(str)
        .str.replace("\\", "/", regex=False)
        .str.rstrip("/")
        .str.split("/")
        .str[-1]
    )
    return station.replace({"": np.nan, "nan": np.nan})


def compute_station_year(events: pd.DataFrame, station_col: str, year_col: str, peak_col: str) -> pd.DataFrame:
    station = station_id_from_series(events[station_col]) if station_col == "folder" else events[station_col]
    frame = pd.DataFrame(
        {
            "station_id": station,
            "year": pd.to_numeric(events[year_col], errors="coerce"),
            "peak_value": pd.to_numeric(events[peak_col], errors="coerce"),
        }
    ).dropna(subset=["station_id", "year", "peak_value"])
    frame = frame[frame["peak_value"] > 0].copy()
    frame["year"] = frame["year"].astype(int)
    return (
        frame.groupby(["station_id", "year"], as_index=False)
        .agg(annual_qpeak=("peak_value", "max"), events_in_year=("peak_value", "size"))
    )


def normalize_by_baseline(
    station_year: pd.DataFrame,
    *,
    year_start: int,
    year_end: int,
    baseline_start: int,
    baseline_end: int,
    min_years_total: int,
    min_years_baseline: int,
    min_years_recent: int,
) -> pd.DataFrame:
    station_year = station_year[station_year["year"].between(year_start, year_end)].copy()
    baseline = (
        station_year[station_year["year"].between(baseline_start, baseline_end)]
        .groupby("station_id")
        .agg(baseline_qpeak_mean=("annual_qpeak", "mean"), baseline_years=("year", "nunique"))
    )
    coverage = (
        station_year.groupby("station_id")
        .agg(total_years=("year", "nunique"))
        .join(baseline, how="left")
    )
    recent = (
        station_year[station_year["year"].between(baseline_end + 1, year_end)]
        .groupby("station_id")["year"]
        .nunique()
        .rename("recent_years")
    )
    coverage = coverage.join(recent, how="left").fillna({"recent_years": 0})
    eligible = coverage[
        (coverage["total_years"] >= min_years_total)
        & (coverage["baseline_years"] >= min_years_baseline)
        & (coverage["recent_years"] >= min_years_recent)
        & (coverage["baseline_qpeak_mean"] > 0)
    ]
    out = station_year[station_year["station_id"].isin(eligible.index)].merge(
        eligible[["baseline_qpeak_mean"]],
        left_on="station_id",
        right_index=True,
        how="left",
    )
    out["qpeak_index_pct"] = out["annual_qpeak"] / out["baseline_qpeak_mean"] * 100.0
    return out


def trend_for_xy(x: np.ndarray, y: np.ndarray) -> dict:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return {
            "slope_pct_per_year": np.nan,
            "slope_pct_per_decade": np.nan,
            "slope_low_pct_per_decade": np.nan,
            "slope_high_pct_per_decade": np.nan,
            "kendall_tau": np.nan,
            "kendall_p": np.nan,
        }
    slope, _intercept, low, high = stats.theilslopes(y, x, alpha=0.95)
    tau, p_value = stats.kendalltau(x, y)
    return {
        "slope_pct_per_year": float(slope),
        "slope_pct_per_decade": float(slope * 10.0),
        "slope_low_pct_per_decade": float(low * 10.0),
        "slope_high_pct_per_decade": float(high * 10.0),
        "kendall_tau": float(tau),
        "kendall_p": float(p_value),
    }


def compute_trends(normalized: pd.DataFrame, min_annual_stations: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    annual = (
        normalized.groupby("year", as_index=False)
        .agg(
            median_index_pct=("qpeak_index_pct", "median"),
            mean_index_pct=("qpeak_index_pct", "mean"),
            p25_index_pct=("qpeak_index_pct", lambda x: np.nanpercentile(x, 25)),
            p75_index_pct=("qpeak_index_pct", lambda x: np.nanpercentile(x, 75)),
            n_stations=("station_id", "nunique"),
        )
    )

    station_rows = []
    for station_id, group in normalized.groupby("station_id"):
        trend = trend_for_xy(group["year"].to_numpy(float), group["qpeak_index_pct"].to_numpy(float))
        trend["station_id"] = station_id
        trend["n_years"] = int(group["year"].nunique())
        if np.isfinite(trend["slope_pct_per_year"]):
            trend["trend_direction"] = "increase" if trend["slope_pct_per_year"] >= 0 else "decrease"
        else:
            trend["trend_direction"] = "unknown"
        trend["trend_significant"] = bool(np.isfinite(trend["kendall_p"]) and trend["kendall_p"] < 0.05)
        station_rows.append(trend)
    station_trends = pd.DataFrame(station_rows)

    annual_for_trend = annual[annual["n_stations"] >= min_annual_stations]
    global_trend = trend_for_xy(
        annual_for_trend["year"].to_numpy(float),
        annual_for_trend["median_index_pct"].to_numpy(float),
    )
    summary = pd.DataFrame(
        [
            {"metric": "n_stations", "value": normalized["station_id"].nunique()},
            {"metric": "n_station_years", "value": len(normalized)},
            {"metric": "annual_median_slope_pct_per_decade", "value": global_trend["slope_pct_per_decade"]},
            {"metric": "annual_median_kendall_tau", "value": global_trend["kendall_tau"]},
            {"metric": "annual_median_kendall_p", "value": global_trend["kendall_p"]},
            {"metric": "station_decrease_count", "value": int((station_trends["slope_pct_per_year"] < 0).sum())},
            {"metric": "station_increase_count", "value": int((station_trends["slope_pct_per_year"] >= 0).sum())},
        ]
    )
    return annual, station_trends, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-csv", type=Path, required=True, help="Flood event summary CSV.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory for trend tables.")
    parser.add_argument("--station-col", help="Station column. Auto-detected if omitted.")
    parser.add_argument("--year-col", default="year", help="Year column.")
    parser.add_argument("--peak-col", help="Peak discharge/runoff column. Auto-detected if omitted.")
    parser.add_argument("--year-start", type=int, default=1960)
    parser.add_argument("--year-end", type=int, default=2015)
    parser.add_argument("--baseline-start", type=int, default=1960)
    parser.add_argument("--baseline-end", type=int, default=1989)
    parser.add_argument("--min-years-total", type=int, default=30)
    parser.add_argument("--min-years-baseline", type=int, default=10)
    parser.add_argument("--min-years-recent", type=int, default=10)
    parser.add_argument("--min-annual-stations", type=int, default=30)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(args.events_csv)
    station_col = find_column(list(events.columns), args.station_col, STATION_CANDIDATES)
    peak_col = find_column(list(events.columns), args.peak_col, PEAK_CANDIDATES)
    if args.year_col not in events.columns:
        raise KeyError(f"Year column not found: {args.year_col}")

    station_year = compute_station_year(events, station_col, args.year_col, peak_col)
    normalized = normalize_by_baseline(
        station_year,
        year_start=args.year_start,
        year_end=args.year_end,
        baseline_start=args.baseline_start,
        baseline_end=args.baseline_end,
        min_years_total=args.min_years_total,
        min_years_baseline=args.min_years_baseline,
        min_years_recent=args.min_years_recent,
    )
    annual, station_trends, summary = compute_trends(normalized, args.min_annual_stations)

    normalized.to_csv(args.out_dir / "station_year_qpeak_normalized.csv", index=False)
    annual.to_csv(args.out_dir / "annual_qpeak_index.csv", index=False)
    station_trends.to_csv(args.out_dir / "station_qpeak_trends.csv", index=False)
    summary.to_csv(args.out_dir / "qpeak_trend_summary.csv", index=False)


if __name__ == "__main__":
    main()
