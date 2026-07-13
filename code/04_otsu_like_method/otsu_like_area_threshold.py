#!/usr/bin/env python3
"""Compute an Otsu-like drainage-area threshold for trend separation.

The method scans candidate thresholds along log10(area) and splits stations
into small/large basins. For each threshold it computes the fraction of
stations with increasing trends on both sides and selects the threshold that
maximizes the squared separation:

    J = (p_small_increase - p_large_increase)^2

This is Otsu-like because it selects the split with the strongest between-class
separation, but the target variable is a binary trend direction rather than a
grayscale histogram.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


AREA_CANDIDATES = (
    "area_km2",
    "basin_area_km2",
    "drainage_area_km2",
    "catchment_area_km2",
    "area",
    "basin_area",
)
TREND_CANDIDATES = (
    "trend_direction",
    "Slope_Trend_MK",
    "Intercept_Trend_MK",
    "Slope_Trend_MK(1980-2015)",
    "Intercept_Trend_MK(1980-2015)",
)


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported table type: {path}")


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


def normalize_area_km2(series: pd.Series, column_name: str) -> pd.Series:
    area = pd.to_numeric(series, errors="coerce")
    finite = area[np.isfinite(area)]
    if len(finite) and "km" not in column_name.lower() and float(finite.max()) > 1e7:
        area = area / 1e6
    return area.where(area > 0)


def trend_to_binary(series: pd.Series, *, invert: bool, require_significant: bool) -> pd.Series:
    text = series.astype(str).str.lower()
    numeric = pd.to_numeric(series, errors="coerce")
    significant = text.str.contains("*", regex=False) | text.str.contains("significant", regex=False)
    increasing = text.str.contains("increas", regex=False) | text.str.contains("positive", regex=False)
    decreasing = text.str.contains("decreas", regex=False) | text.str.contains("negative", regex=False)

    if require_significant:
        increasing = increasing & significant
        decreasing = decreasing & significant

    out = pd.Series(np.nan, index=series.index, dtype="float64")
    out.loc[numeric > 0] = 1.0
    out.loc[numeric < 0] = 0.0
    out.loc[increasing] = 1.0
    out.loc[decreasing] = 0.0
    if invert:
        out = 1.0 - out
    return out


def prepare_input(
    table: pd.DataFrame,
    *,
    area_col: str,
    trend_col: str,
    group_col: str | None,
    group_value: str | None,
    invert_trend: bool,
    require_significant: bool,
) -> pd.DataFrame:
    frame = table.copy()
    if group_col and group_value is not None:
        frame = frame[frame[group_col].astype(str) == str(group_value)].copy()

    out = pd.DataFrame(
        {
            "area_km2": normalize_area_km2(frame[area_col], area_col),
            "is_increase": trend_to_binary(
                frame[trend_col],
                invert=invert_trend,
                require_significant=require_significant,
            ),
        }
    ).dropna()
    out = out[out["area_km2"] > 0].copy()
    out["log10_area"] = np.log10(out["area_km2"])
    return out


def scan_thresholds(
    frame: pd.DataFrame,
    *,
    q_min: float,
    q_max: float,
    q_count: int,
    min_side_n: int,
) -> pd.DataFrame:
    quantiles = np.linspace(q_min, q_max, q_count)
    candidates = np.unique(frame["log10_area"].quantile(quantiles).to_numpy(float))
    rows = []
    for threshold in candidates:
        small = frame["log10_area"] <= threshold
        large = frame["log10_area"] > threshold
        n_small = int(small.sum())
        n_large = int(large.sum())
        if n_small < min_side_n or n_large < min_side_n:
            continue
        p_small = float(frame.loc[small, "is_increase"].mean())
        p_large = float(frame.loc[large, "is_increase"].mean())
        rows.append(
            {
                "threshold_log10_area": float(threshold),
                "threshold_area_km2": float(10 ** threshold),
                "n_small": n_small,
                "n_large": n_large,
                "p_small_increase": p_small,
                "p_large_increase": p_large,
                "delta_large_minus_small": p_large - p_small,
                "J_separation": float((p_small - p_large) ** 2),
            }
        )
    scan = pd.DataFrame(rows)
    if scan.empty:
        return scan
    scan["is_best_threshold"] = False
    scan.loc[scan["J_separation"].idxmax(), "is_best_threshold"] = True
    return scan


def bootstrap_thresholds(
    frame: pd.DataFrame,
    *,
    n_boot: int,
    sample_fraction: float,
    seed: int,
    q_min: float,
    q_max: float,
    q_count: int,
    min_side_n: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    sample_n = max(2 * min_side_n, int(round(len(frame) * sample_fraction)))
    for i in range(n_boot):
        sample = frame.iloc[rng.choice(len(frame), size=sample_n, replace=False)]
        scan = scan_thresholds(sample, q_min=q_min, q_max=q_max, q_count=q_count, min_side_n=min_side_n)
        if scan.empty:
            rows.append({"bootstrap_id": i, "status": "no_valid_threshold"})
            continue
        best = scan.loc[scan["J_separation"].idxmax()].to_dict()
        best["bootstrap_id"] = i
        best["status"] = "ok"
        rows.append(best)
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Input CSV/XLSX station table.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--area-col", help="Area column. Auto-detected if omitted.")
    parser.add_argument("--trend-col", help="Trend-direction column. Auto-detected if omitted.")
    parser.add_argument("--group-col", help="Optional group column.")
    parser.add_argument("--group-value", help="Optional group value to keep.")
    parser.add_argument("--invert-trend", action="store_true", help="Invert increase/decrease labels before scanning.")
    parser.add_argument("--require-significant", action="store_true", help="Treat non-significant trends as missing.")
    parser.add_argument("--q-min", type=float, default=0.10)
    parser.add_argument("--q-max", type=float, default=0.90)
    parser.add_argument("--q-count", type=int, default=101)
    parser.add_argument("--min-side-n", type=int, default=25)
    parser.add_argument("--bootstrap", type=int, default=0, help="Optional bootstrap repetitions.")
    parser.add_argument("--sample-fraction", type=float, default=0.90)
    parser.add_argument("--seed", type=int, default=3098)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    table = read_table(args.input)
    area_col = find_column(list(table.columns), args.area_col, AREA_CANDIDATES)
    trend_col = find_column(list(table.columns), args.trend_col, TREND_CANDIDATES)
    if args.group_col and args.group_col not in table.columns:
        raise KeyError(f"Group column not found: {args.group_col}")

    prepared = prepare_input(
        table,
        area_col=area_col,
        trend_col=trend_col,
        group_col=args.group_col,
        group_value=args.group_value,
        invert_trend=args.invert_trend,
        require_significant=args.require_significant,
    )
    scan = scan_thresholds(
        prepared,
        q_min=args.q_min,
        q_max=args.q_max,
        q_count=args.q_count,
        min_side_n=args.min_side_n,
    )
    if scan.empty:
        raise RuntimeError("No valid threshold. Lower --min-side-n or check input data.")

    best = scan.loc[scan["J_separation"].idxmax()].copy()
    summary = pd.DataFrame(
        [
            {
                "input": str(args.input),
                "area_col": area_col,
                "trend_col": trend_col,
                "group_col": args.group_col or "",
                "group_value": args.group_value or "",
                "n_stations": len(prepared),
                **best.to_dict(),
            }
        ]
    )

    prepared.to_csv(args.out_dir / "otsu_like_prepared_input.csv", index=False)
    scan.to_csv(args.out_dir / "otsu_like_threshold_scan.csv", index=False)
    summary.to_csv(args.out_dir / "otsu_like_threshold_summary.csv", index=False)

    if args.bootstrap > 0:
        boot = bootstrap_thresholds(
            prepared,
            n_boot=args.bootstrap,
            sample_fraction=args.sample_fraction,
            seed=args.seed,
            q_min=args.q_min,
            q_max=args.q_max,
            q_count=args.q_count,
            min_side_n=args.min_side_n,
        )
        boot.to_csv(args.out_dir / "otsu_like_threshold_bootstrap.csv", index=False)


if __name__ == "__main__":
    main()
