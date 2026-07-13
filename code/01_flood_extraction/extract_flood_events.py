#!/usr/bin/env python3
"""Extract annual flood events from daily streamflow CSV files.

The script scans a folder for daily time-series CSV files, detects date and
flow columns, and extracts one flood event per station-year around the annual
peak. For each event it records the rising point, peak, recession point, and
the full event time series.

Outputs:
    ALL_flood_events_summary.csv
    ALL_flood_events_timeseries.csv
    batch_log.csv
"""

from __future__ import annotations

import argparse
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


DATE_CANDIDATES = ("Date", "date", "DATE", "time", "Time", "datetime")
FLOW_CANDIDATES = (
    "Flow (ML)",
    "Flow (ML/day)",
    "flow",
    "FLOW",
    "discharge_m3s",
    "discharge",
    "streamflow",
    "runoff",
    "Runoff",
    "q",
    "Q",
    "discharge",
)


def derive_station_name_from_filename(path: Path) -> str:
    """Infer a readable station name from a daily time-series filename."""
    stem = path.stem
    parts = stem.split("_")
    core = parts[: parts.index("daily")] if "daily" in parts else parts
    if len(core) > 1 and re.match(r"^[A-Za-z]*\d+[A-Za-z]*$", core[0]):
        core = core[1:]
    name = " ".join(core).strip()
    return name if name else stem


def detect_column(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    exact = {str(c): c for c in columns}
    lower = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand in exact:
            return exact[cand]
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def read_daily_ts_csv(path: Path) -> pd.DataFrame:
    """Read a daily streamflow CSV into columns time, runoff, station_name."""
    raw = pd.read_csv(path, comment="#", engine="python", on_bad_lines="skip")
    if raw.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")

    columns = list(raw.columns)
    time_col = detect_column(columns, DATE_CANDIDATES) or columns[0]
    flow_col = detect_column(columns, FLOW_CANDIDATES) or columns[1]

    out = pd.DataFrame(
        {
            "time": pd.to_datetime(raw[time_col], errors="coerce"),
            "runoff": pd.to_numeric(raw[flow_col], errors="coerce"),
        }
    ).dropna(subset=["time"])
    out["station_name"] = derive_station_name_from_filename(path)
    return out.sort_values("time").reset_index(drop=True)


def looks_like_daily_ts(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() == ".csv" and "daily" in name and "ts" in name


def find_event_bounds(runoff: np.ndarray, peak_pos: int, tol: float) -> tuple[int, int]:
    """Find rising and recession bounds around the annual peak."""
    left = int(peak_pos)
    while left > 0 and runoff[left - 1] <= runoff[left]:
        left -= 1

    threshold = runoff[left] * (1.0 + tol)
    right = int(peak_pos)
    while right < runoff.size - 1 and runoff[right] > threshold:
        right += 1
    return left, right


def extract_flood_events_timeseries(
    df: pd.DataFrame,
    *,
    tol: float = 0.05,
    min_valid_ratio: float = 0.5,
) -> list[dict]:
    """Extract annual flood events from one station time series."""
    required = {"time", "runoff", "station_name"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {', '.join(sorted(missing))}")

    frame = df.copy()
    frame["year"] = frame["time"].dt.year
    events: list[dict] = []
    event_counter = 1

    for year, group_raw in frame.groupby("year", sort=True):
        n_total = len(group_raw)
        n_valid = int(group_raw["runoff"].notna().sum())
        if n_total == 0 or n_valid / n_total < min_valid_ratio:
            continue

        group = (
            group_raw.dropna(subset=["runoff"])
            .sort_values("time")
            .reset_index(drop=True)
        )
        if group.empty or not np.isfinite(group["runoff"]).any():
            continue

        runoff = group["runoff"].to_numpy(dtype=float)
        peak_pos = int(np.nanargmax(runoff))
        left, right = find_event_bounds(runoff, peak_pos, tol)

        station_name = str(group.loc[peak_pos, "station_name"])
        event_label = f"{station_name}_{int(year)}_{event_counter}"

        process = group.loc[left:right, ["time", "runoff"]].copy()
        process["station_name"] = station_name
        process["event_id"] = event_label
        process["flag"] = "process"
        process.loc[process.index[0], "flag"] = "rising"
        process.loc[process.index[-1], "flag"] = "recession"
        process.loc[process.index[peak_pos - left], "flag"] = "peak"

        events.append(
            {
                "event_id": event_label,
                "station_name": station_name,
                "year": int(year),
                "rising_time": group.loc[left, "time"],
                "rising_value": float(group.loc[left, "runoff"]),
                "peak_time": group.loc[peak_pos, "time"],
                "peak_value": float(group.loc[peak_pos, "runoff"]),
                "recession_time": group.loc[right, "time"],
                "recession_value": float(group.loc[right, "runoff"]),
                "duration_days": int((group.loc[right, "time"] - group.loc[left, "time"]).days),
                "n_points": int(len(process)),
                "data": process,
            }
        )
        event_counter += 1

    return events


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "station"


def process_one_file(path: Path, out_dir: Path, tol: float, min_valid_ratio: float) -> tuple[pd.DataFrame | None, pd.DataFrame | None, dict]:
    try:
        df = read_daily_ts_csv(path)
        events = extract_flood_events_timeseries(df, tol=tol, min_valid_ratio=min_valid_ratio)
        if not events:
            return None, None, {"file": str(path), "status": "no_events"}

        summary = pd.DataFrame(
            [
                {
                    key: value
                    for key, value in event.items()
                    if key != "data"
                }
                for event in events
            ]
        )
        timeseries = pd.concat([event["data"] for event in events], ignore_index=True)

        station_dir = out_dir / safe_name(str(summary["station_name"].iloc[0]))
        station_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(station_dir / "flood_events_summary.csv", index=False)
        timeseries.to_csv(station_dir / "flood_events_timeseries.csv", index=False)

        return summary, timeseries, {"file": str(path), "status": "ok", "n_events": len(summary)}
    except Exception as exc:  # noqa: BLE001 - batch processing should log and continue.
        return None, None, {"file": str(path), "status": "error", "error": str(exc)}


def run_batch(
    input_root: Path,
    out_dir: Path,
    *,
    tol: float,
    min_valid_ratio: float,
    workers: int,
    include_all_csv: bool,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates = sorted(input_root.rglob("*.csv"))
    if not include_all_csv:
        candidates = [path for path in candidates if looks_like_daily_ts(path)]

    if not candidates:
        raise FileNotFoundError(f"No candidate CSV files found under {input_root}")

    summaries: list[pd.DataFrame] = []
    timeseries: list[pd.DataFrame] = []
    logs: list[dict] = []

    if workers <= 1:
        for path in candidates:
            summary, ts, log = process_one_file(path, out_dir, tol, min_valid_ratio)
            logs.append(log)
            if summary is not None:
                summaries.append(summary)
            if ts is not None:
                timeseries.append(ts)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_one_file, path, out_dir, tol, min_valid_ratio): path
                for path in candidates
            }
            for future in as_completed(futures):
                summary, ts, log = future.result()
                logs.append(log)
                if summary is not None:
                    summaries.append(summary)
                if ts is not None:
                    timeseries.append(ts)

    pd.DataFrame(logs).to_csv(out_dir / "batch_log.csv", index=False)
    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(
            out_dir / "ALL_flood_events_summary.csv",
            index=False,
        )
    if timeseries:
        pd.concat(timeseries, ignore_index=True).to_csv(
            out_dir / "ALL_flood_events_timeseries.csv",
            index=False,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="Folder containing daily streamflow CSV files.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output folder for event CSV files.")
    parser.add_argument("--tol", type=float, default=0.05, help="Recession tolerance relative to rising flow.")
    parser.add_argument("--min-valid-ratio", type=float, default=0.5, help="Minimum valid runoff fraction per station-year.")
    parser.add_argument("--workers", type=int, default=1, help="Parallel worker count.")
    parser.add_argument("--include-all-csv", action="store_true", help="Process every CSV, not only names containing daily and ts.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_batch(
        args.input_root,
        args.out_dir,
        tol=args.tol,
        min_valid_ratio=args.min_valid_ratio,
        workers=args.workers,
        include_all_csv=args.include_all_csv,
    )


if __name__ == "__main__":
    main()
