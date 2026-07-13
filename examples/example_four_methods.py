#!/usr/bin/env python3
"""Example command file for the four GitHub methods.

Edit the paths near the top, run once without --execute to review commands, then
run with --execute when the input files exist.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# Edit these paths for a real case.
DAILY_STREAMFLOW_DIR = ROOT / "Data" / "daily_streamflow_csv"
FLOOD_EVENT_NC_DIR = ROOT / "Data" / "flood_event_nc"
AREA_TREND_TABLE = ROOT / "Data" / "station_area_trend_table.csv"
OUTPUT_DIR = ROOT / "outputs" / "example_case"


def commands() -> list[list[str]]:
    flood_csv_dir = OUTPUT_DIR / "01_flood_events_csv"
    rc_dir = OUTPUT_DIR / "02_event_rc"
    trend_dir = OUTPUT_DIR / "03_qpeak_trends"
    otsu_dir = OUTPUT_DIR / "04_otsu_like_threshold"

    return [
        [
            sys.executable,
            str(ROOT / "code" / "01_flood_extraction" / "extract_flood_events.py"),
            "--input-root",
            str(DAILY_STREAMFLOW_DIR),
            "--out-dir",
            str(flood_csv_dir),
            "--workers",
            "4",
        ],
        [
            sys.executable,
            str(ROOT / "code" / "02_rc_calculation" / "compute_event_rc_from_flood_events.py"),
            "--in-dir",
            str(FLOOD_EVENT_NC_DIR),
            "--out-dir",
            str(rc_dir),
        ],
        [
            sys.executable,
            str(ROOT / "code" / "03_trend_calculation" / "compute_qpeak_trends.py"),
            "--events-csv",
            str(flood_csv_dir / "ALL_flood_events_summary.csv"),
            "--out-dir",
            str(trend_dir),
        ],
        [
            sys.executable,
            str(ROOT / "code" / "04_otsu_like_method" / "otsu_like_area_threshold.py"),
            "--input",
            str(AREA_TREND_TABLE),
            "--out-dir",
            str(otsu_dir),
            "--area-col",
            "area_km2",
            "--trend-col",
            "trend_direction",
        ],
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Run the example commands.")
    args = parser.parse_args()

    for command in commands():
        print("\n" + " ".join(command))
        if args.execute:
            subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
