# -*- coding: utf-8 -*-
"""Run the RC figure pipeline.

The pipeline has two steps:
1. Export figure-ready CSV tables from the station summary table.
2. Plot the two publication figures from those CSV tables.

Use command-line options to point at a different station table, climate
shapefile, or output folders without editing this file.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent

DEFAULT_MAIN_CSV = PROJECT_ROOT / "Data" / "filtered_station_summary_.xlsx"
DEFAULT_CLIMATE_SHP = PROJECT_ROOT / "\u6c14\u5019\u533a\u5212" / "\u6c14\u5019\u533a\u5212.shp"
DEFAULT_CSV_DIR = HERE / "RC_figure_csvs_1minus_beta"
DEFAULT_FIG_DIR = HERE / "RC_final_figures_1minus_beta"


def load_module(path: Path, name: str):
    """Load a Python file whose filename is not import-friendly."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export CSVs and plot the two RC figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--main-csv",
        type=Path,
        default=DEFAULT_MAIN_CSV,
        help="Station summary table. CSV, XLSX, and XLS are supported.",
    )
    parser.add_argument(
        "--climate-shp",
        type=Path,
        default=DEFAULT_CLIMATE_SHP,
        help="Koppen climate-zone shapefile used for Figure 1 panel d.",
    )
    parser.add_argument(
        "--no-climate",
        action="store_true",
        help="Skip climate-zone matching and draw panel d as missing.",
    )
    parser.add_argument(
        "--csv-dir",
        type=Path,
        default=DEFAULT_CSV_DIR,
        help="Folder for intermediate CSV outputs.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=DEFAULT_FIG_DIR,
        help="Folder for final PNG/PDF/SVG figures.",
    )
    parser.add_argument("--dpi", type=int, default=900, help="Raster output DPI.")
    parser.add_argument("--area-bins", type=int, default=60, help="Number of area bins.")
    parser.add_argument("--min-n-per-bin", type=int, default=20, help="Minimum stations per area bin.")
    parser.add_argument("--moving-window-n", type=int, default=150, help="Moving-window station count.")
    parser.add_argument(
        "--climate-net-denominator",
        choices=["all", "directional"],
        default="all",
        help="Denominator used for climate-zone net increase percentages.",
    )
    parser.add_argument("--skip-export", action="store_true", help="Only plot from existing CSV files.")
    parser.add_argument("--skip-plot", action="store_true", help="Only export CSV files.")
    return parser


def resolve_optional_climate_path(path: Path | None, no_climate: bool) -> str | None:
    if no_climate:
        return None
    if path is None:
        return None
    if not path.exists():
        print(f"[WARN] Climate shapefile not found, skip panel d data: {path}")
        return None
    return str(path)


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if args.skip_export and args.skip_plot:
        raise ValueError("Nothing to do: --skip-export and --skip-plot were both set.")
    if not args.skip_export and not args.main_csv.exists():
        raise FileNotFoundError(f"Station summary table not found: {args.main_csv}")

    export_py = HERE / "01_export_csvs_for_two_figures.py"
    plot_py = HERE / "02_plot_two_figures.py"
    if not export_py.exists():
        raise FileNotFoundError(f"Missing: {export_py}")
    if not plot_py.exists():
        raise FileNotFoundError(f"Missing: {plot_py}")

    csv_dir = args.csv_dir
    fig_dir = args.fig_dir
    csv_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    climate_shp = resolve_optional_climate_path(args.climate_shp, args.no_climate)

    print("[CONFIG]")
    print("  main_csv:", args.main_csv)
    print("  climate_shp:", climate_shp or "None")
    print("  csv_dir:", csv_dir)
    print("  fig_dir:", fig_dir)

    if not args.skip_export:
        print("=" * 72)
        print("STEP 1/2: Export figure CSVs")
        print("=" * 72)
        export_mod = load_module(export_py, "export_rc_csvs_1minus_beta")
        export_mod.MAIN_CSV = str(args.main_csv)
        export_mod.CLIMATE_SHP = climate_shp
        export_mod.OUT_DIR = csv_dir
        export_mod.N_AREA_BINS = args.area_bins
        export_mod.MIN_N_PER_BIN = args.min_n_per_bin
        export_mod.MOVING_WINDOW_N = args.moving_window_n
        export_mod.CLIMATE_NET_DENOMINATOR = args.climate_net_denominator
        export_mod.OUT_DIR.mkdir(parents=True, exist_ok=True)
        export_mod.main()

    if not args.skip_plot:
        print("\n" + "=" * 72)
        print("STEP 2/2: Plot two figures")
        print("=" * 72)
        plot_mod = load_module(plot_py, "plot_rc_two_figures")
        plot_mod.CSV_DIR = csv_dir
        plot_mod.FIG_DIR = fig_dir
        plot_mod.DPI = args.dpi
        plot_mod.FIG_DIR.mkdir(parents=True, exist_ok=True)
        plot_mod.main()

    print("\n[DONE]")
    print("CSV output:", csv_dir.resolve())
    print("Figure output:", fig_dir.resolve())


if __name__ == "__main__":
    main()
