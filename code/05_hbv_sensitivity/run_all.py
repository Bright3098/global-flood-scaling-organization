"""Run the self-contained HBV experiment and check the published input settings.

Usage: python run_all.py [--output-dir PATH] [--skip-figures]
No observational input, calibration, network access, or parent modules are needed.
Do not use python -O: the numerical validation includes assertions.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs",
                        help="Output directory; existing files with matching names are replaced.")
    parser.add_argument("--skip-figures", action="store_true",
                        help="Run simulation, fitting and numerical checks using NumPy only.")
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Run without -O so that numerical checks are enabled.")
    started = perf_counter()
    out = args.output_dir.resolve()
    analysis = out / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    step_times = {}
    with (out / "run.log").open("w", encoding="utf-8") as log:
        with redirect_stdout(log):
            step_start = perf_counter()
            model = load_module("hbv_sensitivity")
            model.OUT = analysis
            model.main()
            step_times["simulation_and_model_checks"] = perf_counter() - step_start

            step_start = perf_counter()
            fit = load_module("fit_indices")
            fit.ROOT = analysis
            fit.SOURCE = analysis / "hydrographs.csv"
            fit.main()
            step_times["index_fitting_and_checks"] = perf_counter() - step_start

            if not args.skip_figures:
                step_start = perf_counter()
                plot = load_module("plot_hbv")
                plot.SOURCE = analysis / "hydrographs.csv"
                plot.INDICES = analysis / "fitted_indices.csv"
                plot.OUT = out / "figures"
                plot.main()
                step_times["figure_generation_and_interpolation_checks"] = perf_counter() - step_start

    reference = json.loads((HERE / "reference_values.json").read_text(encoding="utf-8"))
    indices = {r["scenario_id"]: r for r in read_csv(analysis / "fitted_indices.csv")}
    metrics = {r["scenario_id"]: r for r in read_csv(analysis / "scenario_metrics.csv")}
    if set(indices) != set(reference["scenarios"]) or len(indices) != 14:
        raise AssertionError("The 14 scenario labels do not match the reference experiment.")
    maximum_errors = {}
    comparison_rows = []
    for scenario_id, expected in reference["scenarios"].items():
        combined = {**metrics[scenario_id], **indices[scenario_id]}
        for key, expected_value in expected.items():
            error = abs(float(combined[key]) - expected_value)
            maximum_errors[key] = max(maximum_errors.get(key, 0.0), error)
            if error > reference["comparison_tolerance"]:
                raise AssertionError((scenario_id, key, error))
        comparison_rows.append({
            "scenario_id": scenario_id, "panel": combined["panel"],
            **{k: float(combined[k]) for k in
               ("FAI", "FCI", "R_squared", "peak_Q_mm_day",
                "peak_time_from_first_rainfall_days", "runoff_volume_mm",
                "whole_balance_error_mm", "maximum_daily_balance_error_mm")},
        })
    with (out / "scenario_verification.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison_rows[0]))
        writer.writeheader()
        writer.writerows(comparison_rows)

    model_checks = json.loads((analysis / "validation.json").read_text(encoding="utf-8"))
    index_checks = json.loads((analysis / "index_fit_validation.json").read_text(encoding="utf-8"))
    expected_config = (
        model.N == 730 and model.START == 5 and model.TRIBUTARY_LAGS == (0, 2, 3)
        and model_checks["plotted_hydrographs"] == 14
        and index_checks["fit_durations_days"] == list(range(1, 21))
        and index_checks["displayed_calendar_interval_days"] == [0, 20]
        and model_checks["all_checks_passed"] and index_checks["checks_passed"]
    )
    if not expected_config:
        raise AssertionError("Simulation settings differ from the documented experiment.")
    dependencies = {}
    for name in ("numpy", "reportlab", "pypdfium2", "Pillow"):
        try:
            dependencies[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            dependencies[name] = None
    report = {
        "status": "PASS",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "dependencies": dependencies,
        "figures_generated": not args.skip_figures,
        "step_seconds": step_times,
        "total_seconds": perf_counter() - started,
        "timing_scope": "Measured local run including imports, calculation, checks and requested output generation; installation excluded.",
        "reference_version": reference["source_version"],
        "reference_comparison_tolerance": reference["comparison_tolerance"],
        "maximum_reference_errors": maximum_errors,
        "scenario_count": len(indices),
        "simulation_days": model.N,
        "routing_support_widths_days": [1, 3, 7],
        "tributary_lags_days": list(model.TRIBUTARY_LAGS),
        "max_daily_balance_error_mm": model_checks["max_daily_balance_error_mm"],
        "max_whole_balance_error_mm": model_checks["max_whole_balance_error_mm"],
        "panel_d_single_peak_checks": model_checks["panel_d_single_peak_checks"],
        "source_hashes_sha256": {
            name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in ("hbv_sensitivity.py", "fit_indices.py", "plot_hbv.py",
                         "run_all.py", "reference_values.json")
        },
        "all_fourteen_scenario_results": comparison_rows,
        "limitations": [
            "Independent uncalibrated rainfall-only HBV implementation; not official HBV-light.",
            "The 14 labels include the same baseline repeated in four panels.",
            "FAI/FCI use a common 1-20-day fitting range, not an extracted whole-event domain.",
            "Synthetic forcing and states are specified in code; no observations or external data are required.",
            "Installation time and operating systems other than the one recorded above were not tested by this run.",
        ],
    }
    (out / "run_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "scenario_count": report["scenario_count"],
        "figures_generated": report["figures_generated"], "seconds": report["total_seconds"],
        "max_reference_errors": maximum_errors,
        "report": str(out / "run_report.json"),
    }, indent=2))


if __name__ == "__main__":
    main()
