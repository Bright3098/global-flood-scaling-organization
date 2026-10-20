from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


CODE = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = CODE / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


rc = load_module("rc", "02_rc_calculation/compute_event_rc_from_flood_events.py")
events = load_module("events", "01_flood_extraction/extract_flood_events.py")
otsu = load_module("otsu", "04_otsu_like_method/otsu_like_area_threshold.py")


class CanonicalDefinitionTests(unittest.TestCase):
    def test_rc_parameterization(self) -> None:
        t_days = np.arange(1.0, 41.0)
        w = 7.0 * t_days**0.6
        r_log, beta, c_value, n_points, r_squared = rc.fit_rc(w, t_days)
        self.assertAlmostEqual(r_log, np.log(7.0), places=12)
        self.assertAlmostEqual(beta, 0.6, places=12)
        self.assertAlmostEqual(c_value, 0.4, places=12)
        self.assertEqual(n_points, 40)
        self.assertAlmostEqual(r_squared, 1.0, places=12)

    def test_amplitude_invariance_of_c(self) -> None:
        t_days = np.arange(1.0, 31.0)
        w = 5.0 * t_days**0.72
        base = rc.fit_rc(w, t_days)
        scaled = rc.fit_rc(11.0 * w, t_days)
        self.assertAlmostEqual(scaled[0] - base[0], np.log(11.0), places=12)
        self.assertAlmostEqual(scaled[1], base[1], places=12)
        self.assertAlmostEqual(scaled[2], base[2], places=12)

    def test_full_event_window_domain(self) -> None:
        flow = np.asarray([1.0, 3.0, 8.0, 4.0, 2.0])
        vcum = np.cumsum(flow)
        steps = np.arange(1, len(flow) + 1)
        maxima = rc.event_window_maxima(vcum, steps)
        np.testing.assert_allclose(maxima, [8.0, 12.0, 15.0, 17.0, 18.0])

    def test_annual_coverage_and_event_boundary(self) -> None:
        dates = pd.date_range("1980-01-01", "1980-12-31", freq="D")
        flow = np.ones(len(dates), dtype=float)
        flow[148] = 1.2
        flow[150:156] = [2, 3, 5, 9, 4, 1.04]
        frame = pd.DataFrame({"time": dates, "runoff": flow, "station_name": "test"})
        extracted = events.extract_flood_events_timeseries(frame)
        self.assertEqual(len(extracted), 1)
        self.assertEqual(extracted[0]["valid_days"], 366)
        self.assertEqual(extracted[0]["rising_time"], dates[149])
        self.assertEqual(extracted[0]["recession_time"], dates[155])

    def test_permutation_null_is_deterministic(self) -> None:
        frame = pd.DataFrame(
            {
                "area_km2": np.geomspace(10, 100_000, 120),
                "is_increase": np.r_[np.ones(60), np.zeros(60)],
            }
        )
        frame["log10_area"] = np.log10(frame["area_km2"])
        scan = otsu.scan_thresholds(frame, q_min=0.1, q_max=0.9, q_count=51, min_side_n=25)
        _, first = otsu.permutation_null(frame, scan, n_perm=200, seed=3098, min_side_n=25)
        _, second = otsu.permutation_null(frame, scan, n_perm=200, seed=3098, min_side_n=25)
        pd.testing.assert_frame_equal(first, second)


if __name__ == "__main__":
    unittest.main()
