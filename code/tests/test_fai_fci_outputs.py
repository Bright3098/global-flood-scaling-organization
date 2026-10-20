"""End-to-end checks of the canonical indices written to NetCDF.

Each synthetic event has daily volumes A * (d**beta - (d-1)**beta).
For 0 < beta < 1 these decrease with d, so the largest T-day volume is
exactly A * T**beta. The expected FAI and FCI are therefore known without
using the fitting implementation: log10(A) and 1-beta.
"""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import netCDF4 as nc
import numpy as np


CODE = Path(__file__).resolve().parents[1]
SCRIPT = CODE / "02_rc_calculation" / "compute_event_rc_from_flood_events.py"
SPEC = importlib.util.spec_from_file_location("rc_output_tests", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import {SCRIPT}")
rc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rc)


class FaiFciOutputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.workspace = tempfile.TemporaryDirectory(prefix="fai_fci_outputs_")
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)
        self.input_path = self.root / "flood_events__known_power_laws.nc"
        self.amplitudes = np.array([7.0, 77.0, 0.125])
        self.betas = np.array([0.6, 0.6, 0.2])
        self.counts = np.array([31, 37, 23], dtype=np.int32)
        self.event_ids = np.array([101, 203, 305], dtype=np.int64)
        starts = np.r_[0, np.cumsum(self.counts)[:-1]]
        day_ns = np.int64(86_400_000_000_000)

        with nc.Dataset(self.input_path, "w") as dataset:
            dataset.createDimension("event", len(self.counts))
            dataset.createDimension("obs", int(self.counts.sum()))
            dataset.save_obs = "true"
            for name, dtype, values in (
                ("event_id", "i8", self.event_ids),
                ("year", "i4", [1980, 1981, 1982]),
                ("lat", "f4", [10.0, 20.0, 30.0]),
                ("lon", "f4", [100.0, 110.0, 120.0]),
                ("obs_start", "i8", starts),
                ("obs_count", "i4", self.counts),
            ):
                dataset.createVariable(name, dtype, ("event",))[:] = values
            volume = dataset.createVariable("Vcum", "f8", ("obs",))
            volume.units = "mm"
            times = dataset.createVariable("time_ns", "i8", ("obs",))
            times.units = "nanoseconds since 1970-01-01 00:00:00"
            for i, (start, count, amplitude, beta) in enumerate(
                zip(starts, self.counts, self.amplitudes, self.betas)
            ):
                sl = slice(int(start), int(start + count))
                days = np.arange(1, count + 1, dtype=np.float64)
                # Cumulative volume resets at the start of every event and
                # has no leading zero, as required by the input contract.
                volume[sl] = amplitude * days**beta
                times[sl] = (np.arange(count, dtype=np.int64) + 1000 * i) * day_ns

    def process(self, name: str, t_hours=None, cap: int = 0) -> Path:
        output_path = self.root / name
        rc.process_file(
            in_path=self.input_path,
            out_path=output_path,
            t_hours_input=t_hours,
            event_batch=2,
            max_event_len=cap,
            compression=0,
            overwrite=False,
            require_save_obs=True,
        )
        return output_path

    def assert_indices(self, dataset) -> None:
        np.testing.assert_allclose(
            dataset.variables["FAI"][:], np.log10(self.amplitudes), atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            dataset.variables["FCI"][:], 1 - self.betas, atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            dataset.variables["R"][:], np.log(self.amplitudes), atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            dataset.variables["beta"][:], self.betas, atol=1e-12, rtol=0
        )
        np.testing.assert_allclose(
            dataset.variables["R_amplitude"][:], self.amplitudes, atol=1e-12, rtol=0
        )
        np.testing.assert_array_equal(dataset.variables["C"][:], dataset.variables["FCI"][:])
        np.testing.assert_allclose(dataset.variables["r2"][:], 1.0, atol=1e-7, rtol=0)
        np.testing.assert_array_equal(dataset.variables["event_id"][:], self.event_ids)
        self.assertEqual(dataset.complete, "true")

    def test_full_event_netcdf_exports_known_fai_fci(self) -> None:
        with nc.Dataset(self.process("full.nc")) as dataset:
            self.assert_indices(dataset)
            np.testing.assert_array_equal(dataset.variables["n_pts"][:], self.counts)
            np.testing.assert_array_equal(dataset.variables["Tmax_days"][:], self.counts)
            self.assertEqual(dataset.dt_hours, 24.0)
            self.assertEqual(dataset.reference_duration_days, 1.0)
            self.assertEqual(dataset.reference_volume, "one input cumulative-volume unit")
            self.assertEqual(dataset.input_volume_units, "mm")
            self.assertEqual(dataset.max_event_len_steps, 0)
            self.assertIn("complete extracted event", dataset.fit_domain)
            self.assertAlmostEqual(
                float(dataset.variables["FAI"][1] - dataset.variables["FAI"][0]),
                np.log10(11.0),
                places=12,
            )

    def test_fixed_windows_use_actual_integer_sample_durations(self) -> None:
        # The small representational perturbation is accepted, but both the
        # regression and output metadata must use the actual sample duration.
        requested = np.array([24.0 + 1e-10, 48.0, 48.0, 96.0])
        with nc.Dataset(self.process("fixed.nc", t_hours=requested)) as dataset:
            self.assert_indices(dataset)
            np.testing.assert_array_equal(dataset.variables["n_pts"][:], [3, 3, 3])
            np.testing.assert_array_equal(dataset.variables["Tmax_days"][:], [4, 4, 4])
            self.assertEqual(dataset.T_hours_used, "24,48,96")
            self.assertEqual(dataset.fit_domain, "explicit fixed-window sensitivity")

    def test_noninteger_sample_duration_is_rejected_before_output(self) -> None:
        with self.assertRaisesRegex(ValueError, r"integer multiples.*24 h.*36 h"):
            self.process("invalid.nc", t_hours=np.array([24.0, 36.0, 72.0]))
        self.assertFalse((self.root / "invalid.nc").exists())
        self.assertFalse((self.root / "invalid.nc.tmp").exists())

    def test_capped_domain_is_explicit_in_netcdf(self) -> None:
        for name, requested, expected_points, expected_tmax in (
            ("capped_full.nc", None, 10, 10),
            ("capped_fixed.nc", np.array([24.0, 72.0, 480.0]), 2, 3),
        ):
            with self.subTest(name=name):
                with nc.Dataset(self.process(name, t_hours=requested, cap=10)) as dataset:
                    self.assert_indices(dataset)
                    self.assertEqual(dataset.max_event_len_steps, 10)
                    self.assertIn("first 10 observation steps", dataset.fit_domain)
                    self.assertIn("event-length cap", dataset.fit_domain)
                    self.assertNotIn("complete extracted event", dataset.fit_domain)
                    np.testing.assert_array_equal(
                        dataset.variables["n_pts"][:], [expected_points] * 3
                    )
                    np.testing.assert_array_equal(
                        dataset.variables["Tmax_days"][:], [expected_tmax] * 3
                    )


if __name__ == "__main__":
    unittest.main()
