# Software accompanying "Scale-dependent changes in global flood amplitude and concentration"

Repository update: 21 October 2026. Verification records: 20 and 21 October 2026.

This repository contains four core analysis demonstrations and the HBV sensitivity
experiments used in Extended Data Fig. 5. The included examples are small,
independent demonstrations, not the complete global observational and ISIMIP
datasets or a one-command reconstruction of every manuscript result.

## Repository contents

Repository: https://github.com/Bright3098/global-flood-scaling-organization

The current version includes explicit FAI/FCI output variables, the updated
core workflows and the HBV sensitivity module. Definitions and fitting domains
differ from the earlier public commit `fecae272d182068244d3298b83f56d6e9e417704`;
see [algorithm and version notes](docs/ALGORITHMS_AND_VERSION_NOTES.md).

| Directory | Function |
| --- | --- |
| `code/01_flood_extraction/` | Extract annual-maximum flood events from daily streamflow CSV files. |
| `code/02_rc_calculation/` | Fit maximum flood-volume/duration relationships over complete extracted events; export FAI and FCI. |
| `code/03_trend_calculation/` | Demonstrate annual flood-peak trend calculation. |
| `code/04_otsu_like_method/` | Scan basin-area splits for trend-direction separation. |
| `code/05_hbv_sensitivity/` | Run the soil, groundwater and triangular-routing experiments and reproduce the four-panel HBV figure. |
| `examples/case_data/` | Small real-data CSV/NetCDF subsets for the four core demonstrations. |
| `expected_outputs/core/` | Outputs generated during verification of the four core demonstrations. |
| `code/05_hbv_sensitivity/expected_outputs/` | Verified HBV outputs, indices and figure. |
| `docs/` | Software checklist, algorithm notes, version provenance and verification records. |

## Validation on 21 October 2026

The rerun passed nine unit and end-to-end tests, all four core examples and
the full HBV workflow. Core outputs matched the supplied references, and all
14 HBV scenario labels matched their reference results. See the
[validation record](docs/sync_validation_20261021.json). The timings reported
below remain the historical measurements from 20 October 2026.

## System requirements and tested versions

The code uses a CPU. No GPU, non-standard hardware, proprietary hydrological
software or observation downloads are needed for the included examples.
A normal desktop or laptop is sufficient for these small demonstrations.
Requirements for a full global production analysis depend on input size.

On 20 October 2026, two existing environments were tested on Windows 11, build 26200:

| Component | Core workflows and numerical HBV | Complete HBV, including figure |
| --- | --- | --- |
| Python | 3.13.5 | 3.12.14 |
| NumPy | 2.1.3 | 2.3.5 |
| pandas | 2.2.3 | not required |
| SciPy | 1.15.3 | not required |
| netCDF4 | 1.7.2 | not required |
| openpyxl | 3.1.5 | not required |
| ReportLab | not required | 4.4.9 |
| pypdfium2 | not required | 5.14.0 |
| Pillow | not required | 12.3.0 |

The full dependency specifications are in `requirements-core.txt` and
`code/05_hbv_sensitivity/requirements.txt`. The standard-library modules
require no separate installation. Linux and macOS have not been tested.
Use separate environments to reproduce the tested dependency combinations.

## Installation

Install Python first, then create a virtual environment. On Windows, run:

```text
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements-core.txt
```

For the full HBV module, create a separate environment using Python 3.12 and run:

```text
python -m venv .venv-hbv
.venv-hbv\Scripts\activate
python -m pip install -r code/05_hbv_sensitivity/requirements.txt
```

Allow approximately 5-15 minutes for dependency installation on a typical
desktop with a working internet connection. This is a planning estimate;
installation was not timed. The verification below used existing environments.
After installation, the demonstrations run offline.

## Core demonstrations

Run from the repository root:

```text
python examples/example_four_methods.py
python examples/example_four_methods.py --execute
python -m unittest discover -s code/tests -p "test_*.py" -v
```

The first command prints the four commands. The second executes them; the third
runs the existing core-definition tests plus the additional output/window tests. The four demonstrations use
separate input subsets. They are not one continuous global-analysis pipeline.

Outputs are placed in `outputs/case_data_run/`:

```text
01_flood_events_csv/ALL_flood_events_summary.csv
01_flood_events_csv/ALL_flood_events_timeseries.csv
02_event_rc/scaling_RC__barcoo_blackall_case.nc
03_qpeak_trends/qpeak_trend_summary.csv
04_otsu_like_threshold/otsu_like_threshold_summary.csv
```

Reference outputs are included under `expected_outputs/core/`. The four core
demonstrations completed in 5.69 seconds in the tested environment; the five
core-definition tests present on 20 October took 1.20 seconds and passed. The example NetCDF contains four events.
Exported FAI/FCI values were also checked against the fitted intercept/exponent.
See [the core verification record](docs/core_verification.json). These are
measured single-run times from 20 October 2026; runtime varies by computer.

## HBV demonstration

In the HBV environment:

```text
python code/05_hbv_sensitivity/run_all.py
```

For numerical outputs without figure generation, only NumPy is required:

```text
python code/05_hbv_sensitivity/run_all.py --skip-figures
```

The code generates its own rainfall forcings. Four comparisons vary rainfall
amount/duration, initial soil water, triangular-routing duration, and tributary
timing. It retains 14 plotted scenario labels, including repeated baselines.
It simulates 730 daily cells, fits FAI/FCI using every integer duration from
1 to 20 days over the full routed simulation, and displays days 0-20.
The 0-, 2- and 3-day tributary lags retain single-peaked responses.

Outputs include full daily hydrographs and states, parameter tables,
FAI/FCI values, maximum-window volumes, conservation checks and PDF/SVG/PNG
figures. [The HBV documentation](docs/HBV_README.md) provides equations, exact
settings, expected results, measured runtimes and reference comparisons.
The numerical-only run took 0.659 seconds; the complete HBV workflow, including
high-resolution figure generation, took 1.948 seconds in the second tested
environment on 20 October 2026. Allow up to one minute on other desktop systems.
The [complete-run report](code/05_hbv_sensitivity/expected_outputs/run_report.json)
records its measured times and software versions.

This is an independent, uncalibrated implementation of the stated HBV routines,
not the official HBV-light application. Snow is inactive and potential
evapotranspiration is zero. It demonstrates process sensitivity and does not
attribute observed global flood changes.

## Use with other data

```text
python code/01_flood_extraction/extract_flood_events.py --input-root YOUR_DAILY_CSV_FOLDER --out-dir outputs/my_events --workers 4
python code/02_rc_calculation/compute_event_rc_from_flood_events.py --in-dir YOUR_FLOOD_EVENT_NETCDF_FOLDER --out-dir outputs/my_indices
python code/03_trend_calculation/compute_qpeak_trends.py --events-csv YOUR_EVENT_SUMMARY.csv --out-dir outputs/my_peak_trends
python code/04_otsu_like_method/otsu_like_area_threshold.py --input YOUR_AREA_TREND_TABLE.csv --out-dir outputs/my_threshold --area-col area_km2 --trend-col trend_direction
```

Inspect each script's `--help` and the included example input files before
substituting new data. The fitting input is a ragged event NetCDF with
`event_id`, `year`, `lat`, `lon`, `obs_start`, `obs_count` and `Vcum`, optionally
`time_ns`. `Vcum` is cumulative event volume without a leading zero.
Keep its physical units consistent and record them in the input metadata.
Extraction produces CSV; the fitting demonstration uses its separately
provided NetCDF input. No automatic CSV-to-NetCDF bridge is claimed here.

For daily inputs, the current fitting code uses every integer duration from
one day to the complete event by default. It retains the legacy local-code
fields `R` (natural-log intercept) and `C` (concentration), and exports:

```text
ln W(T) = R + beta * ln(T / 1 day)
FAI = R / ln(10)
FCI = C = 1 - beta
```

Here the numeric FAI uses one input cumulative-volume unit as reference volume
at one day. If a different reference volume is required, subtract its base-10
logarithm in the same units. The HBV module uses runoff depth in mm and
Qref = 1 mm/day, Tref = 1 day. Do not interpret absolute intercepts across
different unit/reference choices as physically identical.

The older public version used `R = 10**intercept`, `C = slope` and selected
hour-based durations. Those fields are not interchangeable with the current
`R` and `C`; see [the version comparison](docs/ALGORITHMS_AND_VERSION_NOTES.md).
The core complete-event domain also differs from HBV's fixed 1-20-day domain.

To explore other HBV settings, edit the parameter/scenario declarations in
`code/05_hbv_sensitivity/hbv_sensitivity.py`. Such a change creates a different
experiment and will intentionally fail checks against the supplied fixed
reference values.

## Scope, manuscript locations and license

Algorithm descriptions appear in the manuscript Methods and Supplementary
Notes 1 and 3-6; HBV equations/settings are in Supplementary Note 6, with
results in Extended Data Fig. 5. This repository reproduces the HBV figure
and provides core demonstrations, but does not contain the complete global
datasets or every model-processing and figure pipeline. The optional checklist
item for reproducing all quantitative manuscript results is therefore not met.
See [the software checklist](docs/NATURE_CODE_CHECKLIST.md).

[MIT License](LICENSE) covers the code. Data-provider terms remain applicable
to original third-party observational datasets. HBV forcing and hydrographs
are generated internally. Independent-colleague installation testing and
archival DOI registration have not been completed.
