# HBV sensitivity experiments

This module reproduces the controlled hydrographs and fitted indices in Extended Data Fig. 5 of *Scale-dependent changes in global flood amplitude and concentration*. It uses the final single-peak experiment dated 14 October 2026. All forcing is generated within the code; no observational data or external inputs are required.

## Software and execution

The independent Python implementation is in `code/05_hbv_sensitivity/`. It contains:

- `hbv_sensitivity.py`: daily soil, groundwater and triangular-routing simulation and numerical checks.
- `fit_indices.py`: maximum-window volumes and FAI/FCI fitting.
- `plot_hbv.py`: the four-panel figure and interpolation checks.
- `run_all.py`: the reproducible entry point, reference-result checks and measured runtime report.
- `reference_values.json`: the 14 labelled outputs expected from the final experiment.
- `requirements.txt`: dependency versions available in the bundled verification environment.

Simulation and fitting require only NumPy. Figure generation additionally uses ReportLab, pypdfium2 and Pillow. No GPU, system fonts, parent-folder modules or network connection are used during execution. Built-in Helvetica replaces the original Windows Arial font, changing typography only.

From the repository root:

```text
python -m pip install -r code/05_hbv_sensitivity/requirements.txt
python code/05_hbv_sensitivity/run_all.py
```

For numerical outputs only, in an environment with NumPy:

```text
python code/05_hbv_sensitivity/run_all.py --skip-figures
```

To choose the output location:

```text
python code/05_hbv_sensitivity/run_all.py --output-dir outputs/hbv
```

The default output directory is `code/05_hbv_sensitivity/outputs/`. Existing generated files with matching names are replaced. Python must run without `-O`, so that assertion-based validation remains enabled. Parameter values and the fixed scenario list are in `hbv_sensitivity.py`; modifying them creates a different experiment and is expected to fail the supplied reference-result checks.

On 20 October 2026, the complete workflow was verified on Windows 11, build 26200, with Python 3.12.14, NumPy 2.3.5, ReportLab 4.4.9, pypdfium2 5.14.0 and Pillow 12.3.0. A separate numerical-only run was verified on the same Windows machine with Python 3.13.5 and NumPy 2.1.3; plotting was not tested in that environment. The [complete-run report](../code/05_hbv_sensitivity/expected_outputs/run_report.json) records its versions and measured times. Installation time was not measured, and Linux and macOS were not tested. The [root README](../README.md) provides separate environment instructions and a labelled installation-time estimate.

## Model and scenarios

This is an independent, uncalibrated implementation of the HBV soil, two-store groundwater-response and triangular-routing routines described by Seibert (2005, pp. 4–5), not a run of the official HBV-light software. Snow processes are inactive and potential evapotranspiration is zero.

The model runs 730 one-day cells. The first rainfall begins at zero-based cell index 5, defining time zero. Soil capacity is 200 mm, BETA is 2, LP is 0.7, and both groundwater stores start empty. The baseline initial soil fraction is 0.5. With zero evapotranspiration, LP has no effect.

The soil relation `dSM/dP = 1 − (SM/FC)^2` is integrated analytically over each daily rainfall depth. For `r = SMstart/FC` and `u = tanh(P/FC)`, the updated store is `FC(r+u)/(1+ru)`; recharge equals rainfall minus the increase in soil storage. Recharge enters the upper store, after which `min(1 mm, SUZ)` percolates to the lower store. The same post-percolation upper-store state is used to calculate `Q0 = 0.3 max(SUZ−10 mm, 0)` and `Q1 = 0.1 SUZ`; they are removed together. `Q2 = 0.02 SLZ` is calculated after percolation and removed from the lower store. Rate coefficients have units day−1.

Generated runoff is routed by full convolution with non-negative unit-sum weights. Each weight is the integral of a unit-area triangular density over a successive one-day lag bin. The triangle has support `[0, MAXBAS]` and apex `MAXBAS/2`; its bin-j weight is applied at daily lag j. MAXBAS denotes support width, not mean lag. All convolution tails and residual groundwater are included in water-balance checks.

The baseline is 120 mm of rain over six days, initial soil fraction 0.5 and MAXBAS three days. Rainfall weights are proportional to `min(j+1, D−j)` for `j=0,…,D−1`. Four comparisons retain the same baseline:

- Rainfall totals of 60/120/180 mm over six days, or durations of 3/6/12 days at 120 mm.
- Initial soil fractions of 0.2/0.5/0.8.
- MAXBAS widths of 1/3/7 days.
- Two identical, independently simulated, equal-area tributaries with rainfall-onset lags of 0/2/3 days, combined with weights 0.5/0.5.

There are 14 plotted labels: five in panel a and three in each other panel. The identical baseline is deliberately repeated in all four panels; these are not 14 independent experimental replicates. The selected lag responses are checked to be non-decreasing before their daily maximum and non-increasing thereafter, with a tolerance of 10−10 mm/day.

## Index calculation and figure

For every integer duration `T=1,…,20` days, the maximum accumulated runoff depth `W(T)` is searched over every contiguous T-day window in the complete routed simulation. Daily means are multiplied by the one-day bin width. Unweighted ordinary least squares fits

```text
log10[W(T)/(Qref Tref)] = FAI + (1−FCI) log10(T/Tref),
Qref = 1 mm/day, Tref = 1 day.
```

FAI is a fitted intercept, not an imposed peak. No event-end threshold, quality exclusion, discharge normalization or peak alignment is used. All 14 labels and all 20 durations are retained. These indices describe the common 1–20-day fitting domain; they are not complete extracted-event estimates.

Only days 0–20 after first rainfall onset are displayed. Shape-preserving cubic interpolation connects daily-bin centres for presentation. All simulation, volume, index and single-peak checks use the original daily values.

## Outputs and verification

The `analysis/` subdirectory contains:

- `hydrographs.csv`: every routed daily discharge and rainfall value.
- `states_and_fluxes.csv`: daily soil/groundwater states and water-balance components.
- `scenario_metrics.csv` and `routing_weights.csv`: scenario settings, peaks, totals and kernels.
- `fitted_indices.csv`: all 14 FAI/FCI pairs, R² and fitting metadata.
- `maximum_window_volumes.csv`: all 280 maximizing windows and regression residuals.
- `validation.json`: water balance, routing, soil integration, identical-baseline and single-peak checks.
- `index_fit_validation.json`: full-series window searches, independent summation agreement, unchanged source data and identical baseline indices.
- `methods.txt`: detailed numerical method and parameter documentation.

The `figures/` subdirectory contains PDF, SVG, 600-dpi PNG, preview PNG and `figure_data_audit.json`. The top-level output folder contains `run.log`, `run_report.json` and `scenario_verification.csv`, including all 14 index pairs and balance errors.

Model checks include non-negative states and discharge, daily and whole-record mass conservation, conservative routing, equal baseline responses, identical available recharge for equal-total rainfall, independent Runge–Kutta verification of the analytical soil update, and independent triangular-kernel quadrature. Index checks compare direct window sums with cumulative-sum differences. Figure checks cover interpolation knots, first-derivative continuity, interval overshoot and preservation of daily peak bounds.

Reference CSVs, the four-panel figure in all exported formats and the complete run report are included under [code/05_hbv_sensitivity/expected_outputs/](../code/05_hbv_sensitivity/expected_outputs/). The historical numerical-only verification is summarized below; duplicate numerical-only output files are not included.

The complete workflow took 1.948 seconds in the Python 3.12.14 environment: 0.294 seconds for simulation and model checks, 0.055 seconds for fitting, and 1.534 seconds for figure generation and interpolation checks. These are elapsed times measured inside the runner; interpreter startup and installation are excluded. All 14 FAI/FCI pairs, R² values, peaks, peak times and volumes matched the original final-version results exactly, and the complete hydrograph CSV was byte-identical. The generated figure was visually checked for labels, curves, panel layout and clipping.

The numerical-only verification took 0.659 seconds in the Python 3.13.5 environment. Its maximum daily and whole-record water-balance errors were 2.85×10−14 mm or smaller. Peak magnitude, peak time and total volume matched the original final-version outputs. The maximum FAI and FCI differences were approximately 1.11×10−15 and 1.00×10−15, respectively. All 10,250 daily output rows were separately compared against the original final-version hydrographs, with maximum discharge difference below 1.4×10−23 mm/day and no rainfall difference. These are measured local runs; runtime varies by computer. See the [portability audit](../code/05_hbv_sensitivity/portability_audit.json) and [complete-run report](../code/05_hbv_sensitivity/expected_outputs/run_report.json) for the published verification records.

These checks establish numerical and reference-result consistency for this uncalibrated sensitivity experiment. They do not establish observational validation, independence of FAI and FCI, or universal responses to the parameters.

Reference: Seibert, J. (2005), *HBV light version 2 User’s Manual*, pp. 4–5. https://www.geo.uzh.ch/dam/jcr:c8afa73c-ac90-478e-a8c7-929eed7b1b62/hbv_manual_2005.pdf
