# Algorithm and version notes

## Source provenance and version definitions

The 21 October 2026 repository update incorporates the core implementation
from the local GPLW checkout and the verified HBV sensitivity experiment.
Original core source hashes are recorded in
[core_source_manifest.json](core_source_manifest.json). Explicit FAI/FCI
variables and their metadata were added to the copied fitting script. The
default complete-event regression and window calculations are retained. The
current update also validates custom duration grids and records event-length
caps explicitly, as described below.

The earlier public version is preserved at commit
[`fecae272d182068244d3298b83f56d6e9e417704`](https://github.com/Bright3098/global-flood-scaling-organization/tree/fecae272d182068244d3298b83f56d6e9e417704).
Its `R` and `C` names have different meanings from the legacy local-code fields
retained in the current implementation:

| Feature | Earlier public commit | Current implementation |
| --- | --- | --- |
| Regression | `log10(W) = intercept + slope * log10(T_hours)` | `ln(W) = R + beta * ln(T / 1 day)` |
| `R` | `10**intercept`, a power-law coefficient | Natural-log intercept at one day |
| `C` | Regression slope | `1 - beta`, the concentration index |
| Explicit indices | No `FAI` or `FCI` output | `FAI = R / ln(10)`; `FCI = C = 1 - beta` |
| Default durations | 24, 48, 72, 96, 120, 168, 240, 336, 504 and 720 hours, subject to event support | Every integer observation-step duration through the complete event; for daily input, 1 through N days |

Changing the logarithm base alone preserves the slope and converts the
intercept by division by `ln(10)` when time and volume references are unchanged.
That statement does not make the two versions' outputs interchangeable:
the earlier public version also used hours and a different default fitting
grid. The current code still accepts an explicit `--t-hours` grid for other
analyses; requested durations must be integer multiples of the observation
step. This avoids fitting a nominal duration against a window of a different
length. The default complete-event domain is used in the included core demo.
If a non-default `--max-event-len` cap truncates an event, output metadata records
that restriction; the capped calculation is not a complete-event fit.

In the core fitting script, `W` is numerical cumulative volume in the input
units. FAI refers to one input volume unit at one day. To use a reference volume
`Wref` expressed in those same units, use `FAI - log10(Wref)`. The provided
example NetCDF lacks a physical volume-unit label, which is recorded explicitly
in its output metadata; absolute intercepts should not be compared across
unmatched units or references.

The HBV source comes from the final single-peak experiment of 14 October 2026.
[Its portability audit](../code/05_hbv_sensitivity/portability_audit.json)
records source comparisons and numerical checks. The update makes output
paths portable and uses built-in plotting fonts; the physical calculations
and scenario parameters are unchanged.

## Core algorithms

1. Extract an annual-maximum event from each eligible station-year of daily
   discharge. Defaults are 1980-2020 and at least 300 valid daily observations
   per year. Identify the peak, scan backward over its non-decreasing rising
   limb, and terminate recession at the first value no greater than 1.05 times
   the onset discharge. The extraction script specifies boundary and
   missing-value rules.
2. For an already extracted event, compute maximum accumulated volume over
   all contiguous integer-duration windows, from one observation step through
   the complete event. Fit log volume against log duration. Export the
   intercept, exponent, fit information and explicit FAI/FCI variables,
   referencing time to one day.
3. For the flood-peak demonstration, aggregate annual station peak discharges,
   normalize by the station baseline and estimate Theil-Sen and Kendall
   trends. This workflow demonstrates peak-flow trends, not all FAI/FCI trend
   analyses from the manuscript.
4. Convert station trend direction to an increase indicator, scan candidate
   log-area splits and maximize between-group separation. The implementation
   also includes a permutation-null option and station subsampling; inclusion
   of an option does not imply that it was used in the manuscript.

## HBV algorithm and its distinct fitting domain

Generate synthetic daily rainfall, integrate soil recharge, update upper and
lower groundwater stores, and route generated runoff with a unit-sum triangular
kernel. For tributary timing, combine equal-area responses with prescribed
onset lags. Compute maximum accumulated runoff depth for each integer duration
from 1 to 20 days, searching windows over the full routed simulation. Fit the
base-10 intercept and slope, with `FAI = intercept`, `FCI = 1 - slope`,
`Qref = 1 mm/day` and `Tref = 1 day`.

The figure displays days 0-20 after rainfall onset. This display window does
not restrict the maximum-window search. HBV uses a fixed 1-20-day fitting
domain; the observational core uses complete extracted events by default.
All computations use daily values; interpolation is for display only.

[HBV_README.md](HBV_README.md) describes equations, flux ordering, scenario
parameters and assumptions. The generated method description is also included
in [the reference outputs](../code/05_hbv_sensitivity/expected_outputs/analysis/methods.txt).

## Validation scope

Verification on 20 October 2026 covered the core tests and small demonstrations,
explicit FAI/FCI exports, HBV conservation and numerical checks, reference-result
agreement, and the four-panel HBV figure. See
[core_verification.json](core_verification.json) and the
[HBV run report](../code/05_hbv_sensitivity/expected_outputs/run_report.json).

These checks do not certify a complete rerun of the manuscript's 3,674 stations,
138,937 events or global ISIMIP experiments. This repository does not contain
the complete global datasets or every manuscript figure pipeline.
