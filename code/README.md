# Analysis code

The five workflow directories accompany *Scale-dependent changes in global
flood amplitude and concentration*. Installation commands, tested dependency
versions and measured runtimes are in the [root README](../README.md).

| Directory | Purpose and inputs |
| --- | --- |
| `01_flood_extraction/` | Read daily streamflow CSVs and extract annual-maximum events; write event summary and time-series CSVs. |
| `02_rc_calculation/` | Read ragged `flood_events__*.nc` event files and fit maximum accumulated volume against duration; export FAI, FCI and fit metadata. |
| `03_trend_calculation/` | Aggregate annual peak discharges, normalize by the station baseline and estimate Theil-Sen/Kendall trends. |
| `04_otsu_like_method/` | Scan log-area splits to maximize separation of station trend directions. |
| `05_hbv_sensitivity/` | Run rainfall, soil-water, triangular-routing and tributary-timing experiments and reproduce Extended Data Fig. 5. |
| `tests/` | Core-definition tests plus additional output/window tests. |

## FAI and FCI definitions

For daily observations, the core fitting script uses every integer window
duration from one day through the complete extracted event by default:

```text
ln W(T) = R + beta * ln(T / 1 day)
FAI = R / ln(10)
FCI = C = 1 - beta
```

Here `W` is expressed numerically in the input cumulative-volume units. The
FAI reference volume is one such unit at a reference duration of one day.
Changing volume units or the reference volume changes the intercept; record
the input units before interpreting it. The included event NetCDF does not
specify its physical volume unit, so its FAI demonstrates the calculation
rather than an absolute cross-dataset comparison.

The current `R` is a natural-log intercept and `C` is concentration. The earlier
public commit `fecae272d182068244d3298b83f56d6e9e417704` instead exported
`R = 10**intercept` and `C = slope` from a base-10 fit using hours. These are
different definitions, not only different names. See
[algorithm and version notes](../docs/ALGORITHMS_AND_VERSION_NOTES.md).

HBV uses runoff depth in mm with `Qref = 1 mm/day` and `Tref = 1 day`. Its
FAI/FCI fit uses the common integer-duration grid 1-20 days, searching all
windows in the full routed simulation. The displayed days 0-20 do not limit
the window search. See [HBV methods and verification](../docs/HBV_README.md).

## Demonstrations

From the repository root:

```text
python examples/example_four_methods.py --execute
python -m unittest discover -s code/tests -p "test_*.py" -v
python code/05_hbv_sensitivity/run_all.py
```

Use the separate core and HBV environments described in the root README.
The four core demonstrations have separate [example inputs](../examples/case_data/);
extraction CSVs are not automatically converted to the fitting NetCDF format.
Reference outputs are under [expected_outputs/core](../expected_outputs/core/)
and [05_hbv_sensitivity/expected_outputs](05_hbv_sensitivity/expected_outputs/).
These demonstrations and the HBV figure do not reconstruct every global
observational or ISIMIP result in the manuscript.

The 21 October 2026 rerun passed nine unit and end-to-end tests, four core
examples and the full HBV workflow. All 14 HBV scenario labels matched their
reference results; see the [validation record](../docs/sync_validation_20261021.json).
