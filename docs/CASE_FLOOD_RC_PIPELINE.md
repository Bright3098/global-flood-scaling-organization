# Case: Flood Events To R/C Outputs

## Goal

Reproduce the upstream workflow before figure plotting:

1. Daily discharge data are converted to `flood_events__*.nc`.
2. Event observations are converted to event-level `R` and `C`.
3. Event-level `scaling_RC__*.nc` files are aggregated at basin outlets.
4. Aggregated R/C outputs can feed the existing RC figure workflow.

## Current Repository Boundary

The current committed, runnable boundary starts at `flood_events__*.nc`.

The local workspace contains flood-event extraction notebooks and generated
event files, but the daily-discharge-to-event extraction has not yet been
promoted into a clean script. For a new environment, treat `flood_events__*.nc`
as the required input for this case.

## Flood Event File Contract

Each input NetCDF should be named:

```text
flood_events__<ghm>__<gcm-or-forcing>__<scenario-or-experiment>__<soc-or-variant>.nc
```

Required dimensions:

- `event`
- `obs`

Required variables:

- `event_id(event)`
- `year(event)`
- `lat(event)`
- `lon(event)`
- `obs_start(event)`
- `obs_count(event)`
- `Vcum(obs)`

Optional but recommended:

- `time_ns(obs)`, used to infer the observation timestep
- global attribute `save_obs=true`
- `obs_event_id(obs)`, useful for diagnostics in the full GPU workflow

## Compute Event-Level R/C

The repository includes a CPU case script:

```powershell
python code\compute_event_rc_from_flood_events.py ^
  --in-dir  D:\workroom\GPLW_furture\case_flooding_event ^
  --out-dir D:\workroom\GPLW_furture\case_flooding_event\RC_outputs ^
  --require-save-obs
```

For each event:

```text
W(T) = max_s[Vcum(s + T) - Vcum(s)]
log10(W) = log10(R) + C * log10(T_hours)
```

Default `T_hours`:

```text
24,48,72,96,120,168,240,336,504,720
```

Expected output:

```text
scaling_RC__<same_tag>.nc
```

Output variables:

- `event_id`
- `year`
- `lat`
- `lon`
- `R`
- `C`
- `n_pts`
- `r2`

## Full-Scale GPU Workflow

The local full-production version is preserved in:

```text
D:\workroom\GPLW_furture\R_c.txt
```

It uses CuPy, Numba, NetCDF4, and multiple GPUs. It implements the same
definition of `R` and `C`, but processes large event files in GPU chunks.

Additional dependencies for the GPU workflow:

- `cupy-cuda11x` or `cupy-cuda12x`, matching the installed CUDA runtime
- working NVIDIA drivers

## Aggregate R/C At Basin Outlets

After event-level `scaling_RC__*.nc` files exist, use the ISIMIP outlet
aggregation scripts in the local workspace:

```powershell
python D:\workroom\GPLW_furture\isimip_outlet_scaling_pipeline_keep_year.py ^
  --root D:\workroom\GPLW_furture\ISIMIP\ISIMIP ^
  --out  D:\workroom\GPLW_furture\ISIMIP_RESULTS ^
  --chunk_events 800000
```

For time-window comparisons:

```powershell
python D:\workroom\GPLW_furture\isimip_outlet_scaling_pipeline_windows.py ^
  --root D:\workroom\GPLW_furture\ISIMIP\ISIMIP ^
  --out  D:\workroom\GPLW_furture\ISIMIP_RESULTS ^
  --base_window 1981 2010 ^
  --future_windows "2036-2065,2071-2100" ^
  --chunk_events 800000
```

The outlet scripts:

- build a manifest of `scaling_RC__*.nc`
- identify one outlet per DDM30 basin from basin mask and flow direction files
- extract and aggregate outlet `R`, `C`, `r2`, and `n_pts`
- compute controlled experiment deltas such as `dR` and `dC`
- save manifest, values, experiment tables, and quick maps

## Observed Station Daily-Series Route

For station-style daily Excel input, the local script
`makeUP/china.py` computes yearly flood scaling summaries:

```powershell
python makeUP\china.py ^
  --xlsx 汇总.xlsx ^
  --sheet Sheet1 ^
  --period_start 1980 ^
  --period_end 2020 ^
  --r2_threshold 0.6 ^
  --dump_scaling outputs\station_year_scaling.csv ^
  --out outputs\station_scaling_summary.csv
```

This route computes annual slope/intercept and trend summaries from daily
series. The figure workflow then treats:

- `R` as the intercept/runoff-generation metric
- `beta` as the slope
- `C = 1 - beta` for confluence trend plotting

## Pass Criteria

This case passes when:

- `flood_events__*.nc` inputs have nonzero `event` and `obs` dimensions.
- `compute_event_rc_from_flood_events.py` writes at least one
  `scaling_RC__*.nc`.
- Output files contain finite `R` and `C` values for valid events.
- Outlet aggregation creates manifest and stage-value tables under the output
  directory.
- Generated NetCDF/CSV/figure outputs remain ignored by Git.
