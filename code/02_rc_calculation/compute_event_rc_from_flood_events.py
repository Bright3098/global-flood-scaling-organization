#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compute event-level flood scaling parameters R and C.

Input files are `flood_events__*.nc` NetCDFs with an `event` dimension, an
`obs` dimension, and observed event traces saved as cumulative volume (`Vcum`).

For each flood event, this script computes:

    W(T) = max_s[Vcum(s + T) - Vcum(s)]
    log10(W) = log10(R) + C * log10(T_hours)

The output file mirrors the input name:

    flood_events__tag.nc -> scaling_RC__tag.nc

This CPU implementation is intended for reproducible cases and small or medium
test batches. The original full-production GPU version used in the local
workspace is documented in `docs/CASE_FLOOD_RC_PIPELINE.md`.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable

import numpy as np


DEFAULT_T_HOURS = "24,48,72,96,120,168,240,336,504,720"
FILL_I8 = np.int64(-9223372036854775806)


def parse_t_hours(text: str) -> np.ndarray:
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError("--t-hours must contain at least two positive values")
    arr = np.asarray(vals, dtype=np.float64)
    if np.any(~np.isfinite(arr)) or np.any(arr <= 0):
        raise ValueError("--t-hours values must be finite and positive")
    return arr


def iter_input_files(in_dir: Path) -> list[Path]:
    return sorted(in_dir.glob("**/flood_events__*.nc"))


def output_path_for(in_path: Path, in_dir: Path, out_dir: Path, keep_structure: bool) -> Path:
    rel = in_path.relative_to(in_dir) if keep_structure else Path(in_path.name)
    name = rel.name
    if name.startswith("flood_events__"):
        name = "scaling_RC__" + name[len("flood_events__"):]
    else:
        name = "scaling_RC__" + name
    return out_dir / rel.parent / name


def require_vars(fin, names: Iterable[str]) -> None:
    missing = [name for name in names if name not in fin.variables]
    if missing:
        raise KeyError(f"Missing required variables: {', '.join(missing)}")


def infer_dt_hours(fin) -> float:
    if "time_ns" not in fin.variables:
        return 24.0

    starts = fin.variables["obs_start"][:]
    counts = fin.variables["obs_count"][:]
    valid = np.where(counts >= 2)[0]
    if valid.size == 0:
        return 24.0

    idx = int(valid[0])
    start = int(starts[idx])
    count = int(min(counts[idx], 200))
    t = fin.variables["time_ns"][start:start + count].astype(np.int64)
    dt = np.diff(t).astype(np.float64) / 1e9 / 3600.0
    dt = dt[np.isfinite(dt) & (dt > 0)]
    if dt.size == 0:
        return 24.0
    return float(np.median(dt))


def unique_t_steps(t_hours: np.ndarray, dt_hours: float) -> tuple[np.ndarray, np.ndarray]:
    steps = np.rint(t_hours / dt_hours).astype(np.int32)
    steps[steps < 1] = 1

    keep_hours: list[float] = []
    keep_steps: list[int] = []
    seen: set[int] = set()
    for hour, step in zip(t_hours, steps):
        s = int(step)
        if s in seen:
            continue
        seen.add(s)
        keep_hours.append(float(hour))
        keep_steps.append(s)
    return np.asarray(keep_hours, dtype=np.float64), np.asarray(keep_steps, dtype=np.int32)


def event_window_maxima(vcum: np.ndarray, t_steps: np.ndarray) -> np.ndarray:
    n = int(vcum.size)
    out = np.full(t_steps.size, np.nan, dtype=np.float64)
    if n < 2:
        return out

    v = np.asarray(vcum, dtype=np.float64)
    if np.any(~np.isfinite(v)):
        return out

    # Input Vcum has no leading zero; prepend one so window sums are differences.
    v0 = np.empty(n + 1, dtype=np.float64)
    v0[0] = 0.0
    v0[1:] = v

    for i, step in enumerate(t_steps):
        L = int(step)
        if L <= 0 or L > n:
            continue
        diff = v0[L:] - v0[:-L]
        if diff.size:
            m = np.nanmax(diff)
            if np.isfinite(m) and m > 0:
                out[i] = float(m)
    return out


def fit_rc(w: np.ndarray, t_hours: np.ndarray) -> tuple[float, float, int, float]:
    ok = np.isfinite(w) & (w > 0)
    n = int(ok.sum())
    if n < 2:
        return np.nan, np.nan, n, np.nan

    x = np.log10(t_hours[ok])
    y = np.log10(w[ok])
    slope, intercept = np.polyfit(x, y, 1)
    yhat = intercept + slope * x
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_res = float(np.sum((y - yhat) ** 2))
    r2 = np.nan if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    R = 10.0 ** intercept
    C = slope
    return float(R), float(C), n, float(r2)


def create_output_file(nc, out_path: Path, ne: int, compression: int):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    zopts = {"zlib": compression > 0, "complevel": compression} if compression > 0 else {}
    fout = nc.Dataset(out_path, "w", format="NETCDF4")
    fout.createDimension("event", ne)
    for name, dtype in [("event_id", "i8"), ("year", "i4"), ("lat", "f4"), ("lon", "f4")]:
        fout.createVariable(name, dtype, ("event",), **zopts)
    fout.createVariable("R", "f4", ("event",), fill_value=np.float32(np.nan), **zopts)
    fout.createVariable("C", "f4", ("event",), fill_value=np.float32(np.nan), **zopts)
    fout.createVariable("n_pts", "i2", ("event",), **zopts)
    fout.createVariable("r2", "f4", ("event",), fill_value=np.float32(np.nan), **zopts)
    return fout


def process_file(
    in_path: Path,
    out_path: Path,
    t_hours_input: np.ndarray,
    event_batch: int,
    max_event_len: int,
    compression: int,
    overwrite: bool,
    require_save_obs: bool,
) -> str:
    import netCDF4 as nc

    if out_path.exists() and not overwrite:
        return f"[SKIP] exists: {out_path.name}"

    with nc.Dataset(in_path, "r") as fin:
        if "event" not in fin.dimensions:
            raise KeyError(f"{in_path.name}: missing event dimension")
        if "obs" not in fin.dimensions:
            raise KeyError(f"{in_path.name}: missing obs dimension")
        if require_save_obs and str(getattr(fin, "save_obs", "false")).lower() != "true":
            raise ValueError(f"{in_path.name}: save_obs is not true")

        require_vars(fin, ["event_id", "year", "lat", "lon", "obs_start", "obs_count", "Vcum"])
        ne = int(fin.dimensions["event"].size)
        if ne <= 0:
            return f"[SKIP] no events: {in_path.name}"

        dt_hours = infer_dt_hours(fin)
        t_hours, t_steps = unique_t_steps(t_hours_input, dt_hours)
        if t_steps.size < 2:
            raise ValueError(f"{in_path.name}: fewer than two unique T steps after dt conversion")

        tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
        if tmp_path.exists():
            tmp_path.unlink()

        fout = create_output_file(nc, tmp_path, ne, compression)
        try:
            fout.title = "Event-level flood volume scaling parameters"
            fout.source_file = str(in_path)
            fout.dt_hours = float(dt_hours)
            fout.T_hours_input = ",".join(f"{x:g}" for x in t_hours_input)
            fout.T_hours_used = ",".join(f"{x:g}" for x in t_hours)
            fout.method = "log10(W)=log10(R)+C*log10(T_hours)"
            fout.complete = "false"

            v_event_id = fin.variables["event_id"]
            v_year = fin.variables["year"]
            v_lat = fin.variables["lat"]
            v_lon = fin.variables["lon"]
            v_obs_start = fin.variables["obs_start"]
            v_obs_count = fin.variables["obs_count"]
            v_vcum = fin.variables["Vcum"]

            for e0 in range(0, ne, event_batch):
                e1 = min(ne, e0 + event_batch)
                sl = slice(e0, e1)

                event_id = v_event_id[sl].astype(np.int64)
                year = v_year[sl].astype(np.int32)
                lat = v_lat[sl].astype(np.float32)
                lon = v_lon[sl].astype(np.float32)
                starts = v_obs_start[sl].astype(np.int64)
                counts = v_obs_count[sl].astype(np.int32)

                R = np.full(e1 - e0, np.nan, dtype=np.float32)
                C = np.full(e1 - e0, np.nan, dtype=np.float32)
                npts = np.zeros(e1 - e0, dtype=np.int16)
                r2 = np.full(e1 - e0, np.nan, dtype=np.float32)

                valid = counts > 1
                if np.any(valid):
                    capped = np.minimum(counts, max_event_len).astype(np.int32)
                    ob0 = int(starts[valid].min())
                    ob1 = int((starts[valid] + capped[valid]).max())
                    vbuf = v_vcum[ob0:ob1].astype(np.float64)

                    for i in np.where(valid)[0]:
                        n = int(capped[i])
                        rel = int(starts[i] - ob0)
                        vcum = vbuf[rel:rel + n]
                        w = event_window_maxima(vcum, t_steps)
                        Ri, Ci, ni, r2i = fit_rc(w, t_hours)
                        R[i] = np.float32(Ri)
                        C[i] = np.float32(Ci)
                        npts[i] = np.int16(ni)
                        r2[i] = np.float32(r2i)

                fout.variables["event_id"][sl] = event_id
                fout.variables["year"][sl] = year
                fout.variables["lat"][sl] = lat
                fout.variables["lon"][sl] = lon
                fout.variables["R"][sl] = R
                fout.variables["C"][sl] = C
                fout.variables["n_pts"][sl] = npts
                fout.variables["r2"][sl] = r2
                fout.sync()

            fout.complete = "true"
        finally:
            fout.close()

        os.replace(tmp_path, out_path)
        return f"[OK] {in_path.name} -> {out_path.name} events={ne:,}"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute R/C scaling parameters from flood_events__*.nc files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--in-dir", type=Path, required=True, help="Folder containing flood_events__*.nc files.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Folder for scaling_RC__*.nc outputs.")
    parser.add_argument("--t-hours", default=DEFAULT_T_HOURS, help="Comma-separated T values in hours.")
    parser.add_argument("--event-batch", type=int, default=10_000, help="Events per CPU batch.")
    parser.add_argument("--max-event-len", type=int, default=512, help="Maximum observations per event to use.")
    parser.add_argument("--compression", type=int, default=2, help="NetCDF compression level; 0 disables compression.")
    parser.add_argument("--flat", action="store_true", help="Do not preserve input subfolder structure.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs.")
    parser.add_argument("--require-save-obs", action="store_true", help="Require input attr save_obs=true.")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    t_hours = parse_t_hours(args.t_hours)
    files = iter_input_files(args.in_dir)
    if not files:
        raise FileNotFoundError(f"No flood_events__*.nc files found under {args.in_dir}")

    for in_path in files:
        out_path = output_path_for(in_path, args.in_dir, args.out_dir, keep_structure=not args.flat)
        msg = process_file(
            in_path=in_path,
            out_path=out_path,
            t_hours_input=t_hours,
            event_batch=args.event_batch,
            max_event_len=args.max_event_len,
            compression=args.compression,
            overwrite=args.overwrite,
            require_save_obs=args.require_save_obs,
        )
        print(msg, flush=True)


if __name__ == "__main__":
    main()
