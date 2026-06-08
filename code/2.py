# -*- coding: utf-8 -*-
"""
07_climate_balance_robustness_tests.py

Purpose
-------
Additional robustness tests for climate-zone imbalance.

This script implements three tests:

1. Equal-climate resampling
   - ABCDE_n100: sample 100 stations from each A/B/C/D/E climate zone.
   - ABCD_n250: sample 250 stations from each A/B/C/D zone, excluding E.

2. Leave-one-climate-out
   - Remove one climate zone at a time and recompute global thresholds.

3. Climate-weighted threshold scan
   - Assign station weight w_i = 1 / n_climate.
   - Each climate zone contributes equal total weight.
   - Compute weighted p_small, p_large, and J.

Threshold algorithm
-------------------
Station-level threshold scan along log10(area):

    small = log10(area) <= threshold
    large = log10(area) > threshold

    J = (p_small - p_large)^2

where p_small and p_large are fractions of increasing stations.

Inputs
------
RC_figure_csvs_1minus_beta/
    01_station_map_runoff_generation_R_trend.csv
    02_station_map_confluence_C_trend_from_1minus_beta.csv

Outputs
-------
RC_robustness_climate_balance_tests_1minus_beta/
"""

from __future__ import annotations

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=UserWarning)

try:
    import geopandas as gpd
except Exception:
    gpd = None


# ============================================================
# 0. CONFIG
# ============================================================

REQUIRED_INPUT_FILES = [
    "01_station_map_runoff_generation_R_trend.csv",
    "02_station_map_confluence_C_trend_from_1minus_beta.csv",
]

# 如果自动找不到，手动填绝对路径：
# MANUAL_CSV_DIR = r"D:\workroom\GPLW\code\RC_figure_csvs_1minus_beta"
MANUAL_CSV_DIR = None

CLIMATE_SHP = Path(r"D:\workroom\GPLW\气候区划\气候区划.shp")
CLIMATE_FIELD = "Cli_Zone"

CLIMATE_ZONE_ORDER = ["A", "B", "C", "D", "E"]

ROBUST_FOLDER_NAME = "RC_robustness_climate_balance_tests_1minus_beta"

# 原始阈值扫描设置
SCAN_QMIN = 0.05
SCAN_QMAX = 0.90
SCAN_N_Q = 101
SCAN_MIN_N_SIDE = 25

# Equal-climate resampling
N_RESAMPLING = 1000
RANDOM_SEED = 3098

EQUAL_CLIMATE_SCENARIOS = {
    "ABCDE_n100": {
        "zones": ["A", "B", "C", "D", "E"],
        "n_per_zone": 100,
    },
    "ABCD_n250": {
        "zones": ["A", "B", "C", "D"],
        "n_per_zone": 250,
    },
}

# Weighted threshold scan
# "unweighted": candidate thresholds from ordinary quantiles of log10(area)
# "weighted": candidate thresholds from climate-weighted quantiles of log10(area)
WEIGHTED_CANDIDATE_MODE = "unweighted"

# 如果没有 R_dir/C_dir，则从 R_class/C_class 推断方向
TREAT_NONSIG_AS_NONE = True

DPI = 700

plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.8,
})

COLORS = {
    "R": "#B6783F",
    "C": "#2B8C7D",
    "axis": "#333333",
    "zero": "#9A9A9A",
    "gray": "#777777",
}


# ============================================================
# 1. Path helpers
# ============================================================

def is_valid_csv_dir(folder: Path) -> bool:
    folder = Path(folder)
    return folder.exists() and all((folder / f).exists() for f in REQUIRED_INPUT_FILES)


def auto_find_csv_dir() -> Path:
    if MANUAL_CSV_DIR is not None:
        p = Path(MANUAL_CSV_DIR)
        if is_valid_csv_dir(p):
            print(f"[FOUND] CSV_DIR from MANUAL_CSV_DIR: {p}")
            return p
        print(f"[WARNING] MANUAL_CSV_DIR invalid: {p}")

    cwd = Path.cwd()
    print("[INFO] Current working directory:", cwd)

    candidates = [
        cwd / "RC_figure_csvs_1minus_beta",
        cwd / "RC_figure_csvs",
        cwd,
    ]

    for p in candidates:
        if is_valid_csv_dir(p):
            print("[FOUND] CSV_DIR:", p)
            return p

    for parent in [cwd] + list(cwd.parents):
        p = parent / "RC_figure_csvs_1minus_beta"
        if is_valid_csv_dir(p):
            print("[FOUND] CSV_DIR in parent:", p)
            return p

    print("[INFO] Searching recursively under current directory...")
    for p in cwd.rglob("RC_figure_csvs_1minus_beta"):
        if is_valid_csv_dir(p):
            print("[FOUND] CSV_DIR by recursive search:", p)
            return p

    raise FileNotFoundError(
        "Cannot find required CSV files. "
        "Please set MANUAL_CSV_DIR to the absolute path of RC_figure_csvs_1minus_beta."
    )


CSV_DIR = auto_find_csv_dir()
BASE_DIR = CSV_DIR.parent
ROBUST_DIR = BASE_DIR / ROBUST_FOLDER_NAME
ROBUST_DIR.mkdir(parents=True, exist_ok=True)

print("[INPUT CSV_DIR]", CSV_DIR.resolve())
print("[ROBUST OUT]", ROBUST_DIR.resolve())


# ============================================================
# 2. Utilities
# ============================================================

def read_csv(name: str) -> pd.DataFrame:
    path = CSV_DIR / name
    if not path.exists():
        print("[ERROR] Missing:", path)
        print("[DEBUG] Existing files in CSV_DIR:")
        for p in CSV_DIR.iterdir():
            print("  -", p.name)
        raise FileNotFoundError(path)

    print("[READ]", path)
    return pd.read_csv(path)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        out = ROBUST_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        print("[SAVE]", out)


def clean_axis(ax):
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color(COLORS["axis"])
    ax.tick_params(length=3, width=0.7, color=COLORS["axis"], pad=2)


def area_label(area_km2: float) -> str:
    if not np.isfinite(area_km2):
        return ""
    if area_km2 >= 1e6:
        return f"{area_km2/1e6:.1f}M"
    if area_km2 >= 1e3:
        return f"{area_km2/1e3:.0f}k"
    return f"{area_km2:.0f}"


# ============================================================
# 3. Data preparation
# ============================================================

def detect_area_km2(df: pd.DataFrame) -> pd.Series:
    area_candidates = [
        "area_km2",
        "basin_area_km2",
        "drainage_area_km2",
        "catchment_area_km2",
        "Area_km2",
        "AREA_KM2",
        "area",
        "basin_area",
        "drainage_area",
        "catchment_area",
    ]

    for col in area_candidates:
        if col in df.columns:
            area = pd.to_numeric(df[col], errors="coerce")
            finite = area[np.isfinite(area)]

            if len(finite) > 0:
                if "km" not in col.lower() and np.nanmax(finite) > 1e7:
                    print(f"[INFO] Area column {col} seems to be m2. Convert to km2.")
                    area = area / 1e6

            area = area.where(area > 0)
            print(f"[INFO] Use area column for area_km2: {col}")
            return area

    if "log10_area" in df.columns:
        return 10 ** pd.to_numeric(df["log10_area"], errors="coerce")

    raise KeyError("No basin area column found.")


def detect_log10_area(df: pd.DataFrame) -> pd.Series:
    log_candidates = [
        "log10_area",
        "log10A",
        "log10_area_km2",
        "log10_basin_area",
        "log10_basin_area_km2",
        "log_area",
    ]

    for col in log_candidates:
        if col in df.columns:
            print(f"[INFO] Use log10 area column: {col}")
            return pd.to_numeric(df[col], errors="coerce")

    area = detect_area_km2(df)
    return np.log10(area.where(area > 0))


def detect_climate_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "climate_zone",
        "Cli_Zone",
        "koppen_zone",
        "Koppen_zone",
        "climate",
        "climate_group",
        "CliType",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    return None


def attach_climate_zone(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    climate_col = detect_climate_column(out)
    if climate_col is not None:
        print(f"[INFO] Use existing climate column: {climate_col}")
        out["climate_zone"] = out[climate_col].astype(str)
        return out

    if gpd is None:
        raise ImportError("geopandas is required because station CSV has no climate column.")

    if not CLIMATE_SHP.exists():
        raise FileNotFoundError(CLIMATE_SHP)

    if "lon" not in out.columns or "lat" not in out.columns:
        raise KeyError("Need lon/lat columns for climate-zone spatial join.")

    print("[INFO] No climate column found. Spatial joining to climate shapefile...")

    clim = gpd.read_file(CLIMATE_SHP)

    if clim.crs is None:
        clim = clim.set_crs("EPSG:4326")
    else:
        clim = clim.to_crs("EPSG:4326")

    if CLIMATE_FIELD not in clim.columns:
        raise KeyError(f"{CLIMATE_FIELD} not in climate shapefile.")

    points = gpd.GeoDataFrame(
        out,
        geometry=gpd.points_from_xy(out["lon"], out["lat"]),
        crs="EPSG:4326",
    )

    clim_small = clim[[CLIMATE_FIELD, "geometry"]].rename(
        columns={CLIMATE_FIELD: "climate_zone"}
    )

    joined = gpd.sjoin(points, clim_small, how="left", predicate="within")
    joined = joined.drop(columns=["index_right", "geometry"], errors="ignore")

    print("[INFO] Climate missing:", joined["climate_zone"].isna().sum())

    return pd.DataFrame(joined)


def class_to_direction(s: pd.Series) -> pd.Series:
    text = s.astype(str).str.lower()

    is_up = (
        text.str.contains("↑", regex=False) |
        text.str.contains("up", regex=False) |
        text.str.contains("increase", regex=False) |
        text.str.contains("increasing", regex=False)
    )

    is_down = (
        text.str.contains("↓", regex=False) |
        text.str.contains("down", regex=False) |
        text.str.contains("decrease", regex=False) |
        text.str.contains("decreasing", regex=False)
    )

    is_nsig = (
        text.str.contains("nsig", regex=False) |
        text.str.contains("non", regex=False) |
        text.str.contains("not sig", regex=False)
    )

    out = pd.Series("none", index=s.index, dtype="object")

    if TREAT_NONSIG_AS_NONE:
        out.loc[is_up & (~is_nsig)] = "increase"
        out.loc[is_down & (~is_nsig)] = "decrease"
    else:
        out.loc[is_up] = "increase"
        out.loc[is_down] = "decrease"

    return out


def detect_existing_direction_col(df: pd.DataFrame, metric: str) -> str | None:
    candidates = [
        f"{metric}_dir",
        f"{metric}_direction",
        f"{metric}_trend_dir",
        f"{metric}_trend_direction",
        f"{metric}_direction_class",
        f"{metric}_trend",
        "direction",
        "trend_direction",
    ]

    for col in candidates:
        if col in df.columns:
            vals = df[col].astype(str).str.lower().unique()
            if any(v in ["increase", "decrease", "none"] for v in vals):
                return col

    return None


def prepare_station_df(df: pd.DataFrame, class_col: str, metric: str) -> pd.DataFrame:
    d = attach_climate_zone(df)

    d["area_km2"] = detect_area_km2(d)
    d["log10_area"] = detect_log10_area(d)

    dir_col = detect_existing_direction_col(d, metric)

    if dir_col is not None:
        print(f"[INFO] Use existing direction column for {metric}: {dir_col}")
        d["direction_for_scan"] = d[dir_col].astype(str).str.lower()
        d.loc[~d["direction_for_scan"].isin(["increase", "decrease", "none"]), "direction_for_scan"] = "none"
    else:
        if class_col not in d.columns:
            print("[ERROR] Missing class column:", class_col)
            print("Available columns:")
            for c in d.columns:
                print("  -", c)
            raise KeyError(class_col)

        print(
            f"[INFO] Derive direction_for_scan from {class_col}; "
            f"TREAT_NONSIG_AS_NONE={TREAT_NONSIG_AS_NONE}"
        )
        d["direction_for_scan"] = class_to_direction(d[class_col])

    d = d.replace([np.inf, -np.inf], np.nan)
    d = d.dropna(subset=["area_km2", "log10_area"]).copy()
    d = d[d["area_km2"] > 0].copy()

    d["climate_zone"] = d["climate_zone"].astype(str)
    d.loc[d["climate_zone"].isin(["nan", "None", "NaN"]), "climate_zone"] = np.nan

    d["metric"] = metric

    print(f"\n[{metric}] prepared station data")
    print("  n total:", len(d))
    print("  direction counts:")
    print(d["direction_for_scan"].value_counts(dropna=False))
    print("  climate counts:")
    print(d["climate_zone"].value_counts(dropna=False).sort_index())

    out = ROBUST_DIR / f"prepared_station_data_{metric}.csv"
    d.to_csv(out, index=False, encoding="utf-8-sig")
    print("[SAVE]", out)

    return d


# ============================================================
# 4. Threshold scan: unweighted and weighted
# ============================================================

def weighted_quantile(values, quantiles, sample_weight=None):
    values = np.asarray(values, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)

    if sample_weight is None:
        return np.quantile(values, quantiles)

    sample_weight = np.asarray(sample_weight, dtype=float)

    mask = np.isfinite(values) & np.isfinite(sample_weight) & (sample_weight > 0)
    values = values[mask]
    sample_weight = sample_weight[mask]

    sorter = np.argsort(values)
    values = values[sorter]
    sample_weight = sample_weight[sorter]

    cumulative = np.cumsum(sample_weight)
    cumulative = cumulative / cumulative[-1]

    return np.interp(quantiles, cumulative, values)


def compute_threshold_scan(
    df: pd.DataFrame,
    var_name: str,
    group_name: str,
    qmin: float = SCAN_QMIN,
    qmax: float = SCAN_QMAX,
    n_q: int = SCAN_N_Q,
    min_n_side: int = SCAN_MIN_N_SIDE,
    weight_col: str | None = None,
    weighted_candidates: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """
    Station-level threshold scan.

    If weight_col is None:
        p_small and p_large are ordinary fractions.

    If weight_col is provided:
        p_small and p_large are weighted fractions:
            p = sum(w * I_increase) / sum(w)
    """

    d = df.dropna(subset=["area_km2", "log10_area"]).copy()
    d = d[d["area_km2"] > 0].copy()
    d = d[d["direction_for_scan"].isin(["increase", "decrease", "none"])].copy()

    if weight_col is not None:
        d = d.dropna(subset=[weight_col]).copy()
        d = d[d[weight_col] > 0].copy()

    if len(d) == 0:
        return pd.DataFrame(), {}

    trend = (d["direction_for_scan"] == "increase").to_numpy(float)
    logA = d["log10_area"].to_numpy(float)

    if weight_col is None:
        w = np.ones(len(d), dtype=float)
    else:
        w = d[weight_col].to_numpy(float)

    qs = np.linspace(qmin, qmax, n_q)

    if weight_col is not None and weighted_candidates:
        candidate_logA = weighted_quantile(logA, qs, sample_weight=w)
    else:
        candidate_logA = np.quantile(logA, qs)

    rows = []

    for i, th in enumerate(candidate_logA, start=1):
        small = logA <= th
        large = logA > th

        n_small = int(small.sum())
        n_large = int(large.sum())

        if n_small < min_n_side or n_large < min_n_side:
            p_small = np.nan
            p_large = np.nan
            J = np.nan
            w_small_sum = np.nan
            w_large_sum = np.nan
        else:
            w_small = w[small]
            w_large = w[large]

            w_small_sum = float(np.sum(w_small))
            w_large_sum = float(np.sum(w_large))

            p_small = float(np.sum(w_small * trend[small]) / w_small_sum)
            p_large = float(np.sum(w_large * trend[large]) / w_large_sum)

            J = float((p_small - p_large) ** 2)

        rows.append({
            "var": var_name,
            "group": group_name,
            "candidate_id": i,
            "candidate_log10A": float(th),
            "candidate_A_km2": float(10 ** th),
            "n_total": int(len(d)),
            "n_small": n_small,
            "n_large": n_large,
            "weight_small_sum": w_small_sum,
            "weight_large_sum": w_large_sum,
            "p_small_increase": p_small,
            "p_large_increase": p_large,
            "small_percent_increase": p_small * 100.0 if np.isfinite(p_small) else np.nan,
            "large_percent_increase": p_large * 100.0 if np.isfinite(p_large) else np.nan,
            "delta_large_minus_small_pct": (p_large - p_small) * 100.0 if np.isfinite(p_small) and np.isfinite(p_large) else np.nan,
            "J_separation": J,
            "weighted": weight_col is not None,
        })

    out = pd.DataFrame(rows)

    if out["J_separation"].notna().any():
        best_idx = out["J_separation"].idxmax()
        out["is_best_threshold"] = False
        out.loc[best_idx, "is_best_threshold"] = True
        best = out.loc[best_idx].to_dict()
    else:
        out["is_best_threshold"] = False
        best = {}

    return out, best


def best_to_summary(best: dict, var_name: str, group_name: str, n_group: int) -> dict:
    if not best:
        return {
            "var": var_name,
            "group": group_name,
            "n_total": n_group,
            "best_log10A": np.nan,
            "best_area_km2": np.nan,
            "n_small": np.nan,
            "n_large": np.nan,
            "small_percent_increase": np.nan,
            "large_percent_increase": np.nan,
            "delta_large_minus_small_pct": np.nan,
            "J_separation": np.nan,
            "direction": "insufficient samples",
        }

    delta = float(best["delta_large_minus_small_pct"])
    direction = (
        "larger basins show more increases"
        if delta > 0
        else "smaller basins show more increases"
    )

    return {
        "var": var_name,
        "group": group_name,
        "n_total": int(best["n_total"]),
        "best_log10A": float(best["candidate_log10A"]),
        "best_area_km2": float(best["candidate_A_km2"]),
        "n_small": int(best["n_small"]),
        "n_large": int(best["n_large"]),
        "weight_small_sum": best.get("weight_small_sum", np.nan),
        "weight_large_sum": best.get("weight_large_sum", np.nan),
        "small_percent_increase": float(best["small_percent_increase"]),
        "large_percent_increase": float(best["large_percent_increase"]),
        "delta_large_minus_small_pct": delta,
        "J_separation": float(best["J_separation"]),
        "direction": direction,
        "weighted": bool(best.get("weighted", False)),
    }


def full_threshold(station_df: pd.DataFrame, var_name: str, group_name: str):
    scan, best = compute_threshold_scan(
        station_df,
        var_name=var_name,
        group_name=group_name,
        weight_col=None,
        weighted_candidates=False,
    )
    summary = pd.DataFrame([
        best_to_summary(best, var_name, group_name, len(station_df))
    ])
    return scan, summary


# ============================================================
# 5. Test 1: Equal-climate resampling
# ============================================================

def summarize_resampling(samples: pd.DataFrame, var_name: str, scenario: str, global_full_summary: pd.DataFrame) -> pd.DataFrame:
    valid = samples[samples["valid"] == True].copy()
    fs = global_full_summary.iloc[0]

    if len(valid) == 0:
        return pd.DataFrame([{
            "var": var_name,
            "scenario": scenario,
            "n_valid_samples": 0,
        }])

    full_direction_positive = fs["delta_large_minus_small_pct"] > 0
    sample_direction_positive = valid["delta_large_minus_small_pct"] > 0
    same_direction = sample_direction_positive == full_direction_positive

    return pd.DataFrame([{
        "var": var_name,
        "scenario": scenario,
        "n_valid_samples": len(valid),
        "sample_size": int(valid["sample_size"].iloc[0]),

        "global_full_best_area_km2": fs["best_area_km2"],
        "global_full_delta_large_minus_small_pct": fs["delta_large_minus_small_pct"],

        "threshold_log10A_median": valid["best_log10A"].median(),
        "threshold_log10A_q05": valid["best_log10A"].quantile(0.05),
        "threshold_log10A_q25": valid["best_log10A"].quantile(0.25),
        "threshold_log10A_q75": valid["best_log10A"].quantile(0.75),
        "threshold_log10A_q95": valid["best_log10A"].quantile(0.95),

        "threshold_area_km2_median": valid["best_area_km2"].median(),
        "threshold_area_km2_q05": valid["best_area_km2"].quantile(0.05),
        "threshold_area_km2_q25": valid["best_area_km2"].quantile(0.25),
        "threshold_area_km2_q75": valid["best_area_km2"].quantile(0.75),
        "threshold_area_km2_q95": valid["best_area_km2"].quantile(0.95),

        "delta_large_minus_small_pct_median": valid["delta_large_minus_small_pct"].median(),
        "delta_large_minus_small_pct_q05": valid["delta_large_minus_small_pct"].quantile(0.05),
        "delta_large_minus_small_pct_q25": valid["delta_large_minus_small_pct"].quantile(0.25),
        "delta_large_minus_small_pct_q75": valid["delta_large_minus_small_pct"].quantile(0.75),
        "delta_large_minus_small_pct_q95": valid["delta_large_minus_small_pct"].quantile(0.95),

        "J_separation_median": valid["J_separation"].median(),
        "J_separation_q05": valid["J_separation"].quantile(0.05),
        "J_separation_q95": valid["J_separation"].quantile(0.95),

        "large_more_increase_share": sample_direction_positive.mean(),
        "same_direction_as_global_full_share": same_direction.mean(),
    }])


def equal_climate_resampling(
    station_df: pd.DataFrame,
    var_name: str,
    scenario_name: str,
    zones: list[str],
    n_per_zone: int,
    global_full_summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    d = station_df[station_df["climate_zone"].isin(zones)].copy()

    counts = d["climate_zone"].value_counts()
    print(f"\n[EQUAL-CLIMATE] {var_name}, {scenario_name}")
    print("Counts:")
    print(counts.sort_index())

    for z in zones:
        if counts.get(z, 0) < n_per_zone:
            raise ValueError(
                f"Scenario {scenario_name} requires {n_per_zone} stations in zone {z}, "
                f"but only {counts.get(z, 0)} are available."
            )

    rng = np.random.default_rng(RANDOM_SEED + (1 if var_name == "R" else 2) + abs(hash(scenario_name)) % 10000)

    records = []

    for i in range(N_RESAMPLING):
        sampled_parts = []

        for z in zones:
            sub = d[d["climate_zone"] == z]
            sampled = sub.sample(n=n_per_zone, replace=False, random_state=int(rng.integers(0, 2**31 - 1)))
            sampled_parts.append(sampled)

        sample_df = pd.concat(sampled_parts, ignore_index=True)

        scan, best = compute_threshold_scan(
            sample_df,
            var_name=var_name,
            group_name=scenario_name,
            weight_col=None,
            weighted_candidates=False,
        )

        if not best:
            records.append({
                "var": var_name,
                "scenario": scenario_name,
                "sample_id": i,
                "valid": False,
                "sample_size": len(sample_df),
                "best_log10A": np.nan,
                "best_area_km2": np.nan,
                "delta_large_minus_small_pct": np.nan,
                "J_separation": np.nan,
            })
        else:
            records.append({
                "var": var_name,
                "scenario": scenario_name,
                "sample_id": i,
                "valid": True,
                "sample_size": len(sample_df),
                "best_log10A": float(best["candidate_log10A"]),
                "best_area_km2": float(best["candidate_A_km2"]),
                "n_small": int(best["n_small"]),
                "n_large": int(best["n_large"]),
                "small_percent_increase": float(best["small_percent_increase"]),
                "large_percent_increase": float(best["large_percent_increase"]),
                "delta_large_minus_small_pct": float(best["delta_large_minus_small_pct"]),
                "J_separation": float(best["J_separation"]),
            })

    samples = pd.DataFrame(records)
    summary = summarize_resampling(samples, var_name, scenario_name, global_full_summary)

    return samples, summary


# ============================================================
# 6. Test 2: Leave-one-climate-out
# ============================================================

def leave_one_climate_out(station_df: pd.DataFrame, var_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    scans = []
    summaries = []

    d0 = station_df[station_df["climate_zone"].isin(CLIMATE_ZONE_ORDER)].copy()

    for omitted in CLIMATE_ZONE_ORDER:
        sub = d0[d0["climate_zone"] != omitted].copy()
        group_name = f"Leave_out_{omitted}"

        scan, best = compute_threshold_scan(
            sub,
            var_name=var_name,
            group_name=group_name,
            weight_col=None,
            weighted_candidates=False,
        )

        if len(scan) > 0:
            scan.insert(0, "omitted_climate", omitted)
            scans.append(scan)

        summary = best_to_summary(best, var_name, group_name, len(sub))
        summary["omitted_climate"] = omitted
        summaries.append(summary)

    scan_all = pd.concat(scans, ignore_index=True) if scans else pd.DataFrame()
    summary_all = pd.DataFrame(summaries)

    return scan_all, summary_all


# ============================================================
# 7. Test 3: Climate-weighted threshold scan
# ============================================================

def add_equal_climate_weights(station_df: pd.DataFrame) -> pd.DataFrame:
    d = station_df[station_df["climate_zone"].isin(CLIMATE_ZONE_ORDER)].copy()

    counts = d["climate_zone"].value_counts().to_dict()

    d["climate_weight"] = d["climate_zone"].map(lambda z: 1.0 / counts[z])

    # 每个气候区总权重约为 1
    print("\n[CLIMATE WEIGHTS]")
    print(d.groupby("climate_zone")["climate_weight"].sum())

    return d


def climate_weighted_threshold(station_df: pd.DataFrame, var_name: str):
    d = add_equal_climate_weights(station_df)

    weighted_candidates = WEIGHTED_CANDIDATE_MODE.lower() == "weighted"

    scan, best = compute_threshold_scan(
        d,
        var_name=var_name,
        group_name="Climate_weighted",
        weight_col="climate_weight",
        weighted_candidates=weighted_candidates,
    )

    summary = pd.DataFrame([
        best_to_summary(best, var_name, "Climate_weighted", len(d))
    ])

    summary["candidate_mode"] = WEIGHTED_CANDIDATE_MODE

    return scan, summary


# ============================================================
# 8. Figures
# ============================================================

def plot_equal_climate_summary(equal_summary: pd.DataFrame):
    if equal_summary.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))

    for ax, metric, color in [
        (axes[0], "R", COLORS["R"]),
        (axes[1], "C", COLORS["C"]),
    ]:
        d = equal_summary[equal_summary["var"] == metric].copy()
        if d.empty:
            continue

        x = np.arange(len(d))
        med = d["threshold_log10A_median"].to_numpy(float)
        q05 = d["threshold_log10A_q05"].to_numpy(float)
        q95 = d["threshold_log10A_q95"].to_numpy(float)
        yerr = np.vstack([med - q05, q95 - med])

        ax.errorbar(
            x,
            med,
            yerr=yerr,
            fmt="o",
            color=color,
            ecolor=color,
            capsize=3,
            lw=1.0,
            ms=5,
        )

        labels = [
            f"{row['scenario']}\nn={int(row['sample_size'])}"
            for _, row in d.iterrows()
        ]

        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Best log10(area km$^2$)")
        ax.set_title(f"{metric}: equal-climate resampling")
        clean_axis(ax)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.22, wspace=0.30)
    save_figure(fig, "Fig_equal_climate_resampling_thresholds")
    plt.close(fig)


def plot_leave_one_out(loo_summary: pd.DataFrame):
    if loo_summary.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))

    for ax, metric, color in [
        (axes[0], "R", COLORS["R"]),
        (axes[1], "C", COLORS["C"]),
    ]:
        d = loo_summary[loo_summary["var"] == metric].copy()
        if d.empty:
            continue

        d = d.set_index("omitted_climate").reindex(CLIMATE_ZONE_ORDER).reset_index()

        x = np.arange(len(d))

        ax.plot(
            x,
            d["best_log10A"],
            marker="o",
            color=color,
            lw=1.3,
        )

        ax.set_xticks(x)
        ax.set_xticklabels([f"omit {z}" for z in d["omitted_climate"]])
        ax.set_ylabel("Best log10(area km$^2$)")
        ax.set_title(f"{metric}: leave-one-climate-out")
        clean_axis(ax)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.18, wspace=0.30)
    save_figure(fig, "Fig_leave_one_climate_out_thresholds")
    plt.close(fig)


def plot_weighted_scan(scan_r: pd.DataFrame, scan_c: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))

    for ax, scan, metric, color in [
        (axes[0], scan_r, "R", COLORS["R"]),
        (axes[1], scan_c, "C", COLORS["C"]),
    ]:
        if scan.empty:
            continue

        ax.plot(
            scan["candidate_log10A"],
            scan["J_separation"],
            color=color,
            lw=1.4,
        )

        best = scan[scan["is_best_threshold"] == True]
        if len(best) > 0:
            bx = float(best["candidate_log10A"].iloc[0])
            by = float(best["J_separation"].iloc[0])
            area = float(best["candidate_A_km2"].iloc[0])
            ns = int(best["n_small"].iloc[0])
            nl = int(best["n_large"].iloc[0])

            ax.axvline(bx, color=color, lw=1.0, ls="--")
            ax.scatter(bx, by, s=30, color=color, edgecolor="white", linewidth=0.6, zorder=3)
            ax.text(
                bx,
                by,
                f" {area_label(area)} km²\nn={ns}/{nl}",
                ha="left",
                va="center",
                fontsize=8,
                color=color,
            )

        ax.set_title(f"{metric}: climate-weighted scan")
        ax.set_xlabel("Candidate log10(area km$^2$)")
        ax.set_ylabel("Weighted J")
        clean_axis(ax)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.16, wspace=0.30)
    save_figure(fig, "Fig_climate_weighted_threshold_scan")
    plt.close(fig)


# ============================================================
# 9. Main
# ============================================================

def main():
    print("\n" + "=" * 100)
    print("Climate-balance robustness tests")
    print("=" * 100)
    print("[SCAN_QMIN, SCAN_QMAX, SCAN_N_Q]", SCAN_QMIN, SCAN_QMAX, SCAN_N_Q)
    print("[N_RESAMPLING]", N_RESAMPLING)
    print("[WEIGHTED_CANDIDATE_MODE]", WEIGHTED_CANDIDATE_MODE)

    r_raw = read_csv("01_station_map_runoff_generation_R_trend.csv")
    c_raw = read_csv("02_station_map_confluence_C_trend_from_1minus_beta.csv")

    r_station = prepare_station_df(r_raw, class_col="R_class", metric="R")
    c_station = prepare_station_df(c_raw, class_col="C_class", metric="C")

    # --------------------------------------------------------
    # Baseline global full thresholds
    # --------------------------------------------------------
    scan_global_r, summary_global_r = full_threshold(r_station, "R", "Global")
    scan_global_c, summary_global_c = full_threshold(c_station, "C", "Global")

    global_full = pd.concat([summary_global_r, summary_global_c], ignore_index=True)
    global_full.to_csv(
        ROBUST_DIR / "baseline_global_threshold_summary_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n[Baseline global full thresholds]")
    print(global_full.to_string(index=False))

    # --------------------------------------------------------
    # 1. Equal-climate resampling
    # --------------------------------------------------------
    equal_samples_all = []
    equal_summary_all = []

    for metric, station, global_summary in [
        ("R", r_station, summary_global_r),
        ("C", c_station, summary_global_c),
    ]:
        for scenario_name, cfg in EQUAL_CLIMATE_SCENARIOS.items():
            samples, summary = equal_climate_resampling(
                station_df=station,
                var_name=metric,
                scenario_name=scenario_name,
                zones=cfg["zones"],
                n_per_zone=cfg["n_per_zone"],
                global_full_summary=global_summary,
            )

            samples.to_csv(
                ROBUST_DIR / f"equal_climate_samples_{metric}_{scenario_name}.csv",
                index=False,
                encoding="utf-8-sig",
            )
            summary.to_csv(
                ROBUST_DIR / f"equal_climate_summary_{metric}_{scenario_name}.csv",
                index=False,
                encoding="utf-8-sig",
            )

            equal_samples_all.append(samples)
            equal_summary_all.append(summary)

    equal_samples_combined = pd.concat(equal_samples_all, ignore_index=True)
    equal_summary_combined = pd.concat(equal_summary_all, ignore_index=True)

    equal_samples_combined.to_csv(
        ROBUST_DIR / "equal_climate_samples_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )
    equal_summary_combined.to_csv(
        ROBUST_DIR / "equal_climate_summary_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n[Equal-climate resampling summary]")
    print(equal_summary_combined.to_string(index=False))

    # --------------------------------------------------------
    # 2. Leave-one-climate-out
    # --------------------------------------------------------
    scan_loo_r, summary_loo_r = leave_one_climate_out(r_station, "R")
    scan_loo_c, summary_loo_c = leave_one_climate_out(c_station, "C")

    scan_loo_r.to_csv(
        ROBUST_DIR / "leave_one_climate_out_scan_R.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary_loo_r.to_csv(
        ROBUST_DIR / "leave_one_climate_out_summary_R.csv",
        index=False,
        encoding="utf-8-sig",
    )

    scan_loo_c.to_csv(
        ROBUST_DIR / "leave_one_climate_out_scan_C.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary_loo_c.to_csv(
        ROBUST_DIR / "leave_one_climate_out_summary_C.csv",
        index=False,
        encoding="utf-8-sig",
    )

    loo_summary_combined = pd.concat([summary_loo_r, summary_loo_c], ignore_index=True)
    loo_summary_combined.to_csv(
        ROBUST_DIR / "leave_one_climate_out_summary_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n[Leave-one-climate-out summary]")
    print(loo_summary_combined.to_string(index=False))

    # --------------------------------------------------------
    # 3. Climate-weighted threshold scan
    # --------------------------------------------------------
    scan_weighted_r, summary_weighted_r = climate_weighted_threshold(r_station, "R")
    scan_weighted_c, summary_weighted_c = climate_weighted_threshold(c_station, "C")

    scan_weighted_r.to_csv(
        ROBUST_DIR / "climate_weighted_threshold_scan_R.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary_weighted_r.to_csv(
        ROBUST_DIR / "climate_weighted_threshold_summary_R.csv",
        index=False,
        encoding="utf-8-sig",
    )

    scan_weighted_c.to_csv(
        ROBUST_DIR / "climate_weighted_threshold_scan_C.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary_weighted_c.to_csv(
        ROBUST_DIR / "climate_weighted_threshold_summary_C.csv",
        index=False,
        encoding="utf-8-sig",
    )

    weighted_summary_combined = pd.concat([summary_weighted_r, summary_weighted_c], ignore_index=True)
    weighted_summary_combined.to_csv(
        ROBUST_DIR / "climate_weighted_threshold_summary_combined.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n[Climate-weighted threshold summary]")
    print(weighted_summary_combined.to_string(index=False))

    # --------------------------------------------------------
    # Figures
    # --------------------------------------------------------
    plot_equal_climate_summary(equal_summary_combined)
    plot_leave_one_out(loo_summary_combined)
    plot_weighted_scan(scan_weighted_r, scan_weighted_c)

    print("\n[DONE]")
    print("All outputs saved to:")
    print(ROBUST_DIR.resolve())


if __name__ == "__main__":
    main()