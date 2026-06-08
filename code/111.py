# -*- coding: utf-8 -*-
"""
05_original_algorithm_threshold_robustness.py

Purpose
-------
Use the original threshold algorithm:

1. Equal-sample area bins:
   - sort stations by log10(area)
   - split into equal-size bins
   - compute percent_increase per bin

2. Threshold scan:
   - candidate thresholds from quantiles of log10(area)
   - small = log10A <= threshold
   - large = log10A > threshold
   - p_small = fraction of increase stations in small group
   - p_large = fraction of increase stations in large group
   - J = (p_small - p_large)^2

3. Robustness:
   - one full-sample threshold per climate zone
   - global full-sample threshold
   - global repeated subsampling robustness

Input:
    RC_figure_csvs_1minus_beta/
        01_station_map_runoff_generation_R_trend.csv
        02_station_map_confluence_C_trend_from_1minus_beta.csv

Output:
    RC_robustness_original_algorithm_1minus_beta/
"""

from __future__ import annotations

from pathlib import Path
import warnings
from typing import Dict

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
# MANUAL_CSV_DIR = r"D:\workroom\GPLW\RC_figure_csvs_1minus_beta"
MANUAL_CSV_DIR = None

CLIMATE_SHP = Path(r"D:\workroom\GPLW\气候区划\气候区划.shp")
CLIMATE_FIELD = "Cli_Zone"

CLIMATE_ZONE_ORDER = ["A", "B", "C", "D", "E"]

ROBUST_FOLDER_NAME = "RC_robustness_original_algorithm_1minus_beta"

DPI = 900

# 面积分箱设置
N_AREA_BINS = 12
MIN_N_PER_BIN = 25

# 阈值扫描设置：与原算法一致
SCAN_QMIN = 0.10
SCAN_QMAX = 0.90
SCAN_N_Q = 101
SCAN_MIN_N_SIDE = 25

# 抽样 robustness
N_SAMPLING = 1000
RANDOM_SEED = 3098
SAMPLING_MODE = "subsample"      # "subsample" or "bootstrap"
SUBSAMPLE_FRACTION = 0.90

# 重要：
# True  = Nsig ↑ / Nsig ↓ 统一作为 none，更接近 increase/decrease/none 算法
# False = Nsig ↑ 算 increase，Nsig ↓ 算 decrease
TREAT_NONSIG_AS_NONE = True


# ============================================================
# 1. Style
# ============================================================

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
    "R_light": "#E6CFB4",
    "C_light": "#A9D6CE",
    "axis": "#333333",
    "zero": "#9A9A9A",
}

CLIMATE_ZONE_COLORS = {
    "A": "#3B8EA5",
    "B": "#C9822B",
    "C": "#4C9F70",
    "D": "#7A6BB7",
    "E": "#777777",
}


# ============================================================
# 2. Path helpers
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
# 3. Utilities
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


def panel_title(ax, text: str):
    ax.set_title(text, loc="left", pad=5, fontsize=10)


def area_label(area_km2: float) -> str:
    if not np.isfinite(area_km2):
        return ""
    if area_km2 >= 1e6:
        return f"{area_km2/1e6:.1f}M"
    if area_km2 >= 1e3:
        return f"{area_km2/1e3:.0f}k"
    return f"{area_km2:.0f}"


# ============================================================
# 4. Prepare station data
# ============================================================

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
            print(f"[INFO] Use area column: {col}")
            return np.log10(area)

    print("[ERROR] Cannot detect area column. Available columns:")
    for c in df.columns:
        print("  -", c)
    raise KeyError("No basin area column found.")


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

    # 如果只有 log10_area，就反推 area_km2
    log_area = detect_log10_area(df)
    return 10 ** log_area


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
    """
    Convert R_class / C_class into:
        increase / decrease / none

    If TREAT_NONSIG_AS_NONE=True:
        Sig ↑   -> increase
        Sig ↓   -> decrease
        Nsig ↑  -> none
        Nsig ↓  -> none

    If TREAT_NONSIG_AS_NONE=False:
        Sig ↑ and Nsig ↑ -> increase
        Sig ↓ and Nsig ↓ -> decrease
    """

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
# 5. Original algorithm functions
# ============================================================

def fit_fraction_vs_logarea(tab: pd.DataFrame, x_col: str, y_col: str) -> Dict[str, float]:
    x = pd.to_numeric(tab[x_col], errors="coerce").to_numpy(float)
    y = pd.to_numeric(tab[y_col], errors="coerce").to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)

    if m.sum() < 3:
        return {
            "n_points": int(m.sum()),
            "slope": np.nan,
            "intercept": np.nan,
            "r": np.nan,
            "p": np.nan,
            "r2": np.nan,
        }

    try:
        from scipy.stats import linregress
        res = linregress(x[m], y[m])
        return {
            "n_points": int(m.sum()),
            "slope": float(res.slope),
            "intercept": float(res.intercept),
            "r": float(res.rvalue),
            "p": float(res.pvalue),
            "r2": float(res.rvalue ** 2),
        }
    except Exception:
        k, b = np.polyfit(x[m], y[m], 1)
        r = np.corrcoef(x[m], y[m])[0, 1]
        return {
            "n_points": int(m.sum()),
            "slope": float(k),
            "intercept": float(b),
            "r": float(r),
            "p": np.nan,
            "r2": float(r ** 2),
        }


def compute_area_bin_fraction(
    df: pd.DataFrame,
    dir_col: str,
    var_name: str,
    group_name: str,
    n_bins: int = N_AREA_BINS,
    min_n: int = MIN_N_PER_BIN,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Original algorithm:
    equal-sample area bins, then percent increase per bin.
    """

    d = df.dropna(subset=["area_km2", "log10_area"]).copy()
    d = d[d[dir_col].isin(["increase", "decrease", "none"])].copy()
    d = d.sort_values("log10_area").reset_index(drop=True)

    n = len(d)

    if n == 0:
        return pd.DataFrame(), pd.DataFrame()

    max_bins_by_min_n = max(1, n // min_n)
    n_bins_eff = max(1, min(n_bins, max_bins_by_min_n))

    bin_ids = np.array_split(np.arange(n), n_bins_eff)

    rows = []

    for i, idx in enumerate(bin_ids, start=1):
        sub = d.iloc[idx]

        n_all = len(sub)
        n_inc = int((sub[dir_col] == "increase").sum())
        n_dec = int((sub[dir_col] == "decrease").sum())
        n_none = int((sub[dir_col] == "none").sum())

        rows.append({
            "var": var_name,
            "group": group_name,
            "bin_id": i,
            "n": n_all,
            "n_increase": n_inc,
            "n_decrease": n_dec,
            "n_none": n_none,
            "fraction_increase": n_inc / n_all,
            "percent_increase": n_inc / n_all * 100.0,
            "log10_area_min": float(sub["log10_area"].min()),
            "log10_area_median": float(sub["log10_area"].median()),
            "log10_area_max": float(sub["log10_area"].max()),
            "area_km2_median": float(sub["area_km2"].median()),
        })

    out = pd.DataFrame(rows)

    fit = fit_fraction_vs_logarea(
        out,
        x_col="log10_area_median",
        y_col="percent_increase",
    )
    fit["var"] = var_name
    fit["group"] = group_name

    fit_df = pd.DataFrame([fit])

    return out, fit_df


def compute_threshold_scan_original(
    df: pd.DataFrame,
    dir_col: str,
    var_name: str,
    group_name: str,
    qmin: float = SCAN_QMIN,
    qmax: float = SCAN_QMAX,
    n_q: int = SCAN_N_Q,
    min_n_side: int = SCAN_MIN_N_SIDE,
) -> tuple[pd.DataFrame, dict]:
    """
    Original threshold scan:
        candidate_logA = quantile(logA, qmin..qmax)
        p_small = mean(increase in small)
        p_large = mean(increase in large)
        J = (p_small - p_large)^2
    """

    d = df.dropna(subset=["area_km2", "log10_area"]).copy()
    d = d[d["area_km2"] > 0].copy()
    d = d[d[dir_col].isin(["increase", "decrease", "none"])].copy()

    if len(d) == 0:
        out = pd.DataFrame()
        best = {}
        return out, best

    trend = (d[dir_col] == "increase").to_numpy(bool)
    logA = d["log10_area"].to_numpy(float)

    candidate_logA = np.quantile(logA, np.linspace(qmin, qmax, n_q))

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
        else:
            p_small = float(trend[small].mean())
            p_large = float(trend[large].mean())
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
            "p_small_increase": p_small,
            "p_large_increase": p_large,
            "small_percent_increase": p_small * 100.0 if np.isfinite(p_small) else np.nan,
            "large_percent_increase": p_large * 100.0 if np.isfinite(p_large) else np.nan,
            "delta_large_minus_small_pct": (p_large - p_small) * 100.0 if np.isfinite(p_small) and np.isfinite(p_large) else np.nan,
            "J_separation": J,
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
        "small_percent_increase": float(best["small_percent_increase"]),
        "large_percent_increase": float(best["large_percent_increase"]),
        "delta_large_minus_small_pct": delta,
        "J_separation": float(best["J_separation"]),
        "direction": direction,
    }


# ============================================================
# 6. Climate-zone full thresholds
# ============================================================

def climate_zone_original_algorithm(station_df: pd.DataFrame, var_name: str):
    bin_tables = []
    fit_tables = []
    scan_tables = []
    summaries = []

    for zone in CLIMATE_ZONE_ORDER:
        sub = station_df[station_df["climate_zone"] == zone].copy()
        n_group = len(sub)

        bin_df, fit_df = compute_area_bin_fraction(
            sub,
            dir_col="direction_for_scan",
            var_name=var_name,
            group_name=zone,
            n_bins=N_AREA_BINS,
            min_n=MIN_N_PER_BIN,
        )

        if len(bin_df) > 0:
            bin_tables.append(bin_df)
        if len(fit_df) > 0:
            fit_tables.append(fit_df)

        scan_df, best = compute_threshold_scan_original(
            sub,
            dir_col="direction_for_scan",
            var_name=var_name,
            group_name=zone,
            qmin=SCAN_QMIN,
            qmax=SCAN_QMAX,
            n_q=SCAN_N_Q,
            min_n_side=SCAN_MIN_N_SIDE,
        )

        if len(scan_df) > 0:
            scan_tables.append(scan_df)

        summaries.append(best_to_summary(best, var_name, zone, n_group))

    bins_all = pd.concat(bin_tables, ignore_index=True) if bin_tables else pd.DataFrame()
    fits_all = pd.concat(fit_tables, ignore_index=True) if fit_tables else pd.DataFrame()
    scans_all = pd.concat(scan_tables, ignore_index=True) if scan_tables else pd.DataFrame()
    summary = pd.DataFrame(summaries)

    return bins_all, fits_all, scans_all, summary


# ============================================================
# 7. Global full + global sampling
# ============================================================

def global_original_algorithm(station_df: pd.DataFrame, var_name: str):
    bin_df, fit_df = compute_area_bin_fraction(
        station_df,
        dir_col="direction_for_scan",
        var_name=var_name,
        group_name="Global",
        n_bins=N_AREA_BINS,
        min_n=MIN_N_PER_BIN,
    )

    scan_df, best = compute_threshold_scan_original(
        station_df,
        dir_col="direction_for_scan",
        var_name=var_name,
        group_name="Global",
        qmin=SCAN_QMIN,
        qmax=SCAN_QMAX,
        n_q=SCAN_N_Q,
        min_n_side=SCAN_MIN_N_SIDE,
    )

    summary = pd.DataFrame([
        best_to_summary(best, var_name, "Global", len(station_df))
    ])

    return bin_df, fit_df, scan_df, summary


def global_sampling_original_algorithm(
    station_df: pd.DataFrame,
    var_name: str,
    full_summary: pd.DataFrame,
):
    rng = np.random.default_rng(RANDOM_SEED + (1 if var_name == "R" else 2))

    d = station_df.dropna(subset=["area_km2", "log10_area"]).copy()
    d = d[d["direction_for_scan"].isin(["increase", "decrease", "none"])].copy()
    d = d[d["area_km2"] > 0].copy()

    n = len(d)

    if SAMPLING_MODE.lower() == "bootstrap":
        sample_size = n
        replace = True
    else:
        sample_size = int(np.floor(n * SUBSAMPLE_FRACTION))
        replace = False

    sample_size = max(sample_size, SCAN_MIN_N_SIDE * 2)
    sample_size = min(sample_size, n)

    idx_all = np.arange(n)

    records = []

    for i in range(N_SAMPLING):
        idx = rng.choice(idx_all, size=sample_size, replace=replace)
        sub = d.iloc[idx].copy()

        scan_df, best = compute_threshold_scan_original(
            sub,
            dir_col="direction_for_scan",
            var_name=var_name,
            group_name="Global_sample",
            qmin=SCAN_QMIN,
            qmax=SCAN_QMAX,
            n_q=SCAN_N_Q,
            min_n_side=SCAN_MIN_N_SIDE,
        )

        if not best:
            records.append({
                "var": var_name,
                "sample_id": i,
                "valid": False,
                "sample_mode": SAMPLING_MODE,
                "sample_size": sample_size,
                "best_log10A": np.nan,
                "best_area_km2": np.nan,
                "delta_large_minus_small_pct": np.nan,
                "J_separation": np.nan,
            })
        else:
            records.append({
                "var": var_name,
                "sample_id": i,
                "valid": True,
                "sample_mode": SAMPLING_MODE,
                "sample_size": sample_size,
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
    valid = samples[samples["valid"] == True].copy()

    fs = full_summary.iloc[0]

    if len(valid) == 0:
        summary = pd.DataFrame([{
            "var": var_name,
            "n_full": n,
            "n_valid_samples": 0,
            "full_best_log10A": fs["best_log10A"],
            "full_best_area_km2": fs["best_area_km2"],
        }])
        return samples, summary

    full_direction_positive = fs["delta_large_minus_small_pct"] > 0
    sample_direction_positive = valid["delta_large_minus_small_pct"] > 0
    same_direction = sample_direction_positive == full_direction_positive

    summary = pd.DataFrame([{
        "var": var_name,
        "n_full": n,
        "sample_mode": SAMPLING_MODE,
        "sample_size": sample_size,
        "n_valid_samples": len(valid),

        "full_best_log10A": fs["best_log10A"],
        "full_best_area_km2": fs["best_area_km2"],
        "full_delta_large_minus_small_pct": fs["delta_large_minus_small_pct"],
        "full_J_separation": fs["J_separation"],

        "sampling_threshold_log10A_median": valid["best_log10A"].median(),
        "sampling_threshold_log10A_q05": valid["best_log10A"].quantile(0.05),
        "sampling_threshold_log10A_q95": valid["best_log10A"].quantile(0.95),

        "sampling_threshold_area_km2_median": valid["best_area_km2"].median(),
        "sampling_threshold_area_km2_q05": valid["best_area_km2"].quantile(0.05),
        "sampling_threshold_area_km2_q95": valid["best_area_km2"].quantile(0.95),

        "sampling_delta_large_minus_small_pct_median": valid["delta_large_minus_small_pct"].median(),
        "sampling_delta_large_minus_small_pct_q05": valid["delta_large_minus_small_pct"].quantile(0.05),
        "sampling_delta_large_minus_small_pct_q95": valid["delta_large_minus_small_pct"].quantile(0.95),

        "sampling_J_separation_median": valid["J_separation"].median(),
        "sampling_J_separation_q05": valid["J_separation"].quantile(0.05),
        "sampling_J_separation_q95": valid["J_separation"].quantile(0.95),

        "large_more_increase_share": sample_direction_positive.mean(),
        "same_direction_as_full_share": same_direction.mean(),
    }])

    return samples, summary


# ============================================================
# 8. Figures
# ============================================================

def plot_global_area_bin_fraction(bin_r, bin_c, summary_r, summary_c):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))

    for ax, tab, summary, color, title in [
        (axes[0], bin_r, summary_r, COLORS["R"], "(a) R area-bin increase fraction"),
        (axes[1], bin_c, summary_c, COLORS["C"], "(b) C area-bin increase fraction"),
    ]:
        panel_title(ax, title)

        ax.scatter(
            tab["log10_area_median"],
            tab["percent_increase"],
            s=28,
            color=color,
            alpha=0.55,
            edgecolors="none",
        )

        ax.plot(
            tab["log10_area_median"],
            tab["percent_increase"],
            color=color,
            lw=1.2,
            alpha=0.9,
        )

        fs = summary.iloc[0]
        if np.isfinite(fs["best_log10A"]):
            ax.axvline(fs["best_log10A"], color=color, lw=1.0, ls="--")
            ax.text(
                fs["best_log10A"],
                5,
                f"{area_label(fs['best_area_km2'])} km²",
                rotation=90,
                ha="right",
                va="bottom",
                fontsize=8,
                color=color,
            )

        ax.set_xlabel("log10(Basin area km$^2$)")
        ax.set_ylabel("Increasing stations (%)")
        ax.set_ylim(0, 100)
        clean_axis(ax)

    fig.subplots_adjust(left=0.08, right=0.985, top=0.91, bottom=0.16, wspace=0.30)
    save_figure(fig, "Fig_A_global_area_bin_fraction")
    plt.close(fig)


def plot_global_threshold_scan(scan_r, scan_c):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8))

    for ax, scan, color, title in [
        (axes[0], scan_r, COLORS["R"], "(a) R original threshold scan"),
        (axes[1], scan_c, COLORS["C"], "(b) C original threshold scan"),
    ]:
        panel_title(ax, title)

        ax.plot(scan["candidate_log10A"], scan["J_separation"], color=color, lw=1.4)

        best = scan[scan["is_best_threshold"] == True]
        if len(best) > 0:
            bx = float(best["candidate_log10A"].iloc[0])
            by = float(best["J_separation"].iloc[0])
            area = float(best["candidate_A_km2"].iloc[0])

            ax.axvline(bx, color=color, lw=1.0, ls="--")
            ax.scatter(bx, by, s=30, color=color, edgecolor="white", linewidth=0.6, zorder=3)
            ax.text(bx, by, f" {area_label(area)} km²", ha="left", va="center", fontsize=8, color=color)

        ax.set_xlabel("log10(Basin area km$^2$)")
        ax.set_ylabel("J = (p$_{small}$ - p$_{large}$)$^2$")
        clean_axis(ax)

    fig.subplots_adjust(left=0.08, right=0.985, top=0.91, bottom=0.16, wspace=0.30)
    save_figure(fig, "Fig_B_global_threshold_scan_original")
    plt.close(fig)


def plot_climate_thresholds(summary_r, summary_c):
    r = summary_r.set_index("group").reindex(CLIMATE_ZONE_ORDER)
    c = summary_c.set_index("group").reindex(CLIMATE_ZONE_ORDER)

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.9))

    x = np.arange(len(CLIMATE_ZONE_ORDER))
    width = 0.36

    ax = axes[0]
    panel_title(ax, "(a) Climate-zone thresholds")

    ax.bar(x - width / 2, r["best_log10A"], width=width, color=COLORS["R"], label="R")
    ax.bar(x + width / 2, c["best_log10A"], width=width, color=COLORS["C"], label="C")

    for xi, val, area in zip(x - width / 2, r["best_log10A"], r["best_area_km2"]):
        if np.isfinite(val):
            ax.text(xi, val + 0.04, area_label(area), ha="center", va="bottom", rotation=90, fontsize=7)

    for xi, val, area in zip(x + width / 2, c["best_log10A"], c["best_area_km2"]):
        if np.isfinite(val):
            ax.text(xi, val + 0.04, area_label(area), ha="center", va="bottom", rotation=90, fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(CLIMATE_ZONE_ORDER)
    ax.set_ylabel("Best log10(area km$^2$)")
    ax.set_xlabel("Climate zone")
    ax.legend(frameon=False)
    clean_axis(ax)

    ax = axes[1]
    panel_title(ax, "(b) Large - small increase fraction")

    ax.axhline(0, color=COLORS["zero"], lw=0.8, ls="--")
    ax.bar(x - width / 2, r["delta_large_minus_small_pct"], width=width, color=COLORS["R"], label="R")
    ax.bar(x + width / 2, c["delta_large_minus_small_pct"], width=width, color=COLORS["C"], label="C")

    ax.set_xticks(x)
    ax.set_xticklabels(CLIMATE_ZONE_ORDER)
    ax.set_ylabel("Large - small increase fraction (%)")
    ax.set_xlabel("Climate zone")
    ax.legend(frameon=False)
    clean_axis(ax)

    fig.subplots_adjust(left=0.075, right=0.985, top=0.91, bottom=0.16, wspace=0.32)
    save_figure(fig, "Fig_C_climate_zone_thresholds_original")
    plt.close(fig)


def plot_sampling_threshold(samples_r, samples_c):
    fig, ax = plt.subplots(figsize=(6.4, 4.1))

    data = [
        samples_r.loc[samples_r["valid"] == True, "best_log10A"].dropna().values,
        samples_c.loc[samples_c["valid"] == True, "best_log10A"].dropna().values,
    ]

    bp = ax.boxplot(
        data,
        labels=["R", "C"],
        showfliers=False,
        patch_artist=True,
        widths=0.55,
    )

    for patch, color in zip(bp["boxes"], [COLORS["R"], COLORS["C"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
        patch.set_edgecolor(color)

    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.1)

    ax.set_ylabel("Best log10(area km$^2$)")
    ax.set_title("Global threshold robustness under original algorithm")
    clean_axis(ax)

    save_figure(fig, "Fig_D_global_sampling_threshold_original")
    plt.close(fig)


def plot_sampling_delta(samples_r, samples_c):
    fig, ax = plt.subplots(figsize=(6.4, 4.1))

    data = [
        samples_r.loc[samples_r["valid"] == True, "delta_large_minus_small_pct"].dropna().values,
        samples_c.loc[samples_c["valid"] == True, "delta_large_minus_small_pct"].dropna().values,
    ]

    bp = ax.boxplot(
        data,
        labels=["R", "C"],
        showfliers=False,
        patch_artist=True,
        widths=0.55,
    )

    for patch, color in zip(bp["boxes"], [COLORS["R"], COLORS["C"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
        patch.set_edgecolor(color)

    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.1)

    ax.axhline(0, color=COLORS["zero"], lw=0.8, ls="--")
    ax.set_ylabel("Large - small increase fraction (%)")
    ax.set_title("Direction robustness under original algorithm")
    clean_axis(ax)

    save_figure(fig, "Fig_E_global_sampling_delta_original")
    plt.close(fig)


# ============================================================
# 9. Main
# ============================================================

def main():
    print("\n" + "=" * 100)
    print("Original algorithm threshold robustness")
    print("=" * 100)
    print("[TREAT_NONSIG_AS_NONE]", TREAT_NONSIG_AS_NONE)
    print("[SCAN_QMIN, SCAN_QMAX, SCAN_N_Q]", SCAN_QMIN, SCAN_QMAX, SCAN_N_Q)

    r_raw = read_csv("01_station_map_runoff_generation_R_trend.csv")
    c_raw = read_csv("02_station_map_confluence_C_trend_from_1minus_beta.csv")

    r_station = prepare_station_df(r_raw, class_col="R_class", metric="R")
    c_station = prepare_station_df(c_raw, class_col="C_class", metric="C")

    # --------------------------------------------------------
    # 1. Climate-zone full analysis
    # --------------------------------------------------------
    bins_zone_r, fits_zone_r, scans_zone_r, summary_zone_r = climate_zone_original_algorithm(r_station, "R")
    bins_zone_c, fits_zone_c, scans_zone_c, summary_zone_c = climate_zone_original_algorithm(c_station, "C")

    bins_zone_r.to_csv(ROBUST_DIR / "original_climate_zone_area_bins_R.csv", index=False, encoding="utf-8-sig")
    fits_zone_r.to_csv(ROBUST_DIR / "original_climate_zone_area_bin_linear_fit_R.csv", index=False, encoding="utf-8-sig")
    scans_zone_r.to_csv(ROBUST_DIR / "original_climate_zone_threshold_scan_R.csv", index=False, encoding="utf-8-sig")
    summary_zone_r.to_csv(ROBUST_DIR / "original_climate_zone_threshold_summary_R.csv", index=False, encoding="utf-8-sig")

    bins_zone_c.to_csv(ROBUST_DIR / "original_climate_zone_area_bins_C.csv", index=False, encoding="utf-8-sig")
    fits_zone_c.to_csv(ROBUST_DIR / "original_climate_zone_area_bin_linear_fit_C.csv", index=False, encoding="utf-8-sig")
    scans_zone_c.to_csv(ROBUST_DIR / "original_climate_zone_threshold_scan_C.csv", index=False, encoding="utf-8-sig")
    summary_zone_c.to_csv(ROBUST_DIR / "original_climate_zone_threshold_summary_C.csv", index=False, encoding="utf-8-sig")

    climate_combined = pd.concat([summary_zone_r, summary_zone_c], ignore_index=True)
    climate_combined.to_csv(ROBUST_DIR / "original_climate_zone_threshold_summary_combined.csv", index=False, encoding="utf-8-sig")

    print("\n[Original climate-zone R threshold]")
    print(summary_zone_r.to_string(index=False))

    print("\n[Original climate-zone C threshold]")
    print(summary_zone_c.to_string(index=False))

    # --------------------------------------------------------
    # 2. Global full analysis
    # --------------------------------------------------------
    bin_global_r, fit_global_r, scan_global_r, summary_global_r = global_original_algorithm(r_station, "R")
    bin_global_c, fit_global_c, scan_global_c, summary_global_c = global_original_algorithm(c_station, "C")

    bin_global_r.to_csv(ROBUST_DIR / "original_global_area_bins_R.csv", index=False, encoding="utf-8-sig")
    fit_global_r.to_csv(ROBUST_DIR / "original_global_area_bin_linear_fit_R.csv", index=False, encoding="utf-8-sig")
    scan_global_r.to_csv(ROBUST_DIR / "original_global_threshold_scan_R.csv", index=False, encoding="utf-8-sig")
    summary_global_r.to_csv(ROBUST_DIR / "original_global_threshold_summary_R.csv", index=False, encoding="utf-8-sig")

    bin_global_c.to_csv(ROBUST_DIR / "original_global_area_bins_C.csv", index=False, encoding="utf-8-sig")
    fit_global_c.to_csv(ROBUST_DIR / "original_global_area_bin_linear_fit_C.csv", index=False, encoding="utf-8-sig")
    scan_global_c.to_csv(ROBUST_DIR / "original_global_threshold_scan_C.csv", index=False, encoding="utf-8-sig")
    summary_global_c.to_csv(ROBUST_DIR / "original_global_threshold_summary_C.csv", index=False, encoding="utf-8-sig")

    global_combined = pd.concat([summary_global_r, summary_global_c], ignore_index=True)
    global_combined.to_csv(ROBUST_DIR / "original_global_threshold_summary_combined.csv", index=False, encoding="utf-8-sig")

    print("\n[Original global R threshold]")
    print(summary_global_r.to_string(index=False))

    print("\n[Original global C threshold]")
    print(summary_global_c.to_string(index=False))

    # --------------------------------------------------------
    # 3. Global sampling robustness
    # --------------------------------------------------------
    samples_r, sampling_summary_r = global_sampling_original_algorithm(
        r_station,
        var_name="R",
        full_summary=summary_global_r,
    )

    samples_c, sampling_summary_c = global_sampling_original_algorithm(
        c_station,
        var_name="C",
        full_summary=summary_global_c,
    )

    samples_r.to_csv(ROBUST_DIR / "original_global_sampling_samples_R.csv", index=False, encoding="utf-8-sig")
    sampling_summary_r.to_csv(ROBUST_DIR / "original_global_sampling_summary_R.csv", index=False, encoding="utf-8-sig")

    samples_c.to_csv(ROBUST_DIR / "original_global_sampling_samples_C.csv", index=False, encoding="utf-8-sig")
    sampling_summary_c.to_csv(ROBUST_DIR / "original_global_sampling_summary_C.csv", index=False, encoding="utf-8-sig")

    sampling_combined = pd.concat([sampling_summary_r, sampling_summary_c], ignore_index=True)
    sampling_combined.to_csv(ROBUST_DIR / "original_global_sampling_summary_combined.csv", index=False, encoding="utf-8-sig")

    print("\n[Original global sampling summary]")
    print(sampling_combined.to_string(index=False))

    # --------------------------------------------------------
    # 4. Figures
    # --------------------------------------------------------
    plot_global_area_bin_fraction(bin_global_r, bin_global_c, summary_global_r, summary_global_c)
    plot_global_threshold_scan(scan_global_r, scan_global_c)
    plot_climate_thresholds(summary_zone_r, summary_zone_c)
    plot_sampling_threshold(samples_r, samples_c)
    plot_sampling_delta(samples_r, samples_c)

    print("\n[DONE]")
    print("All outputs saved to:")
    print(ROBUST_DIR.resolve())


if __name__ == "__main__":
    main()