# -*- coding: utf-8 -*-
"""
02_plot_two_figures.py

Purpose
-------
Read the CSV files exported by `01_export_csvs_for_two_figures.py` and draw two
publication-style figures:
1) Figure 1:
   (a) Runoff generation trend map
   (b) Confluence trend map, where C = 1 - beta
   (c) Combinations of R/C trend directions
   (d) Trends by Köppen climate zone

2) Figure 2:
   Area-dependence panels for R and C = 1 - beta

This version is adjusted according to the figure-layout request:
    - remove country borders from maps;
    - remove grid lines from all panels;
    - keep coastlines only for geographic reference;
    - align panel titles and text boxes;
    - match the reference map extent (-180–180°, about -80–85°);
    - slightly reduce station marker size to avoid overplotting;
    - enlarge longitude/latitude tick labels on map panels;
    - use reference-style station markers: translucent colored halo + solid colored core with thin white outline;
    - use a consistent teal-brown color palette;
    - export raster figures at DPI = 900 (>800 dpi).

Key definition
--------------
The confluence metric is:
    C = 1 - beta
Therefore, C trends are opposite to beta trends. The CSV exporter has already
converted the trend direction, and this plotting script directly uses C_class.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================
# 0. CONFIG: only modify this part
# ============================================================

CSV_DIR = Path(r"./RC_figure_csvs_1minus_beta")
FIG_DIR = Path(r"./RC_final_figures_1minus_beta")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Figure sizes. Fig.1 is close to A4 portrait; Fig.2 is a compact 2 x 3 layout.
FIG1_SIZE = (8.2, 11.15)
FIG2_SIZE = (9.0, 6.25)
DPI = 900

# Map style switches.
# COUNTRY_BORDERS is intentionally False to remove internal political boundaries.
SHOW_COUNTRY_BORDERS = False
SHOW_PANEL_GRIDS = False

# Font. If Times New Roman is unavailable, matplotlib will fall back automatically.
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

# ============================================================
# Unified color palette
# ============================================================
# Direction colors on maps:
#     increase = teal, decrease = brown
# Process colors in Fig.2:
#     R / runoff generation = brown
#     C / confluence = teal
COLORS = {
    "sig_up": "#2B8C7D",       # significant increase: dark teal
    "nsig_up": "#A9D6CE",      # non-significant increase: light teal
    "sig_down": "#B6783F",     # significant decrease: dark brown
    "nsig_down": "#E6CFB4",    # non-significant decrease: light brown

    "R": "#B6783F",            # runoff generation panels: brown
    "C": "#2B8C7D",            # confluence panels: teal

    "bar_light": "#9CCFC6",
    "bar_dark": "#20796C",

    "coast": "#555555",
    "axis": "#333333",
    "zero": "#9A9A9A",
}

CLASS_ORDER = ["Sig ↑", "Nsig ↑", "Sig ↓", "Nsig ↓"]
CLASS_COLOR = {
    "Sig ↑": COLORS["sig_up"],
    "Nsig ↑": COLORS["nsig_up"],
    "Sig ↓": COLORS["sig_down"],
    "Nsig ↓": COLORS["nsig_down"],
}

# Station marker style.
# Scatter `s` is marker area in points^2. This version follows the reference
# screenshot: a faint colored halo behind each station plus a larger solid
# colored circle with a thin white outline. This avoids the previous problem
# where a white halo made pale points look blurry or washed out.
MAP_MARKER_SIZE = 8.8
MAP_HALO_SIZE = 19.5
MAP_ALPHA = 0.90
MAP_HALO_ALPHA = 0.22
MAP_EDGE_COLOR = "white"
MAP_EDGE_LINEWIDTH = 0.34
MAP_TICK_LABEL_SIZE = 9.8

# Draw non-significant points first, then significant points, so important
# markers are not covered by pale markers.
CLASS_DRAW_ORDER = ["Nsig ↑", "Nsig ↓", "Sig ↑", "Sig ↓"]

# Match the reference figure: full longitude range, with Antarctica visible but no ±90° tick labels.
MAP_EXTENT = [-180, 180, -80, 85]

PERCENT_FMT = FuncFormatter(lambda v, pos: f"{int(v)}%")


# ============================================================
# 1. Utilities
# ============================================================

def read_csv(name: str, required: bool = True) -> pd.DataFrame | None:
    path = CSV_DIR / name
    if not path.exists():
        if required:
            raise FileNotFoundError(f"Required CSV not found: {path}")
        print(f"[SKIP] Optional CSV not found: {path}")
        return None
    return pd.read_csv(path)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        out = FIG_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        print("[SAVE]", out)


def fmt_lon(x: float) -> str:
    x = int(x)
    if x < 0:
        return f"{abs(x)}°W"
    if x > 0:
        return f"{x}°E"
    return "0°"


def fmt_lat(y: float) -> str:
    y = int(y)
    if y < 0:
        return f"{abs(y)}°S"
    if y > 0:
        return f"{y}°N"
    return "0°"


def panel_title(ax, text: str) -> None:
    """Use the same left-aligned title style for every panel."""
    ax.set_title(text, loc="left", pad=5, fontsize=10)


def clean_axis(ax, keep_left: bool = True, keep_bottom: bool = True) -> None:
    """Remove panel grid lines and standardize spines/ticks."""
    ax.grid(False)
    ax.set_axisbelow(False)
    for side, spine in ax.spines.items():
        spine.set_linewidth(0.8)
        spine.set_color(COLORS["axis"])
    if not keep_left:
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", left=False, labelleft=False)
    if not keep_bottom:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(length=3, width=0.7, color=COLORS["axis"], pad=2)


def try_import_cartopy():
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        return ccrs, cfeature
    except Exception:
        return None, None


# ============================================================
# 2. Figure 1 panels
# ============================================================

def setup_map_axis(ax, ccrs, cfeature):
    ax.set_global()
    ax.set_extent(MAP_EXTENT, crs=ccrs.PlateCarree())

    # Clean background. No country borders, no gridlines.
    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="white", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="white", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.52, edgecolor=COLORS["coast"], zorder=2)
    if SHOW_COUNTRY_BORDERS:
        ax.add_feature(cfeature.BORDERS.with_scale("110m"), linewidth=0.25, edgecolor="0.55", zorder=2)

    xticks = [-180, -120, -60, 0, 60, 120, 180]
    yticks = [-60, -30, 0, 30, 60]
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.set_xticklabels([fmt_lon(x) for x in xticks])
    ax.set_yticklabels([fmt_lat(y) for y in yticks])
    ax.tick_params(length=2.5, width=0.7, pad=2, labelsize=MAP_TICK_LABEL_SIZE)

    # No map grid lines. Keep the map frame for alignment.
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("black")


def plot_station_trend_map(ax, df: pd.DataFrame, class_col: str, title: str, ccrs, cfeature):
    setup_map_axis(ax, ccrs, cfeature)

    for cls in CLASS_DRAW_ORDER:
        sub = df[df[class_col] == cls]
        if len(sub) == 0:
            continue

        # Reference-style translucent colored halo.
        # It gives each point a visible footprint, especially in dense regions,
        # but keeps the map background clean.
        ax.scatter(
            sub["lon"], sub["lat"],
            s=MAP_HALO_SIZE,
            color=CLASS_COLOR[cls],
            alpha=MAP_HALO_ALPHA,
            edgecolors="none",
            linewidths=0,
            transform=ccrs.PlateCarree(),
            zorder=3,
            rasterized=True,
        )

        # Solid colored circle with thin white outline.
        ax.scatter(
            sub["lon"], sub["lat"],
            s=MAP_MARKER_SIZE,
            color=CLASS_COLOR[cls],
            alpha=MAP_ALPHA,
            edgecolors=MAP_EDGE_COLOR,
            linewidths=MAP_EDGE_LINEWIDTH,
            transform=ccrs.PlateCarree(),
            zorder=4,
            rasterized=True,
        )

    panel_title(ax, title)

    # Compact and aligned legend. The labels are ordered consistently in Fig.1a/b.
    handles = []
    for cls in CLASS_ORDER:
        n = int((df[class_col] == cls).sum())
        handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="none",
                markerfacecolor=CLASS_COLOR[cls], markeredgecolor=MAP_EDGE_COLOR,
                markeredgewidth=0.55, markersize=5.1, label=f"{cls}   (n = {n})"
            )
        )
    leg = ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.012, 0.022),
        frameon=True,
        framealpha=0.96,
        borderpad=0.75,
        handletextpad=0.55,
        labelspacing=0.45,
    )
    leg.get_frame().set_linewidth(0.5)


def plot_combination_bar(ax, combo_df: pd.DataFrame):
    panel_title(ax, "(c) Combinations of Trends")

    order = ["R↑ & C↑", "R↑ & C↓", "R↓ & C↑", "R↓ & C↓"]
    combo_df = combo_df.set_index("combo").reindex(order).reset_index()

    x = np.arange(len(combo_df))
    pct_total = combo_df["pct_total"].to_numpy(float)
    pct_sig = combo_df["pct_both_significant"].to_numpy(float)

    ax.bar(x, pct_total, width=0.70, color=COLORS["bar_light"], edgecolor="none", label="Total")
    ax.bar(x, pct_sig, width=0.70, color=COLORS["bar_dark"], edgecolor="none", label="Significant")

    for xi, val in zip(x, pct_total):
        if np.isfinite(val):
            ax.text(xi, val + 0.9, f"{val:.1f}%", ha="center", va="bottom", fontsize=9)
    for xi, val in zip(x, pct_sig):
        if np.isfinite(val) and val > 0.4:
            ax.text(xi, max(val * 0.50, 0.45), f"{val:.1f}%", ha="center", va="center", fontsize=8, color="white")

    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("Percentage of stations (%)")
    ax.set_ylim(0, max(42, np.nanmax(pct_total) * 1.22))
    ax.legend(loc="upper right", frameon=False, handlelength=1.4, borderaxespad=0.2)
    clean_axis(ax)


def plot_climate_scatter(ax, climate_df: pd.DataFrame | None):
    panel_title(ax, "(d) Trends by Climate Zone")

    if climate_df is None or len(climate_df) == 0:
        ax.text(
            0.5, 0.5,
            "Climate-zone CSV not found\nSet CLIMATE_SHP and rerun export script",
            ha="center", va="center", transform=ax.transAxes
        )
        ax.set_axis_off()
        return

    d = climate_df.copy()
    x = d["x_for_panel_d_C_net_increase_pct"].to_numpy(float)
    y = d["y_for_panel_d_R_net_increase_pct"].to_numpy(float)

    xerr_low = d.get("C_xerr_low", pd.Series(np.zeros(len(d)))).to_numpy(float)
    xerr_high = d.get("C_xerr_high", pd.Series(np.zeros(len(d)))).to_numpy(float)
    yerr_low = d.get("R_yerr_low", pd.Series(np.zeros(len(d)))).to_numpy(float)
    yerr_high = d.get("R_yerr_high", pd.Series(np.zeros(len(d)))).to_numpy(float)

    # Keep only zero reference lines; remove gridlines.
    ax.axhline(0, color=COLORS["zero"], lw=0.75, ls="--", zorder=0)
    ax.axvline(0, color=COLORS["zero"], lw=0.75, ls="--", zorder=0)

    ax.errorbar(
        x, y,
        xerr=[xerr_low, xerr_high],
        yerr=[yerr_low, yerr_high],
        fmt="o",
        ms=5,
        mfc=COLORS["C"],
        mec="0.25",
        mew=0.5,
        ecolor="0.2",
        elinewidth=0.7,
        capsize=2,
        alpha=0.92,
        zorder=3,
    )

    roman = ["I", "II", "III", "IV", "V"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        label = roman[i] if i < len(roman) else str(i + 1)
        ax.text(xi + 2.2, yi + 2.2, label, fontsize=9, ha="left", va="bottom")

    ax.set_xlabel("Confluence: Net Increase (%)")
    ax.set_ylabel("Runoff generation: Net Increase (%)")

    xmin = min(-70, np.nanmin(x - xerr_low) - 8)
    xmax = max(70, np.nanmax(x + xerr_high) + 8)
    ymin = min(-70, np.nanmin(y - yerr_low) - 8)
    ymax = max(70, np.nanmax(y + yerr_high) + 8)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    # Align explanatory text boxes to left-bottom and right-top corners.
    note = (
        "Köppen 5 groups:\n"
        "I = Tropical · C\n"
        "II = Arid · R · C\n"
        "III = Temperate · R\n"
        "IV = Cold/Continental · C\n"
        "V = Polar/Highland · R · C"
    )
    ax.text(
        0.03, 0.05, note,
        transform=ax.transAxes,
        ha="left", va="bottom", fontsize=7.8,
        linespacing=1.15,
    )

    sig_note = "C: confluence net trend ≠ 0, P<0.05\nR: runoff net trend ≠ 0, P<0.05"
    ax.text(
        0.98, 0.95, sig_note,
        transform=ax.transAxes,
        ha="right", va="top", fontsize=7.8,
        linespacing=1.20,
        bbox=dict(boxstyle="square,pad=0.30", fc="white", ec=COLORS["sig_down"], lw=0.7),
    )

    clean_axis(ax)


def make_figure1():
    ccrs, cfeature = try_import_cartopy()
    if ccrs is None:
        raise ImportError(
            "cartopy is required for the map panels. Install with: conda install -c conda-forge cartopy"
        )

    r_map = read_csv("01_station_map_runoff_generation_R_trend.csv")
    c_map = read_csv("02_station_map_confluence_C_trend_from_1minus_beta.csv")
    combo = read_csv("03_bar_combinations_R_C_trends.csv")
    climate = read_csv("04_climate_zone_net_trend_summary.csv", required=False)

    fig = plt.figure(figsize=FIG1_SIZE)
    gs = fig.add_gridspec(
        nrows=4, ncols=2,
        height_ratios=[1.18, 1.18, 0.030, 0.82],
        width_ratios=[1.04, 1.00],
        hspace=0.12,
        wspace=0.14,
    )

    ax1 = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(gs[1, :], projection=ccrs.PlateCarree())
    ax3 = fig.add_subplot(gs[3, 0])
    ax4 = fig.add_subplot(gs[3, 1])

    plot_station_trend_map(ax1, r_map, "R_class", "(a) Runoff generation trend", ccrs, cfeature)
    plot_station_trend_map(ax2, c_map, "C_class", "(b) Confluence trend", ccrs, cfeature)
    plot_combination_bar(ax3, combo)
    plot_climate_scatter(ax4, climate)

    fig.subplots_adjust(left=0.065, right=0.985, top=0.985, bottom=0.055)
    save_figure(fig, "Fig1_RC_trends_map_bar_climate")
    plt.close(fig)


# ============================================================
# 3. Figure 2 panels
# ============================================================

def read_best_threshold(scan_df: pd.DataFrame) -> float | None:
    if "is_best_threshold" in scan_df.columns:
        m = scan_df["is_best_threshold"].astype(str).str.lower().isin(["true", "1", "yes"])
        if m.any():
            return float(scan_df.loc[m, "candidate_log10A"].iloc[0])
    if "J_separation" in scan_df.columns and scan_df["J_separation"].notna().any():
        idx = scan_df["J_separation"].idxmax()
        return float(scan_df.loc[idx, "candidate_log10A"])
    return None


def plot_area_bin(ax, df: pd.DataFrame, color: str, panel_label: str):
    x = df["log10_area_median"].to_numpy(float)
    y = df["percent_increase"].to_numpy(float)

    ax.scatter(x, y, s=23, color=color, alpha=0.42, edgecolors="none")

    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() >= 3:
        coef = np.polyfit(x[mask], y[mask], 1)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        yy = coef[0] * xx + coef[1]
        ax.plot(xx, yy, color=color, lw=1.4)

    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(PERCENT_FMT)
    ax.set_ylabel("Stations (%)")
    panel_title(ax, panel_label)
    clean_axis(ax)


def plot_threshold_scan(ax, df: pd.DataFrame, color: str, panel_label: str):
    x = df["candidate_log10A"].to_numpy(float)
    y = df["J_separation"].to_numpy(float)
    best = read_best_threshold(df)

    ax.plot(x, y, color=color, lw=1.5)
    if best is not None and np.isfinite(best):
        ax.axvline(best, color=color, lw=1.0, ls="--")
        ymax = np.nanmax(y) if np.isfinite(np.nanmax(y)) else 0.0
        y_text = ymax * 0.12
        ax.text(best - 0.04, y_text, f"{best:.2g}", color=color, ha="right", va="bottom", fontsize=9)

    ax.set_ylabel("Separation (J)")
    panel_title(ax, panel_label)
    clean_axis(ax)


def plot_moving_window(ax, df: pd.DataFrame, color: str, panel_label: str, threshold: float | None = None):
    x = df["log10_area_median"].to_numpy(float)
    y = df["percent_increase"].to_numpy(float)
    ax.plot(x, y, color=color, lw=1.4)

    if threshold is not None and np.isfinite(threshold):
        ax.axvline(threshold, color=color, lw=1.0, ls="--")

    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(PERCENT_FMT)
    ax.set_ylabel("Local fraction (%)")
    panel_title(ax, panel_label)
    clean_axis(ax)


def make_figure2():
    r_bin = read_csv("05_area_bin_fraction_R_increase.csv")
    r_scan = read_csv("06_threshold_scan_R_increase.csv")
    r_mov = read_csv("07_moving_window_R_increase.csv")

    c_bin = read_csv("08_area_bin_fraction_C_increase.csv")
    c_scan = read_csv("09_threshold_scan_C_increase.csv")
    c_mov = read_csv("10_moving_window_C_increase.csv")

    r_best = read_best_threshold(r_scan)
    c_best = read_best_threshold(c_scan)

    fig, axes = plt.subplots(2, 3, figsize=FIG2_SIZE)
    ax = axes.ravel()

    plot_area_bin(ax[0], r_bin, COLORS["R"], "(a) Runoff generation: % of increasing stations")
    plot_threshold_scan(ax[1], r_scan, COLORS["R"], "(b) Runoff generation: threshold scan")
    plot_moving_window(ax[2], r_mov, COLORS["R"], "(c) Runoff generation: moving window", threshold=r_best)

    plot_area_bin(ax[3], c_bin, COLORS["C"], "(d) Confluence: % of increasing stations")
    plot_threshold_scan(ax[4], c_scan, COLORS["C"], "(e) Confluence: threshold scan")
    plot_moving_window(ax[5], c_mov, COLORS["C"], "(f) Confluence: moving window", threshold=c_best)

    # Consistent x-label layout: only bottom row has x-labels.
    for i in [3, 4, 5]:
        ax[i].set_xlabel("log10(Basin area)")

    fig.subplots_adjust(left=0.075, right=0.985, top=0.96, bottom=0.095, wspace=0.30, hspace=0.36)
    save_figure(fig, "Fig2_area_threshold_moving_window")
    plt.close(fig)


# ============================================================
# 4. Main
# ============================================================

def main():
    print("[INPUT CSV_DIR]", CSV_DIR.resolve())
    print("[OUTPUT FIG_DIR]", FIG_DIR.resolve())
    make_figure1()
    make_figure2()
    print("[DONE] Figures generated.")


if __name__ == "__main__":
    main()
