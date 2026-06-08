# ============================================================
# Fig. 3 final version
# A4 portrait + 3x2 symmetric layout
#
# Panels:
#   a. Global R–C fingerprint, no arrows
#   b. Global effect size
#   c. Warming scale dependence
#   d. Human-activity scale dependence
#   e. Global joint R–C response modes
#   f. Scale-stratified dominant joint R–C response mode
#
# Notes:
#   ΔR is expressed as percentage change.
#   ΔC is expressed as absolute change in scaling exponent C.
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import warnings

warnings.filterwarnings("ignore")

try:
    from IPython.display import display
except Exception:
    display = print


# ============================================================
# 0. Paths
# ============================================================

ROOT_DIR = Path(r"D:\workroom\GPLW_furture\ISIMIP")

CANDIDATE_FILES = [
    ROOT_DIR / "_figures_fig3_warming_human_scale" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_warming_human_scale_adjusted" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_warming_human_scale_v2" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_no_design_v3" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_no_design_v4" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_A4_final" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_A4_final_all_climate" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_A4_final_climate_hist" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_A4_final_climate_grouped_hist" / "fig3_all_drivers_with_ddm30_basin_info.csv",
    ROOT_DIR / "_figures_fig3_A4_final_no_arrow_scale_mode" / "fig3_all_drivers_with_ddm30_basin_info.csv",
]

DF_BASIN_FILE = None
for f in CANDIDATE_FILES:
    if f.exists():
        DF_BASIN_FILE = f
        break

if DF_BASIN_FILE is None:
    matches = list(ROOT_DIR.rglob("fig3_all_drivers_with_ddm30_basin_info.csv"))
    if len(matches) > 0:
        DF_BASIN_FILE = matches[0]

if DF_BASIN_FILE is None:
    raise FileNotFoundError("Cannot find fig3_all_drivers_with_ddm30_basin_info.csv")

FIG_DIR = ROOT_DIR / "_figures_fig3_A4_final_clean_1minus_beta_layout_fixed"
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("Using data file:", DF_BASIN_FILE)
print("Output folder:", FIG_DIR)


# ============================================================
# 1. Style
# ============================================================

mpl.rcParams["font.family"] = "Times New Roman"
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.linewidth"] = 0.6
mpl.rcParams["xtick.major.width"] = 0.6
mpl.rcParams["ytick.major.width"] = 0.6
mpl.rcParams["xtick.major.size"] = 2.5
mpl.rcParams["ytick.major.size"] = 2.5
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

PNG_DPI = 900

COLORS = {
    "1.5°C warming": "#4C9A9A",
    "2.0°C warming": "#1F6F78",
    "Irrigation": "#C08A4B",
    "Water management": "#7A4E2D",
}

ORDER = [
    "1.5°C warming",
    "2.0°C warming",
    "Irrigation",
    "Water management",
]

DISPLAY_LABELS = {
    "1.5°C warming": "1.5°C\nwarming",
    "2.0°C warming": "2.0°C\nwarming",
    "Irrigation": "Irrigation",
    "Water management": "Water\nmanagement",
}

WARMING_ORDER = ["1.5°C warming", "2.0°C warming"]
HUMAN_ORDER = ["Irrigation", "Water management"]

QUAD_ORDER = ["R↑ & C↑", "R↑ & C↓", "R↓ & C↑", "R↓ & C↓"]

QUAD_COLORS = {
    "R↑ & C↑": "#4C9A9A",
    "R↑ & C↓": "#9CCBCB",
    "R↓ & C↑": "#C08A4B",
    "R↓ & C↓": "#D8B58A",
}

DARK_SEGMENTS = {"R↑ & C↑", "R↓ & C↑"}

BASIN_CLASS_ORDER = [
    "small model basins",
    "medium model basins",
    "large model basins",
]

BASIN_XLABELS = ["S", "M", "L"]

R_COL = "delta_R_pct_median"

# ============================================================
# IMPORTANT: redefine confluence as C = 1 - beta
# ------------------------------------------------------------
# In the old figure code, "delta_C_median" was plotted directly.
# Here we treat the old "delta_C_median" as Δbeta / Δslope,
# and convert it to:
#
#     ΔC = Δ(1 - beta) = -Δbeta
#
# Therefore:
#     beta increase  -> C decrease
#     beta decrease  -> C increase
# ============================================================
OLD_BETA_DELTA_COL = "delta_C_median"
C_COL = "delta_C_1minus_beta_median"

BASIN_CLASS_COL = "area_class_tercile"

R_LIM = (-60, 60)
C_LIM = (-0.08, 0.08)


# ============================================================
# 2. Load data
# ============================================================

df = pd.read_csv(DF_BASIN_FILE)

df["lon"] = (((df["lon"] + 180) % 360) - 180).round(4)
df["lat"] = df["lat"].round(4)

# ------------------------------------------------------------
# Convert beta/slope response to C response.
# Priority:
#   1) If an explicit beta/slope delta column exists, use it.
#   2) Otherwise use the old plotted column "delta_C_median" as Δbeta.
# ------------------------------------------------------------
candidate_beta_delta_cols = [
    "delta_beta_median",
    "delta_beta",
    "dBeta",
    "dbeta",
    "Delta_beta",
    "delta_slope_median",
    "delta_slope",
    "dSlope",
    "dslope",
    "Slope_delta",
    OLD_BETA_DELTA_COL,
]

BETA_DELTA_COL = None
for col in candidate_beta_delta_cols:
    if col in df.columns:
        BETA_DELTA_COL = col
        break

if BETA_DELTA_COL is None:
    raise KeyError(
        "Cannot find beta/slope delta column. "
        f"Tried: {candidate_beta_delta_cols}\n"
        f"Available columns: {list(df.columns)}"
    )

df[C_COL] = -df[BETA_DELTA_COL]

print("Using R column:", R_COL)
print("Using beta/slope delta column:", BETA_DELTA_COL)
print("Converted C column:", C_COL, " = -", BETA_DELTA_COL)

df = df.dropna(subset=["driver", "lat", "lon", R_COL, C_COL]).copy()

df["R_plot"] = df[R_COL].clip(R_LIM[0], R_LIM[1])
df["C_plot"] = df[C_COL].clip(C_LIM[0], C_LIM[1])

# Save standardized table with converted C for checking and reuse.
df.to_csv(FIG_DIR / "00_input_with_C_1minus_beta.csv", index=False, encoding="utf-8-sig")

print("Loaded data:", df.shape)

display(
    df.groupby("driver")
    .agg(
        n_grid=("lat", "count"),
        median_R=(R_COL, "median"),
        q25_R=(R_COL, lambda x: np.nanpercentile(x, 25)),
        q75_R=(R_COL, lambda x: np.nanpercentile(x, 75)),
        median_C=(C_COL, "median"),
        q25_C=(C_COL, lambda x: np.nanpercentile(x, 25)),
        q75_C=(C_COL, lambda x: np.nanpercentile(x, 75)),
        median_nGHM=("n_ghms", "median"),
    )
    .reindex(ORDER)
)


# ============================================================
# 3. Summary functions
# ============================================================

def q25(x):
    return np.nanpercentile(x, 25)


def q75(x):
    return np.nanpercentile(x, 75)


def classify_quadrant(r, c):
    if pd.isna(r) or pd.isna(c):
        return np.nan
    if r >= 0 and c >= 0:
        return "R↑ & C↑"
    elif r >= 0 and c < 0:
        return "R↑ & C↓"
    elif r < 0 and c >= 0:
        return "R↓ & C↑"
    else:
        return "R↓ & C↓"


def short_mode_label(mode):
    mapping = {
        "R↑ & C↑": "R↑C↑",
        "R↑ & C↓": "R↑C↓",
        "R↓ & C↑": "R↓C↑",
        "R↓ & C↓": "R↓C↓",
    }
    return mapping.get(mode, str(mode))


def global_effect_summary(df):
    rows = []

    for d in ORDER:
        sub = df[df["driver"] == d].dropna(subset=[R_COL, C_COL]).copy()

        rows.append({
            "driver": d,
            "n_grid": len(sub),

            "R_median_raw": sub[R_COL].median(),
            "R_q25_raw": sub[R_COL].quantile(0.25),
            "R_q75_raw": sub[R_COL].quantile(0.75),

            "C_median_raw": sub[C_COL].median(),
            "C_q25_raw": sub[C_COL].quantile(0.25),
            "C_q75_raw": sub[C_COL].quantile(0.75),

            "median_n_ghms": sub["n_ghms"].median(),
        })

    out = pd.DataFrame(rows)

    for col in ["R_median_raw", "R_q25_raw", "R_q75_raw"]:
        out[col.replace("_raw", "_plot")] = out[col].clip(R_LIM[0], R_LIM[1])

    for col in ["C_median_raw", "C_q25_raw", "C_q75_raw"]:
        out[col.replace("_raw", "_plot")] = out[col].clip(C_LIM[0], C_LIM[1])

    return out


def basin_scale_summary(df, drivers):
    use = df[df["driver"].isin(drivers)].copy()
    use = use.dropna(subset=[BASIN_CLASS_COL, R_COL, C_COL]).copy()

    out = (
        use.groupby(["driver", BASIN_CLASS_COL], dropna=False)
        .agg(
            n_grid=("lat", "count"),

            R_median_raw=(R_COL, "median"),
            R_q25_raw=(R_COL, q25),
            R_q75_raw=(R_COL, q75),

            C_median_raw=(C_COL, "median"),
            C_q25_raw=(C_COL, q25),
            C_q75_raw=(C_COL, q75),

            area_median=("area_km2", "median"),
            log10_area_median=("log10_area", "median"),
            n_ghms_median=("n_ghms", "median"),
        )
        .reset_index()
    )

    for col in ["R_median_raw", "R_q25_raw", "R_q75_raw"]:
        out[col.replace("_raw", "_plot")] = out[col].clip(R_LIM[0], R_LIM[1])

    for col in ["C_median_raw", "C_q25_raw", "C_q75_raw"]:
        out[col.replace("_raw", "_plot")] = out[col].clip(C_LIM[0], C_LIM[1])

    out[BASIN_CLASS_COL] = pd.Categorical(
        out[BASIN_CLASS_COL],
        categories=BASIN_CLASS_ORDER,
        ordered=True
    )

    return out.sort_values(["driver", BASIN_CLASS_COL])


def quadrant_composition(df):
    use = df.copy()
    use["quadrant"] = [
        classify_quadrant(r, c)
        for r, c in zip(use[R_COL], use[C_COL])
    ]

    out = (
        use.groupby(["driver", "quadrant"], dropna=False)
        .size()
        .reset_index(name="n")
    )

    total = use.groupby("driver").size().reset_index(name="total")
    out = out.merge(total, on="driver", how="left")
    out["prop"] = out["n"] / out["total"]

    return out


def scale_quadrant_composition(df):
    use = df.copy()
    use = use.dropna(subset=["driver", BASIN_CLASS_COL, R_COL, C_COL]).copy()

    use["quadrant"] = [
        classify_quadrant(r, c)
        for r, c in zip(use[R_COL], use[C_COL])
    ]

    out = (
        use.groupby(["driver", BASIN_CLASS_COL, "quadrant"], dropna=False)
        .size()
        .reset_index(name="n")
    )

    total = (
        use.groupby(["driver", BASIN_CLASS_COL], dropna=False)
        .size()
        .reset_index(name="total")
    )

    out = out.merge(total, on=["driver", BASIN_CLASS_COL], how="left")
    out["prop"] = out["n"] / out["total"]

    out["driver"] = pd.Categorical(out["driver"], categories=ORDER, ordered=True)
    out[BASIN_CLASS_COL] = pd.Categorical(
        out[BASIN_CLASS_COL],
        categories=BASIN_CLASS_ORDER,
        ordered=True
    )

    return out.sort_values(["driver", BASIN_CLASS_COL, "quadrant"])


def dominant_mode_by_scale(comp_df):
    rows = []

    for d in ORDER:
        for s in BASIN_CLASS_ORDER:
            sub = comp_df[
                (comp_df["driver"] == d) &
                (comp_df[BASIN_CLASS_COL] == s)
            ].copy()

            if len(sub) == 0:
                rows.append({
                    "driver": d,
                    BASIN_CLASS_COL: s,
                    "dominant_mode": np.nan,
                    "dominant_prop": np.nan
                })
                continue

            sub = sub.sort_values(["prop", "quadrant"], ascending=[False, True]).reset_index(drop=True)

            rows.append({
                "driver": d,
                BASIN_CLASS_COL: s,
                "dominant_mode": sub.loc[0, "quadrant"],
                "dominant_prop": sub.loc[0, "prop"]
            })

    out = pd.DataFrame(rows)
    out["driver"] = pd.Categorical(out["driver"], categories=ORDER, ordered=True)
    out[BASIN_CLASS_COL] = pd.Categorical(
        out[BASIN_CLASS_COL],
        categories=BASIN_CLASS_ORDER,
        ordered=True
    )

    return out.sort_values(["driver", BASIN_CLASS_COL])


effect_sum = global_effect_summary(df)
warming_scale = basin_scale_summary(df, WARMING_ORDER)
human_scale = basin_scale_summary(df, HUMAN_ORDER)
quad_sum = quadrant_composition(df)
scale_quad_sum = scale_quadrant_composition(df)
scale_dom = dominant_mode_by_scale(scale_quad_sum)

effect_sum.to_csv(FIG_DIR / "final_global_effect_summary.csv", index=False, encoding="utf-8-sig")
warming_scale.to_csv(FIG_DIR / "final_warming_scale_dependence.csv", index=False, encoding="utf-8-sig")
human_scale.to_csv(FIG_DIR / "final_human_scale_dependence.csv", index=False, encoding="utf-8-sig")
quad_sum.to_csv(FIG_DIR / "final_joint_composition.csv", index=False, encoding="utf-8-sig")
scale_quad_sum.to_csv(FIG_DIR / "final_scale_quadrant_composition.csv", index=False, encoding="utf-8-sig")
scale_dom.to_csv(FIG_DIR / "final_scale_dominant_mode.csv", index=False, encoding="utf-8-sig")

# Metadata for the 1-beta conversion.
metadata_text = f"""Fig. 3 plotting code with C = 1 - beta conversion.

Input file:
{DF_BASIN_FILE}

R column:
{R_COL}

Original beta/slope delta column:
{BETA_DELTA_COL}

Converted C column:
{C_COL} = -{BETA_DELTA_COL}

Meaning:
C = 1 - beta
Delta_C = Delta(1 - beta) = -Delta_beta

Therefore:
beta increase -> C decrease
beta decrease -> C increase
"""

(FIG_DIR / "method_note_C_1minus_beta.txt").write_text(metadata_text, encoding="utf-8")
print(metadata_text)


# ============================================================
# 4. Plotting functions
# ============================================================

def plot_fingerprint_summary(ax, effect_sum):
    """
    Panel a:
    Global R-C fingerprint without arrows.
    Show median points and IQR crosses only.
    """
    for _, row in effect_sum.iterrows():
        d = row["driver"]

        x_plot = row["R_median_plot"]
        y_plot = row["C_median_plot"]

        ax.errorbar(
            x_plot,
            y_plot,
            xerr=[
                [x_plot - row["R_q25_plot"]],
                [row["R_q75_plot"] - x_plot]
            ],
            yerr=[
                [y_plot - row["C_q25_plot"]],
                [row["C_q75_plot"] - y_plot]
            ],
            fmt="o",
            markersize=6,
            color=COLORS[d],
            ecolor=COLORS[d],
            elinewidth=1.0,
            capsize=2.2,
            label=d,
            zorder=5
        )

    ax.axvline(0, color="0.35", lw=0.6, ls="--")
    ax.axhline(0, color="0.35", lw=0.6, ls="--")

    ax.set_xlim(R_LIM)
    ax.set_ylim(C_LIM)

    ax.set_xlabel("ΔR (%)")
    ax.set_ylabel("ΔC")
    ax.set_title("Global R–C fingerprint", fontsize=9)

    ax.grid(lw=0.3, alpha=0.35)
    ax.legend(frameon=False, fontsize=7, loc="best", handlelength=1.0)


def plot_effect_size_panel(ax_r, ax_c, effect_sum):
    ypos = np.arange(len(ORDER))[::-1]

    for y, d in zip(ypos, ORDER):
        row = effect_sum[effect_sum["driver"] == d].iloc[0]

        ax_r.plot(
            [row["R_q25_plot"], row["R_q75_plot"]],
            [y, y],
            color=COLORS[d],
            lw=3.0
        )
        ax_r.scatter(row["R_median_plot"], y, s=28, color=COLORS[d], zorder=5)

        ax_c.plot(
            [row["C_q25_plot"], row["C_q75_plot"]],
            [y, y],
            color=COLORS[d],
            lw=3.0
        )
        ax_c.scatter(row["C_median_plot"], y, s=28, color=COLORS[d], zorder=5)

    ax_r.axvline(0, color="0.35", lw=0.6, ls="--")
    ax_c.axvline(0, color="0.35", lw=0.6, ls="--")

    ax_r.set_xlim(R_LIM)
    ax_c.set_xlim(C_LIM)

    ax_r.set_yticks(ypos)
    ax_r.set_yticklabels([DISPLAY_LABELS[d] for d in ORDER])
    ax_c.set_yticks(ypos)
    ax_c.set_yticklabels([])

    ax_r.set_xlabel("ΔR (%)")
    ax_c.set_xlabel("ΔC")

    ax_r.set_title("Global ΔR", fontsize=9)
    ax_c.set_title("Global ΔC", fontsize=9)

    ax_r.grid(axis="x", lw=0.3, alpha=0.35)
    ax_c.grid(axis="x", lw=0.3, alpha=0.35)


def plot_scale_by_class(ax_r, ax_c, scale_df, drivers, title):
    x = np.arange(len(BASIN_CLASS_ORDER))
    offsets = np.linspace(-0.06, 0.06, len(drivers))

    for off, d in zip(offsets, drivers):
        sub = scale_df[scale_df["driver"] == d].copy()
        sub = sub.sort_values(BASIN_CLASS_COL)
        xpos = x + off

        ax_r.plot(
            xpos,
            sub["R_median_plot"],
            color=COLORS[d],
            lw=1.1,
            marker="o",
            markersize=4.5,
            label=d
        )
        ax_r.errorbar(
            xpos,
            sub["R_median_plot"],
            yerr=[
                sub["R_median_plot"] - sub["R_q25_plot"],
                sub["R_q75_plot"] - sub["R_median_plot"]
            ],
            fmt="none",
            ecolor=COLORS[d],
            elinewidth=0.8,
            capsize=2.0
        )

        ax_c.plot(
            xpos,
            sub["C_median_plot"],
            color=COLORS[d],
            lw=1.1,
            marker="o",
            markersize=4.5,
            label=d
        )
        ax_c.errorbar(
            xpos,
            sub["C_median_plot"],
            yerr=[
                sub["C_median_plot"] - sub["C_q25_plot"],
                sub["C_q75_plot"] - sub["C_median_plot"]
            ],
            fmt="none",
            ecolor=COLORS[d],
            elinewidth=0.8,
            capsize=2.0
        )

    for ax in [ax_r, ax_c]:
        ax.axhline(0, color="0.35", lw=0.6, ls="--")
        ax.set_xticks(x)
        ax.set_xticklabels(BASIN_XLABELS)
        ax.grid(axis="y", lw=0.3, alpha=0.35)

    ax_r.set_ylim(R_LIM)
    ax_c.set_ylim(C_LIM)

    ax_r.set_ylabel("ΔR (%)")
    ax_c.set_ylabel("ΔC")

    ax_r.set_title(title + ": ΔR", fontsize=9)
    ax_c.set_title(title + ": ΔC", fontsize=9)

    ax_r.set_xlabel("Model-basin size")
    ax_c.set_xlabel("Model-basin size")

    ax_c.legend(frameon=False, fontsize=7, loc="best")


def add_segment_labels(ax, y_positions, segment_lefts, segment_widths, segment_names, min_prop=0.08, fontsize=7):
    for i, y in enumerate(y_positions):
        for left, width, name in zip(segment_lefts[i], segment_widths[i], segment_names):
            if width >= min_prop:
                x = left + width / 2
                txt = f"{width*100:.0f}"
                color = "white" if name in DARK_SEGMENTS else "black"
                ax.text(
                    x, y, txt,
                    ha="center", va="center",
                    fontsize=fontsize, color=color
                )


def plot_stacked_composition(ax, summary_df, row_key, row_order, title, add_labels=True, min_label_prop=0.08, fontsize=7):
    ypos = np.arange(len(row_order))[::-1]
    lefts = np.zeros(len(row_order))

    segment_lefts_all = [[] for _ in row_order]
    segment_widths_all = [[] for _ in row_order]

    for q in QUAD_ORDER:
        vals = []

        for r in row_order:
            tmp = summary_df[
                (summary_df[row_key] == r) &
                (summary_df["quadrant"] == q)
            ]
            vals.append(float(tmp["prop"].iloc[0]) if len(tmp) else 0.0)

        ax.barh(
            ypos,
            vals,
            left=lefts,
            height=0.58,
            color=QUAD_COLORS[q],
            edgecolor="white",
            linewidth=0.35,
            label=q
        )

        for i in range(len(row_order)):
            segment_lefts_all[i].append(lefts[i])
            segment_widths_all[i].append(vals[i])

        lefts += np.array(vals)

    if add_labels:
        add_segment_labels(
            ax=ax,
            y_positions=ypos,
            segment_lefts=segment_lefts_all,
            segment_widths=segment_widths_all,
            segment_names=QUAD_ORDER,
            min_prop=min_label_prop,
            fontsize=fontsize
        )

    ax.set_xlim(0, 1)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Proportion of grid cells")

    return ypos


def plot_scale_dominant_mode_heatmap(ax, scale_dom):
    nrow = len(ORDER)
    ncol = len(BASIN_CLASS_ORDER)

    ax.set_xlim(0, ncol)
    ax.set_ylim(0, nrow)
    ax.invert_yaxis()

    for i, d in enumerate(ORDER):
        for j, s in enumerate(BASIN_CLASS_ORDER):
            sub = scale_dom[
                (scale_dom["driver"] == d) &
                (scale_dom[BASIN_CLASS_COL] == s)
            ]

            if len(sub) == 0 or pd.isna(sub["dominant_mode"].iloc[0]):
                facecolor = "white"
                txt = "NA"
                text_color = "black"
            else:
                mode = sub["dominant_mode"].iloc[0]
                prop = sub["dominant_prop"].iloc[0]

                facecolor = QUAD_COLORS.get(mode, "white")
                txt = f"{short_mode_label(mode)}\n{prop*100:.0f}%"
                text_color = "white" if mode in DARK_SEGMENTS else "black"

            rect = Rectangle(
                (j, i), 1, 1,
                facecolor=facecolor,
                edgecolor="white",
                linewidth=1.0
            )
            ax.add_patch(rect)

            ax.text(
                j + 0.5,
                i + 0.5,
                txt,
                ha="center",
                va="center",
                fontsize=7.2,
                color=text_color
            )

    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels(BASIN_XLABELS)

    ax.set_yticks(np.arange(nrow) + 0.5)
    ax.set_yticklabels([DISPLAY_LABELS[d] for d in ORDER])

    ax.tick_params(length=0)

    ax.set_xlabel("Model-basin size")
    ax.set_title("Dominant R–C mode by scale", fontsize=9)

    for spine in ax.spines.values():
        spine.set_visible(False)

    for x in range(ncol + 1):
        ax.axvline(x, color="0.85", lw=0.6, zorder=0)
    for y in range(nrow + 1):
        ax.axhline(y, color="0.85", lw=0.6, zorder=0)


# ============================================================
# 5. Draw main figure
# ============================================================

fig = plt.figure(figsize=(8.27, 11.69))

# More generous margins to avoid clipped labels, titles, panel letters and bottom legend.
gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    height_ratios=[1.03, 1.12, 1.02],
    width_ratios=[1, 1],
    left=0.115,
    right=0.965,
    top=0.955,
    bottom=0.185,
    wspace=0.34,
    hspace=0.50
)

# ------------------------------------------------------------
# a. Global fingerprint
# ------------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])
ax_a.text(-0.105, 1.045, "a", fontsize=11, fontweight="bold", transform=ax_a.transAxes)
plot_fingerprint_summary(ax_a, effect_sum)

# ------------------------------------------------------------
# b. Global effect size
# ------------------------------------------------------------
gs_b = gs[0, 1].subgridspec(1, 2, wspace=0.40)
ax_b1 = fig.add_subplot(gs_b[0, 0])
ax_b2 = fig.add_subplot(gs_b[0, 1])
ax_b1.text(-0.34, 1.045, "b", fontsize=11, fontweight="bold", transform=ax_b1.transAxes)
plot_effect_size_panel(ax_b1, ax_b2, effect_sum)

# ------------------------------------------------------------
# c. Warming scale dependence
# ------------------------------------------------------------
gs_c = gs[1, 0].subgridspec(1, 2, wspace=0.32)
ax_c1 = fig.add_subplot(gs_c[0, 0])
ax_c2 = fig.add_subplot(gs_c[0, 1])
ax_c1.text(-0.32, 1.045, "c", fontsize=11, fontweight="bold", transform=ax_c1.transAxes)
plot_scale_by_class(ax_c1, ax_c2, warming_scale, WARMING_ORDER, "Warming")

# ------------------------------------------------------------
# d. Human scale dependence
# ------------------------------------------------------------
gs_d = gs[1, 1].subgridspec(1, 2, wspace=0.32)
ax_d1 = fig.add_subplot(gs_d[0, 0])
ax_d2 = fig.add_subplot(gs_d[0, 1])
ax_d1.text(-0.32, 1.045, "d", fontsize=11, fontweight="bold", transform=ax_d1.transAxes)
plot_scale_by_class(ax_d1, ax_d2, human_scale, HUMAN_ORDER, "Human activity")

# ------------------------------------------------------------
# e. Global joint composition
# ------------------------------------------------------------
ax_e = fig.add_subplot(gs[2, 0])
ax_e.text(-0.105, 1.045, "e", fontsize=11, fontweight="bold", transform=ax_e.transAxes)

ypos_e = plot_stacked_composition(
    ax=ax_e,
    summary_df=quad_sum,
    row_key="driver",
    row_order=ORDER,
    title="Global joint R–C response modes",
    add_labels=True,
    min_label_prop=0.08,
    fontsize=7
)

ax_e.set_yticks(ypos_e)
ax_e.set_yticklabels([DISPLAY_LABELS[d] for d in ORDER])

# ------------------------------------------------------------
# f. Scale-dominant mode heatmap
# ------------------------------------------------------------
ax_f = fig.add_subplot(gs[2, 1])
ax_f.text(-0.105, 1.045, "f", fontsize=11, fontweight="bold", transform=ax_f.transAxes)

plot_scale_dominant_mode_heatmap(ax_f, scale_dom)

# ------------------------------------------------------------
# Shared legend
# ------------------------------------------------------------
handles = [
    mpl.patches.Patch(facecolor=QUAD_COLORS[q], edgecolor="none", label=q)
    for q in QUAD_ORDER
]

fig.legend(
    handles=handles,
    labels=QUAD_ORDER,
    frameon=False,
    fontsize=7,
    ncol=4,
    loc="lower center",
    bbox_to_anchor=(0.53, 0.118),
    handlelength=1.6,
    columnspacing=1.8
)


for ext in ["pdf", "png", "svg"]:
    out = FIG_DIR / f"Fig3_A4_final_clean_1minus_beta_layout_fixed.{ext}"
    if ext == "png":
        fig.savefig(out, dpi=PNG_DPI, bbox_inches="tight", pad_inches=0.08)
    else:
        fig.savefig(out, bbox_inches="tight", pad_inches=0.08)
    print("Saved:", out)

plt.show()

print("\nDone. Figure saved in:")
print(FIG_DIR)