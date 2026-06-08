# Make two RC figures with consistent colors

This package generates the two target composite figures:

1. `Fig1_RC_trends_map_bar_climate.*`
2. `Fig2_area_threshold_moving_window.*`

The confluence metric is treated as:

```python
C = 1 - beta
```

So the exported confluence trend directions are already converted from beta to C.

## Color convention

The palette is fixed in `02_plot_two_figures.py`:

- Significant increase: dark teal `#2B8C7D`
- Non-significant increase: light teal `#A9D6CE`
- Significant decrease: dark brown `#B6783F`
- Non-significant decrease: light brown `#E6CFB4`
- Runoff generation panels: brown `#B6783F`
- Confluence panels: teal `#2B8C7D`
- Combination bars: light/dark teal

This keeps map signs, R/C process panels, and the bar chart visually consistent.

## Usage

Edit the paths at the top of `run_all_make_two_figures.py` or the two component scripts:

```python
MAIN_CSV = r"C:\\Users\\huawei\\Downloads\\filtered_station_summary_1960_2015_cleaned_coords.csv"
CLIMATE_SHP = r"D:\\workroom\\GPLW\\气候区划\\气候区划.shp"
# CLIMATE_SHP = None

CSV_DIR = Path(r"./RC_figure_csvs_1minus_beta")
FIG_DIR = Path(r"./RC_final_figures_1minus_beta")
```

Then run:

```bash
python run_all_make_two_figures.py
```

Outputs are saved as PNG, PDF, and SVG.


## Latest adjustment
- Map longitude/latitude tick labels are enlarged via `MAP_TICK_LABEL_SIZE = 9.4`.
- Station points now use a very thin boundary: `MAP_EDGE_COLOR = "#3A3A3A"`, `MAP_EDGE_LINEWIDTH = 0.08`.
