import math
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from typing import Any
snakemake: Any


plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "mathtext.rm": "serif",
    "mathtext.it": "serif:italic",
    "mathtext.bf": "serif:bold",
})


############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
regions = list(snakemake.params["regions"])
CF_threshold = float(snakemake.params["CF_threshold"])
##### input
file_gdf_nuts = snakemake.input["gdf_NUTS"]
file_df_summary = snakemake.input["df_summary"]
##### output
file_plot = snakemake.output["plot_venn_all"]
##### wildcards
nuts = snakemake.wildcards["nuts"]


############################## Helpers

def intersection_area(radius_1, radius_2, distance):
    if distance >= radius_1 + radius_2:
        return 0.0

    if distance <= abs(radius_1 - radius_2):
        return math.pi * min(radius_1, radius_2) ** 2

    part_1 = radius_1 ** 2 * math.acos((distance ** 2 + radius_1 ** 2 - radius_2 ** 2) / (2.0 * distance * radius_1))
    part_2 = radius_2 ** 2 * math.acos((distance ** 2 + radius_2 ** 2 - radius_1 ** 2) / (2.0 * distance * radius_2))
    part_3 = 0.5 * math.sqrt(
        (-distance + radius_1 + radius_2)
        * (distance + radius_1 - radius_2)
        * (distance - radius_1 + radius_2)
        * (distance + radius_1 + radius_2)
    )

    return part_1 + part_2 - part_3


def find_distance_for_overlap(capacity_1, capacity_2, overlap_capacity):
    radius_1 = math.sqrt(capacity_1 / math.pi)
    radius_2 = math.sqrt(capacity_2 / math.pi)

    if overlap_capacity <= 0:
        return radius_1 + radius_2

    if overlap_capacity >= min(capacity_1, capacity_2):
        return abs(radius_1 - radius_2)

    lower = abs(radius_1 - radius_2)
    upper = radius_1 + radius_2

    for _ in range(80):
        middle = (lower + upper) / 2.0
        area_middle = intersection_area(radius_1, radius_2, middle)

        if area_middle > overlap_capacity:
            lower = middle
        else:
            upper = middle

    return (lower + upper) / 2.0


def get_first_existing_column(df, candidates):
    for column in candidates:
        if column in df.columns:
            return column
    raise ValueError(f"None of expected columns found: {candidates}")


def compute_grid_layout(n_regions, nuts):
    # Reserve two slots for the legend.
    n_slots = n_regions + 2

    if nuts == "NUTS2":
        return 3, 6

    n_cols = max(2, math.ceil(math.sqrt(1.5 * n_slots)))
    n_rows = math.ceil(n_slots / n_cols)
    return n_rows, n_cols


############################## Load data

df_summary = pd.read_csv(file_df_summary, index_col=0)

gdf_nuts = gpd.read_file(file_gdf_nuts).set_index("NUTS_ID")
name_col = "NUTS_NAME" if "NUTS_NAME" in gdf_nuts.columns else "NAME_LATN"

column_capacity_isa4 = get_first_existing_column(df_summary, ["CAPACITY_ISA4", "CAPACITY_ISA_4"])
column_capacity_cfth = get_first_existing_column(df_summary, ["CAPACITY_CFth", "CAPACITY_CFTH"])
column_capacity_overlap = get_first_existing_column(
    df_summary, ["CAPACITY_CFth_ISA4", "CAPACITY_CFth_ISA_4", "CAPACITY_CFTH_ISA4"]
)

regions_present = [r for r in regions if r in df_summary.index]
if not regions_present:
    raise ValueError(
        f"None of the configured regions for {nuts} are present in {file_df_summary}"
    )

region_capacities = {}
max_capacity = 0.0
for region in regions_present:
    capacity_isa4 = max(float(df_summary.loc[region, column_capacity_isa4]), 0.0)
    capacity_cfth = max(float(df_summary.loc[region, column_capacity_cfth]), 0.0)
    capacity_overlap = max(float(df_summary.loc[region, column_capacity_overlap]), 0.0)
    capacity_overlap = min(capacity_overlap, capacity_isa4, capacity_cfth)
    region_capacities[region] = (capacity_isa4, capacity_cfth, capacity_overlap)
    max_capacity = max(max_capacity, capacity_isa4, capacity_cfth)

if max_capacity <= 0.0:
    raise ValueError(f"All capacities are zero for {nuts}; Venn diagrams are undefined.")


############################## Build plot

color_isa4 = "#1a9850"
color_cfth = "#E7E409"
color_overlap = "#b9d541"
alpha = 0.6

n_rows, n_cols = compute_grid_layout(len(regions_present), nuts)

fontsize = fig_params["sizes"]["LR"]["fontsize"]
dpi = fig_params["sizes"]["LR"]["dpi"]

cell_size = 3.2  # inches per subplot
fig = plt.figure(figsize=(n_cols * cell_size, n_rows * cell_size))
gs = fig.add_gridspec(n_rows, n_cols)

title_fontsize = fontsize * 0.65
legend_fontsize = fontsize * 0.6


def draw_venn(ax, capacity_1, capacity_2, capacity_overlap, title):
    ax.axis("off")
    ax.set_aspect("equal", "box")

    if capacity_1 <= 0.0 and capacity_2 <= 0.0:
        ax.set_title(f"{title}\n(No data)", fontsize=title_fontsize)
        return

    radius_1 = math.sqrt(capacity_1 / math.pi) if capacity_1 > 0 else 0.0
    radius_2 = math.sqrt(capacity_2 / math.pi) if capacity_2 > 0 else 0.0
    distance = find_distance_for_overlap(capacity_1, capacity_2, capacity_overlap)

    center_1_x = -distance / 2.0
    center_2_x = distance / 2.0

    if radius_1 > 0:
        ax.add_artist(plt.Circle(
            (center_1_x, 0.0), radius_1, color=color_isa4, alpha=alpha, ec="black", lw=0.5
        ))
    if radius_2 > 0:
        ax.add_artist(plt.Circle(
            (center_2_x, 0.0), radius_2, color=color_cfth, alpha=alpha, ec="black", lw=0.5
        ))

    x_min = center_1_x - radius_1
    x_max = center_2_x + radius_2
    y_max = max(radius_1, radius_2)
    half_span = max(x_max - x_min, 2.0 * y_max) / 2.0
    cx = (x_min + x_max) / 2.0
    pad = 0.08 * half_span
    ax.set_xlim(cx - half_span - pad, cx + half_span + pad)
    ax.set_ylim(-half_span - pad, half_span + pad)

    ax.set_title(title, fontsize=title_fontsize)


# Legend spans two slots so the labels have room to breathe.
legend_ax = fig.add_subplot(gs[0, 0:2])
legend_ax.axis("off")
legend_isa4 = mpatches.Patch(
    color=color_isa4, alpha=alpha,
    label=r"$P_{_{\mathrm{ISA4}}}$: Technical potential (ISA 4)",
)
legend_cfth = mpatches.Patch(
    color=color_cfth, alpha=alpha,
    label=fr"$P^{{^{{\mathit{{CF}}^*}}}}$: Economic potential ($CF \geq {CF_threshold:g}$)",
)
legend_overlap = mpatches.Patch(
    color=color_overlap, alpha=1.0,
    label=(
        fr"$P_{{_{{\mathrm{{ISA4}}}}}}^{{^{{\mathit{{CF}}^*}}}}$: "
        fr"Techno-economic potential (ISA 4 and $CF \geq {CF_threshold:g}$)"
    ),
)
legend_ax.legend(
    handles=[legend_isa4, legend_cfth, legend_overlap],
    loc="center",
    fontsize=legend_fontsize,
    frameon=False,
    labelspacing=1.4,
    handlelength=1.5,
)


# Regions fill the remaining slots in row-major order, starting after the legend.
slot_idx = 2
for region in regions_present:
    row, col = divmod(slot_idx, n_cols)
    if row >= n_rows:
        break

    ax = fig.add_subplot(gs[row, col])
    capacity_isa4, capacity_cfth, capacity_overlap = region_capacities[region]
    region_name = gdf_nuts.loc[region, name_col] if region in gdf_nuts.index else region
    title = f"{region_name} ({region})"
    draw_venn(ax, capacity_isa4, capacity_cfth, capacity_overlap, title)
    slot_idx += 1


############################## Save output

Path(file_plot).parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout()
fig.savefig(file_plot, bbox_inches="tight", pad_inches=0.15, dpi=dpi)
plt.close(fig)
