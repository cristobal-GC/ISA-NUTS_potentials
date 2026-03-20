import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
##### input
file_df_summary = snakemake.input["df_summary"]
##### output
file_plot = snakemake.output["plot_venn"]
##### wildcards
region = snakemake.wildcards["region"]


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


############################## Load data

df_summary = pd.read_csv(file_df_summary, index_col=0)

if region not in df_summary.index:
    raise ValueError(f"Region {region} not found in summary file {file_df_summary}")

column_capacity_isa4 = get_first_existing_column(df_summary, ["CAPACITY_ISA4", "CAPACITY_ISA_4"])
column_capacity_cfth = get_first_existing_column(df_summary, ["CAPACITY_CFth", "CAPACITY_CFTH"])
column_capacity_overlap = get_first_existing_column(df_summary, ["CAPACITY_CFth_ISA4", "CAPACITY_CFth_ISA_4", "CAPACITY_CFTH_ISA4"])

capacity_isa4 = max(float(df_summary.loc[region, column_capacity_isa4]), 0.0)
capacity_cfth = max(float(df_summary.loc[region, column_capacity_cfth]), 0.0)
capacity_overlap = max(float(df_summary.loc[region, column_capacity_overlap]), 0.0)

capacity_overlap = min(capacity_overlap, capacity_isa4, capacity_cfth)

if capacity_isa4 == 0.0 and capacity_cfth == 0.0:
    raise ValueError(f"Both capacities are zero for region {region}; Venn diagram is undefined.")


############################## Build plot

color_isa4 = "#1a9850"
color_cfth = "#E7E409"
color_overlap = "#b9d541"

alpha = 0.6

radius_1 = math.sqrt(capacity_isa4 / math.pi)
radius_2 = math.sqrt(capacity_cfth / math.pi)
distance = find_distance_for_overlap(capacity_isa4, capacity_cfth, capacity_overlap)

fig, ax = plt.subplots(figsize=(7, 5))

circle_isa4 = plt.Circle((0.0, 0.0), radius_1, color=color_isa4, alpha=alpha, ec="black", lw=0.8)
circle_cfth = plt.Circle((distance, 0.0), radius_2, color=color_cfth, alpha=alpha, ec="black", lw=0.8)

ax.add_artist(circle_isa4)
ax.add_artist(circle_cfth)

x_min = min(-radius_1, distance - radius_2)
x_max = max(radius_1, distance + radius_2)
y_max = max(radius_1, radius_2)

x_pad = max((x_max - x_min) * 0.12, 1e-6)
y_pad = max(y_max * 0.18, 1e-6)

ax.set_xlim(x_min - x_pad, x_max + x_pad)
ax.set_ylim(-y_max - y_pad * 0.2, y_max + y_pad)
ax.set_aspect("equal", "box")
ax.axis("off")

fontsize = fig_params["sizes"]["LR"]["fontsize"]

legend_isa4 = mpatches.Patch(color=color_isa4, alpha=alpha, label=r"$P_{_{\mathrm{ISA4}}}$")
legend_cfth = mpatches.Patch(color=color_cfth, alpha=alpha, label=r"$P^{^{\mathit{CF}^*}}$")
legend_overlap = mpatches.Patch(color=color_overlap, alpha=1.0, label=r"$P_{_{\mathrm{ISA4}}}^{^{\mathit{CF}^*}}$")

fig.legend(
    handles=[legend_isa4, legend_cfth, legend_overlap],
    loc="upper right",
    fontsize=fontsize,
    frameon=False,
    bbox_to_anchor=(1.02, 0.96),
    labelspacing=1.0,
)

plt.tight_layout(rect=[0, 0, 1, 0.95])


############################## Save output

fig.savefig(file_plot, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)
