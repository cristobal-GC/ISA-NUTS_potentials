from pathlib import Path
from collections import defaultdict

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from typing import Any
snakemake: Any


############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
file_gdf_nuts = snakemake.input["gdf_NUTS"]
##### output
file_map_nuts = snakemake.output["map_NUTS"]
##### wildcards
nuts = snakemake.wildcards["nuts"]


############################## Operations

level_map = {"NUTS2": 2, "NUTS3": 3}
if nuts not in level_map:
    raise ValueError(f"Invalid nuts wildcard for this rule: {nuts}. Expected NUTS2 or NUTS3.")

gdf_nuts = gpd.read_file(file_gdf_nuts)
required_cols = {"NUTS_ID", "LEVL_CODE", "geometry"}
missing = required_cols.difference(gdf_nuts.columns)
if missing:
    raise ValueError(f"Missing required columns in {file_gdf_nuts}: {sorted(missing)}")

gdf_level = gdf_nuts.loc[gdf_nuts["LEVL_CODE"] == level_map[nuts]].copy()

# Exclude Ceuta, Melilla, Canarias and associated NUTS3 regions.
excluded_prefixes = ("ES63", "ES64", "ES70")
gdf_level = gdf_level.loc[
    ~gdf_level["NUTS_ID"].astype(str).str.startswith(excluded_prefixes)
].copy()

# Additional name-based safety filter if name columns are present.
name_column = "NUTS_NAME" if "NUTS_NAME" in gdf_level.columns else None
if name_column is None and "NAME_LATN" in gdf_level.columns:
    name_column = "NAME_LATN"

if name_column is not None:
    gdf_level = gdf_level.loc[
        ~gdf_level[name_column].astype(str).str.contains("Ceuta|Melilla|Canarias", case=False, regex=True)
    ].copy()

if gdf_level.empty:
    raise ValueError(f"No geometries found for {nuts} in {file_gdf_nuts}.")

gdf_level = gdf_level.set_index("NUTS_ID", drop=False)
gdf_level["xy"] = gdf_level.index.str[2:].astype(int)
gdf_level = gdf_level.sort_values("xy")

def lighten_color(color, amount):
    # Blend base color with white. amount=0 keeps the original color.
    r, g, b = mcolors.to_rgb(color)
    return (
        r + (1.0 - r) * amount,
        g + (1.0 - g) * amount,
        b + (1.0 - b) * amount,
        1.0,
    )


if nuts == "NUTS2":
    n_regions = len(gdf_level)
    cmap = cm.get_cmap("tab20", n_regions)
    gdf_level["color"] = [cmap(i) for i in range(n_regions)]
else:
    gdf_level["parent_nuts2"] = gdf_level.index.str[:4]
    parent_codes = sorted(gdf_level["parent_nuts2"].unique())
    parent_cmap = cm.get_cmap("tab20", len(parent_codes))
    parent_base_color = {code: parent_cmap(i) for i, code in enumerate(parent_codes)}

    children_by_parent = defaultdict(list)
    for nuts3_code in gdf_level.index:
        children_by_parent[nuts3_code[:4]].append(nuts3_code)

    color_map = {}
    for parent_code, children in children_by_parent.items():
        children_sorted = sorted(children)
        n_children = len(children_sorted)

        for i, child_code in enumerate(children_sorted):
            if n_children == 1:
                amount = 0.18
            else:
                # Slight gradation to keep all children visually tied to the parent color.
                amount = 0.08 + 0.40 * (i / (n_children - 1))

            color_map[child_code] = lighten_color(parent_base_color[parent_code], amount)

    gdf_level["color"] = gdf_level.index.map(color_map)


############################## Create outputs

plt.rcParams.update(
    {
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "legend.title_fontsize": 14,
    }
)

crs = ccrs.PlateCarree()
fig, ax = plt.subplots(1, 1, figsize=(8, 10), subplot_kw={"projection": crs})

ax.add_feature(cfeature.COASTLINE, color="black", linewidth=0.2)

gdf_level.plot(
    ax=ax,
    color=gdf_level["color"],
    edgecolor="black",
    linewidth=0.5,
)

# Add small numeric labels inside each region.
label_fontsize = 10 if nuts == "NUTS2" else 6
for _, row in gdf_level.iterrows():
    point = row["geometry"].representative_point()
    ax.text(
        point.x,
        point.y,
        str(int(row["xy"])),
        ha="center",
        va="center",
        fontsize=label_fontsize,
        color="black",
        transform=crs,
    )

# Zoom to selected regions with a small margin.
xmin, ymin, xmax, ymax = gdf_level.total_bounds
dx = xmax - xmin
dy = ymax - ymin
pad_x = max(dx * 0.03, 1e-6)
pad_y = max(dy * 0.03, 1e-6)
ax.set_extent([xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y], crs=crs)

legend_name_col = "NUTS_NAME" if "NUTS_NAME" in gdf_level.columns else "NUTS_ID"
legend_elements = [
    Patch(
        facecolor=row["color"],
        edgecolor="black",
        label=f"{idx} ({row[legend_name_col]})",
    )
    for idx, row in gdf_level.iterrows()
]

legend_title = "NUTS 2 regions" if nuts == "NUTS2" else "NUTS 3 regions"
ncol_legend = 3 if nuts == "NUTS3" else 2
ax.legend(
    handles=legend_elements,
    title=legend_title,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.08),
    ncol=ncol_legend,
)

#ax.set_title(f"{legend_title} (mainland Spain)")
ax.set_axis_off()

Path(file_map_nuts).parent.mkdir(parents=True, exist_ok=True)
plt.tight_layout()
fig.savefig(file_map_nuts, bbox_inches="tight", pad_inches=0.2)
plt.close(fig)
