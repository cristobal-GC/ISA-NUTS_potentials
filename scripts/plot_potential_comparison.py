from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

from typing import Any
snakemake: Any


############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
CF_threshold = float(snakemake.params["CF_threshold"])
##### input
file_gdf_nuts = snakemake.input["gdf_NUTS"]
file_df_summary = snakemake.input["df_summary"]
##### output
file_map_potential_comparison = snakemake.output["map_potential_comparison"]
##### wildcards
nuts = snakemake.wildcards["nuts"]
resource = snakemake.wildcards["resource"]
year = snakemake.wildcards["year"]


############################## Operations

df_summary = pd.read_csv(file_df_summary, index_col=0)
required_cols = [
    "CAPACITY_ISA4",
    "CAPACITY_CFth",
    "CAPACITY_CFth_ISA4",
    "area_ISA0",
    "area_ISA1",
    "area_ISA2",
    "area_ISA3",
    "area_ISA4",
]
missing_cols = [col for col in required_cols if col not in df_summary.columns]
if missing_cols:
    raise ValueError(
        f"Missing required columns in {file_df_summary}: {missing_cols}"
    )

# Potential density [MW/km2] = CAPACITY [MW] / total ISA area [km2].
area_total = df_summary[["area_ISA0", "area_ISA1", "area_ISA2", "area_ISA3", "area_ISA4"]].sum(axis=1)
denominator = area_total.where(area_total > 0)

df_summary["density_technical"] = df_summary["CAPACITY_ISA4"].div(denominator)
df_summary["density_economic"] = df_summary["CAPACITY_CFth"].div(denominator)
df_summary["density_techno_economic"] = df_summary["CAPACITY_CFth_ISA4"].div(denominator)

gdf_nuts = gpd.read_file(file_gdf_nuts).set_index("NUTS_ID")
if "LEVL_CODE" not in gdf_nuts.columns:
    raise ValueError(f"Missing required column 'LEVL_CODE' in {file_gdf_nuts}")

level_map = {"NUTS2": 2, "NUTS3": 3}
if nuts not in level_map:
    raise ValueError(f"Invalid nuts wildcard for this rule: {nuts}. Expected NUTS2 or NUTS3.")

gdf_level = gdf_nuts.loc[gdf_nuts["LEVL_CODE"] == level_map[nuts]].copy()
for density_col in ["density_technical", "density_economic", "density_techno_economic"]:
    gdf_level[density_col] = df_summary[density_col]


############################## Create outputs

plt.rcParams.update(
    {
        "axes.titlesize": 14,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
    }
)
fontsize = 14

crs = ccrs.PlateCarree()
fig, axes = plt.subplots(1, 3, figsize=(13, 6), subplot_kw={"projection": crs})

plots = [
    ("density_technical", "Technical potential\n(ISA 4)"),
    ("density_economic", f"Economic potential\n(CF>={CF_threshold:g})"),
    ("density_techno_economic", f"Techno-economic potential\n(ISA 4 and CF>={CF_threshold:g})"),
]

has_any_data = False

level_xmin, level_ymin, level_xmax, level_ymax = gdf_level.total_bounds

for ax, (density_col, title) in zip(axes, plots):
    gdf_plot = gdf_level.loc[gdf_level[density_col].notna()].copy()

    if gdf_plot.empty:
        ax.set_extent([level_xmin, level_xmax, level_ymin, level_ymax], crs=crs)
        ax.set_title(f"{title}\n(No data)", fontsize=fontsize)
        ax.coastlines(resolution="10m", color="black", linewidth=0.2)
        ax.set_xticks([])
        ax.set_yticks([])
        continue

    has_any_data = True
    xmin, ymin, xmax, ymax = gdf_plot.total_bounds
    dx = xmax - xmin
    dy = ymax - ymin
    pad_x = max(dx * 0.03, 1e-6)
    pad_y = max(dy * 0.03, 1e-6)

    ax.set_extent([xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y], crs=crs)
    gdf_plot.plot(
        column=density_col,
        cmap="viridis",
        linewidth=0.2,
        edgecolor="black",
        ax=ax,
    )
    mappable = ax.collections[0]
    cbar = fig.colorbar(mappable, ax=ax, orientation="horizontal", pad=0.05)
    cbar.set_label("MW/km2")
    ax.coastlines(resolution="10m", color="black", linewidth=0.2)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xticks([])
    ax.set_yticks([])

if not has_any_data:
    raise ValueError(
        f"No regions with available potential density data found for {nuts} in {file_df_summary}."
    )

plt.subplots_adjust(wspace=0.3)
plt.tight_layout()

Path(file_map_potential_comparison).parent.mkdir(parents=True, exist_ok=True)
dpi = fig_params["sizes"]["LR"]["dpi"]
fig.savefig(file_map_potential_comparison, bbox_inches="tight", pad_inches=0.2, dpi=dpi)
plt.close(fig)
