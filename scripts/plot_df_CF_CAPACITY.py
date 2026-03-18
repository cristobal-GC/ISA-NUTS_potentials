import pandas as pd
import matplotlib.pyplot as plt

from typing import Any
snakemake: Any


############################## Unwrap relevant variables

##### params
fig_params = snakemake.params["fig_params"]
##### input
files_df_CF_CAPACITY = snakemake.input["dfs_CF_CAPACITY"]
##### output
file_plot = snakemake.output["plot_CF_CAPACITY"]
##### wildcards
year = snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]


############################## Operations

##### Load dataframes for all ISA codes
dic_CF_CAPACITY = {}
for i, file_path in enumerate(files_df_CF_CAPACITY):
    dic_CF_CAPACITY[i] = pd.read_csv(file_path)

##### Get ISA colors and labels from config
colors_contraste = fig_params["ISA"]["colors"]
labels = fig_params["ISA"]["labels"]


############################## Create plot

size = fig_params["sizes"]["LR"]["size"]
fontsize = fig_params["sizes"]["LR"]["fontsize"]

fig, ax = plt.subplots(figsize=(size, size/2))

##### Plot stairs for each ISA code
for kk, vv in dic_CF_CAPACITY.items():
    edges = [0] + list(vv["CUM_CAPACITY"].div(1000))
    values = vv['CF']
    
    plt.stairs(
        edges=edges,
        values=values,
        fill=False,
        label=f'{kk}: {labels[kk]}',
        color=colors_contraste[kk]
    )

##### Configure axes and legend
ax.legend(title='ISA code', fontsize=fontsize*0.8)

xmin = 0
xmax = max(d['CAPACITY'].sum() for d in dic_CF_CAPACITY.values()) / 1000 * 1.02
ymin = min(d["CF"].min() for d in dic_CF_CAPACITY.values()) * 0.98
ymax = max(d["CF"].max() for d in dic_CF_CAPACITY.values()) * 1.02

ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xlabel('Cumulative Capacity (GW)', fontsize=fontsize)
ax.set_ylabel('Capacity Factor', fontsize=fontsize)
ax.tick_params(axis='both', labelsize=fontsize)

plt.grid(True, which="both", linestyle="--", linewidth=1., color="gray", alpha=0.3)

##### Save figure
fig.savefig(file_plot, bbox_inches="tight", pad_inches=0.2)
plt.close(fig)
