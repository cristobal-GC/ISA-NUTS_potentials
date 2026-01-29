
import geopandas as gpd
import matplotlib.pyplot as plt



########## Unwrap relevant variables
gdf_file    = snakemake.input["gdf"]
map_file    = snakemake.output["map"]
region      = snakemake.wildcards["region"]



########## Operations

gdf = (gpd.read_file(gdf_file)
          .set_index("id")
          .loc[[region]]
)




########## Make plot

fig, ax = plt.subplots(figsize=(10,10))



gdf.plot(ax=ax, color="none", edgecolor='grey', linewidth=3)

ax.set_title(region)

plt.tight_layout()

fig.savefig(map_file,
            bbox_inches="tight",
            )  # pad=0 para que ocupe todo el cuadrado
        

plt.close(fig)

