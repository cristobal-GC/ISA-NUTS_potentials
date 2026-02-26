import geopandas as gpd
import atlite
import yaml

from typing import Any
snakemake: Any  # This is to avoid my IDE to complain about snakemake variable not being defined, but it is actually defined when running the script with snakemake



############################## Unwrap relevant variables

##### params
cutout_params = snakemake.params["cutout_params"]
CF_params = snakemake.params["CF_params"]
##### input
file_gdf_NUTS = snakemake.input["gdf_NUTS"]
##### output
file_nc_CF = snakemake.output["nc_CF"]
##### wildcards
cutout = snakemake.wildcards["cutout"]
year =snakemake.wildcards["year"]
region = snakemake.wildcards["region"]
resource = snakemake.wildcards["resource"]



############################## Operations

##### Load cutout 
file_cutout = cutout_params[f"{cutout}_{year}"]["path"]
c = atlite.Cutout(file_cutout)

##### Load gdf_NUTS and apply operations
gdf_NUTS = (
    gpd.read_file(file_gdf_NUTS)
    .set_index("NUTS_ID")            # set index 
    .to_crs(c.crs)                   # change crs to that of the cutout
)

##### Filter gdf with only one nuts region    
gdf_NUTS_local = gdf_NUTS.loc[[region]]


##### limit cutout
# Get region bounding box
xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds

# Add margin of one cell
margin_x = c.dx
margin_y = c.dy

xmin -= margin_x
xmax += margin_x
ymin -= margin_y
ymax += margin_y

# Limit cutout
c = c.sel(bounds=(xmin, ymin, xmax, ymax))


##### Obtain CF matrix
if resource=='onwind':
    # Load wind turbine data
    with open(f"data/windturbine/{CF_params[resource]['turbine']}.yaml") as f:
        turbine = yaml.safe_load(f)
    # Compute CF with atlite    
    CF = c.wind(
        turbine=turbine,
        capacity_factor=True  # capacity_factor_timeseries=True is for hCF
    )
elif resource=='solar':
    CF = c.pv(
        panel=CF_params[resource]['panel'],
        orientation=CF_params[resource]['orientation'],
        tracking=None,
        capacity_factor=True
    )
else:
    raise ValueError(f"Resource: {resource} not recognized. Check CF_params for valid resources.")

##### Apply correction factor
CF = CF_params[resource]['correction_factor'] * CF



############################## Create outputs
CF.to_netcdf(file_nc_CF)
