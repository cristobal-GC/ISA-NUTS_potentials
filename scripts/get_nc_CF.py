import yaml
from utils import load_and_limit_cutout, load_gdf_nuts, resolve_user_home_path

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

##### Load gdf_NUTS local (one region)
_, gdf_NUTS_local = load_gdf_nuts(file_gdf_NUTS, region)


##### Load and limit cutout 
file_cutout = resolve_user_home_path(cutout_params[f"{cutout}_{year}"]["path"])
c = load_and_limit_cutout(file_cutout, gdf_NUTS_local)


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
