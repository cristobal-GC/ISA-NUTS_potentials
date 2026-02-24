
from pathlib import Path


########## Set default config file if it exists

if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"



##### Include rules
include: "rules/plotting.smk"
include: "rules/getting.smk",



##### Retrieve relevant information from config file

REGIONS = config["regions"]
RESOURCES = config["resources"]
FORMATS = config["formats"]
RESOLUTIONS = config["resolutions"]



##### Define function to get the NUTS level (e.g. "NUTS2") from the region code (e.g. "ES11")
# Is that useful? current NUTS file ships all the levels

def nuts_from_region(region):

    levels = {
        0: 'NUTS0',
        2: 'NUTS2',
        3: 'NUTS3',
    }

    n_digits = sum(c.isdigit() for c in region)

    try:
        return levels[n_digits]
    except KeyError:
        raise ValueError(f"NUTS level not supported for region: {region}")






rule all:
    input:
#        expand("results/rasters/ISA/ISA_local_{resource}_{region}.tiff", resource=RESOURCES, region=REGIONS),
        expand("results/maps/ISA/{resolution}/{format}/ISA_local_{resource}_{region}_{resolution}.{format}", resource=RESOURCES, region=REGIONS, resolution=RESOLUTIONS, format=FORMATS)




# rule pattern
#     params:
#         param1=config["field1"]
#         param2=config["field2"]
#         param3=lambda wc: config["field"][wc.sample]  <  param that depends on wildcard
#     input:
#         label_input=file_input
#     output:
#         label_output=file_output
#     script:
#         "path_to_script.py"



# To use within a python script:
#
#   snakemake.wildcards["sample"]  <  where {sample} is the wildcard
#   snakemake.params["label"]
#   snakemake.input["label_input"]
#   snakemake.output["label_output"]
#   snakemake.config  <  not recommended, better to use params, except for global params
#   snakemake.log


##### Special functions
#
# workflow.source_path() > to get paths relative to root snakefile
#                        > useful for pointing at scripts in rules/ folder
# workflow.basedir > absolute path to Snakefile 
