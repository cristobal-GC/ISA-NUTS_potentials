
from pathlib import Path


########## Set default config file if it exists

if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"



##### Include rules
include: "rules/getting.smk",
include: "rules/plotting.smk"
include: "rules/retrieving.smk"




##### Retrieve relevant information from config file

REGIONS = config["regions"]
RESOURCES = config["resources"]
CUTOUTS = config["cutouts"]
YEARS = config["years"]
FORMATS = config["formats"]
RESOLUTIONS = config["resolutions"]



rule all:
    input:
        expand("results/maps/ISA/{resolution}/{format}/ISA_local_{resource}_{region}_{resolution}.{format}", resource=RESOURCES, region=REGIONS, resolution=RESOLUTIONS, format=FORMATS),
        expand("results/maps/cutout/{cutout}/{format}/cutout_{resource}_{region}_{cutout}_{year}.{format}", cutout=CUTOUTS, year=YEARS, resource=RESOURCES, region=REGIONS, format=FORMATS) 



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
