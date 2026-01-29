
from pathlib import Path


########## Set default config file if it exists

if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"



##### Include rules
include:
   "rules/plotting.smk",



##### Retrieve relevant information from config file

REGION = config["region"]




##### Define function to get the NUTS level from the region (useful to point at the correct NUTS file)

def nuts_from_region(region):

    levels = {
        #0: 'NUTS0',
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
        expand("results/figs/{region}.pdf", region=REGION)


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
