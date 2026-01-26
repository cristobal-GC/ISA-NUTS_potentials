
from pathlib import Path


########## Set default config file if it exists

if Path("config/config.yaml").exists():

    configfile: "config/config.yaml"



##### Include rules
include:
   "rules/plotting.smk",



##### Retrieve relevant information from config file

REGION = config["region"]




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
#         label_ourput=file_output
#     script:
#         "path_to_script.py"



# To use within a python script:
#
#   snakemake.wildcards["sample"]  <  where {sample} is the wildcard
#   snakemake.params["label"]
#   snakemake.input["label"]
#   snakemake.output["label"]
#   snakemake.config  <  not recommended, better to use params, except for global params
#   snakemake.log


##### Special functions
#
# workflow.source_path() > to get paths relative to root snakefile
#                        > useful for pointing at scripts in rules/ folder
# workflow.basedir > absolute path to Snakefile 
