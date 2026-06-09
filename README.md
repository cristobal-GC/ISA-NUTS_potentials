# ISA-NUTS_potentials

Regional assessment of onshore wind power potential in low environmental
sensitivity areas of Spain.

## Overview

The workflow estimates onshore wind power potential across Spanish NUTS regions
(NUTS0, NUTS2 and NUTS3), restricted to the lowest environmental sensitivity
class (ISA-4) of the ISA index published by MITECO. For each region it:

1. Builds the wind capacity factor (CF) from a reanalysis cutout (e.g. ERA5,
   NEWA) using [atlite](https://github.com/PyPSA/atlite) and a reference wind
   turbine power curve.
2. Selects land with low environmental sensitivity (ISA-4) and CF above a
   configurable threshold.
3. Computes the installable capacity and the associated energy potential, and
   compares them with the actual installed capacity and electricity demand.
4. Produces maps, Venn diagrams and per-region LaTeX summary sheets.

## Requirements

The environment is managed with [pixi](https://pixi.sh):

```bash
pixi install
```

## Usage

The pipeline is orchestrated with Snakemake. To build all outputs:

```bash
pixi run snakemake all --cores <N>
```

Regions, resources, cutouts, years and analysis parameters (turbine, CF
threshold, capacity density, etc.) are configured in
[`config/config.yaml`](config/config.yaml).

## Data

The ISA environmental sensitivity rasters (from MITECO) and the GEBCO bathymetry
are downloaded automatically by dedicated Snakemake rules
([`rules/retrieve.smk`](rules/retrieve.smk)); the downloaded files and all
generated results live outside version control (see `.gitignore`).

The only inputs that must be provided locally are the reanalysis cutouts
referenced in [`config/config.yaml`](config/config.yaml) (`cutout_params`),
which point to machine-specific paths and should be adjusted to your setup.

## License

Released under the MIT License (see [`LICENSE`](LICENSE)).
