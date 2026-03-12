from pathlib import Path
import re

import numpy as np
import pandas as pd
import xarray as xr
from utils import log_xarray_spatial_info

from typing import Any
snakemake: Any


ISA_LEVELS = [0, 1, 2, 3, 4]


############################## Unwrap relevant variables

##### input
files_nc_CAPACITY = list(snakemake.input["ncs_CAPACITY"])
##### params
cap_per_sqkm = float(snakemake.params["cap_per_sqkm"])
##### output
file_df_CAPACITY = snakemake.output["df_CAPACITY"]


############################## Operations

pattern = re.compile(r"CAPACITY_ISA(?P<isa>\d+)_")
capacity_by_isa = {isa: 0.0 for isa in ISA_LEVELS}

for file_path in files_nc_CAPACITY:
    match = pattern.search(Path(file_path).name)
    if not match:
        raise ValueError(f"Could not parse ISA index from file name: {file_path}")

    isa = int(match.group("isa"))
    if isa not in capacity_by_isa:
        continue

    da = xr.open_dataarray(file_path)
    try:
        log_xarray_spatial_info(da, source_label=file_path)
        capacity_by_isa[isa] = round(float(np.nansum(da.values)), 6)
    finally:
        da.close()

rows = []
for isa in ISA_LEVELS:
    capacity = capacity_by_isa[isa]
    area = round(capacity / cap_per_sqkm, 6)
    rows.append({"value": isa, "CAPACITY": capacity, "area": area})

total_capacity = round(sum(row["CAPACITY"] for row in rows), 6)
total_area = round(sum(row["area"] for row in rows), 6)

for row in rows:
    row["porc"] = round((row["area"] / total_area) * 100, 6) if total_area > 0 else 0.0

rows.append({
    "value": "TOTAL",
    "CAPACITY": total_capacity,
    "area": total_area,
    "porc": 100.0 if total_area > 0 else 0.0,
})

df = pd.DataFrame(rows, columns=["value", "CAPACITY", "area", "porc"])


############################## Create outputs

Path(file_df_CAPACITY).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(file_df_CAPACITY, index=False)
