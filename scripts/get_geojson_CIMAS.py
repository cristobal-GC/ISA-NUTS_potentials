from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

from typing import Any
snakemake: Any  # defined by Snakemake at runtime


##############################
# This script builds data/NUTS/CIMAS.geojson from config['CIMAS_domains'].
# Each domain becomes one feature with a rectangular Polygon derived from its
# lat_min/lat_max/lon_min/lon_max bounds. The output follows the same schema as
# the official NUTS GeoJSON (NUTS_ID, LEVL_CODE, geometry, ...) and is written
# in EPSG:4326 so the rest of the workflow can consume it transparently.


# Custom level code for CIMAS domains. It only needs to be distinct and shared
# by all CIMAS features so boundary plotting (which filters siblings by
# LEVL_CODE) draws all domains together.
CIMAS_LEVL_CODE = 9


def _flatten_geometry(raw_geometry):
    """Accept either a flat mapping {lat_min: ...} or a list of single-key
    dicts [{lat_min: ...}, {lat_max: ...}] (the config draft format)."""
    if isinstance(raw_geometry, dict):
        return dict(raw_geometry)

    if isinstance(raw_geometry, list):
        merged = {}
        for item in raw_geometry:
            if not isinstance(item, dict):
                raise ValueError(
                    f"Unexpected geometry entry (expected mapping): {item!r}"
                )
            merged.update(item)
        return merged

    raise ValueError(
        f"Unsupported geometry format: {type(raw_geometry).__name__}. "
        "Expected a mapping or a list of single-key mappings."
    )


def _build_rectangle(bounds, domain_id):
    required = ("lat_min", "lat_max", "lon_min", "lon_max")
    missing = [k for k in required if k not in bounds]
    if missing:
        raise ValueError(
            f"Domain '{domain_id}' is missing geometry bounds {missing}. "
            f"Required: {list(required)}."
        )

    lat_min = float(bounds["lat_min"])
    lat_max = float(bounds["lat_max"])
    lon_min = float(bounds["lon_min"])
    lon_max = float(bounds["lon_max"])

    if lat_min >= lat_max:
        raise ValueError(f"Domain '{domain_id}': lat_min must be < lat_max.")
    if lon_min >= lon_max:
        raise ValueError(f"Domain '{domain_id}': lon_min must be < lon_max.")
    if not (-90 <= lat_min < lat_max <= 90):
        raise ValueError(f"Domain '{domain_id}': latitude out of [-90, 90].")
    if not (-180 <= lon_min < lon_max <= 180):
        raise ValueError(f"Domain '{domain_id}': longitude out of [-180, 180].")

    # box(minx, miny, maxx, maxy) -> (lon_min, lat_min, lon_max, lat_max)
    return box(lon_min, lat_min, lon_max, lat_max)


############################## Unwrap relevant variables

cimas_domains = snakemake.params["cimas_domains"]
file_geojson = snakemake.output["geojson"]


############################## Operations

if not cimas_domains:
    raise ValueError(
        "config['CIMAS_domains'] is empty; cannot build CIMAS GeoJSON."
    )

records = []
for domain_id, spec in cimas_domains.items():
    spec = spec or {}
    nuts_id = str(spec.get("NUTS_ID", domain_id))
    bounds = _flatten_geometry(spec.get("geometry", {}))
    polygon = _build_rectangle(bounds, domain_id)

    records.append(
        {
            "NUTS_ID": nuts_id,
            "LEVL_CODE": int(spec.get("LEVL_CODE", CIMAS_LEVL_CODE)),
            "CNTR_CODE": str(spec.get("CNTR_CODE", "")),
            "NAME_LATN": str(spec.get("NAME_LATN", nuts_id)),
            "NUTS_NAME": str(spec.get("NUTS_NAME", spec.get("NAME_LATN", nuts_id))),
            "geometry": polygon,
        }
    )

gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")

print(
    f"[get_geojson_CIMAS] Built {len(gdf)} CIMAS domain(s): "
    f"{list(gdf['NUTS_ID'])} | CRS: {gdf.crs}"
)


############################## Create outputs

Path(file_geojson).parent.mkdir(parents=True, exist_ok=True)
gdf.to_file(file_geojson, driver="GeoJSON")

print(f"[get_geojson_CIMAS] Wrote {file_geojson}")
