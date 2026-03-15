import atlite
import geopandas as gpd
import logging
import numpy as np
from pathlib import Path
import xarray as xr


logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)


def _select_spatial_dims_xarray(data):
    known_x = {"x", "lon", "longitude", "easting", "eastings"}
    known_y = {"y", "lat", "latitude", "northing", "northings"}

    dims = list(data.dims)
    selected = []
    for dim in dims:
        dim_lower = dim.lower()
        if dim_lower in known_x or dim_lower in known_y:
            selected.append(dim)

    if len(selected) < 2:
        for coord_name in data.coords:
            coord = data.coords[coord_name]
            axis = str(coord.attrs.get("axis", "")).upper()
            if axis in {"X", "Y"} and coord_name in dims and coord_name not in selected:
                selected.append(coord_name)

    if len(selected) < 2:
        non_spatial_guess_exclude = {"time", "month", "year", "step", "band", "variable"}
        for dim in dims:
            if dim.lower() not in non_spatial_guess_exclude and dim not in selected:
                selected.append(dim)
            if len(selected) >= 2:
                break

    if len(selected) < 2:
        for dim in dims:
            if dim not in selected:
                selected.append(dim)
            if len(selected) >= 2:
                break

    return selected[:2]


def _extract_xarray_crs(data):
    for key in ["crs", "spatial_ref"]:
        if key in data.attrs:
            return data.attrs[key]

    grid_mapping_name = data.attrs.get("grid_mapping")
    if grid_mapping_name and grid_mapping_name in data.coords:
        coord = data.coords[grid_mapping_name]
        for key in ["spatial_ref", "crs_wkt", "proj4_params", "crs"]:
            if key in coord.attrs:
                return coord.attrs[key]

    return "unknown"


def log_xarray_spatial_info(data, source_label):
    spatial_dims = _select_spatial_dims_xarray(data)
    crs = _extract_xarray_crs(data)

    details = []
    for dim in spatial_dims:
        n_values = int(data.sizes[dim]) if dim in data.sizes else None
        if dim in data.coords:
            coord_values = data.coords[dim].values
            if coord_values.size > 0 and np.issubdtype(np.asarray(coord_values).dtype, np.number):
                dim_min = float(np.nanmin(coord_values))
                dim_max = float(np.nanmax(coord_values))
            else:
                dim_min = "unknown"
                dim_max = "unknown"
        else:
            dim_min = "unknown"
            dim_max = "unknown"

        details.append((dim, n_values, dim_min, dim_max))

    _log_and_print(f"[spatial-log] source={source_label} | type=xarray | crs={crs}")
    _log_and_print(f"[spatial-log] source={source_label} | spatial_dims={[dim for dim, _, _, _ in details]}")
    for dim, n_values, dim_min, dim_max in details:
        _log_and_print(
            f"[spatial-log] source={source_label} | dim={dim} | n={n_values} | min={dim_min} | max={dim_max}"
        )


def log_raster_spatial_info(raster, source_label):
    if raster.crs is not None and raster.crs.is_geographic:
        x_name, y_name = "lon", "lat"
    else:
        x_name, y_name = "x", "y"

    left, bottom, right, top = raster.bounds

    _log_and_print(f"[spatial-log] source={source_label} | type=raster | crs={raster.crs}")
    _log_and_print(f"[spatial-log] source={source_label} | spatial_dims={[x_name, y_name]}")
    _log_and_print(f"[spatial-log] source={source_label} | dim={x_name} | n={raster.width} | min={left} | max={right}")
    _log_and_print(f"[spatial-log] source={source_label} | dim={y_name} | n={raster.height} | min={bottom} | max={top}")


def _subtract_excluded_geometries(gdf_local, gdf_all):
    """
    Rebuild ES geometry by intersecting with Spanish NUTS1 union excluding ES7.

    This is used to remove non-mainland/islet artifacts that are not fully handled
    by subtracting a fixed set of excluded sub-regions.

    Parameters
    ----------
    gdf_local : geopandas.GeoDataFrame
        Geometry of the requested region (one or more rows).
    gdf_all : geopandas.GeoDataFrame
        Full NUTS GeoDataFrame (indexed by NUTS_ID), used to look up the
        excluded geometries.

    Returns
    -------
    geopandas.GeoDataFrame
        gdf_local with ES geometry intersected with NUTS1 union (excluding ES7).
    """
    from shapely.ops import unary_union

    # Apply only when the requested local geometry corresponds to ES.
    if "ES" not in gdf_local.index:
        return gdf_local

    # Spanish NUTS1 IDs are of the form ES[0-9] (e.g. ES1 ... ES7).
    # Excluding ES7 (Canarias) keeps mainland + remaining continental subregions.
    ids_to_keep = [
        nid
        for nid in gdf_all.index.astype(str)
        if nid.startswith("ES") and len(nid) == 3 and nid[2].isdigit() and nid != "ES7"
    ]
    if not ids_to_keep:
        return gdf_local

    keep_geom = unary_union(gdf_all.loc[ids_to_keep].geometry)

    gdf_local = gdf_local.copy()
    gdf_local["geometry"] = gdf_local.geometry.intersection(keep_geom)

    _log_and_print(
        f"[_subtract_excluded_geometries] Intersected ES with NUTS1 excluding ES7: {ids_to_keep}"
    )

    return gdf_local


def resolve_user_home_path(path_value):
    path = Path(path_value).expanduser()

    if path.is_absolute():
        return str(path)

    return str(Path.home() / path)



def load_gdf_nuts(file_gdf_NUTS, region):
    
    ##### Load full gdf_NUTS and select local region
    gdf_NUTS = (
        gpd.read_file(file_gdf_NUTS)
        .set_index("NUTS_ID")
    )

    _log_and_print(
        f"[load_gdf_nuts] Loaded gdf_NUTS, CRS: {gdf_NUTS.crs}"
    )

    ##### Filter local region. 
    gdf_NUTS_local = gdf_NUTS.loc[[region]]
    
    _log_and_print(
        f"[load_gdf_nuts] Filtered local region: {region}"
    )

    # For ES, also apply the geometry subtraction to remove non-mainland/islet artifacts.
    if region == 'ES':
        gdf_NUTS_local = _subtract_excluded_geometries(gdf_NUTS_local, gdf_NUTS)
        # Log
        _log_and_print(
            f"[load_gdf_nuts] Applied geometry subtraction for ES to remove non-mainland/islet artifacts."
        )

    return gdf_NUTS, gdf_NUTS_local



def load_and_limit_cutout(file_cutout, gdf_local):

    c = atlite.Cutout(file_cutout)

    log_xarray_spatial_info(c.data, source_label=f"cutout:{file_cutout}")

    _log_and_print(
        f"[load_and_limit_cutout] Cutout loaded. CRS is {c.crs}."
    )


    if gdf_local.crs is None:
        raise ValueError("gdf_local must have a defined CRS.")


    if gdf_local.crs != c.crs:
        
        _log_and_print(
            f"[load_and_limit_cutout] Mismatch between gdf_local CRS({gdf_local.crs}) and cutout CRS ({c.crs}). Reprojecting gdf_local to cutout CRS."
        )

        gdf_local = gdf_local.to_crs(c.crs)
     
        

    xmin, ymin, xmax, ymax = gdf_local.total_bounds

    margin_x = c.dx
    margin_y = c.dy

    xmin -= margin_x
    xmax += margin_x
    ymin -= margin_y
    ymax += margin_y

    return c.sel(bounds=(xmin, ymin, xmax, ymax))



def plot_dataarray_on_map(
    data,
    gdf_NUTS,
    gdf_NUTS_local,
    region,
    file_output,
    x_coord="lon",
    y_coord="lat",
    cmap="viridis",
    cbar_label="Value",
    vmin=None,
    vmax=None,
    size=12,
    linewidth=1.5,
    fontsize=16,
    title=None,
    bounds_type="gdf"
):
    """
    Generic function to plot a 2D xarray DataArray on a map with NUTS boundaries.
    
    Parameters
    ----------
    data : xr.DataArray
        2D data array to plot
    gdf_NUTS : geopandas.GeoDataFrame
        Full NUTS boundaries GeoDataFrame
    gdf_NUTS_local : geopandas.GeoDataFrame
        Local region NUTS boundary
    region : str
        Region code (e.g., 'ES11')
    file_output : str
        Output file path
    x_coord : str, optional
        Name of x-coordinate ('lon' or 'x'), by default "lon"
    y_coord : str, optional
        Name of y-coordinate ('lat' or 'y'), by default "lat"
    cmap : str, optional
        Colormap name, by default "viridis"
    cbar_label : str, optional
        Colorbar label, by default "Value"
    vmin : float, optional
        Minimum value for colormap, by default None (auto)
    vmax : float, optional
        Maximum value for colormap, by default None (auto)
    size : int, optional
        Figure size (square), by default 12
    linewidth : float, optional
        NUTS boundary linewidth, by default 1.5
    fontsize : int, optional
        Font size for labels, by default 16
    title : str, optional
        Plot title, by default None
    bounds_type : str, optional
        How to set axis limits: "gdf" (from NUTS boundaries) or "data" (from data coordinates),
        by default "gdf"
    """
    import matplotlib.pyplot as plt
    
    ##### Setup plot parameters
    fig, ax = plt.subplots(figsize=(size, size))
    
    ##### Plot field without automatic colorbar
    plot_kwargs = {
        "ax": ax,
        "x": x_coord,
        "y": y_coord,
        "cmap": cmap,
        "add_colorbar": False,
    }
    
    # Add vmin/vmax if specified
    if vmin is not None:
        plot_kwargs["vmin"] = vmin
    if vmax is not None:
        plot_kwargs["vmax"] = vmax
    
    mappable = data.plot(**plot_kwargs)
    
    ##### Configure axes
    ax.tick_params(axis="both", labelsize=fontsize*0.8)
    ax.set_xlabel(x_coord.capitalize(), fontsize=fontsize*0.8)
    ax.set_ylabel(y_coord.capitalize(), fontsize=fontsize*0.8)
    
    ##### Add colorbar
    cbar = fig.colorbar(
        mappable,
        ax=ax,
        orientation='vertical',
        fraction=0.046,
        pad=0.04
    )
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
    ##### Add NUTS boundaries
    # All NUTS at same level as region
    gdf_NUTS[gdf_NUTS.index.astype(str).str.len() == len(region)].plot(
        ax=ax, color="none", edgecolor='grey', linewidth=linewidth
    )
    
    # Local region with double linewidth
    gdf_NUTS_local.plot(
        ax=ax, color="none", edgecolor='black', linewidth=linewidth*2
    )
    
    ##### Set axis limits
    if bounds_type == "gdf":
        # Use GDF bounds (for geographic coordinates)
        xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds
        km_per_lon = 85
        km_per_lat = 111
        center_x = (xmax + xmin) / 2
        center_y = (ymax + ymin) / 2
        delta_x = xmax - xmin
        delta_y = ymax - ymin
        delta_km = max([km_per_lon*delta_x, km_per_lat*delta_y]) * 1.02
        ax.set_xlim(
            center_x - 0.5*delta_km/km_per_lon,
            center_x + 0.5*delta_km/km_per_lon
        )
        ax.set_ylim(
            center_y - 0.5*delta_km/km_per_lat,
            center_y + 0.5*delta_km/km_per_lat
        )
    elif bounds_type == "data":
        # Use data coordinate bounds
        xmin, xmax = data[x_coord].min().values, data[x_coord].max().values
        ymin, ymax = data[y_coord].min().values, data[y_coord].max().values
        delta_x = xmax - xmin
        delta_y = ymax - ymin
        margin = 0.02
        ax.set_xlim(xmin - margin*delta_x, xmax + margin*delta_x)
        ax.set_ylim(ymin - margin*delta_y, ymax + margin*delta_y)
    
    ##### Add title if provided, otherwise clear any auto-generated title
    if title:
        ax.set_title(title, fontsize=fontsize)
    else:
        ax.set_title("")  # Clear any auto-generated title from xarray
    
    ##### Save figure
    fig.savefig(file_output, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    
    _log_and_print(f"[plot_dataarray_on_map] Saved plot to {file_output}")



