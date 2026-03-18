import atlite
import geopandas as gpd
import logging
import numpy as np
from pathlib import Path
import xarray as xr
from pyproj import CRS



logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)


def _select_spatial_dims_xarray(data):
    """
    Identify spatial dimensions in an xarray Dataset/DataArray.
    Returns two dimensions (x, y) if possible.
    """

    dims = list(data.dims)

    # 1. Prefer standard geographic names
    x_candidates = {"x", "lon", "longitude", "easting"}
    y_candidates = {"y", "lat", "latitude", "northing"}

    x_dim = next((d for d in dims if d.lower() in x_candidates), None)
    y_dim = next((d for d in dims if d.lower() in y_candidates), None)

    if x_dim and y_dim:
        return [x_dim, y_dim]

    # 2. CF axis attribute
    for coord_name in data.coords:
        coord = data.coords[coord_name]
        axis = str(coord.attrs.get("axis", "")).upper()

        if axis == "X" and coord_name in dims:
            x_dim = coord_name
        if axis == "Y" and coord_name in dims:
            y_dim = coord_name

    if x_dim and y_dim:
        return [x_dim, y_dim]

    # 3. Fallback: ignore clearly non-spatial dimensions
    exclude = {"time", "month", "year", "step", "band", "variable"}

    spatial_guess = [d for d in dims if d.lower() not in exclude]

    if len(spatial_guess) >= 2:
        return spatial_guess[:2]

    # 4. Last resort
    return dims[:2]



def _extract_xarray_crs(data):
    """
    Extract CRS from xarray Dataset/DataArray and return a pyproj.CRS object.
    """
    crs_keys = ["crs", "spatial_ref", "crs_wkt", "proj4_params", "epsg"]

    def _try_crs_from_attrs(attrs):
        if not attrs:
            return None

        for key in crs_keys:
            if key not in attrs:
                continue

            value = attrs[key]
            if value is None or value == "":
                continue

            try:
                if key == "epsg":
                    return CRS.from_epsg(int(value))
                return CRS.from_user_input(value)
            except Exception:
                continue

        return None

    def _get_named_mapping_variable(mapping_name):
        if mapping_name in data.coords:
            return data.coords[mapping_name]
        if isinstance(data, xr.Dataset) and mapping_name in data.data_vars:
            return data[mapping_name]
        return None

    crs = _try_crs_from_attrs(getattr(data, "attrs", {}))
    if crs is not None:
        return crs

    grid_mapping_name = getattr(data, "attrs", {}).get("grid_mapping")
    if grid_mapping_name:
        mapping_var = _get_named_mapping_variable(grid_mapping_name)
        if mapping_var is not None:
            crs = _try_crs_from_attrs(getattr(mapping_var, "attrs", {}))
            if crs is not None:
                return crs

    if isinstance(data, xr.Dataset):
        for var in data.data_vars.values():
            crs = _try_crs_from_attrs(getattr(var, "attrs", {}))
            if crs is not None:
                return crs

            grid_mapping_name = var.attrs.get("grid_mapping")
            if not grid_mapping_name:
                continue

            mapping_var = _get_named_mapping_variable(grid_mapping_name)
            if mapping_var is None:
                continue

            crs = _try_crs_from_attrs(getattr(mapping_var, "attrs", {}))
            if crs is not None:
                return crs

    for coord in data.coords.values():
        crs = _try_crs_from_attrs(getattr(coord, "attrs", {}))
        if crs is not None:
            return crs

    return None



def log_xarray_spatial_info(data, source_label):
    """
    Log spatial metadata of an xarray Dataset or DataArray.
    """

    _log_and_print('[log_xarray_spatial_info] Extracting spatial metadata:')


    spatial_dims = _select_spatial_dims_xarray(data)
    crs = _extract_xarray_crs(data)

    # ----------------------------
    # Determine CRS type
    # ----------------------------
    if crs is None:
        crs_type = "undefined"
        units = "unknown"
    elif crs.is_geographic:
        crs_type = "geographic"
        units = "degrees"
    else:
        crs_type = "projected"
        units = "meters"

    _log_and_print(
        f"[log_xarray_spatial_info] source={source_label} | type=xarray | crs={crs} | crs_type={crs_type}"
    )

    _log_and_print(
        f"[log_xarray_spatial_info] source={source_label} | spatial_dims={spatial_dims} | units={units}"
    )

    # ----------------------------
    # Dimension info
    # ----------------------------
    for dim in spatial_dims:

        n_values = int(data.sizes[dim]) if dim in data.sizes else None
        dim_min = "unknown"
        dim_max = "unknown"
        res = "unknown"

        if dim in data.coords:

            coord = data.coords[dim].values

            if coord.size > 0 and np.issubdtype(np.asarray(coord).dtype, np.number):

                dim_min = float(np.nanmin(coord))
                dim_max = float(np.nanmax(coord))

                if coord.size > 1:
                    res = float(np.nanmedian(np.diff(coord)))


        # Format resolution for logging (it could be "unknown")
        res_str = f"{res:.2f}" if isinstance(res, float) else res

        _log_and_print(
            f"[log_xarray_spatial_info] source={source_label} | dim={dim} | n={n_values} | "
            f"min={dim_min} | max={dim_max} | res={res_str}"
        )

    # ----------------------------
    # dtype
    # ----------------------------
    try:
        dtype = str(data.dtype)
        _log_and_print(
            f"[log_xarray_spatial_info] source={source_label} | dtype={dtype}"
        )
    except AttributeError:
        pass



def log_raster_spatial_info(raster, source_label):
    """
    Log spatial metadata of a raster dataset.

    Parameters
    ----------
    raster : rasterio.DatasetReader
        Open raster dataset
    source_label : str
        Identifier of the raster source (e.g., filepath)
    """


    _log_and_print('[log_raster_spatial_info] Extracting spatial metadata:')


    crs = raster.crs

    # Determine CRS type and coordinate naming
    if crs is None:
        crs_type = "undefined"
        x_name, y_name = "x", "y"
        units = "unknown"

    elif crs.is_geographic:
        crs_type = "geographic"
        x_name, y_name = "longitude", "latitude"
        units = "degrees"

    else:
        crs_type = "projected"
        x_name, y_name = "easting", "northing"
        units = "meters"

    # Basic spatial info
    left, bottom, right, top = raster.bounds
    resx, resy = raster.res
    transform = raster.transform

    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | type=raster | crs={crs} | crs_type={crs_type}"
    )

    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | spatial_dims={[x_name, y_name]} | units={units}"
    )

    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | dim={x_name} | n={raster.width} | "
        f"min={left} | max={right} | res={resx}"
    )

    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | dim={y_name} | n={raster.height} | "
        f"min={bottom} | max={top} | res={resy}"
    )

    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | transform={transform}"
    )

    # Optional extra metadata
    _log_and_print(
        f"[log_raster_spatial_info] source={source_label} | dtype={raster.dtypes[0]} | nodata={raster.nodata}"
    )

    # ----------------------------
    # Simple sanity checks
    # ----------------------------

    # Check resolution sign (north-up rasters normally have negative y resolution)
    if resy > 0:
        _log_and_print(
            f"[log_raster_spatial_info] WARNING: source={source_label} | positive y resolution (unexpected orientation)"
        )

    # Check geographic bounds plausibility
    if crs and crs.is_geographic:
        if not (-180 <= left <= 180 and -180 <= right <= 180):
            _log_and_print(
                f"[log_raster_spatial_info] WARNING: source={source_label} | longitude bounds outside [-180,180]"
            )

        if not (-90 <= bottom <= 90 and -90 <= top <= 90):
            _log_and_print(
                f"[log_raster_spatial_info] WARNING: source={source_label} | latitude bounds outside [-90,90]"
            )

    # Check for very large or very small pixel sizes
    if abs(resx) > 10000 or abs(resy) > 10000:
        _log_and_print(
            f"[log_raster_spatial_info] WARNING: source={source_label} | unusually large pixel size"
        )

    if abs(resx) < 1e-6 or abs(resy) < 1e-6:
        _log_and_print(
            f"[log_raster_spatial_info] WARNING: source={source_label} | unusually small pixel size"
        )



def _subtract_excluded_geometries(gdf_local, gdf_all):
    """
    Rebuild ES geometry by subtracting excluded Spanish non-mainland regions.

    This avoids creating a large union geometry in memory.

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
        gdf_local with ES geometry after removing excluded subregions.
    """
    # Apply only when the requested local geometry corresponds to ES.
    if "ES" not in gdf_local.index:
        return gdf_local

    excluded_ids = ["ES7", "ES63", "ES64"]
    available_excluded_ids = [nid for nid in excluded_ids if nid in gdf_all.index]

    if not available_excluded_ids:
        _log_and_print(
            "[_subtract_excluded_geometries] No excluded Spanish subregions found in gdf_all. Returning original ES geometry."
        )
        return gdf_local

    gdf_local = gdf_local.copy()

    for excluded_id in available_excluded_ids:
        excluded_geom = gdf_all.loc[excluded_id, "geometry"]
        gdf_local["geometry"] = gdf_local.geometry.difference(excluded_geom)

    _log_and_print(
        f"[_subtract_excluded_geometries] Subtracted excluded geometries from ES: {available_excluded_ids}"
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
        f"[load_gdf_nuts] gdf_NUTS loaded gdf_NUTS. CRS: {gdf_NUTS.crs}"
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

    """This function loads the cutout and limits it to the bounding box of the local region. It also checks and handles CRS mismatches between the cutout and the local region geometry."""

    ##### Load cutout
    c = atlite.Cutout(file_cutout)

    _log_and_print(
        f"[load_and_limit_cutout] Cutout loaded. CRS: {c.crs}."
    )

    log_xarray_spatial_info(c.data, source_label=f"cutout:{file_cutout}")    


    if gdf_local.crs is None:
        raise ValueError("[load_and_limit_cutout] gdf_local must have a defined CRS.")


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



def load_CF(file_nc_CF):

    """This function loads the CF file."""

    ##### Load CF
    CF = xr.open_dataarray(file_nc_CF)

    _log_and_print(f"[load_CF] CF loaded. CF does not have CRS")

    log_xarray_spatial_info(CF, source_label=file_nc_CF)

    return CF



def load_CAPACITY(file_nc_CAPACITY):

    """This function loads the CAPACITY file."""

    ##### Load CAPACITY
    CAPACITY = xr.open_dataarray(file_nc_CAPACITY)

    _log_and_print(f"[load_CAPACITY] CAPACITY loaded. CAPACITY does not have CRS")

    log_xarray_spatial_info(CAPACITY, source_label=file_nc_CAPACITY)

    return CAPACITY



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
    # Add gdf for regions with the same NUTS code with thin grey lines
    gdf_NUTS[gdf_NUTS['LEVL_CODE'] == gdf_NUTS_local['LEVL_CODE'].iloc[0]].plot(
        ax=ax, color="none", edgecolor='grey', linewidth=linewidth
    )
    
    # Add gdf_local with double linewidth
    gdf_NUTS_local.plot(
        ax=ax, color="none", edgecolor='black', linewidth=linewidth*2
    )
    
    ##### Set axis limits
    if bounds_type == "gdf":
        # Use GDF bounds (for geographic coordinates)
        xmin, ymin, xmax, ymax = gdf_NUTS_local.total_bounds
        center_lat = (ymax + ymin) / 2
        km_per_lat = 111
        km_per_lon = 111 * np.cos(np.deg2rad(center_lat))
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



