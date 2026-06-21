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

    ##### Load full gdf and pick the identifier column.
    # The official NUTS GeoJSON and the CIMAS GeoJSON identify regions by
    # 'NUTS_ID'; the geoBoundaries ADM3 (municipalities) file uses 'shapeName'
    # and lacks the NUTS columns (NUTS_ID, LEVL_CODE, NUTS_NAME, ...). Detect
    # whichever identifier column is present so the rest of the workflow is
    # agnostic to the geometry source.
    gdf = gpd.read_file(file_gdf_NUTS)
    if "NUTS_ID" in gdf.columns:
        id_col = "NUTS_ID"
    elif "shapeName" in gdf.columns:
        id_col = "shapeName"
    else:
        raise ValueError(
            f"[load_gdf_nuts] No known identifier column in {file_gdf_NUTS}. "
            "Expected 'NUTS_ID' or 'shapeName'."
        )

    gdf_NUTS = gdf.set_index(id_col)

    _log_and_print(
        f"[load_gdf_nuts] gdf loaded from {file_gdf_NUTS}. id_col: {id_col}. CRS: {gdf_NUTS.crs}"
    )

    ##### Filter local region.
    # Output filenames use the region id verbatim, so it must not contain spaces;
    # for municipalities (shapeName) the id uses underscores in place of spaces
    # (e.g. "El_Toboso"). Convert back to the real name to match the geometry.
    lookup = region.replace("_", " ") if id_col == "shapeName" else region
    gdf_NUTS_local = gdf_NUTS.loc[[lookup]]

    if len(gdf_NUTS_local) > 1:
        raise ValueError(
            f"[load_gdf_nuts] Region '{region}' matches {len(gdf_NUTS_local)} "
            f"geometries in {file_gdf_NUTS} (ambiguous identifier '{id_col}'). "
            "Use a unique identifier."
        )

    _log_and_print(
        f"[load_gdf_nuts] Filtered local region: {region}"
    )

    # For ES, also apply the geometry subtraction to remove Canarias, Ceuta and Melilla
    if id_col == "NUTS_ID" and region == 'ES':
        gdf_NUTS_local = _subtract_excluded_geometries(gdf_NUTS_local, gdf_NUTS)
        # Log
        _log_and_print(
            "[load_gdf_nuts] Applied geometry subtraction for ES to remove Canarias, Ceuta and Melilla."
        )

    return gdf_NUTS, gdf_NUTS_local



def load_context_boundaries(file_nuts, nuts, clip_to=None):
    """Reference administrative boundaries drawn as thin grey context on a map.

    The region/domain being analysed is always drawn with a thick black outline;
    this returns the lighter background layer drawn behind it:
      - NUTS0 / NUTS2 / NUTS3: the other regions of the SAME NUTS level
        (from the official NUTS GeoJSON).
      - CIMAS: the NUTS3 regions (provinces). A CIMAS domain is a custom
        rectangle with no sibling regions of its own, so NUTS3 boundaries give
        geographic context within the domain.
      - ADM3: the neighbouring municipalities (the whole municipalities GeoJSON;
        that file has no LEVL_CODE, every feature is a municipality).

    Parameters
    ----------
    file_nuts : str
        GeoJSON providing the context geometries (the official NUTS file for
        NUTS/CIMAS, the municipalities file for ADM3). Returned in its CRS
        (EPSG:4326).
    nuts : str
        Level wildcard.
    clip_to : geopandas.GeoDataFrame, optional
        If given, only context geometries intersecting this geometry's bounding
        box (expanded by its own size on each side) are returned. Mainly an
        optimisation for ADM3, where the file holds thousands of municipalities
        but only the local neighbourhood is ever in view.
    """
    levl_by_nuts = {"NUTS0": 0, "NUTS2": 2, "NUTS3": 3, "CIMAS": 3}

    gdf = gpd.read_file(file_nuts)

    if nuts == "ADM3":
        # Municipalities file: no LEVL_CODE; every feature is itself the level.
        context = gdf
    elif nuts in levl_by_nuts:
        context = gdf[gdf["LEVL_CODE"] == levl_by_nuts[nuts]]
    else:
        raise ValueError(
            f"[load_context_boundaries] Unknown nuts level '{nuts}'. "
            f"Expected one of {sorted(levl_by_nuts) + ['ADM3']}."
        )

    if clip_to is not None and len(context):
        xmin, ymin, xmax, ymax = clip_to.to_crs(context.crs).total_bounds
        dx = xmax - xmin
        dy = ymax - ymin
        context = context.cx[xmin - dx:xmax + dx, ymin - dy:ymax + dy]

    return context



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



def load_GEBCO(file_nc_GEBCO):

    """This function loads the GEBCO file as a DataArray."""

    ds_gebco = xr.open_dataset(file_nc_GEBCO)

    if "elevation" in ds_gebco.data_vars:
        gebco = ds_gebco["elevation"]
    else:
        data_vars = list(ds_gebco.data_vars)

        if len(data_vars) != 1:
            raise ValueError(
                f"[load_GEBCO] Could not infer GEBCO variable from {file_nc_GEBCO}. Found variables: {data_vars}"
            )

        gebco = ds_gebco[data_vars[0]]

    _log_and_print(f"[load_GEBCO] GEBCO loaded from {file_nc_GEBCO}")

    log_xarray_spatial_info(gebco, source_label=file_nc_GEBCO)

    return gebco



def set_geographic_square_extent(ax, gdf_local, margin=1.02):
    """Set square, distance-proportionate axis limits for a lon/lat map.

    On a plain (non-projected) matplotlib axis the data are in geographic
    degrees, but one degree of longitude only spans cos(lat) of the ground
    distance of one degree of latitude. Plotting degrees 1:1 would stretch the
    map horizontally. We therefore convert the region's degree extent into
    kilometres (km_per_lon = 111*cos(lat), km_per_lat = 111), equalise the span
    on both axes (taking the larger so the whole region fits), and convert back
    to degrees per axis. The result is a square, undistorted map where 1 km
    looks the same on both axes -- the common convention for every map in this
    workflow, so ISA/CF/cutout/GEBCO are mutually coherent.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Plain (non-cartopy) axis to set limits on.
    gdf_local : geopandas.GeoDataFrame
        Region geometry in geographic CRS (EPSG:4326); its total_bounds set the
        extent.
    margin : float
        Multiplicative padding around the region (1.02 = 2 %).
    """
    xmin, ymin, xmax, ymax = gdf_local.total_bounds
    center_lat = (ymax + ymin) / 2
    km_per_lat = 111
    km_per_lon = 111 * np.cos(np.deg2rad(center_lat))
    center_x = (xmax + xmin) / 2
    center_y = (ymax + ymin) / 2
    delta_x = xmax - xmin
    delta_y = ymax - ymin
    delta_km = max(km_per_lon * delta_x, km_per_lat * delta_y) * margin
    ax.set_xlim(center_x - 0.5 * delta_km / km_per_lon, center_x + 0.5 * delta_km / km_per_lon)
    ax.set_ylim(center_y - 0.5 * delta_km / km_per_lat, center_y + 0.5 * delta_km / km_per_lat)



def geographic_aspect(center_lat):
    """matplotlib aspect for a cartopy PlateCarree GeoAxes so distances are proportionate.

    A PlateCarree GeoAxes defaults to aspect=1, i.e. one degree of longitude is
    drawn as long as one degree of latitude. Because a degree of longitude only
    spans cos(lat) of the ground distance of a degree of latitude, this stretches
    the map horizontally. Returning 1/cos(lat) (the matplotlib aspect = vertical
    units per horizontal unit) makes 1 km look the same on both axes, matching
    the cos(lat) correction used for the non-projected maps
    (see set_geographic_square_extent).
    """
    return 1.0 / np.cos(np.deg2rad(center_lat))



def plot_dataarray_on_map(
    data,
    gdf_context,
    gdf_NUTS_local,
    region,
    file_output,
    x_coord="lon",
    y_coord="lat",
    cmap="viridis",
    cbar_label="Value",
    vmin=None,
    vmax=None,
    vcenter=None,
    size=12,
    linewidth=1.5,
    fontsize=16,
    *,
    dpi,
    title=None,
    bounds_type="gdf"
):
    """
    Generic function to plot a 2D xarray DataArray on a map with NUTS boundaries.
    
    Parameters
    ----------
    data : xr.DataArray
        2D data array to plot
    gdf_context : geopandas.GeoDataFrame
        Reference boundaries drawn as thin grey context (e.g. same-level NUTS
        regions, or NUTS3 provinces for a CIMAS domain). See
        load_context_boundaries.
    gdf_NUTS_local : geopandas.GeoDataFrame
        Region/domain being plotted, drawn as a thick black outline.
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
    vcenter : float, optional
        Threshold value at which the colormap changes color (mapped to the
        midpoint of a diverging colormap via TwoSlopeNorm). The vmin/vmax
        limits are always respected. Only applied when vmin < vcenter < vmax;
        otherwise it is ignored and a regular linear scale is used.
        By default None (no color change, plain linear scale).
    size : int, optional
        Figure size (square), by default 12
    linewidth : float, optional
        NUTS boundary linewidth, by default 1.5
    fontsize : int, optional
        Font size for labels, by default 16
    dpi : int
        Output resolution in dots per inch. Must be provided explicitly.
    title : str, optional
        Plot title, by default None
    bounds_type : str, optional
        How to set axis limits: "gdf" (from NUTS boundaries) or "data" (from data coordinates),
        by default "gdf"
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    ##### Auto-detect spatial dimensions if defaults don't exist
    if x_coord not in data.dims or y_coord not in data.dims:
        detected_dims = _select_spatial_dims_xarray(data)
        if len(detected_dims) >= 2:
            x_coord, y_coord = detected_dims[0], detected_dims[1]
            _log_and_print(f"[plot_dataarray_on_map] Auto-detected spatial dims: x={x_coord}, y={y_coord}")
    
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
    
    # Use a TwoSlopeNorm when a threshold (vcenter) is given so the colormap
    # changes color at that value while still respecting vmin/vmax. This only
    # works if vmin < vcenter < vmax; otherwise fall back to a linear scale.
    use_twoslope = (
        vcenter is not None
        and vmin is not None
        and vmax is not None
        and vmin < vcenter < vmax
    )
    if use_twoslope:
        plot_kwargs["norm"] = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    else:
        if vcenter is not None:
            _log_and_print(
                f"[plot_dataarray_on_map] vcenter={vcenter} is outside "
                f"(vmin={vmin}, vmax={vmax}); ignoring threshold and using a linear scale."
            )
        # Add vmin/vmax if specified
        if vmin is not None:
            plot_kwargs["vmin"] = vmin
        if vmax is not None:
            plot_kwargs["vmax"] = vmax

    mappable = data.plot(**plot_kwargs)
    
    ##### Configure axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")

    ##### Add colorbar
    cbar = fig.colorbar(
        mappable,
        ax=ax,
        orientation='vertical',
        fraction=0.046,
        pad=0.04,
        shrink=0.9
    )
    cbar.set_label(cbar_label, fontsize=fontsize*1.5)
    cbar.ax.tick_params(labelsize=fontsize*1.5)

    ##### Distribute colorbar ticks evenly around the threshold and mark it
    if use_twoslope:
        # By default matplotlib's tick locator picks "nice" round values across
        # [vmin, vmax] with no knowledge of the threshold (vcenter), and the
        # TwoSlopeNorm stretches the two sides of the threshold by different
        # amounts. The result is ticks that look unevenly distributed above vs
        # below the colour change. Place an equal number of evenly spaced ticks
        # on each side of the threshold (each side maps linearly to half the
        # bar, so the ticks come out evenly spaced on the colorbar) plus the
        # threshold itself, so the labels are symmetric about the colour change.
        n_side = 3
        ticks_below = np.linspace(vmin, vcenter, n_side + 1)[:-1]
        ticks_above = np.linspace(vcenter, vmax, n_side + 1)[1:]
        ticks = np.concatenate([ticks_below, [vcenter], ticks_above])
        cbar.set_ticks(ticks)
        cbar.ax.set_yticklabels([f"{t:.2f}" for t in ticks])

        # Threshold line (where the colour changes).
        cbar.ax.axhline(y=vcenter, color="black", linewidth=linewidth, linestyle="--")

    ##### Add boundaries
    # Context regions (thin grey); see load_context_boundaries.
    gdf_context.plot(
        ax=ax, color="none", edgecolor='grey', linewidth=linewidth
    )

    # Region/domain being plotted (thick black).
    gdf_NUTS_local.plot(
        ax=ax, color="none", edgecolor='black', linewidth=linewidth*2
    )
    
    ##### Set axis limits
    if bounds_type == "gdf":
        # Use GDF bounds (geographic coordinates), with the cos(lat) correction
        # so the map is square and undistorted (1 km equal on both axes).
        set_geographic_square_extent(ax, gdf_NUTS_local)
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
        ax.set_title(title, fontsize=fontsize*1.2)
    else:
        ax.set_title("")  # Clear any auto-generated title from xarray
    
    ##### Save figure
    fig.savefig(file_output, bbox_inches="tight", pad_inches=0.2, dpi=dpi)
    plt.close(fig)
    
    _log_and_print(f"[plot_dataarray_on_map] Saved plot to {file_output}")



