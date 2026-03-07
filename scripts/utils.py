import atlite
import geopandas as gpd
import logging


logger = logging.getLogger(__name__)


def _log_and_print(message):
    logger.info(message)
    print(message)



def load_gdf_nuts_local(file_gdf_NUTS, region):
    gdf_NUTS = (
        gpd.read_file(file_gdf_NUTS)
        .set_index("NUTS_ID")
    )

    _log_and_print(
        f"[load_gdf_nuts_local] Loaded gdf_NUTS for region: {region}, CRS: {gdf_NUTS.crs}"
    )

    return gdf_NUTS.loc[[region]]



def load_gdf_nuts_and_local(file_gdf_NUTS, region):
    gdf_NUTS = (
        gpd.read_file(file_gdf_NUTS)
        .set_index("NUTS_ID")
    )

    _log_and_print(
        f"[load_gdf_nuts_and_local] Loaded gdf_NUTS for region: {region}, CRS: {gdf_NUTS.crs}"
    )

    gdf_NUTS_local = gdf_NUTS.loc[[region]]

    return gdf_NUTS, gdf_NUTS_local



def load_and_limit_cutout(file_cutout, gdf_local):

    c = atlite.Cutout(file_cutout)

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



