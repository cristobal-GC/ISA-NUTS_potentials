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
    region = gdf_local.index[0] if len(gdf_local.index) > 0 else "unknown"
    cutout_obj = atlite.Cutout(file_cutout)
    gdf_local_crs_before = gdf_local.crs

    if gdf_local.crs is None:
        raise ValueError("gdf_local must have a defined CRS.")

    if gdf_local.crs != cutout_obj.crs:
        gdf_local = gdf_local.to_crs(cutout_obj.crs)

    _log_and_print(
        f"[load_and_limit_cutout] Loaded cutout from file_cutout: {file_cutout}, region: {region}, cutout CRS: {cutout_obj.crs}, gdf_local CRS(before): {gdf_local_crs_before}, gdf_local CRS(after): {gdf_local.crs}"
    )

    xmin, ymin, xmax, ymax = gdf_local.total_bounds

    margin_x = cutout_obj.dx
    margin_y = cutout_obj.dy

    xmin -= margin_x
    xmax += margin_x
    ymin -= margin_y
    ymax += margin_y

    return cutout_obj.sel(bounds=(xmin, ymin, xmax, ymax))