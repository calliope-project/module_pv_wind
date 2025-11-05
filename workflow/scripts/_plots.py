import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap

cmap_wind = LinearSegmentedColormap.from_list("cmap_wind", ["white", "blue"])
cmap_pv = LinearSegmentedColormap.from_list("cmap_pv", ["white", "orange"])


def average_capacity_factors(
    cf: xr.DataArray, shapes: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Calculate average capacity factors per spatial unit."""
    df_mean_cf = cf.mean(dim="time").to_dataframe(name="average_cf").reset_index()
    df_mean_cf.index.name = "id"

    gdf_mean_cf = gpd.GeoDataFrame(
        df_mean_cf.join(shapes[["shape_id", "geometry"]]),
        geometry="geometry",
        crs=shapes.crs,
    )

    return gdf_mean_cf


def map_capacity_factor(gdf_mean_cf, column, figsize=(4, 4)):
    fig, ax = plt.subplots(figsize=figsize, tight_layout=True)
    print(gdf_mean_cf)
    gdf_mean_cf.plot(ax=ax, column=column, cmap=cmap_wind, legend=True, aspect=None)
    gdf_mean_cf.geometry.boundary.plot(ax=ax, color="black", linewidth=0.5)
    ax.set_title("Average Capacity Factor\nOnshore Wind")
    _blank_axis(ax)

    return fig, ax


def _blank_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
