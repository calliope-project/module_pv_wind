"""Plot capacity factors."""

import geopandas as gpd
import matplotlib.pyplot as plt
import xarray as xr


def plot_timeseries(da: xr.DataArray, name: str) -> plt.Figure:
    """Plot a time series of capacity factors for a given name."""
    fig, ax = plt.subplots(figsize=(6, 3), layout="constrained")
    da.sel(name=name).plot(ax=ax, label=name)
    ax.set_title(name)
    ax.set_xlabel("Time")
    ax.set_ylabel("Capacity factor")
    return fig


def plot_map(gdf: gpd.GeoDataFrame) -> plt.Figure:
    """Plot a map of a GeoDataFrame with capacity factors."""
    fig, ax = plt.subplots(figsize=(4, 3), layout="constrained")

    gdf.plot(ax=ax, column="Capacity factor")

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Annual average of capacity factor")

    return fig


def main(path_capacity_factors, path_shapes, output_plot, output_map):
    """Main function to plot capacity factors and save the plots."""
    da = xr.open_dataarray(path_capacity_factors)
    shapes = gpd.read_file(path_shapes)
    shapes = shapes.set_index(shapes.columns[0])

    plot_timeseries(da, name="AL")
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")

    annual_mean = da.mean("time")
    annual_mean = annual_mean.to_dataframe(name="Capacity factor")
    annual_mean = annual_mean.join(shapes)
    annual_mean = gpd.GeoDataFrame(annual_mean, geometry="geometry", crs=shapes.crs)

    plot_map(annual_mean)
    plt.savefig(output_map, dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main(
        snakemake.input.capacity_factors,
        snakemake.input.shapes,
        snakemake.output.output_plot,
        snakemake.output.output_map,
    )
