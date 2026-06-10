import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from _schemas import Shapes
from matplotlib.colors import LinearSegmentedColormap

cmap_wind = LinearSegmentedColormap.from_list("cmap_wind", ["white", "blue"])
cmap_pv = LinearSegmentedColormap.from_list("cmap_pv", ["white", "orange"])


def create_plot_map(path_capacityfactors, path_shapes, path_map):
    """Plot a map of the mean capacityfactors."""
    cf = xr.open_dataarray(path_capacityfactors)
    shapes = gpd.read_parquet(path_shapes)
    shapes = Shapes.validate(shapes)
    shapes = shapes.set_index("shape_id")

    gdf_mean_cf = shapes.join(
        cf.mean(dim="time").to_dataframe(name="mean_capacityfactor")
    )

    map_capacity_factor(gdf_mean_cf=gdf_mean_cf, column="mean_capacityfactor")
    plt.savefig(path_map)


def create_plot_overview(path_capacityfactors, path_shapes, path_plot):
    """Plot an overview of the mean capacityfactors."""
    cf = xr.open_dataarray(path_capacityfactors)
    shapes = gpd.read_parquet(path_shapes)
    shapes = Shapes.validate(shapes)

    cf_mean = cf.mean(dim="time").to_dataframe(name="mean_capacityfactor").reset_index()

    gdf_mean_cf = pd.merge(shapes[["shape_id", "country_id"]], cf_mean, on="shape_id")

    plot_overview(cf=gdf_mean_cf, column="mean_capacityfactor")
    plt.savefig(path_plot)


def plot_overview(cf, column, ax=None, color="k"):
    """Plots the summary statistics of the capacity factors."""
    COLS = ["shape_id", "country_id", column]

    missing_cols = [col for col in COLS if col not in cf.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in cf: {missing_cols}")

    if ax is None:
        fig, ax = plt.subplots()

    _cf = cf.sort_values(["country_id", column], ascending=[True, False]).reset_index(
        drop=True
    )
    _cf.index.name = "index"

    ax.scatter(
        x=_cf.index,
        y=_cf[column],
        marker=".",
        linestyle="",
        linewidth=1,
        color=color,
        alpha=0.7,
        label="Mean CF",
    )
    ax.grid(which="both", axis="y", linestyle="--", alpha=0.5)

    # set a major xtick where a new country starts, and label it with the country code
    major = _cf.reset_index().groupby("country_id", as_index=False).first()
    ax.set_xticks(major["index"])

    # set minor xticks and labels at the midpoinst between the major ticks
    ticks = ax.get_xticks()
    midpoints = (ticks[:-1] + ticks[1:]) / 2
    midpoints = midpoints.tolist() + [ticks[-1]]  # add a final midpoint

    ax.set_xticks(midpoints, minor=True)
    ax.set_xticklabels(major["country_id"], rotation=90, fontsize=8, minor=True)

    ax.tick_params(axis="x", which="minor", length=0)  # hide minor tick marks
    ax.tick_params(axis="x", which="major", labelbottom=False)  # hide major tick labels
    plt.xticks(fontsize=8)

    # draw vertical lines for wind and pv at the major xticks
    ax.vlines(
        major["index"],
        ymin=0,
        ymax=major[column],
        color=color,
        linestyle="-",
        alpha=0.5,
    )

    ax.set_title("Capacity Factors")
    ax.set_ylabel("Mean CF")
    ax.set_xlabel("Shapes")

    return fig, ax


def map_capacity_factor(gdf_mean_cf, column, figsize=(4, 4)):
    fig, ax = plt.subplots(figsize=figsize, tight_layout=True)

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
