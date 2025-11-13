"""Prepare PV capacityfactors, given a cutout, a layout, spatial units to aggregate to and technology specifications."""

import _plots
import geopandas as gpd
import pandas as pd
import xarray as xr
import yaml
from _schemas import PointLayout, Shapes

import workflow.scripts._backend_atlite as _backend_atlite


def read_yaml(filepath):
    """Open a yaml file as python dictionary."""
    with open(filepath) as file:
        return yaml.safe_load(file)


def prepare_capacityfactors_point_layout(
    path_cutout, path_spatial_units, path_tech_specs, path_layout, path_output
):
    """Prepare capacityfactors aggregated to spatial units weighted by a point layout."""
    # load inputs
    spatial_units = gpd.read_parquet(path_spatial_units)
    spatial_units = Shapes.validate(spatial_units)
    layout = pd.read_csv(path_layout, index_col=0)
    layout = PointLayout.validate(layout)
    tech_specs = read_yaml(path_tech_specs)

    # prepare inputs
    spatial_units = spatial_units.set_index("shape_id")

    # compute capacityfactors
    capacityfactors = _backend_atlite.cf_agg_from_point_layout(
        path_cutout=path_cutout,
        layout=layout,
        spatial_units=spatial_units,
        tech_specs=tech_specs,
    )

    # save output
    capacityfactors.to_netcdf(path_output)


def plot(path_capacityfactors, path_spatial_units, path_map):
    """Plot capacityfactors."""
    # load inputs
    cf = xr.open_dataarray(path_capacityfactors)
    spatial_units = gpd.read_parquet(path_spatial_units)
    spatial_units = Shapes.validate(spatial_units)
    spatial_units = spatial_units.set_index("shape_id")
    gdf_mean_cf = spatial_units.join(
        cf.mean(dim="time").to_dataframe(name="wind_onshore")
    )

    # plot a map of annual capacityfactors
    fig, ax = _plots.map_capacity_factor(gdf_mean_cf=gdf_mean_cf, column="wind_onshore")
    fig.savefig(path_map)


if __name__ == "__main__":
    prepare_capacityfactors_point_layout(
        path_cutout=snakemake.input.cutout,
        path_spatial_units=snakemake.input.spatial_units,
        path_tech_specs=snakemake.input.tech_specs,
        path_layout=snakemake.input.layout,
        path_output=snakemake.output.data,
    )
    plot(
        path_capacityfactors=snakemake.output.data,
        path_spatial_units=snakemake.input.spatial_units,
        path_map=snakemake.output.plot_map,
    )
