"""Prepare PV capacityfactors, given a cutout, a layout, spatial units to aggregate to and technology specifications."""

import backend_atlite
import geopandas as gpd
import pandas as pd
import yaml


def read_yaml(filepath):
    """Open a yaml file as python dictionary."""
    with open(filepath) as file:
        return yaml.safe_load(file)


if __name__ == "__main__":
    path_cutout = snakemake.input.cutout
    spatial_units = gpd.read_file(snakemake.input.spatial_units)
    spatial_units = spatial_units.set_index(spatial_units.columns[0])
    tech_specs = read_yaml(snakemake.input.tech_specs)

    layout = pd.read_csv(snakemake.input.layout)

    capacityfactors = backend_atlite.cf_agg_from_point_layout(
        path_cutout=path_cutout,
        layout=layout,
        spatial_units=spatial_units,
        tech_specs=tech_specs,
    )

    capacityfactors.to_netcdf(snakemake.output[0])
