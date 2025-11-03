"""Prepare PV capacityfactors, given a cutout, a layout, spatial units to aggregate to and technology specifications."""

import backend_atlite
import geopandas as gpd
import rioxarray as rxr
import yaml


def read_yaml(filepath):
    """Open a yaml file as python dictionary."""
    with open(filepath) as file:
        return yaml.safe_load(file)


def prepare_capacityfactors_raster_layout(
    path_cutout, path_spatial_units, path_tech_specs, path_layout, path_output
):
    """Prepare capacityfactors aggregated to spatial units weighted by a raster layout."""
    # load inputs
    spatial_units = gpd.read_file(path_spatial_units)
    tech_specs = read_yaml(path_tech_specs)
    layout = rxr.open_rasterio(path_layout, masked=True)

    # prepare inputs
    spatial_units = spatial_units.set_index(spatial_units.columns[0])
    layout = layout.fillna(0)

    # compute capacityfactors
    capacityfactors = backend_atlite.cf_agg_from_raster_layout(
        path_cutout=path_cutout,
        layout=layout,
        spatial_units=spatial_units,
        tech_specs=tech_specs,
    )

    # save output
    capacityfactors.to_netcdf(path_output)


if __name__ == "__main__":
    prepare_capacityfactors_raster_layout(
        path_cutout=snakemake.input.cutout,
        path_spatial_units=snakemake.input.spatial_units,
        path_tech_specs=snakemake.input.tech_specs,
        path_layout=snakemake.input.layout,
        path_output=snakemake.output[0],
    )
