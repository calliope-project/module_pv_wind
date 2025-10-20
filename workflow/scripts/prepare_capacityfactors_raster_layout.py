"""Prepare PV capacityfactors, given a cutout, a layout, spatial units to aggregate to and technology specifications."""

import atlite
import backend_atlite
import geopandas as gpd
import rioxarray as rxr
import yaml


def read_yaml(filepath):
    """Open a yaml file as python dictionary."""
    with open(filepath) as file:
        return yaml.safe_load(file)


if __name__ == "__main__":
    cutout = atlite.Cutout(snakemake.input.cutout)
    spatial_units = gpd.read_file(snakemake.input.spatial_units)
    spatial_units = spatial_units.set_index(spatial_units.columns[0])
    tech_specs = read_yaml(snakemake.input.tech_specs)

    layout = rxr.open_rasterio(snakemake.input.layout, masked=True)
    layout = layout.fillna(0)

    capacityfactors = backend_atlite.cf_agg_from_raster_layout(
        cutout=cutout, layout=layout, spatial_units=spatial_units, tech_specs=tech_specs
    )

    capacityfactors.to_netcdf(snakemake.output[0])
