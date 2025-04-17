"""Prepare PV capacityfactors, given a cutout, a layout, spatial units to aggregate to and technology specifications."""

import atlite
import geopandas as gpd
import rasterio as rio
import rioxarray as rxr
import yaml


def read_yaml(filepath):
    """Open a yaml file as python dictionary."""
    with open(filepath) as file:
        return yaml.safe_load(file)


if __name__ == "__main__":
    cutout = atlite.Cutout(snakemake.input.cutout)
    layout = rxr.open_rasterio(snakemake.input.layout)
    spatial_units = gpd.read_file(snakemake.input.spatial_units).set_index("id")
    tech_specs = read_yaml(snakemake.input.tech_specs)

    # resample layout to the resolution of the cutout
    match = (
        cutout.uniform_layout()
        .rio.write_crs(cutout.crs)
        .rio.write_transform(cutout.transform)
    )
    layout_matched = layout.squeeze(drop=True)
    layout_matched = layout_matched.rio.reproject_match(
        match, resampling=rio.enums.Resampling.sum, nodata=0
    )

    capacityfactors_pv = cutout.pv(
        shapes=spatial_units, layout=layout_matched, **tech_specs
    )
    capacityfactors_pv.to_netcdf(snakemake.output[0])
