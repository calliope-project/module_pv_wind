import atlite
import rioxarray as rxr
import rasterio as rio
import geopandas as gpd
from pathlib import Path


def read_yaml(filepath):
    with open(filepath, "r") as file:
        return yaml.safe_load(file)


if __name__ == "__main__":
    here = Path(__file__).parent
    path_availability = here / "data" / "availability_NLD.tif"
    path_spatial_units = here / "data" / "regional_NLD.geojson"
    path_cutout = here / "data" / "era5.nc"
    path_cf_pv = here / "results" / "cf_pv_no_layout_raster_to_raster.nc"
    path_cf_wind = here / "results" / "cf_wind_no_layout_raster_to_raster.nc"

    # load data
    cutout = atlite.Cutout(path_cutout)
    availability = rxr.open_rasterio(path_availability)
    spatial_units = gpd.read_file(path_spatial_units).set_index("id")

    # resample availability to the resolution of the cutout
    match = cutout.uniform_layout().rio.write_crs(cutout.crs).rio.write_transform(cutout.transform)
    layout = availability.squeeze(drop=True)
    layout = layout.rio.reproject_match(
        match,
        resampling=rio.enums.Resampling.sum,
        nodata=0
    )

    # experimental/hacky: create a matrix with an extra dimension with a coordinate 
    # for each pixel to "aggregate" to each pixel, which means not aggregating.
    matrix = layout.stack(spatial=["y", "x"])
    matrix = matrix.expand_dims(dim={"name":matrix.spatial.values}, axis=0)
    matrix = matrix.data

    # compute capacity factors
    capacityfactors_pv = cutout.pv(
        panel="CSi", 
        orientation={"slope": 30.0, "azimuth": 180.0}, 
        capacity_factor_timeseries=True
    )
    capacityfactors_pv.to_netcdf(path_cf_pv)

    capacityfactors_wind = cutout.wind(
        turbine="Vestas_V90_3MW", 
        capacity_factor_timeseries=True
    )
    capacityfactors_wind.to_netcdf(path_cf_wind)
