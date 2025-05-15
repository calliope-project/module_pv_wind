import atlite
import rioxarray as rxr
import rasterio as rio
import pandas as pd
import geopandas as gpd
from pathlib import Path


def read_yaml(filepath):
    with open(filepath, "r") as file:
        return yaml.safe_load(file)


if __name__ == "__main__":
    here = Path(__file__).parent
    path_layout = here / "data" / "layout_points_NLD.csv"
    path_spatial_units = here / "data" / "regional_NLD.geojson"
    path_cutout = here / "results" / "cutout_era5.nc"
    path_cf_pv = here / "results" / "cf_pv_layout_point.nc"
    path_cf_wind = here / "results" / "cf_wind_layout_point.nc"

    # load data
    cutout = atlite.Cutout(path_cutout)
    layout_point = pd.read_csv(path_layout)
    spatial_units = gpd.read_file(path_spatial_units).set_index("id")


    # prepare layout from list of points
    layout_point = layout_point.rename(columns={"lon": "x", "lat": "y"})
    layout = cutout.layout_from_capacity_list(layout_point, col="capacity")

    # compute capacity factors
    capacityfactors_pv = cutout.pv(
        panel="CSi", 
        orientation={"slope": 30.0, "azimuth": 180.0}, 
        layout=layout,
        shapes=spatial_units,
        per_unit=True,
    )
    capacityfactors_pv.to_netcdf(path_cf_pv)

    capacityfactors_wind = cutout.wind(
        turbine="Vestas_V90_3MW", 
        layout=layout,
        shapes=spatial_units,
        per_unit=True,
    )
    capacityfactors_wind.to_netcdf(path_cf_wind)