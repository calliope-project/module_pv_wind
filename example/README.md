# Minimal example for computing capacity factors using atlite

This is a minimal example that produces capacity factors for pv and wind
for the Netherlands. The input data has been cropped and downsampled to keep the files small.

To run the example, first create the conda environment

    conda env create -f environment.yaml

To download the ERA5 cutout, you need to provide your credentials in a file
called .cdsapirc. Given that, you can run

    python example/download_cutout_era5.py

, which will save the cutout in example/results/cutout_era5.nc. Next, you can run

    prepare_cf_from_raster_layout.py

and will find two files, example/results/cf_pv.nc and example/results/cf_wind.nc, containing the capacity factors.

## Use cases
  
Provide a raster layout, aggregate results to polygons. See `prepare_cf_from_raster_layout.py`

Provide a point layout, aggregate results to polygons. See `prepare_cf_from_point_layout.py`

Provide no layout, get capacity factors as raster at the resolution of the cutout. Same as one of the above, but do not pass `shapes`, `layout`, `matrix`, `per_unit`. Instead, pass `capacity_factor_timeseries=True`. See `prepare_cf_from_raster_no_layout_to_raster.py`.

Provide a raster or point layout, do not aggregate, get capacity factors as raster at cutout resolution. See `prepare_cf_from_raster_layout_to_raster.py`.
