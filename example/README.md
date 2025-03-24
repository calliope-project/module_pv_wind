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
