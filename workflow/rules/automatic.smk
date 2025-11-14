"""Rules to used to download automatic resource files."""

if config["download_cutout"]:

    path_cutout = ancient("resources/automatic/cutout_era5.nc")

    rule download_cutout:
        output:
            "resources/automatic/cutout_era5.nc",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/download_cutout.py"

else:
    path_cutout = ancient("resources/user/cutout_{name_cutout}.nc")
