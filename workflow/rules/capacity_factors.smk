rule download_cutout:
    output:
        "results/cutout_era5_download.nc",
    conda:
        "../envs/geo.yaml"
    script:
        "../scripts/download_cutout.py"


if config["layout"] == "raster":

    rule prepare_capacityfactors_raster_layout:
        input:
            cutout=ancient("resources/user/cutout_{name_cutout}.nc"),  # TODO: Replace with results/cutout.nc as soon as the download works again
            tech_specs="resources/user/tech_specs_{name_tech}.yaml",
            layout="resources/user/layout_{name_layout}.tif",
            spatial_units="resources/user/spatial_units_{name_spatial_units}.geojson",
        output:
            "results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.nc",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_raster_layout.py"

elif config["layout"] == "point":

    rule prepare_capacityfactors_point_layout:
        input:
            cutout=ancient("resources/user/cutout_{name_cutout}.nc"),  # TODO: Replace with results/cutout.nc as soon as the download works again
            tech_specs="resources/user/tech_specs_{name_tech}.yaml",
            layout="resources/user/layout_{name_layout}.csv",
            spatial_units="resources/user/spatial_units_{name_spatial_units}.geojson",
        output:
            "results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.nc",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_point_layout.py"
