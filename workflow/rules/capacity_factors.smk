if config["layout"] == "raster":

    rule prepare_capacityfactors_raster_layout:
        input:
            cutout=path_cutout,
            tech_specs="resources/user/tech_specs/{name_tech}.yaml",
            layout="resources/user/layout/{name_layout}.tif",
            spatial_units="resources/user/spatial_units/{name_spatial_units}.parquet",
        output:
            data="results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.nc",
            plot_map="results/{name_cutout}/{name_spatial_units}/{name_layout}/annual_capacity_factors_{name_tech}_map.png",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_raster_layout.py"

elif config["layout"] == "point":

    rule prepare_capacityfactors_point_layout:
        input:
            cutout=path_cutout,
            tech_specs="resources/user/tech_specs/{name_tech}.yaml",
            layout="resources/user/layout/{name_layout}.csv",
            spatial_units="resources/user/spatial_units/{name_spatial_units}.parquet",
        output:
            data="results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.nc",
            plot_map="results/{name_cutout}/{name_spatial_units}/{name_layout}/annual_capacity_factors_{name_tech}_map.png",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_point_layout.py"
