if config["layout"] == "raster":

    rule prepare_capacityfactors_raster_layout:
        input:
            cutout=path_cutout,
            tech_specs="<tech_specs>",
            layout="<layout_raster>",
            shapes="<shapes>",
        output:
            data="<capacity_factors>",
            plot_map="<plot_map>",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_raster_layout.py"

elif config["layout"] == "point":

    rule prepare_capacityfactors_point_layout:
        input:
            cutout=path_cutout,
            tech_specs="<tech_specs>",
            layout="<layout_point>",
            shapes="<shapes>",
        output:
            data="<capacity_factors>",
            plot_map="<plot_map>",
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/prepare_capacityfactors_point_layout.py"
