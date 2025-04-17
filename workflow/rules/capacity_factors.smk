rule download_cutout:
    output:
        "results/cutout.nc",
    conda:
        "../envs/geo.yaml"
    script:
        "../scripts/download_cutout.py"


rule prepare_capacityfactors_wind:
    input:
        cutout=ancient("resources/user/era5.nc"),
        layout="resources/user/layout_{tech}.tif",
        spatial_units="resources/user/spatial_units/{resolution}.geojson",
        tech_specs="resources/user/tech_specs/{tech}_{subtech}.yaml",
    output:
        "results/{resolution}/capacityfactors_{tech}_{subtech}.nc",
    wildcard_constraints:
        tech="wind_offshore|wind_onshore",
    conda:
        "../envs/atlite.yaml"
    script:
        "../scripts/prepare_capacityfactors_wind.py"


rule prepare_capacityfactors_pv:
    input:
        cutout=ancient("resources/user/era5.nc"),
        layout="resources/user/layout_{tech}.tif",
        spatial_units="resources/user/spatial_units/{resolution}.geojson",
        tech_specs="resources/user/tech_specs/{tech}_{subtech}.yaml",
    output:
        "results/{resolution}/capacityfactors_{tech}_{subtech}.nc",
    wildcard_constraints:
        tech="pv_rooftop|pv_open_field",
    conda:
        "../envs/atlite.yaml"
    script:
        "../scripts/prepare_capacityfactors_pv.py"
