rule plot:
    input:
        capacity_factors="results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.nc",
        shapes="resources/user/spatial_units/{name_spatial_units}.geojson",
    output:
        output_plot="results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}.png",
        output_map="results/{name_cutout}/{name_spatial_units}/{name_layout}/capacityfactors_{name_tech}_map.png",
    conda:
        "../envs/atlite.yaml"
    script:
        "../scripts/plot.py"
