"""Rules to used to download automatic resource files."""

if config["download_cutout"]:

    path_cutout = ancient("<cutout>")

    rule download_cutout:
        output:
            path_cutout,
        conda:
            "../envs/atlite.yaml"
        script:
            "../scripts/download_cutout.py"

else:
    path_cutout = ancient("<cutout>")
