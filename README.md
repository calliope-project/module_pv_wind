# PV and wind capacity factors

This data module produces capacity factors for PV and wind at arbitrary spatial resolution.
It is a configurable, modular `snakemake` workflow as part of the [`clio`](https://clio.readthedocs.io/) data modules.
The module uses [atlite](https://atlite.readthedocs.io/) as a backend.
The module runs the following steps:

- Download meteorological reanalysis data (optional, can also be provided by the user).
- Given a technology specification, a capacity layout (raster or point data), and regions for aggregation, the module produces time series that are the weighted average of the capacity factor for each region.
- Outputs are saved as `.nc` files and visualised as maps.

## Using this module

This module can be imported directly into any `snakemake` workflow.
Please consult the integration example in `tests/integration/Snakefile` for more information.

To use this module, you need to specify the wind or PV technology, provide a capacity layout that is used as weights, and the regions for aggregation.
Meteorological reanalysis data can be either manually placed in the directory `resources/user`, or can be downloaded internally.
Please refer to `INTERFACE.yaml` for a full documentation of the module's interface.

## Development

We use [`pixi`](https://pixi.sh/) as our package manager for development.
Once installed, run the following to clone this repo and install all dependencies.

```shell
git clone git@github.com:calliope-project/module_pv_wind.git
cd module_pv_wind
pixi install --all
```

For testing, simply run:

```shell
pixi run test
```

To view the documentation locally, use:

```shell
pixi run serve-docs
```

To test a minimal example of a workflow using this module:

```shell
pixi shell    # activate this project's environment
cd tests/integration/  # navigate to the integration example
snakemake --use-conda  # run the workflow!
```
