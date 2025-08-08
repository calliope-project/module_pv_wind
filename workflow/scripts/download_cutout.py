"""Download a cutout using atlite."""

import logging

import atlite

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    logger.info(f"Using atlite version: {atlite.__version__}")
    cutout_params = snakemake.config["cutout_params"]

    cutout_params["x"] = slice(*cutout_params["x"])
    cutout_params["y"] = slice(*cutout_params["y"])
    cutout_params["time"] = slice(*cutout_params["time"])

    features = cutout_params.pop("features", None)
    monthly_requests = cutout_params.pop("monthly_requests", False)
    logger.info(f"Preparing cutout with cutout_params: {cutout_params}")
    logger.info(f"Preparing cutout with features: {features}")
    cutout = atlite.Cutout(snakemake.output[0], **cutout_params)
    cutout.prepare(features=features, monthly_requests=monthly_requests)
