import atlite
import logging
from pathlib import Path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    here = Path(__file__).parent
    path_cutout = here / "results" / "cutout_era5.nc"

    features = None
    cutout = atlite.Cutout(
        path_cutout, 
        module="era5",
        x=slice(3, 7.5),
        y=slice(50.5, 54),
        time=slice("2013-01-01", "2013-01-02"),
    )
    cutout.prepare(features=features)
