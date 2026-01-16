#%% 
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# load data
here = Path(__file__).parent

cf_wind_layout_point = xr.open_dataset(here / "results" / "cf_wind_layout_point.nc")

cf_wind_layout_raster = xr.open_dataset(here / "results" / "cf_wind_layout_raster.nc")

cf_wind_layout_raster_to_raster = xr.open_dataset(here / "results" / "cf_wind_layout_raster_to_raster.nc")
cf_wind_layout_raster_to_raster

#%%
df_cf = cf_wind_layout_point.to_dataframe()

idx=pd.IndexSlice

fig,ax=plt.subplots()
for id in df_cf.index.get_level_values("id").unique():
    timeseries = df_cf.loc[idx[:, id], "__xarray_dataarray_variable__"]
    timeseries.plot(ax=ax, label=id)

plt.legend(
    loc="upper right",
    bbox_to_anchor=(1.2, 1),
    bbox_transform=ax.transAxes,
)

#%%
df_cf = cf_wind_layout_raster.to_dataframe()

fig,ax=plt.subplots()
for id in df_cf.index.get_level_values("id").unique():
    timeseries = df_cf.loc[idx[:, id], "__xarray_dataarray_variable__"]
    timeseries.plot(ax=ax, label=id)

plt.legend(
    loc="upper right",
    bbox_to_anchor=(1.2, 1),
    bbox_transform=ax.transAxes,
)

#%%
cf_wind_layout_raster_to_raster
