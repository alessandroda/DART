from netCDF4 import Dataset
import numpy as np
import xarray as xr
import math

# Open the original NetCDF file
input_file = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/CAMEO/FARM/conc_g1_20201221_SO2.nc"

ds = xr.open_dataset(input_file, decode_times=False)
for i, value in enumerate(ds.time.values):
    days = math.floor(value)
    fractional_seconds = int((value - days) * 86400)
    days += 109207
    ds.isel(time=i).to_netcdf(
        f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/CAMEO/FARM/hourly_nc/model_{int(fractional_seconds)}_{days}.nc"
    )
