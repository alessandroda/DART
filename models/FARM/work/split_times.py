from netCDF4 import Dataset
import numpy as np
import xarray as xr

# Open the original NetCDF file
input_file = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/template_farm.nc"

ds = xr.open_dataset(input_file, use_cftime=True)
for time_ith in range(len(ds.time.values) - 1):
    print(time_ith)
    ds.isel(time=time_ith).to_netcdf(f"farm_time{time_ith}.nc")
    print(f"farm_time{time_ith}.nc")
