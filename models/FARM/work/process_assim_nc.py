import xarray as xr
import glob

files = glob.glob("/mnt/mumbai_n4r5/dausilio/SATDA_OFFLINE/satassim_NO2_202101*")
result = None
for file in files:
    # check result is empty
    if not result:
        result = xr.open_dataset(file)
        result = result[
            ["mdlspace_ys", "mdlspace_yr", "mdlspace_ya", "background", "analysis"]
        ]
        continue
    ds = xr.open_dataset(file)
    ds = ds[["mdlspace_ys", "mdlspace_yr", "mdlspace_ya", "background", "analysis"]]
    result = xr.concat([result, ds], dim="time")

result.sortby("time").resample(time="M").mean().to_netcdf(
    "/mnt/mumbai_n4r5/dausilio/SATDA_OFFLINE/01_assim.nc"
)
