import netCDF4
import xarray as xr

# ds = xr.open_dataset(
#     "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/conc_g2_20210105.nc"
# )
# ds_meteo = xr.open_dataset(
#     "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/GAP_g2_20210105.nc"
# )
# ds["P"] = ds_meteo["P"]
# ds["SP"] = ds_meteo["SP"]
# ds.to_netcdf(
#     "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/conc_g2_20210105_3.nc"
# )

file_name_list = [
    "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/output_file.nc"
]
for file_name in file_name_list:
    nc_file = netCDF4.Dataset(file_name, "r+")

    # Convert time to days since 1900-01-01
    nc_file["time"][:] = nc_file["time"][:] / 24
    nc_file["time"].units = "days since 1900-01-01"
    nc_file["time"].calendar = "gregorian"

    # remove scale_factor and offset from the file
    del nc_file["P"].add_offset
    del nc_file["P"].scale_factor

    del nc_file["c_NO2"].add_offset
    del nc_file["c_NO2"].scale_factor

    nc_file.close()
