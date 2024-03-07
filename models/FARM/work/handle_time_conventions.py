import netCDF4

file_name_list = ["conc_g2_20210105.nc", "conc_g2_20210106.nc", "conc_g2_20210107.nc"]
for file_name in file_name_list:
    nc_file = netCDF4.Dataset(file_name, "r+")

    # Convert time to days since 1900-01-01
    nc_file["time"][:] = nc_file["time"][:] / 24
    nc_file["time"].units = "days since 1900-01-01"
    nc_file["time"].calendar = "gregorian"

    # remove scale_factor and offset from the file
    del nc_file["c_NO2"].add_offset
    del nc_file["c_NO2"].scale_factor

    nc_file.close()
