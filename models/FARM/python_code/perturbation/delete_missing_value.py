from glob import glob
import os
import subprocess
import xarray as xr

path_emissions = "/gporq3/minni/FARM-DART/perturbation_fields/emission_base"
path_emission_changed = "/gporq3/minni/FARM-DART/perturbation_fields/emission_base_mod"
name_netcdfs = "HERMESv3_*.nc"

if not os.path.exists(path_emission_changed):
    os.makedirs(path_emission_changed)

netcdfs = glob(str(path_emissions / name_netcdfs))
for netcdf_emi in netcdfs:
    try:
        ds = xr.open_dataset(netcdf_emi)
    except:
        continue
    file_name = f"{os.path.basename(netcdf_emi)}"
    file_path = path_emission_changed / file_name
    ds.to_netcdf(file_path)
    command_fill_value = f"ncatted -a _FillValue,,d,, {file_name}"
    subprocess.run(command_fill_value, shell=True)
