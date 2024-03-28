import netCDF4
import xarray as xr
import math
import numpy as np
import yaml
import subprocess


# Convert FARM file to DART readable format
def add_meteo_to_farm(input_file_path, meteo_file_path):
    # Open original FARM file
    ds = xr.open_dataset(input_file_path, mask_and_scale=False)

    # Open meteorological file
    ds_meteo = xr.open_dataset(meteo_file_path)

    # Add meteorological variables to FARM dataset
    ds["P"] = ds_meteo["P"]
    ds["SP"] = ds_meteo["SP"]
    ds["T"] = ds_meteo["T"]

    return ds


# Modify time dimension and remove scale_factor/offset
def modify_nc_file(ds, pol):
    # Convert time to days since 1900-01-01
    ds["time"] = (
        (ds["time"] - np.datetime64("1900-01-01T00:00:00")) / np.timedelta64(1, "h")
    ).astype(int)
    ds["time"] = ds["time"] / 24
    ds["time"].attrs["units"] = "days since 1900-01-01"
    ds["time"].attrs["calendar"] = "gregorian"

    # # Remove scale_factor and offset from specific variables
    # del ds["P"].attrs["add_offset"]
    # del ds["P"].attrs["scale_factor"]
    # del ds[pol].attrs["add_offset"]
    # del ds[pol].attrs["scale_factor"]

    return ds


# Convert time to hours and days since a specific date and save as new files
def convert_time_to_hours(ds, output_folder):
    for i, value in enumerate(ds.time.values):
        days = math.floor(value)
        fractional_seconds = int((value - days) * 86400)
        days += 109207
        output_file_path = f"{output_folder}/model_{int(fractional_seconds)}_{days}.nc"
        ds.isel(time=i).to_netcdf(output_file_path)
        command_add_offset = f"ncatted -a add_offset,,d,, {output_file_path}"
        command_scale_factor = f"ncatted -a scale_factor,,d,, {output_file_path}"
        command_fill_value = f"ncatted -a _FillValue,,d,, {output_file_path}"
        # Execute the command
        # add print info
        print(f"Deleting Fill_value attrs in {output_file_path}")
        subprocess.run(command_fill_value, shell=True)
        print(f"Deleting add_offset attrs in {output_file_path}")
        subprocess.run(command_add_offset, shell=True)
        print(f"Deleting scale_factor attrs in {output_file_path}")
        subprocess.run(command_scale_factor, shell=True)
    ds.close()


# Main function to perform all tasks
def main():
    # Read YAML file with settings
    with open(
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/farm_to_dart.yaml", "r"
    ) as yaml_file:
        settings = yaml.safe_load(yaml_file)

    input_farm_file = settings["input_farm_file"]
    input_meteo_file = settings["input_meteo_file"]
    output_folder = settings["output_folder"]
    pol = settings["pol"]

    # Convert FARM file to DART readable format
    ds = add_meteo_to_farm(input_farm_file, input_meteo_file)

    # Modify FARM file
    ds = modify_nc_file(ds, pol)

    # Convert time to hours and days since a specific date
    convert_time_to_hours(ds, output_folder)


if __name__ == "__main__":
    main()
