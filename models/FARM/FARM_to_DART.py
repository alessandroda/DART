import netCDF4
import xarray as xr
import math
import numpy as np
import yaml
import subprocess
from datetime import datetime, timedelta
import os


def round_to_closest_hour(timestamp):
    if timestamp.minute >= 30:
        # Round up to the next hour
        rounded_timestamp = timestamp.replace(minute=0, second=0) + timedelta(hours=1)
    else:
        # Round down to the current hour
        rounded_timestamp = timestamp.replace(minute=0, second=0)
    return rounded_timestamp


# Define a function to check the status of the submitted job
def check_job_status(job_id):
    try:
        output = subprocess.check_output(["squeue", "-j", job_id])
        return True if job_id in output.decode() else False
    except subprocess.CalledProcessError:
        return False


def submit_slurm_job(case, option):
    job_name = f"DART_{case}"
    error_file = f"err_{case}_{option}.log"
    output_file = f"out_{case}_{option}.log"
    match option:
        case "filter":
            command_execute = f"mpirun -np 16 ./filter"

    slurm_script = f"""#!/bin/sh


#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=16
#SBATCH --error=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/errs/{error_file}
#SBATCH --output=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/outs/{output_file}


conda activate /home/dausilio/miniconda3/envs/dartenv
cd /mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work
{command_execute}
"""

    slurm_script_file = f"submit_job_{case}.sh"

    with open(slurm_script_file, "w") as file:
        file.write(slurm_script)

    output = subprocess.run(
        ["sbatch", slurm_script_file], capture_output=True, text=True
    )
    jobid = output.stdout.strip().split()[-1]
    print(f"Job submitted for {case} with job ID : {jobid}")
    return jobid


def set_date_gregorian(year, month, day, hours=0, minutes=0, seconds=0):
    """
    Computes time corresponding to date for Gregorian calendar.
    """

    base_year = 1601

    # Check for valid date and time
    if (
        seconds > 59
        or seconds < 0
        or minutes > 59
        or minutes < 0
        or hours > 23
        or hours < 0
        or day < 1
        or month > 12
        or month < 1
        or year < base_year
    ):

        errstring = f"year,mon,day,hour,min,sec {year} {month} {day} {hours} {minutes} {seconds} not a valid date."
        raise ValueError(errstring)

    if month != 2 and any([day > month_day for month_day in days_per_month]):
        raise ValueError(f"month ({month}) does not have {day} days.")

    # Check for leap year
    leap = is_leap_year(year)

    if month == 2 and (day > 29 or (not leap and day > 28)):
        raise ValueError(
            f"month ({month}) does not have {day} days in a non-leap year."
        )

    # Compute number of leap years fully past since base_year
    nleapyr = (
        (year - base_year) // 4 - (year - base_year) // 100 + (year - base_year) // 400
    )

    # Count up days in this year
    ndays = sum(
        days_per_month[m - 1] + (1 if leap and m == 2 else 0) for m in range(1, month)
    )

    totseconds = seconds + 60 * (minutes + 60 * (hours + 24 * ndays))
    totdays = day - 1 + ndays + 365 * (year - base_year - nleapyr) + 366 * nleapyr

    return totseconds, totdays


def is_leap_year(year):
    """
    Checks if the given year is a leap year.
    """
    if year % 4 != 0:
        return False
    elif year % 100 != 0:
        return True
    elif year % 400 != 0:
        return False
    else:
        return True


days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # Days in each month


def run_command_in_directory(command, directory):
    # Save the current working directory
    original_directory = os.getcwd()

    try:
        # Change to the specified directory
        os.chdir(directory)
        command = directory + "/" + command
        # Execute the command
        subprocess.call(command, shell=True)

    finally:
        # Change back to the original working directory
        os.chdir(original_directory)


def searchFile(t1, t2, listing):
    orbit_filename = listing[["filename", "start_time"]][
        (listing["start_time"] >= t1) & (listing["start_time"] <= t2)
    ]
    return orbit_filename


def replace_nml_template(
    input_nml_path: str, entries_tbr_dict: dict, output_nml_path: str
):
    # Validate input dictionary
    if not isinstance(entries_tbr_dict, dict):
        print("Error: 'entries_tbr_dict' must be a dictionary.")
        return

    # Read input file
    try:
        with open(input_nml_path, "r") as f1:
            input_nml = f1.read()
    except FileNotFoundError:
        print(f"Error: Input file '{input_nml_path}' not found.")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Replace entries
    for key, value in entries_tbr_dict.items():
        input_nml = input_nml.replace(key, value)

    # Write to output file
    try:
        with open(output_nml_path, "w") as f2:
            f2.write(input_nml)
    except Exception as e:
        print(f"Error writing to output file: {e}")
        return

    print("Replacement completed successfully.")


def open_dataset(path: str):
    return xr.open_dataset(path, mask_and_scale=False)


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
def modify_nc_file(ds, pol, time_list):
    is_present = False
    for timestamp in time_list:
        if timestamp in ds["time"].values:
            is_present = True
            break

    if not is_present:
        return

    ds = ds.sel(time=time_list)
    # Convert time to days since 1900-01-01
    ds["time"] = (
        (ds["time"] - np.datetime64("1900-01-01T00:00:00")) / np.timedelta64(1, "h")
    ).astype(int)
    ds["time"] = ds["time"] / 24
    ds["time"].attrs["units"] = "days since 1900-01-01"
    ds["time"].attrs["calendar"] = "gregorian"
    return ds


# Convert time to hours and days since a specific date and save as new files
def convert_time_to_hours(ds, output_folder, seconds, days, **kwargs):
    ens_member = kwargs.get("ens_member", None)
    for i, value in enumerate(ds.time.values):

        output_file_path = f"{output_folder}/model_{seconds}_{days}_{ens_member}.nc"
        ds.isel(time=i).to_netcdf(output_file_path)
        command_add_offset = f"ncatted -a add_offset,,d,, {output_file_path}"
        command_scale_factor = f"ncatted -a scale_factor,,d,, {output_file_path}"
        command_fill_value = f"ncatted -a _FillValue,,d,, {output_file_path}"
        # Execute the command
        # add print info
        print(f"Changing offset, scale_factor, fille_value attrs {output_file_path}")
        subprocess.run(command_fill_value, shell=True)
        # print(f"Deleting add_offset attrs in {output_file_path}")
        subprocess.run(command_add_offset, shell=True)
        # print(f"Deleting scale_factor attrs in {output_file_path}")
        subprocess.run(command_scale_factor, shell=True)
    ds.close()
