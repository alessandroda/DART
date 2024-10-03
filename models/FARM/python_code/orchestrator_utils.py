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


# Define a function to check the status of the submitted job
def check_job_status_cresco(job_id):
    """Check if the job is still running."""
    command = f"bjobs {job_id}"
    try:
        output = subprocess.check_output(command, shell=True).decode("utf-8").strip()
        if "DONE" in output or "EXIT" in output:
            return True  # Job is finished
        return False  # Job is still running
    except subprocess.CalledProcessError as e:
        return False  # Job not found or already completed


def submit_slurm_job(case, option):
    job_name = f"DART{case}"
    error_file = f"err_{case}_{option}.log"
    output_file = f"out_{case}_{option}.log"
    match option:
        case "filter":
            command_execute = f"mpirun -np 17 ./filter"

    slurm_script = f"""#!/bin/sh


#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=17
#SBATCH --error=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/{error_file}
#SBATCH --output=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/{output_file}


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

    totseconds = seconds + 60 * (minutes + 60 * (hours))
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
        command = directory / command
        # Execute the command

        subprocess.call(str(command), shell=True)
        #subprocess.run(command, capture_output=True, text=True)
        # jobid = output.stdout.strip().split()[-1]
        # print(f"Job submitted for {case} with job ID : {jobid}")
    finally:
        # Change back to the original working directory
        os.chdir(original_directory)
    # return jobid

def run_command_in_directory_bsub(command, directory):
    original_directory = os.getcwd()
    try:
        breakpoint()
        os.chdir(directory)
        command = directory / command
        #subprocess.call(command,shell=True)
        output = subprocess.run(command, capture_output=True, text=True)
        jobid = output.stdout.strip().split()[-1]
        print(f"Job submitted for {command} with job ID : {jobid}")
    finally:
        os.chdir(original_directory)
    return jobid

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
        input_nml = input_nml.replace(key, str(value))

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


def prepare_farm_to_dart_nc(timestamp_farm, rounded_timestamp, seconds_model,days_model):
    breakpoint()
    meteo_file = f'/gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmwf_d0_g1_{timestamp_farm.strftime("%Y%m%d")}' 
    
    output_meteo_1hr = path.manager.base_path / 'RUN/data/temp/met_1hr.nc'
    output_meteo_1hr.parent.mkdir(parents = True, exist_ok = True)
    
    output_meteo_1hr_plus1 = path_manager.base_path /'RUN/data/temp/met_1hr+1.nc'
    
    arconv_input = path_manager.base_path / f'RUN/data/OUTPUT/OUT/ic_g1_{rounded_timestamp.strftime(%Y%m%d")}'
    
    arconv_output = path_manager.base_path / f'RUN/data/temp/ic_arconved.nc'

    conc_file = path_manager.base_path / f'RUN/data/to_DART/ic_g1_{seconds}_{days}_00.nc'
    conc_file.parent.mkdir(parents = True, exist_ok = True)

    # step 1: Select P, SP, T from the input file
    subprocess.run(["cdo", "selname,SP,P,T", meteo_file, output_meteo_1hr])
    # step 2: Select timestep from the rounded_timestamp
    subproces.run(["cdo", f"seltimestep,rounded_timesamp.hour",output_meteo_1hr, output_meteo_1hr])
    # step 3: Shift time by 1 hour
    subprocess.run(["cdo", "shifttime,1hour", output_meteo_1hr,output_meteo_1hr_plus1])
    # step 4: Convert FARM concentrations using arconv
    subprocess.run(["/gporq3/minni/FARM-DART/arconv-2.5.10", arconv_input,arconv_output, "1"])
    # step 5: Use ncks to append P, SP, and T variables to the FARM
    # concentration file
    subprocess.run(["ncks", "-A", "-v", "P,SP,T",output_meteo_1hr_plus1,arconv_ouput])

    # step 6
    subprocess.run(["cdo", "-setreftime,1900-01-01,00:00:00,days", arconv_output, arconv_output])

    # step 7 final file
    subprocess.run(["cdo", "-setcalendar,gregorian", arconv_ouput, conc_file])







