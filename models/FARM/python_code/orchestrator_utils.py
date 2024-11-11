import netCDF4
import xarray as xr
import math
import numpy as np
import yaml
import subprocess
from datetime import datetime, timedelta
import os
from pathlib import Path
import time
import logging
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)

def process_member(mem, path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model):
    try:
        meteo_file = f'/gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmwf_d0_g1_{timestamp_farm.strftime("%Y%m%d")}.nc'
        temp_output_meteo = path_manager.base_path / f'RUN/data/temp/output_meteo_{mem}.nc'
        temp_output_meteo_plus1 = path_manager.base_path / f'RUN/data/temp/output_meteo_plus1_{mem}.nc'
        temp_output_meteo_selected = path_manager.base_path / f'RUN/data/temp/output_meteo_selected_{mem}.nc'
        
        arconv_input_file = path_manager.base_path / f'RUN/data/OUTPUT_{mem}/OUT/ic_g1_{rounded_timestamp.strftime("%Y%m%d%H")}.nc'
        arconv_output_file = path_manager.base_path / f'RUN/data/temp/arconv_output_{mem}.nc'
        
        final_concentration_file = path_manager.base_path / f'RUN/data/to_DART/ic_g1_{seconds_model}_{days_model}_{mem}.nc'
        temp_concentration_file = path_manager.base_path / f'RUN/data/to_DART/temp_conc_{mem}.nc'
        temp1_concentration_file = path_manager.base_path / f'RUN/data/to_DART/temp1_conc_{mem}.nc'

        # Ensure the output directory exists
        final_concentration_file.parent.mkdir(parents=True, exist_ok=True)
        with open(f"logs_orchestrator/subprocess_out_{mem}.log", "a") as log_file:  # Append log file
            logging.info(f"Processing member {mem}")

            # Step 1: Select SP, P, and T from the input meteo file
            subprocess.run(
                ["cdo", "selname,SP,P,T", meteo_file, temp_output_meteo],
                stdout=log_file, stderr=log_file, check=True
            )

            # Step 2: Select the timestep from the rounded timestamp
            subprocess.run(
                ["cdo", f"seltimestep,{rounded_timestamp.hour}", temp_output_meteo, temp_output_meteo_selected],
                stdout=log_file, stderr=log_file, check=True
            )

            # Step 3: Shift time by 1 hour
            subprocess.run(
                ["cdo", "shifttime,1hour", temp_output_meteo_selected, temp_output_meteo_plus1],
                stdout=log_file, stderr=log_file, check=True
            )

            # Step 4: Convert FARM concentrations using arconv
            subprocess.run(
                ["/gporq3/minni/FARM-DART/arconv-2.5.10", arconv_input_file, arconv_output_file, "1"],
                stdout=log_file, stderr=log_file, check=True
            )

            # Step 5: Use ncks to append SP, P, and T variables to the FARM concentration file
            subprocess.run(
                ["ncks", "-A", "-v", "P,SP,T", temp_output_meteo_plus1, arconv_output_file],
                stdout=log_file, stderr=log_file, check=True
            )

            # Step 6: Copy the result to the final concentration file
            subprocess.run(["cp", arconv_output_file, final_concentration_file], check=True)

            # Step 7: Set reference time in the concentration file
            subprocess.run(
                ["cdo", "-setreftime,1900-01-01,00:00:00,days", final_concentration_file, temp_concentration_file],
                check=True
            )

            # Step 8: Set calendar to Gregorian
            subprocess.run(
                ["cdo", "-setcalendar,gregorian", temp_concentration_file, temp1_concentration_file],
                check=True
            )

            # Step 9: Remove unnecessary attributes from the concentration file
            subprocess.run(
                ["ncatted", "-a", "add_offset,,d,,", temp1_concentration_file], check=True
            )
            subprocess.run(
                ["ncatted", "-a", "scale_factor,,d,,", temp1_concentration_file], check=True
            )
            subprocess.run(
                ["ncatted", "-a", "_FillValue,,d,,", temp1_concentration_file], check=True
            )
            subprocess.run(
                ["ncatted", "-a", "missing_value,,d,,", temp1_concentration_file], check=True
            )

            # Step 10: Copy the cleaned file to the final concentration file location
            subprocess.run(["cp", temp1_concentration_file, final_concentration_file], check=True)

    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed for member {mem}: {e.cmd}")
        logging.error(f"Error output: {e.output}")
        raise
    finally:
        # Cleanup: Remove temporary files
        temp_files = [
            temp_output_meteo,
            temp_output_meteo_plus1,
            temp_output_meteo_selected,
            arconv_output_file,
            temp_concentration_file,
            temp1_concentration_file,
        ]
        for temp_file in temp_files:
            if temp_file.exists():
                temp_file.unlink(missing_ok=True)

def prepare_farm_to_dart_nc_par(path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model, no_mems):
    # Get number of workers from LSF or default to 1
    max_workers = 10
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_member, mem, path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model): mem for mem in range(no_mems)}
        for future in as_completed(futures):
            try:
                future.result()  # Check if there were exceptions
            except Exception as e:
                logging.error(f"An error occurred during processing of member {futures[future]}: {e}")



class CleanupContext:
    def __init__(self, temp_preproc_farm_folder: Path):
        self.temp_folder = temp_preproc_farm_folder

    def __enter__(self):
        os.makedirs(self.temp_folder, exist_ok=True)
        return self

    def __exit__(self):
        breakpoint()
        if os.path.exists(self.temp_folder):
            for filename in os.listdir(self.temp_folder):
                file_path = self.temp_folder / filename
                logger.info("I WILL DELETE IF NEEEEEDED")
                # file_path.unlink(missing_ok=True)


class TimeManager:
    def __init__(self, start_time: str, end_time: str, dt_seconds: int):
        """
        Initializes the TimeManager with the start, end times and the time delta.
        """
        self.start_time = pd.to_datetime(start_time)
        self.end_time = pd.to_datetime(end_time)
        self.current_time = self.start_time
        self.simulated_time = None
        self.timestamp_farm_run = None
        self.sat_obs = None
        self.dt = pd.Timedelta(dt_seconds, unit="s")

    def increment_time(self):
        """
        Increments the current time by the delta (dt).
        """
        self.current_time += self.dt

    @staticmethod
    def round_to_closest_hour(timestamp):
        if timestamp.minute >= 30:
            # Round up to the next hour
            rounded_timestamp = timestamp.replace(minute=0, second=0) + timedelta(
                hours=1
            )
        else:
            # Round down to the current hour
            rounded_timestamp = timestamp.replace(minute=0, second=0)
        return rounded_timestamp

    def is_within_bounds(self):
        """
        Checks if the current time is within the start and end bounds.
        """
        return self.current_time <= self.end_time

    def get_formatted_time(self, time_format: str = "%Y%m%d_%H%M%S"):
        """
        Returns the current time formatted as a string according to the specified format.
        Default format is YYYYMMDD_HHMMSS.
        """
        return self.current_time.strftime(time_format)

    def update_simulated_time(self, new_time):
        """
        Updates the simulated time to track the last simulated timestamp.
        """
        self.simulated_time = new_time


class PathManager:
    """
    Class to handle the paths and prevent errors along the execution.

    Attributes:
        base_path: The base path used for constructing all other paths.
        env_path: Path to the environment executable.
        listing_path: Path to the CSO listing file.
        run_submit_farm_template: Path to the submit farm template.
        path_submit_bsh: Path to the farm and dart submission directory.
        path_filter: Path to the filter executable.
    """

    def __init__(
        self,
        base_path,
        env_python,
        listing_file,
        run_submit_farm_template,
        path_submit_bsh,
        path_filter,
        path_data,
        log_paths=True,
    ):
        """
        Initializes the PathManager with a base path and relative paths, and checks if they exist.

        """
        self.base_path = Path(base_path).resolve()
        self.env_dir = self.base_path / env_python
        self.listing_file = self.base_path / listing_file
        self.run_submit_farm_template = self.base_path / run_submit_farm_template
        self.path_submit_bsh = self.base_path / path_submit_bsh
        self.path_filter = self.base_path / path_filter
        self.log_paths = log_paths
        self.path_data = self.base_path / path_data

        self.check_paths_exist()

    def check_paths_exist(self):
        """Checks if the base, environment, and listing paths exist."""
        logger.info("Checking paths")
        paths_to_check = {
            "Base path": self.base_path,
            "Environment directory": self.env_dir,
            "Listing file": self.listing_file,
            "Submit farm template": self.run_submit_farm_template,
            "Farm submission path": self.path_submit_bsh,
            "Submit filter path": self.path_filter,
            "Path data": self.path_data,
        }

        for path_name, path_value in paths_to_check.items():
            if not path_value.exists():
                raise FileNotFoundError(f"{path_name} does not exist: {path_value}")
            if self.log_paths:
                logger.info(f"{path_name} exists: {path_value}")

    def get_paths(self):
        """Returns all paths as a tuple."""
        return (
            self.base_path,
            self.env_dir,
            self.listing_file,
            self.run_submit_farm_template,
            self.path_submit_farm,
            self.path_submit_filter,
        )


# def round_to_closest_hour(timestamp):
#     if timestamp.minute >= 30:
#         # Round up to the next hour
#         rounded_timestamp = timestamp.replace(minute=0, second=0) + timedelta(hours=1)
#     else:
#         # Round down to the current hour
#         rounded_timestamp = timestamp.replace(minute=0, second=0)
#     return rounded_timestamp


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

    output_file = directory / f"{command}_output.log"
    original_directory = os.getcwd()

    try:
        os.chdir(directory)
        command = directory / command
        logging.info(f"running command: {command}")
        with open(output_file, "a") as log_file:
            subprocess.call(str(command), shell=True, stdout=log_file, stderr=log_file)
    finally:
        os.chdir(original_directory)


def run_command_in_directory_bsub(command, directory, farm=True):
    #breakpoint()
    original_directory = os.getcwd()
    try:
        os.chdir(directory)
        command = directory / command
        # subprocess.call(command,shell=True)
        subprocess.run(["chmod", "+x", command])
        # 08.10.2024
        # subprocess.run should return a list of ids when running on the
        # members
        output = subprocess.run(command, capture_output=True, text=True)
        # breakpoint()
        if farm:
            lines = output.stdout.strip().splitlines()
            jobid = [line.split('id:')[1] for line in lines if 'id:' in line]
            print(f"Job submitted for {command} with job IDs : {jobid}")

        else:
            jobid = output.stdout.strip().split()[1]
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

    logger.info(
        f"Replacement {input_nml_path} to {output_nml_path} completed successfully."
    )


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


def prepare_dart_to_farm_nc(path_manager, output_sim_folder, time_model, ass_var, no_mems):
    """
    Prepares posterior NetCDF files to a FARM standard by adjusting reference time,
    calendar format, and appending the 'c_SO2' variable to the FARM prior file.
    """
    for mem in range(no_mems):
        try:
            posterior_file = f"{output_sim_folder}/ic_g1_posterior_{time_model}_{mem}.nc"

            tmp_dir = Path(output_sim_folder) / "tmp"

            tmp0_posterior = tmp_dir / "tmp0.nc"
            tmp1_posterior = tmp_dir / "tmp1.nc"
            Path(tmp0_posterior).parent.mkdir(parents=True, exist_ok=True)

            # prior in FARM format
            prior_farm_folder = Path(f"{path_manager.path_data}/OUTPUT_{mem}/OUT/prior/")
            prior_farm_folder.mkdir(parents=True, exist_ok=True)
            #breakpoint()
            prior_from_farm_file = (
                prior_farm_folder / f'ic_g1_{time_model}.nc'
            )

            # first to delete variable c_SO2
            result_tmp = Path(
                f'{path_manager.path_data}/OUTPUT_{mem}/OUT/ic_g1_{time_model}_tmp.nc'
            )

            result = Path(
                f'{path_manager.path_data}/OUTPUT_{mem}/OUT/ic_g1_{time_model}.nc'
            )



            with open("subprocess_output.log", "w") as log_file:
                logger.info("0: Storing prior elsewhere")
                #breakpoint()

                subprocess.run(
                    ["mv", result, prior_from_farm_file],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )

                logger.info("1: Setreftime of posterior to FARM standard")
                subprocess.run(
                    [
                        "cdo",
                        "setreftime,1900-01-01,00:00:00,hours",
                        str(posterior_file),
                        str(tmp0_posterior),
                    ],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )
                #breakpoint() 
                logger.info("2: setcalendar of posterior to FARM standard")
                subprocess.run(
                    [
                        "cdo",
                        "setcalendar,proleptic_gregorian",
                        str(tmp0_posterior),
                        str(tmp1_posterior),
                    ],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )

                logger.info("3: rename x,y to lat/lon")
                subprocess.run(
                    ["ncrename", "-d", "x,lon", "-d", "y,lat", str(tmp1_posterior)],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )
                logger.info(f"4-5: delete {ass_var} and append the assimilated one")
                ds = xr.open_dataset(str(prior_from_farm_file))
                #ds_result = ds.drop_vars(ass_var)
                ds_tmp1_posterior = xr.open_dataset(str(tmp1_posterior))
                ds['c_SO2'].values = ds_tmp1_posterior[ass_var].values 
                ds.to_netcdf(str(result_tmp))
                # Step 4: Convert FARM concentrations using arconv
                logger.info("6: convert FARM with arconv")
                subprocess.run(
                    [
                        "/gporq3/minni/FARM-DART/arconv-2.5.10",
                        result_tmp,
                        result,
                        "1",
                    ],
                    stdout=log_file,
                    stderr=log_file,
                )

                #logger.info("4: delete var from FARM prior")
                #subprocess.run(
                #    ["cdo", f"delname,{ass_var}",str(prior_from_farm_file), str(result)],
                #    stdout=log_file,
                #    stderr=log_file,
                #    check=True,
                #)

                #logger.info("5: append assimilated c_SO2 from posterior DART")
                #subprocess.run(
                #    ["ncks", "-A","-v", str(ass_var), str(tmp1_posterior), str(result)],
                #    stdout=log_file,
                #    stderr=log_file,
                #    check=True,
                #)
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {e.cmd}")
            logging.error(f"Error output: {e.output}")
            logging.info(
                f"something went wrong restore at least the original core in {result} to continue farm run"
            )
            with open("subprocess_output.log","a") as log_file:
                subprocess.run(
                    ["mv", prior_from_farm_file, result],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )
            raise
        finally:
            # Cleanup: Remove temporary files
            temp_files = [tmp0_posterior, tmp1_posterior]
            for temp_file in temp_files:
                temp_file.unlink(missing_ok=True)
            

def prepare_farm_to_dart_nc(
    path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model,
    no_mems
):
    for mem in range(no_mems):
        try:
            logging.info
            meteo_file = f'/gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmwf_d0_g1_{timestamp_farm.strftime("%Y%m%d")}.nc'

            temp_output_meteo = path_manager.base_path / "RUN/data/temp/output_meteo.nc"
            temp_output_meteo_plus1 = (
                path_manager.base_path / "RUN/data/temp/output_meteo_plus1.nc"
            )
            temp_output_meteo_selected = (
                path_manager.base_path / "RUN/data/temp/output_meteo_selected.nc"
            )
            logging.info('Running ARCONV for the members')
            arconv_input_file = (
                path_manager.base_path
                / f'RUN/data/OUTPUT_{mem}/OUT/ic_g1_{rounded_timestamp.strftime("%Y%m%d%H")}.nc'
            )
            arconv_output_file = path_manager.base_path / "RUN/data/temp/arconv_output.nc"

            final_concentration_file = (
                path_manager.base_path
                / f"RUN/data/to_DART/ic_g1_{seconds_model}_{days_model}_{mem}.nc"
            )
            temp_concentration_file = (
                path_manager.base_path / "RUN/data/to_DART/temp_conc.nc"
            )
            temp1_concentration_file = (
                path_manager.base_path / "RUN/data/to_DART/temp1_conc.nc"
            )

            final_concentration_file.parent.mkdir(parents=True, exist_ok=True)
            with open("logs_orchestrator/subprocess_out.log", "w") as log_file:
                logger.info("1: Select SP, P, and T from the input meteo file")
                subprocess.run(
                    ["cdo", "selname,SP,P,T", meteo_file, temp_output_meteo],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )

                logger.info("2: Select the timestep from the rounded timestamp")
                subprocess.run(
                    [
                        "cdo",
                        f"seltimestep,{rounded_timestamp.hour}",
                        temp_output_meteo,
                        temp_output_meteo_selected,
                    ],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )

                # Step 3: Shift time by 1 hour
                logger.info("3: shift time by 1 hour")
                subprocess.run(
                    [
                        "cdo",
                        "shifttime,1hour",
                        temp_output_meteo_selected,
                        temp_output_meteo_plus1,
                    ],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )

                # Step 4: Convert FARM concentrations using arconv
                logger.info("4: convert FARM with arconv")
                subprocess.run(
                    [
                        "/gporq3/minni/FARM-DART/arconv-2.5.10",
                        arconv_input_file,
                        arconv_output_file,
                        "1",
                    ],
                    stdout=log_file,
                    stderr=log_file,
                )

                # Step 5: Use ncks to append SP, P, and T variables to the FARM concentration file
                logger.info("5: use ncks to append SP, P and T variables to the FARM concentration file")
                subprocess.run(
                    [
                        "ncks",
                        "-A",
                        "-v",
                        "P,SP,T",
                        temp_output_meteo_plus1,
                        arconv_output_file,
                    ]
                )

                # Step 6: Copy the result to the final concentration file
                subprocess.run(["cp", arconv_output_file, final_concentration_file])

                # Step 7: Set reference time in the concentration file
                subprocess.run(
                    [
                        "cdo",
                        "-setreftime,1900-01-01,00:00:00,days",
                        final_concentration_file,
                        temp_concentration_file,
                    ]
                )

                # Step 8: Set calendar to Gregorian
                subprocess.run(
                    [
                        "cdo",
                        "-setcalendar,gregorian",
                        temp_concentration_file,
                        temp1_concentration_file,
                    ]
                )

                # Step 9: Remove unnecessary attributes from the concentration file
                subprocess.run(
                    ["ncatted", "-a", "add_offset,,d,,", temp1_concentration_file]
                )
                subprocess.run(
                    ["ncatted", "-a", "scale_factor,,d,,", temp1_concentration_file]
                )
                subprocess.run(
                    ["ncatted", "-a", "_FillValue,,d,,", temp1_concentration_file]
                )
                subprocess.run(
                    ["ncatted", "-a", "missing_value,,d,,", temp1_concentration_file]
                )

                # Step 10: Copy the cleaned file to the final concentration file location
                logger.info("10: copy the cleaned file to the final concentration")
                subprocess.run(["cp", temp1_concentration_file, final_concentration_file])
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {e.cmd}")
            logging.error(f"Error output: {e.output}")
            raise
        finally:
            # Cleanup: Remove temporary files
            temp_files = [
                temp_output_meteo,
                temp_output_meteo_plus1,
                temp_output_meteo_selected,
                arconv_output_file,
                temp_concentration_file,
                temp1_concentration_file,
            ]
            for temp_file in temp_files:
                temp_file.unlink(missing_ok=True)


def are_terminated_correctly(job_ids : list, path_manager : PathManager, timestamp_farm : str, no_mems : int) -> bool:
    # breakpoint()
    while True:
        running_jobs = []
        for jobid in job_ids:
            if not check_job_status_cresco(jobid):
                running_jobs.append(jobid)

        if not running_jobs:
            logger.info(f"{job_ids} have finished")
            return ic_g1_exist(path_manager, timestamp_farm, no_mems)
        else:
            logger.info(
                f"Jobs still running: {running_jobs}. Waiting for them to finish..."
            )
            time.sleep(30)

def submit_and_wait(path_manager : PathManager, commands_with_directories: list, timestamp_farm : str, no_mems : int) -> bool:
    # breakpoint()
    for command, directory in commands_with_directories:
        job_ids = run_command_in_directory_bsub(command, directory)
    if not are_terminated_correctly(job_ids, path_manager, timestamp_farm, no_mems):
        job_ids = run_command_in_directory_bsub(command,directory)

def ic_g1_exist(path_manager : PathManager, timestamp_farm : str, no_mems : int):
    file_name = f'ic_g1_{timestamp_farm}.nc'
    all_exist = True
    for mem in range(no_mems):
        file_path  = Path(path_manager.base_path / f'RUN/data/OUTPUT_{mem}/OUT/{file_name}')
        if os.path.exists(file_path):
            print(f"The core {file_name} for mem {mem} exists in the directory")
        else:
            print(f"The core {file_name} for mem {mem} does not exist in the directory")
            all_exist = False
            break
    return all_exist