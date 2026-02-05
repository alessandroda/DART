from dataclasses import dataclass
from typing import Optional, List, Tuple
import netCDF4
import xarray as xr
import math
import numpy as np
import yaml
import subprocess
from datetime import timedelta
import os
from pathlib import Path
import time
import logging
import pandas as pd
import re
from pipeline_errors import SchedulerError
from mimesi_types import ModelType, Scheduler
from paths import PathManager
import shlex
from scheduler import submit_job, wait_for_slurm_jobs

logger = logging.getLogger(__name__)

@dataclass
class CommandSpec:
    command: str  # executable/script name
    directory: Path  # where to run it
    args: Optional[List[str]] = None

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
        self.last_perturbed_day = None
        self.end_file_date_control_run = self.start_time - timedelta(days=1)
        self.end_file_date = None

        self.check_start_ahead_end()

    def check_start_ahead_end(self):
        """check if end time follows start time"""
        if self.start_time > self.end_time:
            logger.warning(
                f"End time ({self.end_time}) is ahead start time {self.start_time}"
            )

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


# def round_to_closest_hour(timestamp):
#     if timestamp.minute >= 30:
#         # Round up to the next hour
#         rounded_timestamp = timestamp.replace(minute=0, second=0) + timedelta(hours=1)
#     else:
#         # Round down to the current hour
#         rounded_timestamp = timestamp.replace(minute=0, second=0)
#     return rounded_timestamp


# Define a function to check the status of the submitted job
def check_job_status_cresco(job_id, **kwargs):
    """Check if the job is still running."""
    command = f"bjobs {job_id}"
    which_run = kwargs.get("which_run", None)
    try:
        output = subprocess.check_output(command, shell=True).decode("utf-8").strip()
        if "DONE" in output:
            logger.info(f"{which_run} {job_id}: status DONE")
            return True
        elif "EXIT" in output:
            logger.info(f"{which_run} {job_id}: status EXIT")
            return True
        elif "RUN" in output:
            logger.info(f"{which_run} {job_id}: status RUN")
            return False
        elif "PEND" in output:
            logger.info(f"{which_run} {job_id}: status PEND")
            return False
        else:
            logger.warning(f"{job_id}: Unknown status")
            return True  # Job is still running
    except subprocess.CalledProcessError as e:
        return True  # Job not found or already completed


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

    # if month != 2 and any([day > month_day for month_day in days_per_month]):
    #    raise ValueError(f"month ({month}) does not have {day} days.")
    if day > days_per_month[month - 1]:
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

def run_command_in_directory(spec: CommandSpec) -> Tuple[int, Optional[str]]:
    logger = logging.getLogger(__name__)
    original_directory = os.getcwd()

    try:
        logger.info(f"[CMD] Entering directory: {spec.directory}")
        os.chdir(spec.directory)
        
        """
        command_path = Path(spec.command)
        if not command_path.is_absolute():
            command_path = spec.directory / command_path

        if not command_path.exists():
            raise FileNotFoundError(command_path)

        cmd = [str(command_path)]
        subprocess.run(["chmod", "+x", str(command_path)], check=True)
        """     
        cmd = shlex.split(spec.command)
        if not cmd:
            raise ValueError("Empty command")

        # If first token is a path, ensure it exists + executable
        first = Path(cmd[0])
        if first.exists(): 
            first = first.resolve()
            # Make executable if needed
            subprocess.run(
                ["chmod", "+x", str(first)],
                check=True,
            )
            cmd[0] = str(first)
        
        ### end of changes
        if spec.args:
            cmd.extend(map(str, spec.args))

        logger.info(
            "[CMD] Running: %s",
            " ".join(shlex.quote(c) for c in cmd),
        )

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        ) #stdout='Submitted Batch Session 3561287\n'
        stdout = result.stdout.strip()
        stderr = result.stderr.strip()

        stdout = result.stdout
        
        # Try strict format: ID:12345
        job_ids = re.findall(r"^ID:(\d+)$", stdout, re.MULTILINE)

        # Fallback: any number
        if not job_ids:
            job_ids = re.findall(r"\d+", stdout)
        
        if not job_ids:
            raise RuntimeError(
                "No Slurm job IDs found in output.\n" "Expected lines like: ID:<jobid>"
            )
        return result.returncode, job_ids
    finally:
        os.chdir(original_directory)


def run_command_in_directory_bsub(
    command, directory, farm=True, replace_emissions=False
):

    if farm and replace_emissions:
        raise ValueError("Only one of 'farm' or 'replace_emissions' can be True.")

    # breakpoint()
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
            jobid = [line.split("id:")[1] for line in lines if "id:" in line]
            print(f"Job submitted for {command} with job IDs : {jobid}")
        elif replace_emissions:
            lines = output.stdout.strip().splitlines()
            jobid = [
                re.search(r"Job <(\d+)>", line).group(1)
                for line in lines
                if re.search(r"Job <(\d+)>", line)
            ]
            print(f"Job submitted for {command} with job IDs: {jobid}")
        else:
            jobid = output.stdout.strip().split()[1]
            print(f"Job submitted for {command} with job ID : {jobid}")
    finally:
        os.chdir(original_directory)
    return jobid

def submit_irene(spec: CommandSpec) -> str:
    rc, job_id = run_command_in_directory(spec)

    if rc != 0:
        raise SchedulerError(
            f"Submission command failed: {spec.command} " f"(return code {rc})"
        )
    if not job_id:
        raise SchedulerError(f"No job id returned by command {spec.command}")
    
    time.sleep(5)
    logger.info(f"[TGCC-IRENE] Submitted job with ID:{job_id}")

    return job_id

def submit_and_wait_cineca(
    path_manager: PathManager,
    spec: CommandSpec,
    timestamp_chimere: str,
    no_mems: int,
    scheduler: Scheduler,
    model_type: ModelType,
):

    rc, job_ids = run_command_in_directory(spec)

    if rc != 0:
        raise SchedulerError(
            f"Submission command failed: {spec.command} " f"(return code {rc})"
        )

    if not job_ids:
        raise SchedulerError(f"No job id returned by command {spec.command}")

    logger.info(f"[SLURM] Submitted job {job_ids}")
    time.sleep(10)

    return job_ids

def submit_and_wait_slurm(
    model_type: ModelType,
    path_manager: PathManager,
    scheduler: Scheduler,
    commands_with_directories: list[tuple[Path, Path]],
    timestamp_model: str,
    no_mems: int,
    case_dir: str,
    queue: str,
    max_retries: int = 2,
) -> bool:

    attempt = 0

    while attempt <= max_retries:
        attempt += 1
        logger.info(f"SLURM submission attempt {attempt}")

        # --- submit ---
        all_job_ids = []
        for command, directory in commands_with_directories:
            job_ids = submit_job(scheduler, command, directory)
            all_job_ids.extend(job_ids)
            time.sleep(2)

        # --- wait until finished ---
        wait_for_slurm_jobs(all_job_ids)

        # --- inspect outputs ---
        mems_to_rerun = get_list_mems_to_rerun(
            all_job_ids,
            path_manager,
            timestamp_model,
            no_mems,
        )

        if not mems_to_rerun:
            logger.info("All ensemble members completed successfully")
            return True

        logger.warning(f"Members to rerun: {mems_to_rerun}")

        if attempt >= max_retries:
            raise RuntimeError(
                f"SLURM retries exceeded. Failed members: {mems_to_rerun}"
            )

        # --- prepare rerun script ---
        list_mems = [str(mem) for mem in mems_to_rerun]

        replace_nml_template(
            input_nml_path=path_manager.base_path.run_submit_model_template,
            entries_tbr_dict={
                "da_date_start": timestamp_model,
                "da_date_end": timestamp_model,
                "@no_mems_list": str(tuple(list_mems)).replace(",", ""),
                "@case_dir": case_dir,
                "@cresco_queue": queue,
            },
            output_nml_path=commands_with_directories[0][1]
            / commands_with_directories[0][0],
        )

        # Only rerun failed members
        commands_with_directories = [commands_with_directories[0]]

    return False


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


def prepare_dart_to_farm_nc(
    path_manager, output_sim_folder, time_model, ass_var, no_mems
):
    """
    Prepares posterior NetCDF files to a FARM standard by adjusting reference time,
    calendar format, and appending the ass_var variable to the FARM prior file.
    """
    for mem in range(no_mems):
        try:
            posterior_file = (
                f"{output_sim_folder}/ic_g1_posterior_{time_model}_{mem}.nc"
            )

            tmp_dir = Path(output_sim_folder) / "tmp"

            tmp0_posterior = tmp_dir / "tmp0.nc"
            tmp1_posterior = tmp_dir / "tmp1.nc"
            Path(tmp0_posterior).parent.mkdir(parents=True, exist_ok=True)

            # prior in FARM format
            prior_farm_folder = Path(
                f"{path_manager.path_data}/OUTPUT_{mem}/OUT/prior/"
            )
            prior_farm_folder.mkdir(parents=True, exist_ok=True)
            # breakpoint()
            prior_from_farm_file = prior_farm_folder / f"ic_g1_{time_model}.nc"

            # first to delete variable ass_var
            result_tmp = Path(
                f"{path_manager.path_data}/OUTPUT_{mem}/OUT/ic_g1_{time_model}_tmp.nc"
            )

            result = Path(
                f"{path_manager.path_data}/OUTPUT_{mem}/OUT/ic_g1_{time_model}.nc"
            )

            with open("subprocess_output.log", "w") as log_file:
                logger.info("0: Storing prior elsewhere")
                # breakpoint()

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
                # breakpoint()
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
                # ds_result = ds.drop_vars(ass_var)
                ds_tmp1_posterior = xr.open_dataset(str(tmp1_posterior))
                ds[ass_var].values = ds_tmp1_posterior[ass_var].values
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
                # Step 5: put in list files in to_DART dir to be deleted at the
                # end of the execution
                to_dart_path = Path(f"{path_manager.path_data}/to_DART/")
                to_dart_files = list(to_dart_path.rglob("ic_g1*.nc"))

                # logger.info("4: delete var from FARM prior")
                # subprocess.run(
                #    ["cdo", f"delname,{ass_var}",str(prior_from_farm_file), str(result)],
                #    stdout=log_file,
                #    stderr=log_file,
                #    check=True,
                # )

                # logger.info("5: append assimilated c_SO2 from posterior DART")
                # subprocess.run(
                #    ["ncks", "-A","-v", str(ass_var), str(tmp1_posterior), str(result)],
                #    stdout=log_file,
                #    stderr=log_file,
                #    check=True,
                # )
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed: {e.cmd}")
            logging.error(f"Error output: {e.output}")
            logging.info(
                f"something went wrong restore at least the original core in {result} to continue farm run"
            )
            with open("subprocess_output.log", "a") as log_file:
                subprocess.run(
                    ["mv", prior_from_farm_file, result],
                    stdout=log_file,
                    stderr=log_file,
                    check=True,
                )
            raise
        finally:
            # Cleanup: Remove temporary files
            temp_files = [
                tmp0_posterior,
                tmp1_posterior,
                result_tmp,
                prior_from_farm_file,
            ]  # + to_dart_files
            for temp_file in temp_files:
                temp_file.unlink(missing_ok=True)


def prepare_farm_to_dart_nc(
    path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model, no_mems
):

    os.makedirs(path_manager.path_data / "/temp", exist_ok=True)

    for mem in range(no_mems):
        try:
            logging.info
            meteo_file = f'/gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmwf_d0_g1_{timestamp_farm.strftime("%Y%m%d")}.nc'

            temp_output_meteo = path_manager.path_data / "/temp/output_meteo.nc"
            temp_output_meteo_plus1 = (
                path_manager.path_data / "/temp/output_meteo_plus1.nc"
            )
            temp_output_meteo_selected = (
                path_manager.path_data / "/temp/output_meteo_selected.nc"
            )
            logging.info("Running ARCONV for the members")
            arconv_input_file = (
                path_manager.path_data
                / f'/OUTPUT_{mem}/OUT/ic_g1_{rounded_timestamp.strftime("%Y%m%d%H")}.nc'
            )
            arconv_output_file = path_manager.path_data / "/temp/arconv_output.nc"

            final_concentration_file = (
                path_manager.path_data
                / f"/to_DART/ic_g1_{seconds_model}_{days_model}_{mem}.nc"
            )
            temp_concentration_file = path_manager.path_data / "/to_DART/temp_conc.nc"
            temp1_concentration_file = path_manager.path_data / "/to_DART/temp1_conc.nc"

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
                logger.info(
                    "5: use ncks to append SP, P and T variables to the FARM concentration file"
                )
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
                subprocess.run(
                    ["cp", temp1_concentration_file, final_concentration_file]
                )
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


def modify_yaml_date(file_path, new_date):
    try:
        with open(file_path, "r") as file:
            data = yaml.safe_load(file)
        if "time" in data and "start_time" in data["time"]:
            data["time"]["start_time"] = new_date
        else:
            logger.info("Not changed start_time in yaml config")
            return
        with open(file_path, "w") as file:
            yaml.safe_dump(data, file)

        logger.info(f"Next run starts from {new_date}")

    except Exception as e:
        logger.error(f"An error occurred: {e}")


def submit_and_wait(
    path_manager: PathManager,
    commands_with_directories: list,
    timestamp_farm: str,
    no_mems: int,
    case_dir: str,
    cresco_queue: str,
) -> bool:
    for command, directory in commands_with_directories:
        job_ids = run_command_in_directory_bsub(command, directory)
        time.sleep(10)
    mems_to_rerun = get_list_mems_to_rerun(
        job_ids, path_manager, timestamp_farm, no_mems
    )

    if mems_to_rerun:
        list_mems = [str(mem) for mem in mems_to_rerun]
        replace_nml_template(
            input_nml_path=path_manager.run_submit_farm_template,
            entries_tbr_dict={
                "da_date_start": timestamp_farm,  # .strftime('%Y%m%d%H'),
                "da_date_end": timestamp_farm,  # .strftime('%Y%m%d%H'),
                "@no_mems_list": str(tuple(list_mems)).replace(",", ""),
                "@case_dir": case_dir,
                "@cresco_queue": cresco_queue,
            },
            output_nml_path=directory / command,
        )
        return submit_and_wait(
            path_manager,
            [(command, directory)],
            timestamp_farm,
            no_mems,
            case_dir,
            cresco_queue,
        )  # added case_dir and cresco_queue
    return True


def get_list_mems_to_rerun(
    job_ids: list[str],
    scheduler: Scheduler,
    model_type: ModelType,
    path_manager: Optional[PathManager] = None,
    timestamp_model: Optional[str] = None,
    no_mems: Optional[int] = None
) -> list[int]:

    #datetime_model_p1 = pd.to_datetime(timestamp_model, format="%Y%m%d%H") + timedelta(
    #    hours=1
    #)

    while True:
        running_jobs = []

        for jobid in job_ids:
            if scheduler == Scheduler.SLURM:
                finished = check_job_status_slurm(
                    jobid, which_run=model_type.value.upper()
                ) #da cambiare, ritorna finished anche se é fallito o se é pending !
            else:
                finished = check_job_status_cresco(
                    jobid, which_run=model_type.value.upper()
                )

            if not finished:
                running_jobs.append(jobid)

        if not running_jobs:
            logger.info(f"Jobs {job_ids} have finished")
            return False
        
            #return check_ic_g1_existing(
            #    path_manager=path_manager,
            #    model=model_type,
            #    datetime_model=datetime_model_p1,
            #    no_mems=no_mems,
            #)

        logger.info(f"Jobs still running: {running_jobs}. Waiting...")
        time.sleep(30)


def check_ic_g1_existing(
    path_manager: PathManager,
    model: ModelType,
    datetime_farm: pd.Timestamp,
    no_mems: int,
) -> list[int]:

    mems_to_rerun = []

    for mem in range(no_mems):
        ic_path = path_manager.get_ic_g1_path(
            model=model,
            mem=mem,
            timestamp=datetime_farm,
        )

        if ic_path.exists() and ic_path.stat().st_size > 0:
            logger.info(
                f"{model} | ic_g1 exists for mem {mem} "
                f"({ic_path.stat().st_size} bytes)"
            )
        else:
            logger.warning(f"{model} | ic_g1 missing for mem {mem}")
            mems_to_rerun.append(mem)

    return mems_to_rerun


def ic_g1_not_existing(
    path_manager: PathManager, datetime_farm: pd.Timestamp, no_mems: int
):
    timestamp_farm_p1 = datetime_farm.strftime("%Y%m%d%H")
    file_name = f"ic_g1_{timestamp_farm_p1}.nc"
    mems_to_rerun = []
    for mem in range(no_mems):
        file_path = path_manager.chimere_output_runs_dir(mem)
        if os.path.exists(file_path):
            logger.info(
                f"The core {file_name} for mem {mem} exists in the directory {os.path.getsize(file_path)} bytes"
            )
        else:
            logger.info(
                f"The core {file_name} for mem {mem} does not exist in the directory"
            )
            mems_to_rerun.append(mem)
    return mems_to_rerun


def replace_priorinflation(path_manager: PathManager, timestamp_farm: str):
    """
    Replaces the priorinflation files in the FARM model's working directory.

    Args:
        path_manager (PathManager): Object containing the base path.
        timestamp_farm (str): Timestamp for logging purposes.
    """
    # Define file names
    file_mappings = {
        "output_priorinf_mean.nc": "input_priorinf_mean.nc",
        "output_priorinf_sd.nc": "input_priorinf_sd.nc",
    }

    # Define working directory
    work_path = Path(path_manager.base_path / path_manager.path_filter)

    logger.info(
        f"Starting renaming of priorinflation files for next run: {timestamp_farm}"
    )

    for input_file, output_file in file_mappings.items():
        try:
            # Resolve full paths
            src = work_path / input_file  # input_file has output_*.ncs
            dest = work_path / output_file

            # Rename the file
            if src.exists():
                src.rename(dest)
                logger.info(f"Renamed {src} to {dest}")
            else:
                logger.warning(f"Source file not found: {src}")
        except Exception as e:
            logger.error(f"Failed to rename {input_file} to {output_file}: {e}")

    logger.info(f"Finished renaming priorinflation files for run: {timestamp_farm}")


def filter_dates(dates, mode):
    """Filter hourly dates according to backup mode."""
    mode = str(mode).strip().lower()

    if mode == "5daily":
        logger.info("Saving ic_g1 at 00:00 for days: [1, 5, 10, 15, 20, 25, 30]")
        keep_days = [1, 5, 10, 15, 20, 25, 30]
        return dates[(dates.hour != 0) | (~dates.day.isin(keep_days))]

    elif mode == "daily":
        logger.info("Saving ic_g1 at 00:00 for each day")
        return dates[dates.hour != 0]

    elif mode == "hourly":
        logger.info("Saving ic_g1 for every hour (no exclusions)")
        return dates

    else:
        logger.warning(
            f"Unknown days_backup option '{mode}'. Defaulting to 'daily'. "
            f"Accepted values: ['daily', '5daily', 'hourly']"
        )
        return dates[dates.hour != 0]


def check_job_status_slurm(job_id, **kwargs):
    """
    Returns True if the job has finished (COMPLETED, FAILED, CANCELLED, etc.)
    Returns False if the job is still RUNNING or PENDING
    """
    which_run = kwargs.get("which_run", None)

    try:
        result = subprocess.check_output(
            ["squeue", "-j", str(job_id), "-h"],
            stderr=subprocess.DEVNULL,
            text=True,
        )

        if result.strip():
            logger.info(f"{which_run} {job_id}: status RUNNING/PENDING")
            return False  # still in queue
        else:
            logger.info(f"{which_run} {job_id}: status FINISHED")
            return True  # finished

    except subprocess.CalledProcessError:
        # squeue error or job disappeared → finished
        logger.info(f"{which_run} {job_id}: status FINISHED (not in squeue)")
        return True


def get_list_mems_to_rerun_slurm(
    job_ids: list[str],
    path_manager: PathManager,
    timestamp_model: str,
    no_mems: int,
) -> list[int]:

    datetime_farm_p1 = pd.to_datetime(timestamp_model, format="%Y%m%d%H") + timedelta(
        hours=1
    )

    while True:
        running_jobs = []

        for jobid in job_ids:
            # SLURM: job still in squeue → still running
            if not check_job_status_slurm(jobid):
                running_jobs.append(jobid)

        if not running_jobs:
            logger.info(f"SLURM jobs {job_ids} have finished")
            return ic_g1_not_existing(
                path_manager,
                datetime_farm_p1,
                no_mems,
            )

        logger.info(f"SLURM jobs still running: {running_jobs}. Waiting...")
        time.sleep(30)


def safe_symlink(target: Path, link: Path):
    """Create a symlink safely:
    - If it exists and points correctly → do nothing
    - If it is broken → recreate
    - Never overwrite real files/directories
    """
    if not link.parent.exists():
        link.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Try Python 3.11+ approach
        link.symlink_to(target, exist_ok=True)
        logger.info(f"Symlink ensured: {link} -> {target}")
    except TypeError:
        # Older Python (<3.11) doesn't support exist_ok
        try:
            link.symlink_to(target)
            logger.info(f"Symlink created: {link} -> {target}")
        except FileExistsError:
            if link.is_symlink():
                if not link.exists():
                    # Broken symlink → fix it
                    logger.info(f"Broken symlink detected. Recreating: {link}")
                    link.unlink()
                    link.symlink_to(target)
                else:
                    # Symlink exists and is valid → nothing to do
                    logger.info(f"Symlink already exists and is valid: {link}")
            else:
                # Path exists but is a file or directory → skip
                logger.warning(f"{link} exists and is not a symlink. Skipping. \
                               Most likely this run is a continuation of previous runs. If not, \
                               this run then needs to start from the contrul run from this current \
                               time: please remove from the ensemble folders the specific end files \
                               to allow the linking")
                
def check_and_clean_broken_links(run_dir: Path) -> bool:
    """
    Check for broken symlinks in run_dir (depth <= 2), remove them,
    and report if any were found.

    Returns:
        True if broken links were found (submission should be skipped)
        False otherwise
    """

    logger.info(f">> Checking links...")

    broken_found = False

    # Scan up to depth 2
    base_depth = len(run_dir.resolve().parts)

    for path in run_dir.rglob("*"):

        # Limit depth to maxdepth=2
        if len(path.resolve().parts) - base_depth > 2:
            continue

        if path.is_symlink():

            # Broken if target does not exist
            if not path.exists():

                logger.info(f"   [!] BROKEN LINK FOUND: {path}")

                try:
                    target = path.resolve(strict=False)
                    logger.info(f"      Points to: {target}")
                except Exception:
                    logger.info("      Points to: <unresolvable>")

                # Remove broken symlink
                path.unlink()
                logger.info("       >>> Unlinked broken reference!")

                broken_found = True

    return broken_found

def from_liststr_to_listdict(ensemble_list: list[str], labels: list[str]) -> list[dict]:
    """
    Convert a list of colon-separated strings into a list of dictionaries.
    Automatically assigns MemberID.
    Maps parts to provided labels, creates ExtraID_n if needed.

    Args:
        ensemble_list: List like ["00", "00:01", "00:01:02", ...]
        labels: List of labels like ["EmisID", "MeteoID", "ChemID"]

    Returns:
        List of dictionaries with integer values.
    """
    ensemble_dicts = []

    for member_id, item in enumerate(ensemble_list):
        parts = item.split(":")

        entry = {
            "MemberID": member_id
        }

        for i, value in enumerate(parts):

            if i < len(labels):
                key = labels[i]
            else:
                key = f"ExtraID_{i - len(labels) + 1}"
                logger.info("only MeteoID and EmisID are supported for the moment: specific functions (in pipeline and paths) need to be created to allow more")
                logger.info("Extra IDs were given: untracked")

            entry[key] = int(value)

        ensemble_dicts.append(entry)

    return ensemble_dicts



