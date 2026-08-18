from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, List, Tuple, Protocol
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
from pipeline_errors import SchedulerError
from mimesi_types import ModelType, Scheduler
import shlex
from scheduler import submit_job, wait_for_slurm_jobs
from typing import Iterable, Optional, Union
from dataclasses import dataclass
from typing import List, Optional, Tuple
import shlex
from typing import TYPE_CHECKING
from pipeline_time import AssimWindow, TimeManager

if TYPE_CHECKING:
    from pipelines.chimere2017.paths import Chimere2017Paths


class FarmPathManagerLike(Protocol):
    base_path: Path
    path_filter: Path
    run_submit_farm_template: Path


class ChimereOutputPathsLike(Protocol):
    def chimere_output_runs_dir(self, mem: int) -> Path:
        ...

logger = logging.getLogger(__name__)

def format_bytes(num_bytes: int) -> str:
    units = ["bytes", "KB", "MB", "GB", "TB", "PB"]
    size = float(num_bytes)
    for unit in units:
        if size < 1024.0 or unit == units[-1]:
            if unit == "bytes":
                return f"{int(size)} {unit}"
            return f"{size:.2f} {unit}"
        size /= 1024.0

@dataclass
class CommandSpec:
    command: str  # executable/script name
    directory: Path  # where to run it
    args: Optional[List[str]] = None

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
        #if not job_ids:
        #    job_ids = re.findall(r"\d+", stdout)
        
        if not job_ids:
            job_ids = re.findall(r"Submitted Batch Session (\d+)$", stdout, re.MULTILINE)
            logger.info(f"Found: {job_ids}")

        #if not job_ids:
        #    raise RuntimeError(
        #        "No Slurm job IDs found in output.\n" "Expected lines like: ID:<jobid>"
        #    )
        
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

def monitor_job_dart(self, job_id):
        logger.info(f"Monitoring job {job_id}")
        job_id = job_id.strip()[1:-1]

        while True:
            if check_job_status_cresco(job_id, which_run="FARM"):
                print("Job completed successfully.")
                # Handle successful job completion: move files
                self.move_analysis_files()
                replace_priorinflation(
                    self.path_manager,
                    self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                )
                break
            else:
                print("Job is still running. Waiting...")
                time.sleep(10)

def submit_irene(spec: CommandSpec) -> str:
    rc, job_id = run_command_in_directory(spec)

    if rc != 0:
        raise SchedulerError(
            f"Submission command failed: {spec.command} " f"(return code {rc})"
        )
    if not job_id:
        logger.info(f"No job id returned by command {spec.command}")
        logger.info(f"No monitoring will be performed")
        return None
    
    time.sleep(5)
    #breakpoint()
    logger.info(f"[TGCC-IRENE] Submitted job with ID:{job_id}")

    return job_id[0]

def submit_and_wait_cineca(
    spec: CommandSpec,
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


def searchFile(t1, t2, listing):
    orbit_filename = listing[["filename", "start_time"]][
        (listing["start_time"] >= t1) & (listing["start_time"] <= t2)
    ]
    return orbit_filename


def replace_nml_template(
    input_nml_path: str,
    entries_tbr_dict: dict,
    output_nml_path: str,
):
    if not isinstance(entries_tbr_dict, dict):
        raise TypeError("'entries_tbr_dict' must be a dictionary")

    try:
        with open(input_nml_path, "r") as f1:
            input_nml = f1.read()
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Input template not found: {input_nml_path}") from e
    except Exception as e:
        raise RuntimeError(f"Error reading input template {input_nml_path}") from e

    for key, value in entries_tbr_dict.items():
        input_nml = input_nml.replace(key, str(value))

    try:
        with open(output_nml_path, "w") as f2:
            f2.write(input_nml)
        os.chmod(output_nml_path, 0o775)
    except Exception as e:
        raise RuntimeError(f"Error writing output file {output_nml_path}") from e

    logger.info(
        "Replacement %s → %s completed successfully.",
        os.path.basename(input_nml_path),
        os.path.basename(output_nml_path),
    )
    logger.debug(
        "Replacement %s → %s completed successfully.",
        input_nml_path,
        output_nml_path,
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
    path_manager: FarmPathManagerLike,
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
        job_ids=job_ids,
        scheduler=Scheduler.LSF,
        model_type=ModelType.FARM,
        path_manager=path_manager,
        timestamp_model=timestamp_farm,
        no_mems=no_mems,
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
    *,
    job_ids: list[str],
    scheduler: Scheduler,
    model_type: ModelType,
    path_manager: Optional[object] = None,
    timestamp_model: Optional[str] = None,
    no_mems: Optional[int] = None,
    run_hours: int
) -> list[int]:
    """Block until all jobs in job_ids have finished.
    Return list of ensemble members that did not produce valid outputs.
    Empty list means success.
    """

    while True:
        running_jobs = []

        for jobid in job_ids:
            if scheduler == Scheduler.SLURM:
                finished = check_job_status_slurm(
                    jobid, which_run=model_type.value.upper()
                )
            else:
                finished = check_job_status_cresco(
                    jobid, which_run=model_type.value.upper()
                )

            if not finished:
                running_jobs.append(jobid)

        if not running_jobs:
            logger.info(f"Jobs {job_ids} have finished")
            logger.info(f"Checking if runs failed ...")
            if path_manager is None or timestamp_model is None or no_mems is None:
                logger.warning(
                    "Skipping output-file validation because path_manager, "
                    "timestamp_model, or no_mems is missing."
                )
                return []

            timestamp_dt = pd.to_datetime(timestamp_model, format="%Y%m%d%H")

            #datetime_model_p1 = timestamp_dt + timedelta(hours=1)
            return check_ic_files_exist(
                path_manager=path_manager,
                model=model_type,
                datetime_model=timestamp_dt,
                no_mems=no_mems,
                run_hours= run_hours
            )

        logger.info(f"Jobs still running: {running_jobs}. Waiting...")
        time.sleep(30)

def check_restart_files_exist_irene(
    ic_paths: list,
    model: ModelType,
    no_mems: int,
) -> list[int]:

    mems_to_rerun = []

    for mem in range(1, no_mems+1):
        if ic_paths[mem-1].exists() and ic_paths[mem-1].stat().st_size > 0:
            logger.info(
                f"{model} | restart_file exists for mem {mem}: {ic_paths[mem-1]}"
                f"({ic_paths[mem-1].stat().st_size} bytes)"
            )
        else:
            logger.warning(f"{model} | resatrt_file is missing for mem {mem}: {ic_paths[mem-1]}")
            mems_to_rerun.append(mem)

    return mems_to_rerun


def check_restart_files_exist(
    path_manager: Chimere2017Paths,
    model: ModelType,
    no_mems: int,
    datetime_model: Optional[str] = None,
) -> list[int]:

    if datetime_model is None:
        raise ValueError("datetime_model is required when checking restart files")


    mems_to_rerun = []

    for mem in range(no_mems):
        end_file = path_manager.get_chimere_output_path(
            model=model,
            mem=mem,
            timestamp=datetime_model,
            prefix="end",
            offset=1,
        )

        if end_file.exists() and end_file.stat().st_size > 0:
            logger.info(
                f"{model} | restart_file {end_file} exists for mem {mem} "
                f"({format_bytes(end_file.stat().st_size)})"
            )
        else:
            logger.warning(f"{model} | restart_file {end_file} missing for mem {mem}")
            mems_to_rerun.append(mem)

    return mems_to_rerun


def check_ic_files_exist(
    path_manager: Chimere2017Paths,
    model: ModelType,
    no_mems: int,
    datetime_model: pd.Timestamp,
    run_hours: int
) -> list[int]:
    mems_to_rerun = []

    for mem in range(no_mems):
        ic_path = path_manager.get_chimere_output_path(
            model=model,
            mem=mem,
            timestamp=datetime_model,
            prefix='end',
            offset=run_hours
        )

        if ic_path.exists() and ic_path.stat().st_size > 0:
            logger.info(
                f"{model} | ic_g1 exists for mem {mem} "
                f"({format_bytes(ic_path.stat().st_size)})"
            )
        else:
            logger.warning(f"{model} | ic_g1 missing for mem {mem}: {ic_path}")
            mems_to_rerun.append(mem)

    return mems_to_rerun


def ic_g1_not_existing(
    path_manager: ChimereOutputPathsLike, datetime_farm: pd.Timestamp, no_mems: int
):
    timestamp_farm_p1 = datetime_farm.strftime("%Y%m%d%H")
    file_name = f"ic_g1_{timestamp_farm_p1}.nc"
    mems_to_rerun = []
    for mem in range(no_mems):
        file_path = path_manager.chimere_output_runs_dir(mem)
        if os.path.exists(file_path):
            logger.info(
                f"The core {file_name} for mem {mem} exists in the directory {format_bytes(os.path.getsize(file_path))}"
            )
        else:
            logger.info(
                f"The core {file_name} for mem {mem} does not exist in the directory"
            )
            mems_to_rerun.append(mem)
    return mems_to_rerun


def replace_priorinflation(path_filter: Path, timestamp: str):
    """
    Promote prior inflation outputs to inputs for the next DART cycle.

    Args:
        path_filter (Path): DART/filter work directory.
        timestamp (str): Timestamp for logging purposes.
    """
    # Define file names
    file_mappings = {
        "output_priorinf_mean.nc": "input_priorinf_mean.nc",
        "output_priorinf_sd.nc": "input_priorinf_sd.nc",
    }

    # Define working directory
    work_path = Path(path_filter)

    logger.info(
        f"Starting renaming of prior inflation files for next run: {timestamp} (work_dir={work_path})"
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

    logger.info(f"Finished renaming prior inflation files for run: {timestamp}")


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
    path_manager: ChimereOutputPathsLike,
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


def monitor_job_status(
    job_ids: list[str],
    scheduler: Scheduler,
    model_type: Optional[ModelType] = None,
    ):

    while True:
        running_jobs = []

        for jobid in job_ids:
            if scheduler == Scheduler.SLURM:
                finished = check_job_status_slurm(
                    jobid#, which_run=model_type.value.upper()
                )
            else:
                finished = check_job_status_cresco(
                    jobid#, which_run=model_type.value.upper()
                )

            if not finished:
                running_jobs.append(jobid)

        if not running_jobs:
            logger.info(f"Jobs {job_ids} have finished")
            return

        logger.info(f"Jobs still running: {running_jobs}. Waiting...")
        time.sleep(15)



def safe_symlink(target: Path, link: Path, do_copy: Optional[bool]=False):
    """Create a symlink safely:
    - If it exists and points correctly → do nothing
    - If it is broken → recreate
    - If it points to the wrong target → recreate
    - Never overwrite real files/directories
    """
    if not link.parent.exists():
        link.parent.mkdir(parents=True, exist_ok=True)

    try:
        if not do_copy:
            link.symlink_to(target, exist_ok=True)
            logger.info(f"Symlink ensured: {link} -> {target}")
        else:
            import shutil
            if link.is_symlink():
                link.unlink()
            shutil.copy2(target, link)
            logger.info(f"Copied ensured: {link} of {target}")
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
                    # Symlink exists and is valid → verify target
                    try:
                        current_target = link.resolve(strict=False)
                        expected_target = target.resolve()

                        if current_target != expected_target:
                            logger.info(
                                f"Symlink points to wrong target "
                                f"({current_target}). Recreating: {link}"
                            )
                            link.unlink()
                            link.symlink_to(target)
                        else:
                            logger.info(f"Symlink already exists and is correct: {link}")
                    except OSError as e:
                        logger.error(f"Error verifying symlink {link}: {e}")
            else:
                # Path exists but is a file or directory → skip
                logger.warning(f"{link} exists and is not a symlink. Skipping linking.")
                
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
                    logger.info(f"      Points to target: {target}")
                    logger.info(f"   [!] CHECK IF TARGET EXITS")
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
            "MemberID": member_id + 1  # Start MemberID from 1 for clarity
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

def compute_hourly(data: xr.Dataset, time: int, path_saving_data: Path, path_saving_list: Optional[Path]=None) -> Path:
    data_sel = data.sel(Time=slice(time, time+2)) #to keep Time dimension
    if len(data_sel.Times.values) == 0: #when time is saved as float (isel drops Time even with drop=False)
        data_sel = data.sel(Time=slice(data.Time.values[time], data.Time.values[time+1]))
    
    path_saving_data.parent.mkdir(parents=True, exist_ok=True)
    data_sel.to_netcdf(path_saving_data)
    if path_saving_list:
        path_saving_list.write_text("1\n" + str(path_saving_data) + "\n")
        logger.info("Hourly dataset computed and listing created")
    else:
        logger.info("Hourly dataset computed")

def cut_block(data: xr.Dataset, time: int, hours: int, path_saving_data: Path, path_saving_list: Optional[Path]=None) -> Path:
    data_sel = data.sel(Time=slice(time, time+hours+1)) #to keep Time dimension
    if len(data_sel.Times.values) == 0: #when time is saved as float (isel drops Time even with drop=False)
        data_sel = data.sel(Time=slice(data.Time.values[time], data.Time.values[time+hours]))
    
    path_saving_data.parent.mkdir(parents=True, exist_ok=True)
    data_sel.to_netcdf(path_saving_data)
    if path_saving_list:
        path_saving_list.write_text("1\n" + str(path_saving_data) + "\n")
        logger.info("Hourly dataset computed and listing created")
    else:
        logger.info("Hourly dataset computed")


def add_missing_variable_and_select_last_timestep(no_mems: int, var_to_add: str, domain: str, out_file_func: Callable, orig_file_func: Callable, to_dart_func: Callable, date_ymdH: str, NHOURS: int, date_ymdH_chim: str, skip_model_part: bool):
    for mem in range(1, no_mems+1):
        out_file_name=out_file_func(mem, date_ymdH_chim, NHOURS)
        orig_file_name=orig_file_func(mem, domain, date_ymdH_chim, NHOURS)
        to_dart_file_name=to_dart_func(mem, date_ymdH) 
        
        if not os.path.exists(to_dart_file_name) and not skip_model_part:
            logger.info(f'From {out_file_name} created {to_dart_file_name}')
            """#ds = xr.open_dataset(f'/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/first_tests_202002_06-15/ENS{i}/chim_ENS{i}_2020021413_1_out.nc')
            #vcmeteo = xr.open_dataset(f'/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/first_tests_202002_06-15/ENS{i}/exdomout_2020021413_1_EUROCOMEX3.nc')
            ds = xr.open_dataset(out_file_name)
            meteo = xr.open_dataset(orig_file_name)

            meteo = meteo.rename({"Time": "time_counter", "south_north": "y", "west_east": "x" })
            meteo = meteo.isel(time_counter=slice(0,1))
            meteo = meteo.assign_coords(time_counter=ds.time_counter, y=ds.y, x=ds.x)

            ds[var_to_add] = meteo.psfc.astype("float32")
            #ds.to_netcdf(f'chim_ENS{i}_2020021413_1_out_psfc_float.nc')
            ds.to_netcdf(out_file_name, mode="a")

            ds.close()
            meteo.close()"""

            with xr.open_dataset(out_file_name) as ds:
                if var_to_add not in ds.data_vars:
                    with xr.open_dataset(orig_file_name) as meteo:
                        meteo_sub = (
                            meteo.rename({"Time": "time_counter", "south_north": "y", "west_east": "x"})
                                .isel(time_counter=slice(0, -1))
                                .assign_coords(time_counter=ds.time_counter, y=ds.y, x=ds.x)
                        )

                        ds[var_to_add] = meteo_sub[var_to_add].astype("float32")
                        ## start changes ##
                        ds = ds.isel(time_counter=-1) 
                        #tmp = str(out_file_name) + ".tmp"
                        #ds.to_netcdf(tmp)
                        #os.replace(tmp, out_file_name)
                        ds.to_netcdf(to_dart_file_name)
        else:
            logger.info(f'File exists: {to_dart_file_name}')


def write_dart_filter_list(list_file_func: Path, out_file_func: Callable, no_mems: int, **kwargs):
    try:    
        list_file_func.write_text(
            "\n".join(
                str(out_file_func(mem, **kwargs))
                for mem in range(1, no_mems+1)
            ) + "\n"
        )
        logger.info(f"Wrote: {list_file_func}")
    except:
        logger.warning(f"Writing of the following failed: {list_file_func}")

def _get_nc4_encoding(dataset: xr.Dataset) -> dict:
    """Keep only the netCDF4-compression keys present in the source encoding."""
    valid_keys = {
        'shuffle', 'zlib', 'quantize_mode', 'blosc_shuffle',
        'fletcher32', 'least_significant_digit', 'endian',
        'complevel', 'szip_coding', 'dtype',
        'szip_pixels_per_block', 'compression',
        'significant_digits', '_FillValue'
    }
    return {
        var: {key: value for key, value in dataset[var].encoding.items() if key in valid_keys}
        for var in dataset.variables
    }


def _write_dataset_with_encoding(dataset: xr.Dataset, target: Path, *, unlimited_dims=None):
    """Write dataset to disk preserving compression metadata only where supported."""
    encoding = _get_nc4_encoding(dataset)
    tmp_path = Path(str(target) + ".tmp")
    dataset.to_netcdf(
        tmp_path,
        encoding=encoding,
        unlimited_dims=unlimited_dims,
    )
    os.replace(tmp_path, target)


def _compute_emission_ratio(val_correct_da: xr.DataArray, val_prev_da: xr.DataArray, emis_prev_da: xr.DataArray):
    """Compute the per-pixel correction ratio while tolerating vertical-dimension mismatches."""
    logger.info("Computing the ratio...")
    #prepare a matrix of ones
    ratio = np.ones_like(emis_prev_da.values)

    if "bottom_top" in val_correct_da.dims and "bottom_top" in emis_prev_da.dims:
        # Alignment check vertical levels (es. 20 vs 7 levels)
        if val_correct_da.sizes["bottom_top"] != emis_prev_da.sizes["bottom_top"]:
            logger.info(
                " Vertical level mismatch detected "
                "(DART: %s, EMIS: %s). Scaling LEVEL 0 only.",
                val_correct_da.sizes["bottom_top"],
                emis_prev_da.sizes["bottom_top"],
            )
            v_correct_lvl0 = val_correct_da.isel(bottom_top=0).values
            v_prev_lvl0 = val_prev_da.isel(bottom_top=0).values
            ratio_lvl0 = np.divide(
                v_correct_lvl0,
                v_prev_lvl0,
                out=np.zeros_like(v_correct_lvl0), #fallback if v_prev_lvl0 == 0
                where=v_prev_lvl0 != 0,
            )
            #insert lvl0 ratio in the matrix
            ratio[0, :, :] = ratio_lvl0
            #if on the other end, we would like to use lvl0 ratio for all levels
            #ratio[:] = ratio_lvl0[np.newaxis, :, :]
            logger.info(f"Emission ratio computed: ready to apply it. Mean: {ratio.mean()}")
            return ratio
    
    # if there is no vertical mismatch (or bottom_top does not exist), ratios are specific of and computed at each level
    ratio = np.divide(
        val_correct_da.values,
        val_prev_da.values,
        out=np.zeros_like(val_correct_da.values), #fallback if val_prev_da.values == 0
        where=val_prev_da.values != 0,
    )
    logger.info(f"Emission ratio computed: ready to apply it. Mean: {ratio.mean()}")
    return ratio


def _load_or_init_ratio_memory(ratio: np.ndarray, ratio_memory_file: Optional[Path], *, extend_from_prev_slot: bool, persistence_until_next_day: bool, curr_time: pd.Timestamp): 
    #when calling the function, since * is given as an input, the boolean inputs (since they follow *) need to be explicitly called out 
    """Load existing memory or initialize it, preserving the legacy persistence semantics."""
    if extend_from_prev_slot == False and persistence_until_next_day == False:
        logger.info("extend_from_prev_slot=False and persistence_until_next_day=False: avoiding any type of persistency)")
        if ratio_memory_file is not None and ratio_memory_file.exists():
            logger.info("Removing persistent orbit memory that is probably a leftover")
            ratio_memory_file.unlink()
        return (
            np.ones_like(ratio),
            np.zeros_like(ratio, dtype=bool),
            np.full(ratio.shape, np.datetime64("NaT"), dtype="datetime64[ns]"),
        )

    if ratio_memory_file is not None and ratio_memory_file.exists():
        logger.info(f"Loading orbit memory: {ratio_memory_file}")
        with xr.open_dataset(ratio_memory_file) as mem_ds:
            ratio_prev = mem_ds["ratio"].values
            footprint_prev = mem_ds["footprint"].values.astype(bool)
            assim_prev = mem_ds["assim_time"].values.astype("datetime64[ns]")
        valid_times = assim_prev[~np.isnat(assim_prev)]
        # Condition: If there are no valid times (all NaT) OR the latest valid time is <= curr_time
        if valid_times.size == 0 or pd.Timestamp(valid_times.max()) <= curr_time:
            if not valid_times.size == 0:
                logger.info(f"Latest valid time in the loaded file: {pd.Timestamp(valid_times.max())}")
            return ratio_prev, footprint_prev, assim_prev
        else:
            logger.info(f"Ignoring orbit memory {ratio_memory_file} (assim_time > curr_time)")

    logger.info("No previous orbit memory found.")
    return (
        np.ones_like(ratio),
        np.zeros_like(ratio, dtype=bool),
        np.full(ratio.shape, np.datetime64("NaT"), dtype="datetime64[ns]"),
    )


def _get_pixel(obj, h=None, x=35, y=34):
    """
    Legge DIRETTAMENTE lo scalare del pixel (x, y) dall'oggetto Xarray o Numpy,
    senza eseguire alcuna operazione aritmetica.
    """
    if isinstance(obj, (xr.DataArray, xr.Dataset)):
        da = obj if isinstance(obj, xr.DataArray) else obj[list(obj.data_vars)[0]]
        if h is not None and "Time" in da.dims:
            da = da.isel(Time=h)
        if "bottom_top" in da.dims:
            da = da.isel(bottom_top=0)
        
        if "south_north" in da.dims and "west_east" in da.dims:
            return float(da.isel(south_north=y, west_east=x).values)
        elif "y" in da.dims and "x" in da.dims:
            return float(da.isel(y=y, x=x).values)
        else:
            return float(da.values[..., y, x])
            
    elif isinstance(obj, np.ndarray):
        if obj.ndim == 3:
            return float(obj[0, y, x])
        elif obj.ndim == 2:
            return float(obj[y, x])
        else:
            return float(obj[..., y, x])
    else:
        return float(obj)


def _record_paris_tracking(
    tracking_file: Path,
    call_time: pd.Timestamp,
    pollutant: str,
    sim_start_time: pd.Timestamp,
    sim_end_time: pd.Timestamp,
    updated_times: list,
    vals_before: list,
    vals_after: list,
    ratios_applied: list
):
    """
    Registra sul CSV lo stato grezzo estrapolato dai dataset.
    """
    time_cols = pd.date_range(sim_start_time, sim_end_time, freq="1h")
    col_names = ["call_time", "pollutant", "type"] + [t.strftime("%Y-%m-%d %H:00") for t in time_cols]

    if tracking_file and tracking_file.exists():
        df_old = pd.read_csv(tracking_file)
        prev_after = df_old[(df_old["pollutant"] == pollutant) & (df_old["type"] == "AFTER")]
        if not prev_after.empty:
            base_before = prev_after.iloc[-1].to_dict()
            base_after = prev_after.iloc[-1].to_dict()
        else:
            base_before = {c: np.nan for c in col_names}
            base_after = {c: np.nan for c in col_names}
    else:
        df_old = pd.DataFrame(columns=col_names)
        base_before = {c: np.nan for c in col_names}
        base_after = {c: np.nan for c in col_names}

    row_before = base_before.copy()
    row_after = base_after.copy()
    row_ratio = {c: 1.0 for c in col_names}

    call_str = call_time.strftime("%Y-%m-%d %H:%M")
    row_before.update({"call_time": call_str, "pollutant": pollutant, "type": "BEFORE"})
    row_after.update({"call_time": call_str, "pollutant": pollutant, "type": "AFTER"})
    row_ratio.update({"call_time": call_str, "pollutant": pollutant, "type": "RATIO"})

    for t, b, a, r in zip(updated_times, vals_before, vals_after, ratios_applied):
        t_str = t.strftime("%Y-%m-%d %H:00")
        if t_str in col_names:
            row_before[t_str] = b
            row_after[t_str] = a
            row_ratio[t_str] = r

    new_rows = pd.DataFrame([row_before, row_after, row_ratio])
    df_updated = pd.concat([df_old, new_rows], ignore_index=True)
    
    if tracking_file:
        tracking_file.parent.mkdir(parents=True, exist_ok=True)
        df_updated.to_csv(tracking_file, index=False)


def update_pollutant_in_end(
    dart_file: Path,
    end_file: Path,
    out_file: Path,
    emis_file: Path,
    emis_file_preIM: Path,
    dart_in_file: Path,
    next_emis_file: Path,
    ratio_memory_file: Path,
    pollutant: str,
    next_slot_time: pd.Timestamp,
    use_ratios_avg: bool = False,
    persistence_until_next_day: bool = True,
    extend_from_prev_slot: bool = True,
    damping_active: bool = False,
    hours_forward: int = 6,
    is_first_mem: bool = False,
    sim_start_time: pd.Timestamp = None,
    sim_end_time: pd.Timestamp = None,
    tracking_file: Path = None,
    paris_x: int = 35,
    paris_y: int = 34,
):
    """
    Apply DART-derived emission corrections using baseline snapshots to prevent re-correction.
    
    CRITICAL PATTERN: Capture a clean copy of emissions BEFORE any corrections, then apply
    ratios to this baseline for ALL subsequent hours/days. This ensures idempotent corrections
    and prevents cascading amplification (R² instead of R).
    """
    ratio_memory_path = Path(ratio_memory_file) if ratio_memory_file is not None else None

    # === VALIDATION: Ensure all required input files exist ===
    for f in [dart_file, end_file, out_file, emis_file, emis_file_preIM, dart_in_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} is missing")

    tmp_emis = None
    tmp_end = None

    # === LOAD DATASETS ===
    with xr.open_dataset(dart_file, decode_timedelta=True) as dart_ds, \
         xr.open_dataset(end_file, decode_timedelta=True) as end_ds, \
         xr.open_dataset(out_file, decode_timedelta=True) as out_ds, \
         xr.open_dataset(emis_file, decode_timedelta=True) as emis_ds, \
         xr.open_dataset(emis_file_preIM, decode_timedelta=True) as emispreIM_ds, \
         xr.open_dataset(dart_in_file, decode_timedelta=True) as dart_in_ds:

        # === EMISSION CORRECTION LOGIC (for EMISA/EMISB pollutants) ===
        if pollutant in ['EMISA', 'EMISB']:
            # Map EMISA->NO, EMISB->NO2 for tracking in emissions file
            poll_name = 'NO2' if pollutant == 'EMISB' else 'NO'
            logger.info(f"Proceeding updating {poll_name} emissions ...")

            # Extract current timestep boundaries
            t_start = end_ds.Times.isel(Time=-2).item().decode("utf-8").replace("_", " ")
            t_end = end_ds.Times.isel(Time=-1).item().decode("utf-8").replace("_", " ")
            t_start = pd.to_datetime(t_start)
            t_end = pd.to_datetime(t_end)
            idx_start = t_start.hour
            #t_start is basically the time of the emissions corrected by DART
            #t_end is the time of the emissions to be corrected to move forward in time
           
            # Load DART posterior 
            poll_dart = dart_ds[pollutant].load()
            # Load DART prior
            poll_dart_in = dart_in_ds[pollutant].load()
            # Load original emis (perturbed)
            emis_ds[poll_name] = emis_ds[poll_name].load()

            # === COMPUTE CORRECTION RATIO ===
            # Ratio = DART posterior / DART prior (correction factor)
            val_correct_da = xr.where(poll_dart < 0, 0, poll_dart).isel(time_counter=0)
            val_prev_da = xr.where(poll_dart_in < 0, 0, poll_dart_in)
            emis_prev_da = emis_ds[poll_name].isel(Time=idx_start) #used only to inherit shape 
            ratio = _compute_emission_ratio(val_correct_da, val_prev_da, emis_prev_da)
            # === IDENTIFY NEW FOOTPRINT ===
            # Mark pixels where ratio changed significantly (new satellite observation)
            # NB Temporary solution: should be replaced by the real observation mask (only in obs_seq.out/final or in original listed obs)
            footprint_new = np.abs(ratio - 1) > 1e-10
            # === LOAD ORBITAL MEMORY ===
            # Retrieve previous orbit footprints, ratios, and assimilation times to track
            # which pixels were corrected and when if available; elsewhere mimic it with ones on current one
            ratio_prev, footprint_prev, assim_prev = _load_or_init_ratio_memory(
                ratio, #only its shape is used
                ratio_memory_path,
                extend_from_prev_slot=extend_from_prev_slot,
                persistence_until_next_day=persistence_until_next_day,
                curr_time=t_end
            )
            # Update orbit memory with the previous satellite pass
            ratio_memory = ratio_prev.copy()
            footprint_memory = footprint_prev.copy()
            assim_time_memory = assim_prev.copy()
            
            #  === UPDATE MEMORY WITH CURRENT OVERPASS ===
            # Baseline Assignment: Add current overpass on top.
            # Pixels belonging to the current orbit overwrite previous corrections, automatically handling:
            #   - overlap areas: new orbit wins
            #   - new areas: new correction is inserted
            #   - old areas not covered anymore: previous correction is preserved
            ratio_memory[footprint_new] = ratio[footprint_new]
            footprint_memory[footprint_new] = True
            assim_time_memory[footprint_new] = np.datetime64(t_end, "ns")

            # Conditional Refinement: pixels obsevred both now and previously this day get overwritten 
            #   - overlap areas: 
                # If satellite corrected same pixel today: average the ratios
                # If satellite corrected pixel today but corrected it on a different day: use new ratio
            #   - new areas: new correction is inserted
            #   - old areas not covered anymore: previous correction is preserved
            if use_ratios_avg:
                curr_day_np = np.datetime64(t_end, "D")
                prev_days_np = assim_prev.astype("datetime64[D]")
                # valid_prev selects previous pixels with a valid assim time (not NaT) and of the previous orbit (footprint_prev = True)
                valid_prev = ~np.isnat(assim_prev) & footprint_prev #np.isnat Test element-wise for NaT and return result as a boolean array
                # is_same_day has True where the 3 are satisfied: it was previously observed today and it is also part of the current orbit
                is_same_day = valid_prev & footprint_new & (prev_days_np == curr_day_np)
                # for is_same_day we compute the average ratio (this safely overwrites the baseline assignment made above)
                ratio_memory[is_same_day] = (ratio_prev[is_same_day] + ratio[is_same_day]) / 2.0
            
            track_times, track_before, track_after, track_ratios = [], [], [], []
            # Pre IM emissions: necessary to avoid apply percistency ratios over and over again on following hours (only one should be used
            baseline_emis = emispreIM_ds[poll_name].copy(deep=True)  # <-- KEY LINE

            if persistence_until_next_day:
                # === APPLY RATIO TO CURRENT DAY (using baseline snapshot) ===
                # Loop through each hour and apply ratio to baseline emissions.
                # NEVER read from modified array; always use baseline_emis copy.
                current_ratio = np.where(footprint_memory, ratio_memory, 1.0)
                num_hours_curr = emis_ds.sizes.get("Time", 24)
                base_date = t_end.floor("D")

                for h in range(idx_start, num_hours_curr):
                    time_h = base_date + pd.Timedelta(hours=h)
                    # Capture BEFORE from baseline snapshot (not from already-modified array)
                    val_curr_before = _get_pixel(baseline_emis, h=h, x=paris_x, y=paris_y)
                    r_pixel = _get_pixel(current_ratio, x=paris_x, y=paris_y)
                    idx_tuple = tuple(h if dim == "Time" else slice(None) for dim in emis_ds[poll_name].dims)
                    # KEY: Apply ratio to baseline, never to already-corrected values
                    base_values = baseline_emis.isel(Time=h).values
                    emis_ds[poll_name].values[idx_tuple] = base_values * current_ratio
                    val_curr_after = _get_pixel(emis_ds[poll_name], h=h, x=paris_x, y=paris_y)

                    # Record pixel values for CSV tracking
                    if is_first_mem:
                        track_times.append(time_h)
                        track_before.append(val_curr_before)
                        track_after.append(val_curr_after)
                        track_ratios.append(r_pixel)

                # Write modified emissions file (current day)
                tmp_emis = Path(str(emis_file) + ".tmp")
                _write_dataset_with_encoding(emis_ds, tmp_emis, unlimited_dims=['Time'] if 'Time' in emis_ds.dims else None)

                if next_emis_file:
                    logger.info("Applying persistent orbit corrections to next day's emission file: %s", next_emis_file)
                    # === EXTEND CORRECTIONS TO NEXT DAY ===
                    # Carry orbital memory forward: same corrections applied tomorrow using baseline snapshot
                    with xr.open_dataset(next_emis_file, decode_timedelta=True) as next_emis_ds:
                        next_emis_ds[poll_name] = next_emis_ds[poll_name].load()
                        # Capture baseline for next day before any modifications
                        baseline_next_emis = next_emis_ds[poll_name].copy(deep=True)
                        num_hours = next_emis_ds.sizes.get("Time", 24)
                        next_base_date = base_date + pd.Timedelta(days=1)

                        for h in range(0, num_hours):
                            time_h = next_base_date + pd.Timedelta(hours=h)
                            # Capture BEFORE from baseline snapshot (not from already-modified array)
                            val_curr_before = _get_pixel(baseline_next_emis, h=h, x=paris_x, y=paris_y)
                            r_pixel = _get_pixel(current_ratio, x=paris_x, y=paris_y)
                            idx_tuple = tuple(h if dim == "Time" else slice(None) for dim in next_emis_ds[poll_name].dims)
                            # KEY: Apply ratio to baseline of NEXT DAY file
                            base_values = baseline_next_emis.isel(Time=h).values
                            next_emis_ds[poll_name].values[idx_tuple] = base_values * current_ratio
                            val_curr_after = _get_pixel(next_emis_ds[poll_name], h=h, x=paris_x, y=paris_y)

                            # Record pixel values for CSV tracking
                            if is_first_mem:
                                track_times.append(time_h)
                                track_before.append(val_curr_before)
                                track_after.append(val_curr_after)
                                track_ratios.append(r_pixel)

                        # Write modified emissions file (next day)
                        tmp_next = Path(str(next_emis_file) + ".tmp")
                        _write_dataset_with_encoding(next_emis_ds, tmp_next, unlimited_dims=['Time'] if 'Time' in next_emis_ds.dims else None)
                        os.replace(tmp_next, next_emis_file)
            else:
                # === NON-PERSISTENCE MODE (damping or limited extension) ===
                # Apply corrections for a limited time window with optional decay
                if extend_from_prev_slot:
                    end_apply = min(next_slot_time, t_end + pd.Timedelta(hours=hours_forward))
                    time_range = pd.date_range(start=t_end, end=end_apply, freq="1h")
                else:
                    time_range = pd.DatetimeIndex([t_end])

                for t in time_range:
                    # Compute time elapsed since each pixel was corrected
                    elapsed_hours = (np.datetime64(t) - assim_time_memory) / np.timedelta64(1, "h")
                    elapsed_hours = np.where(footprint_memory, elapsed_hours, 999.0)

                    if damping_active:
                        # === DAMPING MODE: ratio decays over time ===
                        # weight decreases linearly from 1 to 0, causing ratio to decay toward 1
                        weight = np.clip(1.0 - elapsed_hours / (hours_forward + 1), 0.0, 1.0)
                        current_ratio = 1.0 + (ratio_memory - 1.0) * weight
                    else:
                        # === EXTENSION MODE: ratio persists, then stops ===
                        # Apply correction within hours_forward window, then revert to no correction
                        current_ratio = np.where(elapsed_hours <= hours_forward, ratio_memory, 1.0)

                    idx_time = t.hour
                    # Capture BEFORE from baseline snapshot (not from already-modified array)
                    val_curr_before = _get_pixel(baseline_emis, h=idx_time, x=paris_x, y=paris_y)
                    r_pixel = _get_pixel(current_ratio, x=paris_x, y=paris_y)
                    idx_tuple = tuple(idx_time if dim == "Time" else slice(None) for dim in emis_ds[poll_name].dims)
                    # KEY: Apply decayed/extended ratio to baseline
                    base_values = baseline_emis.isel(Time=idx_time).values
                    emis_ds[poll_name].values[idx_tuple] = base_values * current_ratio
                    val_curr_after = _get_pixel(emis_ds[poll_name], h=idx_time, x=paris_x, y=paris_y)

                    # Record pixel values for CSV tracking
                    if is_first_mem:
                        track_times.append(t)
                        track_before.append(val_curr_before)
                        track_after.append(val_curr_after)
                        track_ratios.append(r_pixel)

                # Write modified emissions file
                tmp_emis = Path(str(emis_file) + ".tmp")
                _write_dataset_with_encoding(emis_ds, tmp_emis, unlimited_dims=['Time'] if 'Time' in emis_ds.dims else None)

            if is_first_mem and tracking_file and sim_start_time and sim_end_time:
                # === RECORD TRACKING CSV ===
                # Log pixel values (BEFORE/AFTER/RATIO) for validation and debugging
                _record_paris_tracking(
                    tracking_file=tracking_file,
                    call_time=t_end,
                    pollutant=pollutant,
                    sim_start_time=sim_start_time,
                    sim_end_time=sim_end_time,
                    updated_times=track_times,
                    vals_before=track_before,
                    vals_after=track_after,
                    ratios_applied=track_ratios,
                )

            if extend_from_prev_slot or persistence_until_next_day:
                # === PERSIST ORBITAL MEMORY ===
                # Save ratio, footprint, and assimilation time for next update cycle
                if ratio_memory_path is not None:
                    memory_ds = xr.Dataset(
                        {
                            "ratio": (val_correct_da.dims, ratio_memory),
                            "footprint": (val_correct_da.dims, footprint_memory.astype(np.int8)),
                            "assim_time": (val_correct_da.dims, assim_time_memory.astype("datetime64[ns]")),
                        },
                        coords=val_correct_da.coords,
                    )
                    memory_tmp = Path(str(ratio_memory_path) + ".tmp")
                    memory_ds.to_netcdf(memory_tmp)
                    os.replace(memory_tmp, ratio_memory_path)
                else:
                    logger.info("No ratio-memory path provided; skipping orbit memory persistence.")

            logger.info("Emission update completed using pixel-based orbit memory and damping.")
        else:
            # === OTHER POLLUTANTS (non-EMISA/EMISB) ===
            # Direct replacement of end_file values with DART posterior
            poll = dart_ds[pollutant].load()
            airm = out_ds['airm'].sel(time_counter=slice(out_ds['airm'].time_counter.values[-1], out_ds['airm'].time_counter.values[-1])).load()
            poll = xr.where(poll < 0, 0, poll)
            if poll.shape != airm.shape:
                logger.warning("Chimere original out file and dart outputs differ in shape")
            airm = airm.broadcast_like(poll)
            poll_molec = (1e-9 * poll * airm).astype(end_ds[pollutant].dtype)
            poll_molec = poll_molec.rename({'y': 'south_north', 'x': 'west_east', 'time_counter': 'Time'})
            end_ds[pollutant].loc[dict(Time=end_ds.Time[-1])] = poll_molec.isel(Time=-1).values
            logger.info(f"DART's updated {pollutant} successfully replaced into {end_file}")
            tmp_end = Path(str(end_file) + ".tmp")
            end_ds.to_netcdf(tmp_end)

    # === FINALIZE: Replace temporary files with actual output ===
    if pollutant in ['EMISA', 'EMISB']:
        if tmp_emis is not None:
            os.replace(tmp_emis, emis_file)
    else:
        if tmp_end is not None:
            os.replace(tmp_end, end_file)


def save_diff(file_a: Path, file_b: Path, out_path: Path, label: str):
    """Memory-optimized subtraction using Dask lazy-loading."""
    try:
        # 'chunks={}' enables Dask. 
        # You can also specify specific dimensions like chunks={'time': 1, 'lev': 5}
        with xr.open_dataset(file_a, chunks={'time': 1}) as ds_a, \
                xr.open_dataset(file_b, chunks={'time': 1}) as ds_b:
            
            ds_b = ds_b.sel(time_counter=slice(ds_b['airm'].time_counter.values[-1], ds_b['airm'].time_counter.values[-1])) 
            
            # Keep only variables present in BOTH datasets
            common_vars = list(set(ds_a.data_vars) & set(ds_b.data_vars))

            if not common_vars:
                raise ValueError("No common variables between datasets")

            ds_a = ds_a[common_vars]
            ds_b = ds_b[common_vars]

            # This operation is now "lazy" - no math happens yet
            diff = ds_a - ds_b
            #relative_diff = xr.where(ds_b != 0, ((ds_a - ds_b) / ds_b)*100, 0)
            valid_mask = (
                (ds_b != 0)
                & np.isfinite(ds_b)
                & np.isfinite(ds_a)
            )

            safe_num = (ds_a - ds_b).where(valid_mask)
            safe_den = ds_b.where(valid_mask)

            relative_diff = (safe_num / safe_den) * 100
            # Optional cleanup of inf values
            relative_diff = relative_diff.where(np.isfinite(relative_diff))

            # Delete existing files before saving to avoid conflicts
            out_path.unlink(missing_ok=True)
            relative_out_path = out_path.with_suffix('.relative.nc')
            relative_out_path.unlink(missing_ok=True)
            
            # The computation and writing happen chunk-by-chunk to the disk
            diff.to_netcdf(out_path)
            relative_diff.to_netcdf(relative_out_path)
            
            logger.info(f"[{label}] Memory-optimized diff saved to {out_path}")
            logger.info(f"[{label}] Memory-optimized relative diff saved to {relative_out_path}")
    except Exception as e:
        logger.error(f"Failed to compute {label}: {e}")
    
def remove_negative_values(obs_file_path: Path, obs_file_out: Path):
    logger.info(f"Filtering negative values in {obs_file_path} ...")
    try:
        with xr.open_dataset(obs_file_path) as ds:
            ds['vcd'] = ds['vcd'].where(ds['vcd'] >= 0, 0)  # Set negative values to 0
            os.makedirs(obs_file_out.parent, exist_ok=True)
            ds.to_netcdf(obs_file_out) 
        logger.info(f"Negative values filtered successfully in {obs_file_out}.")
    except Exception as e:
        logger.error(f"Error filtering negative values in {obs_file_path}: {e}")
        raise
