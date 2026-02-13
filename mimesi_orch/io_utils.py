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
import re

from paths import PathManager
from scheduler import Scheduler, submit_job, wait_for_slurm_jobs


logger = logging.getLogger(__name__)


def process_member(
    mem, path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model
):
    try:
        meteo_file = f'/gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmwf_d0_g1_{timestamp_farm.strftime("%Y%m%d")}.nc'
        temp_output_meteo = path_manager.path_data / f"temp/output_meteo_{mem}.nc"
        temp_output_meteo_plus1 = (
            path_manager.path_data / f"temp/output_meteo_plus1_{mem}.nc"
        )
        temp_output_meteo_selected = (
            path_manager.path_data / f"temp/output_meteo_selected_{mem}.nc"
        )
        arconv_input_file = (
            path_manager.path_data
            / f'OUTPUT_{mem}/OUT/ic_g1_{rounded_timestamp.strftime("%Y%m%d%H")}.nc'
        )
        arconv_output_file = path_manager.path_data / f"temp/arconv_output_{mem}.nc"

        final_concentration_file = (
            path_manager.path_data
            / f"to_DART/ic_g1_{seconds_model}_{days_model}_{mem}.nc"
        )
        temp_concentration_file = path_manager.path_data / f"to_DART/temp_conc_{mem}.nc"
        temp1_concentration_file = (
            path_manager.path_data / f"to_DART/temp1_conc_{mem}.nc"
        )

        # Ensure the output directory exists
        # final_concentration_file.parent.mkdir(parents=True, exist_ok=True)
        with open(
            f'logs_orchestrator/farm_to_dart_full_logs/subprocess_out_{mem}_{rounded_timestamp.strftime("%Y%m%d%H")}.log',
            "a",
        ) as log_file:  # Append log file
            logging.info(f"Processing member {mem}")

            # Step 1: Select SP, P, and T from the input meteo file
            subprocess.run(
                ["cdo", "selname,SP,P,T", meteo_file, temp_output_meteo],
                stdout=log_file,
                stderr=log_file,
                check=True,
            )

            # Step 2: Select the timestep from the rounded timestamp
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
            subprocess.run(
                [
                    "/gporq3/minni/FARM-DART/arconv-2.5.10",
                    arconv_input_file,
                    arconv_output_file,
                    "1",
                ],
                stdout=log_file,
                stderr=log_file,
                check=True,
            )

            # Step 5: Use ncks to append SP, P, and T variables to the FARM concentration file
            subprocess.run(
                [
                    "ncks",
                    "-A",
                    "-v",
                    "P,SP,T",
                    temp_output_meteo_plus1,
                    arconv_output_file,
                ],
                stdout=log_file,
                stderr=log_file,
                check=True,
            )

            # Step 6: Copy the result to the final concentration file
            subprocess.run(
                ["cp", arconv_output_file, final_concentration_file], check=True
            )

            # Step 7: Set reference time in the concentration file
            subprocess.run(
                [
                    "cdo",
                    "-setreftime,1900-01-01,00:00:00,days",
                    final_concentration_file,
                    temp_concentration_file,
                ],
                check=True,
            )

            # Step 8: Set calendar to Gregorian
            subprocess.run(
                [
                    "cdo",
                    "-setcalendar,gregorian",
                    temp_concentration_file,
                    temp1_concentration_file,
                ],
                check=True,
            )

            # Step 9: Remove unnecessary attributes from the concentration file
            subprocess.run(
                ["ncatted", "-a", "add_offset,,d,,", temp1_concentration_file],
                check=True,
            )
            subprocess.run(
                ["ncatted", "-a", "scale_factor,,d,,", temp1_concentration_file],
                check=True,
            )
            subprocess.run(
                ["ncatted", "-a", "_FillValue,,d,,", temp1_concentration_file],
                check=True,
            )
            subprocess.run(
                ["ncatted", "-a", "missing_value,,d,,", temp1_concentration_file],
                check=True,
            )

            # Step 10: Copy the cleaned file to the final concentration file location
            subprocess.run(
                ["cp", temp1_concentration_file, final_concentration_file], check=True
            )

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


def prepare_farm_to_dart_nc_par(
    path_manager, timestamp_farm, rounded_timestamp, seconds_model, days_model, no_mems
):
    os.makedirs(path_manager.path_data / "temp", exist_ok=True)
    os.makedirs(path_manager.path_data / "to_DART", exist_ok=True)
    max_workers = 48
    logger.info("Starting the orchestration of FARM to DART NetCDF conversion")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_member,
                mem,
                path_manager,
                timestamp_farm,
                rounded_timestamp,
                seconds_model,
                days_model,
            ): mem
            for mem in range(no_mems)
        }
        for future in as_completed(futures):
            mem = futures[future]
            try:
                future.result()  # Check if there were exceptions
                logger.info(f"Member {mem} processed successfully.")
            except Exception as e:
                logging.error(
                    f"An error occurred during processing of member {mem}: {e}"
                )


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
 