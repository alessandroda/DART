import shutil
import time
import pandas as pd
import os
import subprocess
import math
from pathlib import Path
import numpy as np
import yaml
from datetime import datetime, timedelta
import logging

# import pdb; pdb.set_trace()

from orchestrator_utils import (
    check_job_status_cresco,
    modify_nc_file,
    open_dataset,
    replace_nml_template,
    round_to_closest_hour,
    run_command_in_directory,
    run_command_in_directory_bsub,
    searchFile,
    set_date_gregorian,
    submit_slurm_job,
    prepare_farm_to_dart_nc,
    PathManager,
    CleanupContext,
    submit_and_wait,
    prepare_dart_to_farm_nc,
)


logging.basicConfig(
    filename="farm_to_dart.log", format="%(asctime)s %(message)s", level=logging.INFO
)
# logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

logger = logging.getLogger(__name__)

# filter_job = [path_manager.env_dir, "submit.py"]

dir_obs_converter = "DART/observations/obs_converters/S5P_TROPOMI_L3/work/"
obs_converter_command = "convert_s5p_tropomi_l3"

path_manager = PathManager(
    base_path="/gporq3/minni/FARM-DART/",
    env_python="miniconda3/bin/python3.10",
    listing_file="DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03__listing.csv",
    run_submit_farm_template="RUN/script/templates/new_run_submit.template.bsh",
    path_submit_bsh="RUN/script/auto_runs/",
    path_filter="DART/models/FARM/work/",
    path_data="RUN/data/",
    log_paths=True,
)

# TIME SETTINGS
start_time = "2023-03-01 09:00:00"
start_time = pd.to_datetime(start_time)
end_time = "2023-03-01 11:00:00"
end_time = pd.to_datetime(end_time)
dt = 180
dt = pd.Timedelta(dt, unit="s")

# PATH SETTINGS
listing = pd.read_csv(
    path_manager.listing_file,
    sep=";",
)
listing["start_time"] = pd.to_datetime(listing["start_time"])

# ENSEMBLE SETTINGS
ens_members = list(np.arange(0, 20, 1))

logging.info("TIME LOOP BEGINS")
t = start_time
simulated_time = None
while t <= end_time:
    # if not path_manager.check_existing_farm:
    #    run_farm()
    """
    Run FARM in Hourly Mode:
    The timestamp_farm provided for the simulation will produce a file with a timestamp that is one hour ahead of the timestamp_farm.
    Replace timestamp_farm in the template bash file with the timestamp you are using.
    2023030109 -> ore 2023030110 nel file
    timestamp_farm -> the one on the file and use to submit the simulation
    rounded_timestamp -> the one appearing in the farm output file and named with timestamp_farm (+1)
    """

    logger.info(f"Time loop iter: {t}")
    timestamp_farm = round_to_closest_hour(t)
    logger.info(f"1.----------Running FARM for hour {t}")
    if timestamp_farm != simulated_time:
        string_to_replace_template = f'{timestamp_farm.strftime("%Y%m%d%H")}.bsh'
        """
            The run_submit_farm_template should be still the same but to be
            replaced with the members value and use replace_nml_template to set
            the right EMISSION, BC, and core which should end with two digits
            xx.nc that are the memeber's value
        """
        # strings_to_replace_template =[f'{timestamp_farm.strftime("%Y%m%d%H")}_{mem}.bsh' for mem in ens_members]
        # commands_farm_run = [ 'run_submit' + string_to_replace_template for string_to_replace_template in strings_to_replace_template()]
        # path_runs = [path_submit_farm + 'run_submit' + string_to_replace_template for string_to_replace_template in strings_to_replace_template()]
        # for path_run in path_runs:
        #    replace_nml_template(
        #        input_nml_path=path_manager.run_submit_farm_template,
        #        entries_tbr_dict={
        #         "da_date_start": timestamp_farm.strftime("%Y%m%d%H"),
        #         "da_date_end": timestamp_farm.strftime("%Y%m%d%H"),
        #     },
        #     output_nml_path= path_run
        #    )
        # commands_with_directories = [(command_farm_run, path_manager.path_submit.bsh) for command_farm_run in commands_farm_run]
        # submit_and_wait(commands_with_directories)

        command_farm_run = "run_submit" + string_to_replace_template
        path_run = (
            path_manager.path_submit_bsh / f"run_submit{string_to_replace_template}"
        )

        replace_nml_template(
            input_nml_path=path_manager.run_submit_farm_template,
            entries_tbr_dict={
                "da_date_start": timestamp_farm.strftime("%Y%m%d%H"),
                "da_date_end": timestamp_farm.strftime("%Y%m%d%H"),
            },
            output_nml_path=path_run,
        )

        commands_with_directories = [(command_farm_run, path_manager.path_submit_bsh)]

        # submit_and_wait(commands_with_directories)

        simulated_time = timestamp_farm
    # if not path_manager.check_exisiting_obs_sat:
    #    run_obs_converter
    logger.info(f"2. Increment time and search for orbit")
    t += dt
    formatted_t_str = t.strftime("%Y%m%d_%H%M%S")
    logger.info(formatted_t_str)
    rounded_timestamp = round_to_closest_hour(t)
    time_list = [rounded_timestamp, rounded_timestamp + timedelta(hours=-1)]
    orbit_filename = searchFile(t - 0.5 * dt, t + 0.5 * dt, listing)
    if orbit_filename.empty:
        continue
    logging.info(
        f"Orbit file found: {orbit_filename['filename'].values[0]}  corresponding to time: {orbit_filename['start_time']}"
    )
    t_sat_obs = pd.to_datetime(orbit_filename["start_time"].values[0])
    seconds, days = set_date_gregorian(
        t_sat_obs.year,
        t_sat_obs.month,
        t_sat_obs.day,
        t_sat_obs.hour,
        t_sat_obs.minute,
        t_sat_obs.second,
    )
    obs_seq_name = f"obs_seq_{seconds}_{days}.out"
    # OBS CONVERTER INPUT.NML
    replace_nml_template(
        path_manager.base_path
        / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input_template.nml",
        entries_tbr_dict={
            "$file_path_s5p": path_manager.base_path
            / f'DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/{orbit_filename["filename"].values[0]}',
            "$file_out": path_manager.base_path
            / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/{obs_seq_name}",
        },
        output_nml_path=path_manager.base_path
        / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input.nml",
    )

    # run_command_in_directory(obs_converter_command, path_manager.base_path / dir_obs_converter)
    # preprocessing model farm nc move before converter
    time_model = rounded_timestamp
    seconds_model, days_model = set_date_gregorian(
        time_model.year,
        time_model.month,
        time_model.day,
        time_model.hour,
        time_model.minute,
        time_model.second,
    )
    output_sim_folder = path_manager.path_data / f"/posteriors/{formatted_t_str}"
    Path(output_sim_folder).mkdir(parents=True, exist_ok=True)
    analysis_sim_folder = path_manager.path_data / f"/analysis/{formatted_t_str}"
    Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
    preassim_sim_folder = path_manager.path_data / f"/preassim/{formatted_t_str}"
    Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)

    # RUN FILTER: FILTER INPUT NML
    replace_nml_template(
        path_manager.base_path / "DART/models/FARM/work/input_template.nml",
        entries_tbr_dict={
            "$obs_sequence_name": obs_seq_name,
            "$folder_path": output_sim_folder,
            "$date_assim": formatted_t_str,
            "$template_farm": path_manager.base_path
            / f"RUN/data/to_DART/ic_g1_{seconds_model}_{days_model}_00.nc",
            "$init_time_days": str(days_model),
            "$init_time_seconds": str(seconds_model),
            "$first_obs_days": str(days),
            "$first_obs_seconds": str(seconds),
        },
        output_nml_path=path_manager.base_path / "DART/models/FARM/work/input.nml",
    )

    # FILTER_INPUT_LIST.TXT
    replace_nml_template(
        path_manager.base_path / "DART/models/FARM/work/filter_input_list_template.txt",
        entries_tbr_dict={
            "$fme -d x,lon -d y,lat tmp1.ncncrename -d x,lon -d y,lat tmp1.ncolder_path": path_manager.base_path
            / f"RUN/data/to_DART/",
            "$days": str(seconds_model),
            "$seconds": str(days_model),
        },
        output_nml_path=path_manager.base_path
        / "DART/models/FARM/work/filter_input_list.txt",
    )

    # FILTER_OUTPUT_LIST.TXT
    replace_nml_template(
        path_manager.base_path
        / "DART/models/FARM/work/filter_output_list_template.txt",
        entries_tbr_dict={"$folder_path": output_sim_folder, "$date": formatted_t_str},
        output_nml_path=path_manager.base_path
        / "DART/models/FARM/work/filter_output_list.txt",
    )
    # with CleanupContext(path_manager.base_path / 'RUN/data/temp/'):
    # prepare_farm_to_dart_nc(path_manager, timestamp_farm, rounded_timestamp, seconds_model,days_model)

    with open(
        path_manager.base_path / "DART/models/FARM/python_code/orchestrator_info.yaml",
        "r",
    ) as yaml_file:
        settings = yaml.safe_load(yaml_file)

    for ens_member in settings["ens_members"]:
        break
    replace_nml_template(
        path_manager.base_path / "RUN/script/templates/submit_filter.template.bsh",
        entries_tbr_dict={
            "CURRENT_DATE": rounded_timestamp.strftime("%Y%m%d%H"),
            "CORES": str(5),
        },
        output_nml_path=path_manager.path_submit_bsh / "submit_filter.bsh",
    )
    replace_nml_template(
        path_manager.base_path / "RUN/script/templates/run_filter.template.bsh",
        entries_tbr_dict={"CORES": str(5)},
        output_nml_path=path_manager.path_submit_bsh / "run_filter.bsh",
    )
    """
    # RUN THE FILTER
    job_id = run_command_in_directory_bsub("./submit_filter.bsh",path_manager.path_submit_bsh, farm = False)
    # A.D'A 04/10/2024 should pass to orchestrator_utils
    job_id = job_id.strip()[1:-1] 
    while True:
        if check_job_status_cresco(job_id):
            print("Job completed successfully.")
            # Move files if job completed successfully
            for filename in os.listdir(f"{path_manager.path_filter}"):
                if filename.startswith("analysis_"):
                    try:
                        shutil.move(
                            os.path.join(
                                path_manager.path_filter,
                                filename,
                            ),
                            analysis_sim_folder,
                        )
                    except shutil.Error:
                        print(
                            f"Failed to move '{filename}' to '{analysis_sim_folder}' because it already exists."
                        )
                elif filename.startswith("preassim_"):
                    try:
                        shutil.move(
                            os.path.join(
                                path_manager.path_filter,
                                filename,
                            ),
                            os.path.join(preassim_sim_folder, filename),
                        )
                    except shutil.Error:
                        print(
                            f"Failed to move '{filename}' to '{preassim_sim_folder}' because it already exists."
                        )

            break
        else:
            print("Job is still running. Waiting...")
            time.sleep(10)
    """
    # Handle the posterior
    # create a method that picks the SO2 from the posterior files and place the
    # SO2 in the core of FARM in the next run.
    breakpoint()
    prepare_dart_to_farm_nc(
        path_manager, output_sim_folder, formatted_t_str, time_model, ass_var="c_SO2"
    )
    print("------------------")
