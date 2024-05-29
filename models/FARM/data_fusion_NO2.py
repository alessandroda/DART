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
from FARM_to_DART import (
    check_job_status,
    convert_time_to_hours,
    modify_nc_file,
    open_dataset,
    replace_nml_template,
    round_to_closest_hour,
    run_command_in_directory,
    searchFile,
    set_date_gregorian,
    submit_slurm_job,
)

# SETTINGS
cwd = os.getcwd()
start_time = "2021-01-05 00:00:00"
start_time = pd.to_datetime(start_time)
end_time = "2021-01-31 23:00:00"
end_time = pd.to_datetime(end_time)
dt = 180
dt = pd.Timedelta(dt, unit="s")
dir_obs_converter = "/mnt/mumbai_n4r5/dausilio/projects/DART/observations/obs_converters/S5P_TROPOMI_L3/work/"
obs_converter_command = "convert_s5p_tropomi_l3"
dir_work_farm = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/"
filter_job = ["/home/dausilio/miniconda3/envs/dartenv/bin/python", "submit.py"]
listing = pd.read_csv(
    "/mnt/mumbai_n4r5/dausilio/SATDA_OFFLINE/OBS/listing-NO2-QA.csv", sep=";"
)

listing["start_time"] = pd.to_datetime(listing["start_time"])
print("------------------")
print("TIME LOOP BEGINS")
print("------------------")
t = start_time
while t <= end_time:
    print(f"time: {t}")
    t = pd.to_datetime(t) + dt
    formatted_t_str = str(t).replace(" ", "_").replace(":", "").replace("-", "")
    rounded_timestamp = round_to_closest_hour(t)
    time_list = [rounded_timestamp, rounded_timestamp + timedelta(hours=-1)]
    orbit_filename = searchFile(t - 0.5 * dt, t + 0.5 * dt, listing)
    if orbit_filename.empty:
        continue
    print(f"Orbit file found: {orbit_filename['filename'].values[0]}")
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
        "/mnt/mumbai_n4r5/dausilio/projects/DART/observations/obs_converters/S5P_TROPOMI_L3/work/input_template.nml",
        entries_tbr_dict={
            "$file_path_s5p": orbit_filename["filename"].values[0],
            "$file_out": obs_seq_name,
        },
        output_nml_path="/mnt/mumbai_n4r5/dausilio/projects/DART/observations/obs_converters/S5P_TROPOMI_L3/work/input.nml",
    )

    run_command_in_directory(obs_converter_command, dir_obs_converter)
    time_model = rounded_timestamp
    seconds_model, days_model = set_date_gregorian(
        time_model.year,
        time_model.month,
        time_model.day,
        time_model.hour,
        time_model.minute,
        time_model.second,
    )
    output_sim_folder = f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/posteriors/{formatted_t_str}"
    Path(output_sim_folder).mkdir(parents=True, exist_ok=True)
    analysis_sim_folder = f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/analysis/{formatted_t_str}"
    Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
    preassim_sim_folder = f"/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/IRIDE_ens/preassim/{formatted_t_str}"
    Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)

    # FILTER INPUT NML
    replace_nml_template(
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/input_template.nml",
        entries_tbr_dict={
            "$obs_sequence_name": obs_seq_name,
            "$folder_path": output_sim_folder,
            "$date_assim": formatted_t_str,
            "$template_farm": f"model_{seconds_model}_{days_model}_00.nc",
            "$init_time_days": str(days_model),
            "$init_time_seconds": str(seconds_model),
            "$first_obs_days": str(days),
            "$first_obs_seconds": str(seconds),
        },
        output_nml_path="/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/input.nml",
    )
    # FILTER_INPUT_LIST.TXT
    replace_nml_template(
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/filter_input_list_template.txt",
        entries_tbr_dict={
            "$days": str(seconds_model),
            "$seconds": str(days_model),
        },
        output_nml_path="/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/filter_input_list.txt",
    )

    # FILTER_OUTPUT_LIST.TXT
    replace_nml_template(
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/filter_output_list_template.txt",
        entries_tbr_dict={
            "$folder_path": output_sim_folder,
        },
        output_nml_path="/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/filter_output_list.txt",
    )

    with open(
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/farm_to_dart.yaml", "r"
    ) as yaml_file:
        settings = yaml.safe_load(yaml_file)

    for ens_member in settings["ens_members"]:
        # replace {value_member} in settings["input_farm_file"] with ens_member
        input_farm_file = settings["input_farm_file"].format(value_member=ens_member)
        # input_meteo_file = settings["input_meteo_file"]
        output_folder = settings["output_folder"]
        pol = settings["pol"]

        # Convert FARM file to DART readable format
        # ds = add_meteo_to_farm(input_farm_file, input_meteo_file)
        ds = open_dataset(input_farm_file)
        # Modify FARM file
        ds = modify_nc_file(ds, pol, time_list=[rounded_timestamp])
        if ds is None:
            continue

        # Convert time to hours and days since a specific date
        convert_time_to_hours(
            ds, output_folder, seconds_model, days_model, ens_member=str(ens_member)
        )

    job_id = submit_slurm_job(f"{ens_member}", "filter")

    while True:
        if not check_job_status(job_id):
            print("Job completed successfully.")
            # Move files if job completed successfully
            for filename in os.listdir(
                "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work"
            ):
                if filename.startswith("analysis_"):
                    try:
                        shutil.move(
                            os.path.join(
                                "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work",
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
                                "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work",
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

    print("------------------")
