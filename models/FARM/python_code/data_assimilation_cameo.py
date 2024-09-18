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
    check_job_status_cresco,
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

path_cresco = "/gporq3/minni/FARM-DART/"
env_cresco = "/gporq3/minni/FARM-DART/miniconda3/bin/python3.10"
listing_cresco = "/gporq3/minni/FARM-DART/DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03__listing.csv"
run_submit_path_cresco = "/gporq3/minni/FARM-DART/RUN/script/run_submit.template.sh"
path_submit_farm = "/gporq3/minni/FARM-DART/RUN/script/"
path_submit_filter = "/gporq3/minni/FARM-DART/RUN/script"

path_mumbai = "/mnt/mumbai_n4r5/dausilio/projects/"
env_mumbai = "/home/dausilio/miniconda3/envs/dartenv/bin/python"
listing_mumbai = "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/projects/CAMEO/listing-SO2-test_CAMEO.csv"


def get_base_path():
    cwd = os.getcwd()

    # Example: Check if cwd contains 'gporq3' to determine if it's Cresco or Mumbai
    if "gporq3" in cwd:
        print("Detected Cresco environment")
        return path_cresco, env_cresco, listing_cresco
    elif "mumbai" in cwd:
        print("Detected Mumbai environment")
        return path_mumbai, env_mumbai, listing_mumbai
    else:  # Add more conditions for other environments if needed
        raise ValueError("Unknown environment. Please set the base path manually.")


# Use the smartly selected path
base_path, env_path, listing_path = get_base_path()


# TIME SETTINGS
start_time = "2023-03-01 09:00:00"
start_time = pd.to_datetime(start_time)
end_time = "2023-03-01 11:00:00"
end_time = pd.to_datetime(end_time)
dt = 180
dt = pd.Timedelta(dt, unit="s")

# PATH SETTINGS
dir_obs_converter = base_path + "DART/observations/obs_converters/S5P_TROPOMI_L3/work/"
obs_converter_command = "convert_s5p_tropomi_l3"
dir_work_farm = base_path + "/gporq3/minni/FARM-DART/DART/models/FARM/work/"
filter_job = [env_path, "submit.py"]
listing = pd.read_csv(
    listing_path,
    sep=";",
)
listing["start_time"] = pd.to_datetime(listing["start_time"])

print("------------------")
print("TIME LOOP BEGINS")
print("------------------")
t = start_time
simulated_time = None
while t <= end_time:
    '''
        Run FARM in Hourly Mode:
        The timestamp_farm provided for the simulation will produce a file with a timestamp that is one hour ahead of the timestamp_farm.
        Replace timestamp_farm in the template bash file with the timestamp you are using.
        2023030109 -> ore 2023030110 nel file
        timestamp_farm -> the one on the file and use to submit the simulation
        rounded_timestamp -> the one appearing in the farm output file and named with timestamp_farm (+1)
    '''
    
    print(f"time: {t}")
    timestamp_farm = round_to_closest_hour(t)
    if timestamp_farm != simulated_time:
        # path_run = run_submit_path_cresco.replace("template.sh", "test.sh")
        # replace_nml_template(
        #     input_nml_path=run_submit_path_cresco,
        #     entries_tbr_dict={
        #         "da_date_start": timestamp_farm.strftime("%Y%m%d%H"),
        #         "da_date_end": timestamp_farm.strftime("%Y%m%d%H"),
        #     },
        #     output_nml_path=path_run,
        # )

        # stdout, stderr = run_command_in_directory("sh run_submit.bsh", path_run)
        # job_id = stdout.strip()
        # # Check if the job_id is valid
        # if job_id:
        #     print(f"Submitted job with ID: {job_id}")

        #     # Wait for the job to complete
        #     while not check_job_status_cresco(job_id):
        #         print(f"Waiting for job {job_id} to complete...")
        #         time.sleep(30)  # Sleep for 30 seconds before checking again

        #     print(f"Job {job_id} is complete. Continuing...")

            simulated_time = timestamp_farm + timedelta(hours=1)
        # else:
        #     print("No job ID returned, something went wrong with submission.")

    # increment dt to search for the hours to assimilate
    t += dt
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
        base_path
        + "observations/obs_converters/S5P_TROPOMI_L3/work/input_template.nml",
        entries_tbr_dict={
            "$file_path_s5p": orbit_filename["filename"].values[0],
            "$file_out": obs_seq_name,
        },
        output_nml_path=base_path
        + "observations/obs_converters/S5P_TROPOMI_L3/work/input.nml",
    )

    # run_command_in_directory(obs_converter_command, dir_obs_converter)
    time_model = rounded_timestamp
    seconds_model, days_model = set_date_gregorian(
        time_model.year,
        time_model.month,
        time_model.day,
        time_model.hour,
        time_model.minute,
        time_model.second,
    )
    output_sim_folder = (
        f"{base_path}/models/FARM/projects/CAMEO/posteriors/{formatted_t_str}"
    )
    Path(output_sim_folder).mkdir(parents=True, exist_ok=True)
    analysis_sim_folder = (
        f"{base_path}/DART/models/FARM/projects/CAMEO/analysis/{formatted_t_str}"
    )
    Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
    preassim_sim_folder = (
        f"{base_path}/models/FARM/projects/CAMEO/preassim/{formatted_t_str}"
    )
    Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)

    # FILTER INPUT NML
    replace_nml_template(
        f"{base_path}/models/FARM/work/input_template.nml",
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
    # Passaggi da eseguire per portare il file di FARM con P, SP, T necessarie per DART in realta' isolerei anche C_SO2 per peso
    # cdo selname,SP,P,T /gporq3/minni/CAMEO/RUN/data/INPUT/METEO/ifsecmc dest_path
    # cdo seltimestep,10 met.nc met_1hr.nc
    # cdo shifttime,1hour met_1hr.nc met_1hr+1.nc 
    # arconv input output 1 (conc di FARM)
    # ncks -A -v P,SP,T ../met_1hr+1.nc conc_g1_
    # 

    # FILTER_INPUT_LIST.TXT
    replace_nml_template(
        f"{base_path}/models/FARM/work/filter_input_list_template.txt",
        entries_tbr_dict={
            "$days": str(seconds_model),
            "$seconds": str(days_model),
        },
        output_nml_path="{base_path}/models/FARM/work/filter_input_list.txt",
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
        "/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/python_code/farm_to_dart_cameo.yaml",
        "r",
    ) as yaml_file:
        settings = yaml.safe_load(yaml_file)

    for ens_member in settings["ens_members"]:
        # replace {value_member} in settings["input_farm_file"] with ens_member
        input_farm_file = (
            settings["input_farm_file"]
            .format(date=timestamp_farm.strftime("%Y%m%d%H"))
            # .format(value_member=ens_member)
        )
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

    replace_nml_template(
        f"{path_submit_filter}/template/run_submit.template.bsh",
        entries_tbr_dict={
            "CURRENT_DATE": rounded_timestamp.strftime("%Y%m%d%H"),
            "NO_PROC": 20
        },
        output_nml_path=f"{path_submit_filter}/run_submit.bsh",
    )
    replace_nml_template(
        f"{path_submit_filter}/template/run_filter.template.bsh",
        entries_tbr_dict={
            "NO_PROC": 20
        },
        output_nml_path=f"{path_submit_filter}/run_filter.bsh",
    )
    # RUN THE FILTER
    stdout, stderr = run_command_in_directory("bash run_submit.bsh", path_submit_filter)
    job_id = stdout.strip()
    # job_id = submit_slurm_job(f"{ens_member}", "filter")

    while True:
        if not check_job_status_cresco(job_id):
            print("Job completed successfully.")
            # Move files if job completed successfully
            for filename in os.listdir(f"{base_path}/models/FARM/work"):
                if filename.startswith("analysis_"):
                    try:
                        shutil.move(
                            os.path.join(
                                f"{base_path}/models/FARM/work",
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
                                f"{base_path}/models/FARM/work",
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
