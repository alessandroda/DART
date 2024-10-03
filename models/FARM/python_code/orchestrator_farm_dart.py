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
    prepare_farm_to_dart_nc
)

#filter_job = [path_manager.env_dir, "submit.py"]

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

    def __init__(self, base_path, env_python, listing_file, run_submit_farm_template, path_submit_bsh, path_filter, log_paths=True):
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

        self.check_paths_exist()

    def check_paths_exist(self):
        """Checks if the base, environment, and listing paths exist."""
        paths_to_check = {
            "Base path": self.base_path,
            "Environment directory": self.env_dir,
            "Listing file": self.listing_file,
            'Submit farm template': self.run_submit_farm_template,
            'Farm submission path': self.path_submit_bsh,
            'Submit filter path': self.path_filter
        }

        for path_name, path_value in paths_to_check.items():
            if not path_value.exists():
                raise FileNotFoundError(f"{path_name} does not exist: {path_value}")
            if self.log_paths:
                print(f"{path_name} exists: {path_value}")

    def get_paths(self):
        """Returns all paths as a tuple."""
        return (self.base_path, self.env_dir, self.listing_file, self.run_submit_farm_template, self.path_submit_farm, self.path_submit_filter)






dir_obs_converter =  "DART/observations/obs_converters/S5P_TROPOMI_L3/work/"
obs_converter_command = "convert_s5p_tropomi_l3"
    
path_manager = PathManager(base_path="/gporq3/minni/FARM-DART/",
                            env_python= "miniconda3/bin/python3.10",
                            listing_file="DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03__listing.csv",
                            run_submit_farm_template="RUN/script/templates/run_submit.template.bsh",
                            path_submit_bsh="RUN/script/auto_runs/",
                            path_filter= "DART/models/FARM/work/",
                            log_paths=True)

# TIME SETTINGS
start_time = "2023-03-01 09:00:00"
start_time = pd.to_datetime(start_time)
end_time = "2023-03-01 11:00:00"
end_time = pd.to_datetime(end_time)
dt = 3600
dt = pd.Timedelta(dt, unit="s")

# PATH SETTINGS
listing = pd.read_csv(
    path_manager.listing_file,
    sep=";",
)
listing["start_time"] = pd.to_datetime(listing["start_time"])

# ENSEMBLE SETTINGS
ens_members = list(np.arange(0,20,1))


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
        import pdb; pdb.set_trace()
        string_to_replace_template =  f'{timestamp_farm.strftime("%Y%m%d%H")}.bsh'
        '''
            TODO MEMBERS LOOP
            strings_to_replace_template =[f'{timestamp_farm.strftime("%Y%m%d%H")}_{mem}.bsh' for mem in ens_members]
                the members run differ in the EMISSION and BCS path to be specified in
                the submit.sh. Therefore, submit_$mem_value.sh needs to be created. In
                other words we need more placeholders to palce in the submit.sh
            commands_farm_run = [ 'run_submit' + string_to_replace_template for string_to_replace_template in strings_to_replace_template()]
            path_runs = [path_submit_farm + 'run_submit' + string_to_replace_template for string_to_replace_template in strings_to_replace_template()]
                For each string_to_replace in strings_to_replace we neeed to replace
                two templates, one for submit.sh and another for run_submit.sh
                because submit.sh will contain the members placeholders and
                run_submit.sh will use submit_mem_number.sh to run the specific farm
                 member
        '''

        command_farm_run = 'run_submit' + string_to_replace_template 
        path_run = path_manager.path_submit_bsh / f'run_submit{string_to_replace_template}'

        replace_nml_template(
             input_nml_path=path_manager.run_submit_farm_template,
             entries_tbr_dict={
                 "da_date_start": timestamp_farm.strftime("%Y%m%d%H"),
                 "da_date_end": timestamp_farm.strftime("%Y%m%d%H"),
             },
             output_nml_path= path_run 
         )
        
        
        #job_id = run_command_in_directory_bsub(command_farm_run, path_manager.path_submit_bsh)
        # Check if the job_id is valid
        #if job_id:
        #    print(f"Submitted job with ID: {job_id}")
        #    # Wait for the job to complete
        #    while not check_job_status_cresco(job_id):
        #        print(f"Waiting for job {job_id} to complete...")
        #        time.sleep(30)  # Sleep for 30 seconds before checking again
        #    print(f"Job {job_id} is complete. Continuing...")
        #else:
        #     print("No job ID returned, something went wrong with submission.")
        #     break 
        # NO NEED ANYMORE FOR THIS CHECK
        simulated_time = timestamp_farm
 
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
        path_manager.base_path / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input_template.nml",
        entries_tbr_dict={
            "$file_path_s5p":
            path_manager.base_path / f'DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/{orbit_filename["filename"].values[0]}',
            "$file_out":
            path_manager.base_path / f'DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/{obs_seq_name}',
        },
        output_nml_path=path_manager.base_path / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input.nml",
    )
    import pdb; pdb.set_trace()

    run_command_in_directory(obs_converter_command, path_manager.base_path / dir_obs_converter)
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
        path_manager.base_path / f'RUN/data//posteriors/{formatted_t_str}'
    )
    Path(output_sim_folder).mkdir(parents=True, exist_ok=True)
    analysis_sim_folder = (
        path_manager.base_path / f'RUN/data/analysis/{formatted_t_str}'
    )
    Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
    preassim_sim_folder = (
        path_manager.base_path / f'RUN/data/preassim/{formatted_t_str}'
    )
    Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)

    # FILTER INPUT NML
    replace_nml_template(
        path_manager.base_path / "DART/models/FARM/work/input_template.nml",
        entries_tbr_dict={
            "$obs_sequence_name": obs_seq_name,
            "$folder_path": output_sim_folder,
            "$date_assim": formatted_t_str,
            "$template_farm": path_manager.base_path / f"RUN/data/to_DART/ic_g1_{seconds_model}_{days_model}_00.nc",
            "$init_time_days": str(days_model),
            "$init_time_seconds": str(seconds_model),
            "$first_obs_days": str(days),
            "$first_obs_seconds": str(seconds),
        },
        output_nml_path= path_manager.base_path / "DART/models/FARM/work/input.nml",
    )

   
    # FILTER_INPUT_LIST.TXT
    replace_nml_template(
        path_manager.base_path / "DART/models/FARM/work/filter_input_list_template.txt",
        entries_tbr_dict={
            "$folder_path": path_manager.base_path / f"RUN/data/to_DART/",
            "$days": str(seconds_model),
            "$seconds": str(days_model),
        },
        output_nml_path=path_manager.base_path / "DART/models/FARM/work/filter_input_list.txt",
    )

    # FILTER_OUTPUT_LIST.TXT
    replace_nml_template(
        path_manager.base_path / "DART/models/FARM/work/filter_output_list_template.txt",
        entries_tbr_dict={
            "$folder_path": output_sim_folder,
            "$date": formatted_t_str
        },
        output_nml_path= path_manager.base_path / "DART/models/FARM/work/filter_output_list.txt",
    )

    prepare_farm_to_dart_nc(path_manager, timestamp_farm, rounded_timestamp, seconds_model,days_model)

    with open(
        path_manager.base_path /"DART/models/FARM/python_code/orchestrator_info.yaml",
        "r",
    ) as yaml_file:
        settings = yaml.safe_load(yaml_file)

    for ens_member in settings["ens_members"]:
        break
    breakpoint()
    replace_nml_template(
        path_manager.base_path / "RUN/script/templates/submit_filter.template.bsh",
        entries_tbr_dict={
            "CURRENT_DATE": rounded_timestamp.strftime("%Y%m%d%H"),
            "NO_PROC": str(20)
        },
        output_nml_path=path_manager.path_submit_bsh / "submit_filter.bsh",
    )
    replace_nml_template(
        path_manager.base_path / "RUN/script/templates/run_filter.template.bsh",
        entries_tbr_dict={
            "NO_PROC": str(20)
        },
        output_nml_path=path_manager.path_submit_bsh / "run_filter.bsh",
    )
    # RUN THE FILTER
    job_id = run_command_in_directory_bsub("./submit_filter.bsh", path_manager.path_submit_bsh)
    import pdb; pdb.set_trace()
    
    while True:
        if not check_job_status_cresco(job_id):
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

    print("------------------")


