import shutil
import time
import pandas as pd
import os
from pathlib import Path
from datetime import timedelta
import logging
import yaml
from orchestrator_utils import (
    check_job_status_cresco,
    replace_nml_template,
    run_command_in_directory,
    run_command_in_directory_bsub,
    searchFile,
    set_date_gregorian,
    prepare_farm_to_dart_nc,
    prepare_farm_to_dart_nc_par,
    PathManager,
    TimeManager,
    submit_and_wait,
    prepare_dart_to_farm_nc,
    modify_yaml_date
)
 
logging.basicConfig(filename=f'logs_orchestrator/farm_to_dart_{time.strftime("%Y%m%d_%H%M%S")}.log', format="%(asctime)s %(message)s", level=logging.INFO)

logger = logging.getLogger(__name__)


CONFIG_PATH = "config_orchestrator.yaml"


def load_config(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)



class FarmToDartPipeline:
    def __init__(self):
        # Load configuration from YAML file
        self.config = load_config(CONFIG_PATH)
        # Setup path manager and time manager with loaded configurations
        self.path_manager = PathManager(
            base_path=self.config['paths']['base_path'],
            env_python=self.config['paths']['env_python'],
            listing_file=self.config['paths']['listing_file'],
            run_submit_farm_template=self.config['paths']['run_submit_farm_template'],
            path_submit_bsh=self.config['paths']['path_submit_bsh'],
            path_filter=self.config['paths']['path_filter'],
            path_data=self.config['paths']['path_data'],
            log_paths=True,
        )

        self.time_manager = TimeManager(
            start_time=self.config['time']['start_time'],
            end_time=self.config['time']['end_time'],
            dt_seconds=self.config['time']['dt_seconds'],
        )

        self.listing = pd.read_csv(self.path_manager.listing_file, sep=";")
        self.listing['start_time'] = pd.to_datetime(self.listing['start_time'])
        self.days_obs = 0
        self.seconds_obs = 0
        self.days_model = 0
        self.seconds_model = 0
        self.output_sim_folder = None
        self.ass_var = self.config['assimilation']['ass_var']
        self.no_mems = self.config['assimilation']['no_mems']
        

    def run_farm(self):
        logger.info(
            f"1.----------Running FARM for hour {self.time_manager.current_time}"
        )
        timestamp_farm = TimeManager.round_to_closest_hour(
            self.time_manager.current_time
        ).strftime("%Y%m%d%H")
        string_to_replace_template = f'{timestamp_farm}.bsh'

        command_farm_run = "run_submit_ens" + string_to_replace_template
        path_run = (
            self.path_manager.path_submit_bsh
            / f"run_submit_ens{string_to_replace_template}"
        )
        list_mems = [str(mem) for mem in range(self.no_mems)]
        replace_nml_template(
            input_nml_path=self.path_manager.run_submit_farm_template,
            entries_tbr_dict={
                "da_date_start": timestamp_farm,
                "da_date_end": timestamp_farm,
                "@no_mems_list": str(tuple(list_mems)).replace(',','')
            },
            output_nml_path=path_run,
        )
        commands_with_directories = [
            (command_farm_run, self.path_manager.path_submit_bsh)
        ]
        submit_and_wait(self.path_manager, commands_with_directories, timestamp_farm, self.no_mems)
        
    def process_satellite_data(self):
        logger.info(f"2. Increment time and search for orbit")
        self.time_manager.increment_time()

        orbit_filename = searchFile(
            self.time_manager.current_time - 0.5 * self.time_manager.dt,
            self.time_manager.current_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return False

        logger.info(f"Orbit file found: {orbit_filename['filename'].values[0]}")
        self.time_manager.sat_obs = pd.to_datetime(
            orbit_filename["start_time"].values[0]
        )
        return orbit_filename['filename'].values[0]

    def run_obs_converter(self, orbit_filename):
        self.seconds_obs, self.days_obs = set_date_gregorian(
            self.time_manager.sat_obs.year,
            self.time_manager.sat_obs.month,
            self.time_manager.sat_obs.day,
            self.time_manager.sat_obs.hour,
            self.time_manager.sat_obs.minute,
            self.time_manager.sat_obs.second,
        )
        obs_seq_name = f"obs_seq_{self.seconds_obs}_{self.days_obs}.out"

        replace_nml_template(
            self.path_manager.base_path
            / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input_template.nml",
            entries_tbr_dict={
                "$file_path_s5p": self.path_manager.base_path
                / f'DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/{orbit_filename}',
                "$file_out": self.path_manager.base_path
                / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/{obs_seq_name}",
            },
            output_nml_path=self.path_manager.base_path
            / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/input.nml",
        )
        try:
            run_command_in_directory(
                "convert_s5p_tropomi_l3",
                self.path_manager.base_path
                / "DART/observations/obs_converters/S5P_TROPOMI_L3/work/",
            )
        except Exception as e:
            logger.error(f"Error running obs converter: {e}")
            return False
        return obs_seq_name


    def set_days_seconds_model(self):
            self.seconds_model, self.days_model = set_date_gregorian(
            self.time_manager.simulated_time.year,
            self.time_manager.simulated_time.month,
            self.time_manager.simulated_time.day,
            self.time_manager.simulated_time.hour,
            self.time_manager.simulated_time.minute,
            self.time_manager.simulated_time.second,
        )


    def run_dart(self, obs_seq_name):
        logger.info("Running DART")
        
        self.output_sim_folder = (
            self.path_manager.path_data
            / f"posteriors/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(self.output_sim_folder).mkdir(parents=True, exist_ok=True)

        replace_nml_template(
            self.path_manager.base_path / "DART/models/FARM/work/input_template.nml",
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.output_sim_folder,
                "$folder_obs_path": self.path_manager.base_path
                / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/",
                "$date_assim": self.time_manager.current_time.strftime("%Y%m%d_%H%M%S"),
                "$template_farm": self.path_manager.base_path
                / f"RUN/data/to_DART/ic_g1_{self.seconds_model}_{self.days_model}_0.nc",
                "$init_time_days": str(self.days_model),
                "$init_time_seconds": str(self.seconds_model),
                "$first_obs_days": str(self.days_obs),
                "$first_obs_seconds": str(self.seconds_obs),
                "$no_mems": str(self.no_mems)
            },
            output_nml_path=self.path_manager.base_path
            / "DART/models/FARM/work/input.nml",
        )
        # FILTER_INPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path / "DART/models/FARM/work/filter_input_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.path_manager.base_path / f"RUN/data/to_DART/",
                "$days": str(self.seconds_model),
                "$seconds": str(self.days_model),
            },
            output_nml_path=self.path_manager.base_path / "DART/models/FARM/work/filter_input_list.txt",
        )

        # FILTER_OUTPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path
            / "DART/models/FARM/work/filter_output_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.output_sim_folder,
                "$date": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
            },
            output_nml_path=self.path_manager.base_path / "DART/models/FARM/work/filter_output_list.txt",
        )

        replace_nml_template(
            self.path_manager.base_path
            / "RUN/script/templates/submit_filter.template.bsh",
            entries_tbr_dict={
                "CURRENT_DATE": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                "CORES": str(5),
            },
            output_nml_path=self.path_manager.path_submit_bsh / "submit_filter.bsh",
        )

        job_id = run_command_in_directory_bsub(
            "./submit_filter.bsh", self.path_manager.path_submit_bsh, farm =
            False
        )
        time.sleep(10)
        self.monitor_job(job_id)

    def monitor_job(self, job_id):
        logger.info(f"Monitoring job {job_id}")
        job_id = job_id.strip()[1:-1]

        while True:
            if check_job_status_cresco(job_id):
                print("Job completed successfully.")
                # Handle successful job completion: move files
                self.move_analysis_files()
                break
            else:
                print("Job is still running. Waiting...")
                time.sleep(10)

    def move_analysis_files(self):
        analysis_sim_folder = (
            self.path_manager.path_data
            / f"analysis/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(analysis_sim_folder).mkdir(parents =True,exist_ok=True)
        preassim_sim_folder = (
            self.path_manager.path_data
            / f"preassim/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)
        for filename in os.listdir(f"{self.path_manager.path_filter}"):
                if filename.startswith("analysis_"):
                    try:
                        shutil.move(
                            os.path.join(
                                self.path_manager.path_filter,
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
                                self.path_manager.path_filter,
                                filename,
                            ),
                            os.path.join(preassim_sim_folder, filename),
                        )
                    except shutil.Error:
                        print(
                            f"Failed to move '{filename}' to '{preassim_sim_folder}' because it already exists."
                        )


    def run_pipeline(self):
        logger.info("TIME LOOP BEGINS")
        while self.time_manager.current_time <= self.time_manager.end_time:
            self.run_farm()  # Run FARM executable
            self.time_manager.simulated_time = self.time_manager.current_time + timedelta(hours=1)
            modify_yaml_date(CONFIG_PATH, self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00")) 
            self.set_days_seconds_model() 
            orbit_filename = self.process_satellite_data()  # Process satellite data
            if orbit_filename:  # Only proceed if satellite data is found
                obs_seq_name = self.run_obs_converter(orbit_filename)
                
                prepare_farm_to_dart_nc(
                    self.path_manager,
                    self.time_manager.simulated_time,
                    self.time_manager.simulated_time,
                    self.seconds_model,
                    self.days_model,
                    self.no_mems
                )
                self.run_dart(obs_seq_name)  # Run DART assimilation
                prepare_dart_to_farm_nc(
                    self.path_manager,
                    self.output_sim_folder,
                    self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                    self.ass_var,
                    self.no_mems
                )
            logger.info(
                f"Completed processing for {self.time_manager.current_time.strftime('%Y-%m-%d %H:%M')}"
            )
        logger.info("Pipeline execution completed.")


# Instantiate and run the pipeline
pipeline = FarmToDartPipeline()
pipeline.run_pipeline()
