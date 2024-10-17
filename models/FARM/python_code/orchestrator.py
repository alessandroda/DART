import shutil
import time
import pandas as pd
import os
from pathlib import Path
from datetime import timedelta
import logging

from orchestrator_utils import (
    check_job_status_cresco,
    replace_nml_template,
    run_command_in_directory,
    run_command_in_directory_bsub,
    searchFile,
    set_date_gregorian,
    prepare_farm_to_dart_nc,
    PathManager,
    TimeManager,
    submit_and_wait,
    prepare_dart_to_farm_nc,
)

logging.basicConfig(
    filename="farm_to_dart.log", format="%(asctime)s %(message)s", level=logging.INFO
)

logger = logging.getLogger(__name__)


class FarmToDartPipeline:
    def __init__(self):
        self.path_manager = PathManager(
            base_path="/gporq3/minni/FARM-DART/",
            env_python="miniconda3/bin/python3.10",
            listing_file="DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03__listing.csv",
            run_submit_farm_template="RUN/script/templates/new_run_submit.template.bsh",
            path_submit_bsh="RUN/script/auto_runs/",
            path_filter="DART/models/FARM/work/",
            path_data="RUN/data/",
            log_paths=True,
        )

        self.time_manager = TimeManager(
            start_time="2023-03-01 09:00:00",
            end_time="2023-03-01 18:00:00",
            dt_seconds=3600,
        )

        self.listing = pd.read_csv(self.path_manager.listing_file, sep=";")
        self.listing["start_time"] = pd.to_datetime(self.listing["start_time"])
        self.days_obs = 0
        self.seconds_obs = 0
        self.days_model = 0
        self.seconds_model = 0
        self.output_sim_folder = None
        self.ass_var = "c_SO2"

    def run_farm(self):
        logger.info(
            f"1.----------Running FARM for hour {self.time_manager.current_time}"
        )
        timestamp_farm = TimeManager.round_to_closest_hour(
            self.time_manager.current_time
        )
        string_to_replace_template = f'{timestamp_farm.strftime("%Y%m%d%H")}.bsh'

        command_farm_run = "run_submit" + string_to_replace_template
        path_run = (
            self.path_manager.path_submit_bsh
            / f"run_submit{string_to_replace_template}"
        )

        replace_nml_template(
            input_nml_path=self.path_manager.run_submit_farm_template,
            entries_tbr_dict={
                "da_date_start": timestamp_farm.strftime("%Y%m%d%H"),
                "da_date_end": timestamp_farm.strftime("%Y%m%d%H"),
            },
            output_nml_path=path_run,
        )

        commands_with_directories = [
            (command_farm_run, self.path_manager.path_submit_bsh)
        ]
        submit_and_wait(commands_with_directories)
        self.time_manager.simulated_time = timestamp_farm + timedelta(hours=1)

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
        return orbit_filename

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
                / f'DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/{orbit_filename["filename"].values[0]}',
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

    def run_dart(self, obs_seq_name):
        logger.info("Running DART")
        self.seconds_model, self.days_model = set_date_gregorian(
            self.time_manager.simulated_time.year,
            self.time_manager.simulated_time.month,
            self.time_manager.simulated_time.day,
            self.time_manager.simulated_time.hour,
            self.time_manager.simulated_time.minute,
            self.time_manager.simulated_time.second,
        )

        output_sim_folder = (
            self.path_manager.path_data
            / f"posteriors/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(output_sim_folder).mkdir(parents=True, exist_ok=True)

        replace_nml_template(
            self.path_manager.base_path / "DART/models/FARM/work/input_template.nml",
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": output_sim_folder,
                "$folder_obs_path": self.path_manager.base_path
                / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/",
                "$date_assim": self.time_manager.current_time.strftime("%Y%m%d_%H%M%S"),
                "$template_farm": self.path_manager.base_path
                / f"RUN/data/to_DART/ic_g1_{self.seconds_model}_{self.days_model}_00.nc",
                "$init_time_days": str(self.days_model),
                "$init_time_seconds": str(self.seconds_model),
                "$first_obs_days": str(self.days_obs),
                "$first_obs_seconds": str(self.seconds_obs),
            },
            output_nml_path=self.path_manager.base_path
            / "DART/models/FARM/work/input.nml",
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
            "./submit_filter.bsh", self.path_manager.path_submit_bsh
        )
        self.monitor_job(job_id, output_sim_folder)

    def monitor_job(self, job_id, output_sim_folder):
        logger.info(f"Monitoring job {job_id}")
        job_id = job_id.strip()[1:-1]

        while True:
            if check_job_status_cresco(job_id):
                print("Job completed successfully.")
                # Handle successful job completion: move files
                self.move_analysis_files(output_sim_folder)
                break
            else:
                print("Job is still running. Waiting...")
                time.sleep(10)

    def move_analysis_files(self, output_sim_folder):
        analysis_sim_folder = (
            self.path_manager.path_data
            / f"analysis/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        preassim_sim_folder = (
            self.path_manager.path_data
            / f"preassim/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )

        for filename in os.listdir(self.path_manager.path_filter):
            if filename.startswith("analysis_"):
                shutil.move(
                    os.path.join(self.path_manager.path_filter, filename),
                    analysis_sim_folder,
                )
            elif filename.startswith("preassim_"):
                shutil.move(
                    os.path.join(self.path_manager.path_filter, filename),
                    preassim_sim_folder,
                )

    def run_pipeline(self):
        logger.info("TIME LOOP BEGINS")
        while self.time_manager.current_time <= self.time_manager.end_time:
            self.run_farm()  # Run FARM executable

            orbit_filename = self.process_satellite_data()  # Process satellite data
            if orbit_filename:  # Only proceed if satellite data is found

                obs_seq_name = self.run_obs_converter(orbit_filename)

                prepare_farm_to_dart_nc(
                    self.path_manager,
                    self.time_manager.simulated_time,
                    self.time_manager.simulated_time,
                    self.seconds_model,
                    self.days_model,
                )

                self.run_dart(obs_seq_name)  # Run DART assimilation

                prepare_dart_to_farm_nc(
                    self.path_manager,
                    self.output_sim_folder,
                    self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                    self.ass_var,
                )
            logger.info(
                f"Completed processing for {self.time_manager.current_time.strftime('%Y-%m-%d %H:%M')}"
            )
        logger.info("Pipeline execution completed.")


# Instantiate and run the pipeline
pipeline = FarmToDartPipeline()
pipeline.run_pipeline()
