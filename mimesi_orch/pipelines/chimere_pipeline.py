from datetime import timedelta
from pathlib import Path
import shutil
import time

from mimesi_orch.scheduler.scheduler import Scheduler
from mimesi_orch.config_models import AppConfig
from mimesi_orch.paths import PathManager
from pipelines.base_pipeline import BaseAssimilationPipeline
import logging
import os
import pandas as pd
from mimesi_orch.orchestrator_utils import (
    check_job_status_cresco,
    modify_yaml_date,
    filter_dates,
    TimeManager,
    prepare_farm_to_dart_nc_par,
    replace_nml_template,
    replace_priorinflation,
    searchFile,
    set_date_gregorian,
    submit_and_wait,
    run_command_in_directory,
    run_command_in_directory_bsub,
    submit_and_wait_slurm,
)


logger = logging.getLogger(__name__)


class Chimere2017DartPipeline(BaseAssimilationPipeline):
    def __init__(
        self,
        time_manager: TimeManager,
        path_manager: PathManager,
        config: AppConfig,
    ):
        super().__init__(time_manager)
        # Load configuration from YAML file
        self.config = config
        # Setup path manager and time manager with loaded configurations
        self.path_manager = path_manager
        self.listing = pd.read_csv(self.path_manager.listing_file, sep=";")
        self.listing["start_time"] = pd.to_datetime(self.listing["start_time"])
        self.days_obs = 0
        self.seconds_obs = 0
        self.days_model = 0
        self.seconds_model = 0
        self.output_sim_folder = None
        a = self.config.assimilation

        self.model_type = a.model_type
        logger.info(f"Running assimilation with model_type={self.model_type}")

        self.ass_var = a.ass_var
        self.emi_var = a.emi_var
        self.no_mems = a.no_mems
        self.case_dir = a.case_dir
        self.obs_type = a.obs_type
        self.state_variable_conc = a.ass_var
        self.state_variable_qty = a.state_variable_qty
        self.run_assimilation_flag = a.run_assimilation_flag
        self.case_emi_dir = a.case_emi_dir

        self.cineca_queue = self.config.cluster.cluster_queue

        self.backup_perturb_days = self.config.time.backup_perturb_days
        self.backup_ic_hours = self.config.time.backup_ic_hours
        self.backup_ic_option = self.config.time.backup_ic_option

        logger.info(f"Using scheduler={self.scheduler}, queue={self.cineca_queue}")

        self.scheduler = self.config.cluster.scheduler

    def before_step(self):
        pass

    def cleanup_CHIMERE2017(self):
        pass

    def replace_perturb_into_original_emissions(self):
        pass

    def finalize_step(self):
        """
        Cleanup + YAML update.
        """

        modify_yaml_date(
            self.config["_config_path"],
            self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00"),
        )
        self.cleanup_FARM()

    # former run_farm
    def run_model(self):
        """
        Run CHIMERE for the current time step.
        """
        logger.info(f"CHIMERE running model at {self.time_manager.current_time}")
        timestamp_arg_run_chimere = self.time_manager.current_time.strftime(
            "%Y-%m-%d %H:00"
        )
        # to be understood time_emi, date_emi
        timestamp_chimere = TimeManager.round_to_closest_hour(
            self.time_manager.current_time
        ).strftime("%Y%m%d%H")
        commands_with_directories = []
        for mem in self.no_mems:
            string_to_replace_template = f"{timestamp_chimere}_mem_{mem}.sh"
            command_chimere_run = (
                "run_mimesi-ITA7_"
                + string_to_replace_template
                + f" '{timestamp_arg_run_chimere}'"
            )
            path_run = self.path_manager.chimere_name_run_sub_ens_bash(
                string_to_replace_template
            )
            replace_nml_template(
                input_nml_path=self.path_manager.base_path
                / self.path_manager.run_submit_model_template,
                entries_tbr_dict={
                    "@mimesi_dh_inizio": "1",
                    "@mimesi_nhours_list": "1",
                },
                output_nml_path=path_run,
            )
            commands_with_directories.append(
                (command_chimere_run, self.path_manager.chimere_base_dir)
            )
        submit_and_wait_slurm(
            self.model_type,
            self.path_manager,
            self.scheduler,
            commands_with_directories,
            timestamp_chimere,
            self.no_mems,
            self.case_dir,
            self.cineca_queue,
        )

        # ./run_mimesi-ITA7.sh '2026-01-26 0:00'
        # il bash esegue un altro bash lancia_chimere_m_nh.sh
        # dentro a questo bash si esegue chimere.sh
        # chimere .sh esegue finalmente lo slurm
        #         export NP

        # sbatch --wait --job-name=${dom}.${idatestart}.${simclab} \
        #         --account=arpae_aqm \
        #         --output=${job_o_log} --error=${job_e_log}    ${chimere_root}/scripts/run_chimere.job
        # exitstato_run=$?
        # set +x

        # else
        # echo "No such file ${chimere_tmp}/chimere.e ! Bye."
        # exit 1
        # fi

        pass

    def process_satellite_data(self):

        orbit_filename = searchFile(
            self.time_manager.current_time - 0.5 * self.time_manager.dt,
            self.time_manager.current_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return False

        orbit_filename["start_time"] = pd.to_datetime(orbit_filename["start_time"])
        orbit_filename = orbit_filename[orbit_filename["start_time"].dt.hour >= 10]

        if orbit_filename.empty:
            logger.info(
                f'No valid orbit file found after 10 AM.{orbit_filename["start_time"]}'
            )
            return False

        logger.info(f"Orbit file found: {orbit_filename['filename'].values[0]}")
        self.time_manager.sat_obs = pd.to_datetime(
            orbit_filename["start_time"].values[0]
        )
        return orbit_filename["filename"].values[0]

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
            self.path_manager.dart_s5p_input_template(),
            entries_tbr_dict={
                "$file_path_s5p": self.path_manager.dart_file_s5p_orbit(orbit_filename),
                "$file_out": self.path_manager.dart_obs_seq(
                    self.seconds_obs, self.days_obs
                ),
                "$obs_type": self.obs_type,
            },
            output_nml_path=self.path_manager.dart_s5p_input(),
        )
        try:
            run_command_in_directory(
                "convert_s5p_tropomi_l3",
                self.path_manager.dart_s5p_work(),
            )
        except Exception as e:
            logger.error(f"Error running obs converter: {e}")
            return False
        return obs_seq_name

    def after_model(self):
        prepare_farm_to_dart_nc_par(
            self.path_manager,
            self.time_manager.simulated_time,
            self.time_manager.simulated_time,
            self.seconds_model,
            self.days_model,
            self.no_mems,
        )

    def run_dart(self, obs_seq_name):
        logger.info("Running DART")

        self.output_sim_folder = self.path_manager.dart_posteriors_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        Path(self.output_sim_folder).mkdir(parents=True, exist_ok=True)

        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "input_template.nml",
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.output_sim_folder,
                "$folder_obs_path": self.path_manager.dart_s5p_output_dir(),
                "$date_assim": self.time_manager.current_time.strftime("%Y%m%d_%H%M%S"),
                "$template_farm": self.path_manager.path_data
                / f"to_DART/ic_g1_{self.seconds_model}_{self.days_model}_0.nc",
                "$init_time_days": str(self.days_model),
                "$init_time_seconds": str(self.seconds_model),
                "$first_obs_days": str(self.days_obs),
                "$first_obs_seconds": str(self.seconds_obs),
                "$no_mems": str(self.no_mems),
                "$obs_type": str(self.obs_type),
                "$state_variable_conc": str(self.state_variable_conc),
                "$state_variable_qty": str(self.state_variable_qty),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "input.nml",
        )
        # FILTER_INPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_input_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.path_manager.path_data / f"to_DART/",
                "$days": str(self.seconds_model),
                "$seconds": str(self.days_model),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_input_list.txt",
        )

        # FILTER_OUTPUT_LIST.TXT
        replace_nml_template(
            self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_output_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.output_sim_folder,
                "$date": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
            },
            output_nml_path=self.path_manager.base_path
            / self.path_manager.path_filter
            / "filter_output_list.txt",
        )
        # SUBMIT_FILTER.BSH
        replace_nml_template(
            self.path_manager.base_path
            / "RUN/script/templates/submit_filter.template.bsh",
            entries_tbr_dict={
                "CURRENT_DATE": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                "CORES": str(20),
                "QUEUE": self.cresco_queue,
                "DEST_LOG_PATH": self.output_sim_folder,
            },
            output_nml_path=self.path_manager.path_submit_bsh / "submit_filter.bsh",
        )
        # RUN_FILTER.BSH
        replace_nml_template(
            self.path_manager.base_path
            / "RUN/script/templates/run_filter.template.bsh",
            entries_tbr_dict={
                "CORES": str(20),
                "@ABS_FILTER_PATH": self.path_manager.base_path
                / self.path_manager.path_filter,
            },
            output_nml_path=self.path_manager.path_submit_bsh / "run_filter.bsh",
        )

        job_id = run_command_in_directory_bsub(
            "./submit_filter.bsh", self.path_manager.path_submit_bsh, farm=False
        )
        time.sleep(10)
        self.monitor_job_dart(job_id)

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

    def move_analysis_files(self):
        analysis_sim_folder = (
            self.path_manager.path_data
            / f"analysis/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
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
                        os.path.join(analysis_sim_folder, filename),
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

    def run_assimilation_if_needed(self):
        """
        Run DART assimilation only if:
        - assimilation is enabled
        - satellite observations are available
        """
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config")
            return

        orbit_filename = self.process_satellite_data()
        if not orbit_filename:
            logger.info("[DART] No satellite data found, skipping assimilation")
            return

        obs_seq_name = self.run_obs_converter(orbit_filename)

        obs_path = (
            self.path_manager.base_path
            / f"DART/observations/obs_converters/S5P_TROPOMI_L3/data/SO2-COBRA/C03dart/{obs_seq_name}"
        )
        if not obs_path.exists():
            logger.info("[DART] obs_seq not found, skipping assimilation")
            return

        self.run_dart(obs_seq_name)
