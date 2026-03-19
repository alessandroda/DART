from datetime import timedelta
from pathlib import Path
import shutil
import time

from pipeline_errors import FatalPipelineError, ModelRunError
from mimesi_types import Scheduler
from config_models import AppConfig
from pipelines.base_pipeline import BaseAssimilationPipeline
import logging
import os
import pandas as pd
import subprocess
from pipelines.chimere2017.paths import Chimere2017Paths
import xarray as xr
from time_utils import set_date_gregorian
from orchestrator_utils import (
    CommandSpec,
    check_job_status_cresco,
    check_job_status_slurm,
    get_list_mems_to_rerun,
    modify_yaml_date,
    TimeManager,
    replace_nml_template,
    replace_priorinflation,
    searchFile,
    run_command_in_directory,
    run_command_in_directory_bsub,
    submit_and_wait_cineca,
)


logger = logging.getLogger(__name__)


class Chimere2017DartPipeline(BaseAssimilationPipeline):
    def __init__(
        self,
        time_manager: TimeManager,
        paths: Chimere2017Paths,
        config: AppConfig,
    ):
        super().__init__(time_manager)
        # Load configuration from YAML file
        self.config = config
        # Setup path manager and time manager with loaded configurations
        self.paths = paths
        self.listing = pd.read_csv(self.paths.listing_file, sep=";")
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

        self.scheduler = self.config.cluster.scheduler
        logger.info(f"Using scheduler={self.scheduler}, queue={self.cineca_queue}")
        self._pending_orbit_filename = None
        self._pending_orbit_start = None
        self._pending_orbit_time = None

    def before_step(self):
        self.update_ibc_inputs()
        self.update_meteo_input()

    def cleanup_CHIMERE2017(self):
        pass

    def _read_boun_list_target(self, list_path: Path) -> Path:
        try:
            lines = [
                line.strip()
                for line in list_path.read_text().splitlines()
                if line.strip()
            ]
        except FileNotFoundError as e:
            raise FatalPipelineError(f"Missing IBC list file: {list_path}") from e
        except Exception as e:
            raise FatalPipelineError(f"Failed reading IBC list {list_path}: {e}") from e

        if not lines:
            raise FatalPipelineError(f"IBC list is empty: {list_path}")

        if len(lines) > 1 and lines[0].isdigit():
            return Path(lines[1])

        return Path(lines[-1])

    def _resolve_boun_daily_list(self, daily_list_name: str) -> Path:
        ibc_dir = self.paths.path_data / "basecase/IBC"
        candidate = ibc_dir / daily_list_name
        if candidate.exists():
            return candidate
        raise FatalPipelineError(
            f"Daily BOUN list not found in any IBC dir: {daily_list_name}"
        )

    def update_ibc_inputs(self):
        """
        Update INI/BOUN list files and extract hourly BOUN files for CHIMERE runs.
        """
        current_time = self.time_manager.current_time
        start_ts = current_time.strftime("%Y%m%d%H")
        end_ts = (current_time + timedelta(hours=1)).strftime("%Y%m%d%H")
        hour_index = current_time.hour

        daily_start = current_time.replace(hour=0, minute=0, second=0, microsecond=0)
        daily_end = daily_start + timedelta(days=1)

        daily_list_name = f"BOUN_CONCS.{daily_start:%Y%m%d%H}_{daily_end:%Y%m%d%H}_ITA7.list"
        
        daily_boun_list_path = self._resolve_boun_daily_list(daily_list_name)
        daily_boun_path = self._read_boun_list_target(daily_boun_list_path)

        if not daily_boun_path.exists():
            raise FatalPipelineError(
                f"Daily BOUN netcdf not found: {daily_boun_path}"
            )

        templates_dir = self.paths.path_submit_bsh / "templates"
        ini_template = (
            templates_dir / "INI_CONCS.YYYYMMDDHH_YYYYMMDDHH+dh_ITA7_template.list"
        )
        boun_template = (
            templates_dir / "BOUN_CONCS.YYYYMMDDHH_YYYYMMDDHH+dh_ITA7_template.list"
        )
        if not ini_template.exists():
            raise FatalPipelineError(f"Missing INI template: {ini_template}")
        if not boun_template.exists():
            raise FatalPipelineError(f"Missing BOUN template: {boun_template}")

        for mem in range(self.no_mems):
            mem_ibc_dir = self.paths.path_data / f"RUN_{mem}/IBC"
            mem_ibc_dir.mkdir(parents=True, exist_ok=True)

            hourly_boun_nc = (
                mem_ibc_dir / f"BOUN_CONCS.{start_ts}_{end_ts}_ITA7.nc"
            )
            try:
                subprocess.run(
                    [
                        "ncks",
                        "-O",
                        "-d",
                        f"Time,{hour_index},{hour_index+1}",
                        str(daily_boun_path),
                        str(hourly_boun_nc),
                    ],
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                raise FatalPipelineError(
                    f"ncks failed extracting hour {hour_index} from {daily_boun_path}"
                ) from e

            replace_nml_template(
                input_nml_path=str(boun_template),
                entries_tbr_dict={
                    "@mimesi_path_data": str(self.paths.path_data),
                    "@mimesi_ens_memeber": mem,
                    "@mimesi_start_date": start_ts,
                    "@mimesi_end_date": end_ts,
                },
                output_nml_path=str(
                    mem_ibc_dir / f"BOUN_CONCS.{start_ts}_{end_ts}_ITA7.list"
                ),
            )

            replace_nml_template(
                input_nml_path=str(ini_template),
                entries_tbr_dict={
                    "YYYYMMDDHH+dh": end_ts,
                    "YYYYMMDDHH": start_ts,
                    "@mimesi_path_data": str(self.paths.path_data),
                    "@mimesi_ens_member": mem,
                },
                output_nml_path=str(
                    mem_ibc_dir / f"INI_CONCS.{start_ts}_{end_ts}_ITA7.list"
                ),
            )

    def update_meteo_input(self):
        current_time = self.time_manager.current_time
        start_ts = current_time.strftime("%Y%m%d%H")
        end_ts = (current_time + timedelta(hours=1)).strftime("%Y%m%d%H")
        hour_index = current_time.hour

        daily_start = current_time.replace(hour=0, minute=0, second=0, microsecond=0)
        daily_end = daily_start + timedelta(days=1)

        meteo_dir = self.paths.path_data / "ATM"
        meteo_daily_name = (
            f"exdomout.{daily_start:%Y%m%d%H}_{daily_end:%Y%m%d%H}_ITA7.nc"
        )
        meteo_daily_path = meteo_dir / meteo_daily_name
        if not meteo_daily_path.exists():
            raise FatalPipelineError(f"Daily meteo netcdf not found: {meteo_daily_path}")

        meteo_dir.mkdir(parents=True, exist_ok=True)
        hourly_meteo_nc = meteo_dir / f"exdomout.{start_ts}_{end_ts}_ITA7.nc"
        try:
            subprocess.run(
                [
                    "ncks",
                    "-O",
                    "-d",
                    f"Time,{hour_index},{hour_index+1}",
                    str(meteo_daily_path),
                    str(hourly_meteo_nc),
                ],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise FatalPipelineError(
                f"ncks failed extracting hour {hour_index} from {meteo_daily_path}"
            ) from e


    def replace_perturb_into_original_emissions(self):
        pass

    def finalize_step(self):
        """
        Cleanup + YAML update.
        """
    
        modify_yaml_date(
            self.config._config_path,
            self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00"),
        )
    def _prepare_chimere_run_assets(self) -> tuple[Path, Path]:
        
        run_dir = self.paths.path_submit_bsh / "runs"
        pars_dir = self.paths.path_submit_bsh / "pars"
        templates_dir = self.paths.path_submit_bsh / "templates"
        
        lancia_script = (
            self.paths.path_submit_bsh / "lancia" / "lancia_chimere_m_nh.sh"
        )

        for d in (run_dir, pars_dir, templates_dir):
            d.mkdir(parents=True, exist_ok=True)

        
        par_template = self.paths.chimere_par_template
        if par_template is None:
            par_template = (
                self.paths.base_path
                / "catena_aria_test/config/mimesi/chimere.mimesi-ITA7_template.par"
            )
            logger.warning(
                f"chimere_par_template not set; using default: {par_template}"
            )
        if not par_template.exists():
            raise FatalPipelineError(f"CHIMERE par template not found: {par_template}")

        
        par_name = par_template.name
        if par_name.endswith("_template.par"):
            par_prefix = par_name.replace("_template.par", "")
        elif par_name.endswith(".par"):
            par_prefix = Path(par_name).stem
        else:
            par_prefix = par_name

        for mem in range(self.no_mems):
            output_par = pars_dir / f"{par_prefix}_{mem}.par"
            try:
                replace_nml_template(
                    input_nml_path=str(par_template),
                    entries_tbr_dict={"@mimesi_ens_member": mem},
                    output_nml_path=str(output_par),
                )
            except Exception as e:
                raise FatalPipelineError(
                    f"Failed to prepare par file {output_par}: {e}"
                )

        if not lancia_script.exists():
            raise FatalPipelineError(f"Missing lancia script: {lancia_script}")

        return run_dir, lancia_script


    def run_model(self):
        """
        Run CHIMERE for the current time step.
        """
        logger.info(f"[STEP] Running CHIMERE model at {self.time_manager.current_time}")

        timestamp_arg_run_chimere = self.time_manager.current_time.strftime(
            "%Y-%m-%d %H:00"
        )

        timestamp_chimere = TimeManager.round_to_closest_hour(
            self.time_manager.current_time
        ).strftime("%Y%m%d%H")

        run_dir, lancia_script = self._prepare_chimere_run_assets()

        file_run_ens = self.paths.chimere_name_run_sub_ens_bash(
            timestamp_chimere
        )
        try:
            replace_nml_template(
                input_nml_path=self.paths.base_path
                / self.paths.run_submit_model_template,
                entries_tbr_dict={
                    "@mimesi_dh_inizio": "0",  # this becomes variable
                    "@mimesi_nhours_list": "1",  # this becomes variable
                    "@mimesi_ens_size": self.no_mems,
                    "@mimesi_lancia_script": lancia_script
                },
                output_nml_path=self.paths.path_submit_bsh / file_run_ens,
            )
        except Exception as e:
            raise FatalPipelineError(f"Failed to prepare CHIMERE submit script: {e}")

        command = CommandSpec(
            command=file_run_ens,
            args=[timestamp_arg_run_chimere],
            directory=self.paths.path_submit_bsh,
        )

        job_ids = submit_and_wait_cineca(
            command,
        )

        mems_to_rerun = get_list_mems_to_rerun(
            job_ids=job_ids,
            scheduler=self.scheduler,
            model_type=self.model_type,
            path_manager=self.paths,
            timestamp_model=timestamp_chimere,
            no_mems=self.no_mems,
        )

        if mems_to_rerun:
            raise ModelRunError(f"Ensemble members failed: {mems_to_rerun}")

    def process_satellite_data(self):

        orbit_filename = searchFile(
            self.time_manager.current_time - 0.5 * self.time_manager.dt,
            self.time_manager.current_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return False

        orbit_filename["start_time"] = pd.to_datetime(orbit_filename["start_time"])
        
        logger.info(f"Orbit file found: {orbit_filename['filename'].values[0]}")
        self.time_manager.sat_obs = pd.to_datetime(
            orbit_filename["start_time"].values[0]
        )
        return orbit_filename["filename"].values[0]

    def _find_orbit_for_time(self, target_time):
        orbit_filename = searchFile(
            target_time - 0.5 * self.time_manager.dt,
            target_time + 0.5 * self.time_manager.dt,
            self.listing,
        )
        if orbit_filename.empty:
            return None

        orbit_filename["start_time"] = pd.to_datetime(orbit_filename["start_time"])
        start_time = pd.to_datetime(orbit_filename["start_time"].values[0])
        filename = orbit_filename["filename"].values[0]
        return filename, start_time

    def _clear_pending_orbit(self):
        self._pending_orbit_filename = None
        self._pending_orbit_start = None
        self._pending_orbit_time = None

    def run_obs_converter(self, orbit_filename):
        self.seconds_obs, self.days_obs = set_date_gregorian(
            self.time_manager.sat_obs.year,
            self.time_manager.sat_obs.month,
            self.time_manager.sat_obs.day,
            self.time_manager.sat_obs.hour,
            self.time_manager.sat_obs.minute,
            self.time_manager.sat_obs.second,
        )
        existing_obs_seq = {
            path.name for path in self.paths.dart_s5p_output_dir().glob("obs_seq_*.out")
        }
        obs_seq_name = self.time_manager.sat_obs.strftime("obs_seq_%Y%m%dT%H%M%S.out")
        self.paths.dart_s5p_output_dir().mkdir(parents=True, exist_ok=True)

        replace_nml_template(
            self.paths.dart_s5p_input_template(),
            entries_tbr_dict={
                "$file_path_s5p": self.paths.dart_file_s5p_orbit(orbit_filename),
                "$file_out": self.paths.dart_obs_seq(
                    obs_seq_name
                ),
                "$obs_type": self.obs_type,
            },
            output_nml_path=self.paths.dart_s5p_input(),
        )
        try:
            run_command_in_directory(
                CommandSpec(
                    command="convert_s5p_tropomi_l3",
                    directory=self.paths.dart_s5p_work(),
                )
            )
        except Exception as e:
            raise FatalPipelineError(f"Error running obs converter: {e}") from e

        generated_obs_seq = [
            path
            for path in self.paths.dart_s5p_output_dir().glob("obs_seq_*.out")
            if path.name not in existing_obs_seq
        ]
        if generated_obs_seq:
            return max(generated_obs_seq, key=lambda p: p.stat().st_mtime).name

        expected_obs_seq = self.paths.dart_obs_seq(
            self.time_manager.sat_obs.strftime("obs_seq_%Y%m%dT%H%M%S.out")
        )
        if expected_obs_seq.exists():
            return expected_obs_seq.name

        fallback_obs_seq = sorted(
            self.paths.dart_s5p_output_dir().glob("obs_seq_*.out"),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if fallback_obs_seq:
            logger.warning(
                "Could not uniquely identify newly generated obs_seq. Using latest file: %s",
                fallback_obs_seq[0].name,
            )
            return fallback_obs_seq[0].name

        logger.error("No obs_seq output found in %s", self.paths.dart_s5p_output_dir())
        return False

    def after_model(self):
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config; skipping after_model")
            return

        target_time = self.time_manager.current_time + self.time_manager.dt
        orbit_info = self._find_orbit_for_time(target_time)
        if not orbit_info:
            logger.info(
                "[DART] No satellite data for %s, skipping after_model",
                target_time.strftime("%Y%m%d%H"),
            )
            self._clear_pending_orbit()
            return

        self._pending_orbit_filename, self._pending_orbit_start = orbit_info
        self._pending_orbit_time = target_time

        to_dart_dir = self.paths.path_data / "to_DART"
        to_dart_dir.mkdir(parents=True, exist_ok=True)       

        for mem in range(self.no_mems):

            #temporary files
            temp_out_psfc = to_dart_dir / f"tmp_psfc_mem_{mem}.nc"
            tmp_ts = to_dart_dir / f"tmp_ts_mem_{mem}.nc"
            tmp_pres = to_dart_dir / f"tmp_pres_temp_mem{mem}.nc"

            out_file = self.paths.get_chimere_output_path(
                self.model_type,
                mem,
                self.time_manager.current_time,
                "out",
                1,
            )

            if not out_file.exists():
                logger.error("No CHIMERE out.*.nc files found after run_model().")
                return

            logger.info("Post-processing CHIMERE output (mem %s)", mem)

            # ------------------------------------------------
            # Select LAST timestep (IMPORTANT)
            # ------------------------------------------------
            subprocess.run(
                ["ncks", "-d", "Time,-1", out_file, tmp_ts],
                check=True,
            )

            # ------------------------------------------------
            # Compute surface pressure
            # ------------------------------------------------
            with xr.open_dataset(tmp_ts) as ds:
                ds = ds.load()
                k = 0
                ds["psfc"] = (ds.pres[:, k, :, :] - ds.a_vcoord[k] * 1e5) / ds.b_vcoord[k]

            ds.to_netcdf(temp_out_psfc)

            # ------------------------------------------------
            # Extract met fields needed by DART
            # ------------------------------------------------
            subprocess.run(
                ["cdo", "selname,pres,temp,psfc,lat,lon", temp_out_psfc, tmp_pres],
                check=True,
            )

            t1 = self.time_manager.current_time.strftime("%Y%m%d%H")
            tp = (self.time_manager.current_time + pd.Timedelta(hours=1)).strftime(
                "%Y%m%d%H"
            )
            out_file_with_mem = to_dart_dir / f"out.{t1}_{tp}_{mem}.nc"
            subprocess.run(["cdo", f"selname,{self.ass_var}", tmp_ts, out_file_with_mem], check=True)

            subprocess.run(
                ["ncks", "-A", tmp_pres, out_file_with_mem],
                check=True,
            )  
            subprocess.run([
                "ncatted", "-O",
                "-a", "_FillValue,NO2,o,f,NaN",
                out_file_with_mem
            ], check=True)
            # Remove ALL missing_value attributes
            subprocess.run([
            "ncatted", "-O",
            "-a", "missing_value,,d,,",
            out_file_with_mem
            ], check=True)            
            for f in [temp_out_psfc, tmp_ts, tmp_pres]:
                if f.exists():
                    f.unlink()


    def run_dart(self, obs_seq_name):
        logger.info("Running DART")

        self.output_sim_folder = self.paths.dart_posteriors_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        Path(self.output_sim_folder).mkdir(parents=True, exist_ok=True)
        filter_cores = (
            self.config.monitoring.cores
            if self.config.monitoring is not None
            else 20
        )
        input_start_time = self.time_manager.simulated_time - self.time_manager.dt
        input_end_time = self.time_manager.simulated_time
        t1 = input_start_time.strftime("%Y%m%d%H")
        tp = input_end_time.strftime("%Y%m%d%H")
        template_chimere_path = self.paths.path_data / "to_DART" / f"out.{t1}_{tp}_0.nc"

        replace_nml_template(
            self.paths.base_path
            / self.paths.path_filter
            / "input_template.nml",
            entries_tbr_dict={
                "$obs_sequence_name": obs_seq_name,
                "$folder_path": self.output_sim_folder,
                "$folder_obs_path": self.paths.dart_s5p_output_dir(),
                "$date_assim": self.time_manager.simulated_time.strftime("%Y%m%d_%H%M%S"),
                "$template_chimere": template_chimere_path,
                "$init_time_days": str(self.days_model),
                "$init_time_seconds": str(self.seconds_model),
                "$first_obs_days": str(self.days_obs),
                "$first_obs_seconds": str(self.seconds_obs),
                "$no_mems": str(self.no_mems),
                "$obs_type": str(self.obs_type),
                "$state_variable_conc": str(self.state_variable_conc),
                "$state_variable_qty": str(self.state_variable_qty),
            },
            output_nml_path=self.paths.base_path
            / self.paths.path_filter
            / "input.nml",
        )
        # FILTER_INPUT_LIST.TXT
        replace_nml_template(
            self.paths.base_path
            / self.paths.path_filter
            / "filter_input_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.paths.path_data / f"to_DART/",
                "$t1": t1,
                "$tp": tp,
            },
            output_nml_path=self.paths.base_path
            / self.paths.path_filter
            / "filter_input_list.txt",
        )
        logger.info(
            "[DART] filter_input_list interval start=%s end=%s simulated_time=%s",
            input_start_time.strftime("%Y-%m-%d %H:%M:%S"),
            input_end_time.strftime("%Y-%m-%d %H:%M:%S"),
            self.time_manager.simulated_time.strftime("%Y-%m-%d %H:%M:%S"),
        )

        # FILTER_OUTPUT_LIST.TXT
        replace_nml_template(
            self.paths.base_path
            / self.paths.path_filter
            / "filter_output_list_template.txt",
            entries_tbr_dict={
                "$folder_path": self.output_sim_folder,
                "$date": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
            },
            output_nml_path=self.paths.base_path
            / self.paths.path_filter
            / "filter_output_list.txt",
        )
        # SUBMIT_FILTER.BSH
        submit_filter_entries = {
            "CORES": str(filter_cores),
            "CURRENT_DATE": self.time_manager.simulated_time.strftime("%Y%m%d%H"),
        }
        if self.scheduler == Scheduler.SLURM:
            submit_filter_entries.update(
                {
                    "SBATCH_PARTITION": self.cineca_queue,
                    "DEST_LOG_PATH": self.output_sim_folder,
                }
            )
        else:
            submit_filter_entries.update(
                {
                    "QUEUE": self.cineca_queue,
                    "DEST_LOG_PATH": self.output_sim_folder,
                }
            )

        replace_nml_template(
            self.paths.base_path / self.paths.path_submit_bsh / "templates/run_mimesi-filter_template.bsh",
            entries_tbr_dict=submit_filter_entries,
            output_nml_path=self.paths.path_submit_bsh / "submit_filter.bsh",
        )
        

        if self.scheduler == Scheduler.SLURM:
            job_ids = submit_and_wait_cineca(
                CommandSpec(
                    command="./submit_filter.bsh",
                    directory=self.paths.path_submit_bsh,
                )
            )
            time.sleep(10)
            self.monitor_job_dart(job_ids)
        else:
            job_id = run_command_in_directory_bsub(
                "./submit_filter.bsh", self.paths.path_submit_bsh, farm=False
            )
            time.sleep(10)
            self.monitor_job_dart(job_id)

    def monitor_job_dart(self, job_id):
        if self.scheduler == Scheduler.SLURM:
            job_ids = job_id if isinstance(job_id, (list, tuple)) else [job_id]
            logger.info(f"Monitoring SLURM jobs {job_ids}")

            while True:
                running_jobs = []
                for jid in job_ids:
                    if not check_job_status_slurm(jid, which_run="DART"):
                        running_jobs.append(jid)

                if not running_jobs:
                    print("Job completed successfully.")
                    self.move_analysis_files()
                    replace_priorinflation(
                        self.paths,
                        self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                    )
                    break

                print(f"Jobs still running: {running_jobs}. Waiting...")
                time.sleep(10)
        else:
            logger.info(f"Monitoring job {job_id}")
            if isinstance(job_id, (list, tuple)):
                job_id = job_id[0]
            job_id = job_id.strip()[1:-1]

            while True:
                if check_job_status_cresco(job_id, which_run="DART"):
                    print("Job completed successfully.")
                    self.move_analysis_files()
                    replace_priorinflation(
                        self.paths,
                        self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                    )
                    break
                else:
                    print("Job is still running. Waiting...")
                    time.sleep(10)

    def move_analysis_files(self):
        analysis_sim_folder = (
            self.paths.path_data
            / f"analysis/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(analysis_sim_folder).mkdir(parents=True, exist_ok=True)
        preassim_sim_folder = (
            self.paths.path_data
            / f"preassim/{self.time_manager.simulated_time.strftime('%Y%m%d%H')}"
        )
        Path(preassim_sim_folder).mkdir(parents=True, exist_ok=True)
        for filename in os.listdir(f"{self.paths.path_filter}"):
            if filename.startswith("analysis_"):
                try:
                    shutil.move(
                        os.path.join(
                            self.paths.path_filter,
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
                            self.paths.path_filter,
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

        orbit_filename = None
        if (
            self._pending_orbit_time is not None
            and self._pending_orbit_time == self.time_manager.current_time
            and self._pending_orbit_filename
        ):
            orbit_filename = self._pending_orbit_filename
            self.time_manager.sat_obs = self._pending_orbit_start
        else:
            orbit_filename = self.process_satellite_data()
        self._clear_pending_orbit()
        if not orbit_filename:
            logger.info("[DART] No satellite data found, skipping assimilation")
            return

        obs_seq_name = self.run_obs_converter(orbit_filename)

        if not obs_seq_name:
            logger.info("[DART] Observation conversion failed, skipping assimilation")
            return

        obs_path = self.paths.dart_s5p_output_dir() / obs_seq_name
        if not obs_path.exists():
            logger.info("[DART] obs_seq not found, skipping assimilation")
            return

        self.run_dart(obs_seq_name)

from pipelines.registry import register_pipeline
from pipelines.chimere2017.config import Chimere2017PipelineConfig


class Chimere2017Pipeline(Chimere2017DartPipeline):
    """Concrete CHIMERE 2017 pipeline with in-folder implementation."""


def _build(config, time_manager):
    p_cfg = Chimere2017PipelineConfig.from_app_config(config)
    p_paths = Chimere2017Paths(p_cfg)
    return Chimere2017Pipeline(time_manager, p_paths, config)


register_pipeline("chimere2017", _build)
