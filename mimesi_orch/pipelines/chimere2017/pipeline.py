from datetime import timedelta
from pathlib import Path
import shutil
import time
import numpy as np
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
    get_list_mems_to_rerun,
    modify_yaml_date,
    replace_nml_template,
    replace_priorinflation,
    searchFile,
    run_command_in_directory,
    submit_and_wait_cineca,
)
from pipeline_time import AssimWindow, TimeManager


logger = logging.getLogger(__name__)

_CDFCHECK_5D_PREFIX = (
    "Warning (cdfCheckVars): 5 dimensional variables are not supported"
)


def _run_cdo_suppress_5d(args):
    result = subprocess.run(
        ["cdo", *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.stderr:
        lines = [line.strip() for line in result.stderr.splitlines() if line.strip()]
        filtered = [
            line for line in lines if not line.startswith(_CDFCHECK_5D_PREFIX)
        ]
        if filtered:
            logger.warning("cdo warnings: %s", " | ".join(filtered))
        else:
            logger.debug("cdo warnings suppressed: %s", " | ".join(lines))
    return result


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
        self.emi_perturbations = getattr(a, "emi_perturbations", None)
        self.emi_perturbation_dir = a.emi_perturbation_dir
        self.cineca_queue = self.config.cluster.cluster_queue

        self.backup_perturb_days = self.config.time.backup_perturb_days
        self.backup_ic_hours = self.config.time.backup_ic_hours
        self.backup_ic_option = self.config.time.backup_ic_option

        self.scheduler = self.config.cluster.scheduler
        logger.info(f"Using scheduler={self.scheduler}, queue={self.cineca_queue}")
        self._pending_orbit_filename = None
        self._pending_orbit_start = None
        self._pending_orbit_time = None
        self._generated_daily_emission_files: list[Path] = []
        self._generated_daily_emission_stamp: str | None = None

    def build_assim_window(self) -> AssimWindow:
        """
        Build the current cycle window.

        Windows start and end on model hour boundaries. Observations are
        assigned to the nearest cycle end using the +/- 30 minute rule.
        """
        start_time = self.time_manager.current_time
        
        if start_time.hour == 0:
            return AssimWindow(
                start_time=start_time,
                end_time=start_time + self.time_manager.dt,
                run_hours=1,
                has_assimilation=False,   # optional
            )
        half_dt = self.time_manager.dt / 2
        day_end = start_time.replace(
            hour=0,
            minute=0,
            second=0,
            microsecond=0,
        ) + timedelta(days=1)
        run_limit = min(day_end, self.time_manager.end_time + self.time_manager.dt)

        slot_time = start_time + self.time_manager.dt
        while slot_time <= run_limit:
            orbit_matches = searchFile(
                slot_time - half_dt,
                slot_time + half_dt,
                self.listing,
            )
            if not orbit_matches.empty:
                orbit_matches = orbit_matches.copy()
                orbit_matches["start_time"] = pd.to_datetime(
                    orbit_matches["start_time"]
                )
                orbit_matches = orbit_matches.sort_values("start_time")
                orbit_row = orbit_matches.iloc[0]
                run_hours = int((slot_time - start_time) / self.time_manager.dt)
                return AssimWindow(
                    start_time=start_time,
                    end_time=slot_time,
                    run_hours=run_hours,
                    has_assimilation=True,
                    obs_time=pd.to_datetime(orbit_row["start_time"]),
                    orbit_filename=orbit_row["filename"],
                )
            slot_time += self.time_manager.dt

        if run_limit <= start_time:
            raise FatalPipelineError(
                f"Invalid run limit for assimilation window: start={start_time} end={run_limit}"
            )

        return AssimWindow(
            start_time=start_time,
            end_time=run_limit,
            run_hours=int((run_limit - start_time) / self.time_manager.dt),
            has_assimilation=False,
        )

    def before_step(self):
        self.replace_perturb_into_original_emissions()
        self.update_ibc_inputs()
        self.update_meteo_inputs()
        self.update_emission_inputs()

    def _subset_netcdf_inplace(self, nc_path: Path, keep_vars: list[str]) -> None:
        """
        Reduce a NetCDF file to a minimal variable set.

        Uses NCO (ncks) when available; otherwise falls back to xarray.
        """
        if not nc_path.exists():
            logger.warning("[CLEANUP] NetCDF not found (skip subsetting): %s", nc_path)
            return

        keep_vars = [str(v).strip() for v in keep_vars if str(v).strip()]
        if not keep_vars:
            logger.warning("[CLEANUP] Empty keep_vars for %s (skip subsetting)", nc_path)
            return

        try:
            with xr.open_dataset(nc_path, decode_cf=False) as ds:
                available = set(ds.variables)
        except Exception as e:
            logger.warning(
                "[CLEANUP] Failed reading NetCDF header for %s (skip subsetting): %s",
                nc_path,
                e,
            )
            return

        existing = [v for v in keep_vars if v in available]
        if not existing:
            logger.warning(
                "[CLEANUP] None of the requested variables exist in %s. Requested=%s",
                nc_path,
                ",".join(keep_vars),
            )
            return

        tmp_path = nc_path.with_suffix(f"{nc_path.suffix}.tmp")
        ncks = shutil.which("ncks")
        if ncks:
            try:
                subprocess.run(
                    [
                        "ncks",
                        "-O",
                        "-v",
                        ",".join(existing),
                        str(nc_path),
                        str(tmp_path),
                    ],
                    check=True,
                )
                tmp_path.replace(nc_path)
                return
            except subprocess.CalledProcessError as e:
                logger.warning(
                    "[CLEANUP] ncks failed subsetting %s (will try xarray fallback): %s",
                    nc_path,
                    e,
                )
            except Exception as e:
                logger.warning(
                    "[CLEANUP] Failed replacing subset NetCDF for %s (will try xarray fallback): %s",
                    nc_path,
                    e,
                )

        try:
            with xr.open_dataset(nc_path, decode_cf=False) as ds:
                ds[existing].to_netcdf(tmp_path)
            tmp_path.replace(nc_path)
        except Exception as e:
            logger.warning("[CLEANUP] xarray subsetting failed for %s: %s", nc_path, e)
            if tmp_path.exists():
                tmp_path.unlink()

    def _cleanup_window_inputs(self) -> None:
        """
        Remove window-specific inputs created in before_step().

        Safe during the time loop:
        - deletes IBC window files
        - deletes emission *window* slice (keeps daily file)
        - deletes meteo window slice
        """
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        cleanup_cfg = getattr(self.config, "cleanup", None)
        if cleanup_cfg is None or not cleanup_cfg.enabled:
            return

        start_ts = self.current_window.start_time.strftime("%Y%m%d%H")
        end_ts = self.current_window.end_time.strftime("%Y%m%d%H")

        # IBC (per member)
        if cleanup_cfg.delete_window_ibc:
            for mem in range(self.no_mems):
                ibc_dir = self.paths.path_data / f"RUN_{mem}/IBC"
                if not ibc_dir.exists():
                    continue
                for name in (
                    f"BOUN_CONCS.{start_ts}_{end_ts}_ITA7.nc",
                    f"BOUN_CONCS.{start_ts}_{end_ts}_ITA7.list",
                    f"INI_CONCS.{start_ts}_{end_ts}_ITA7.list",
                ):
                    p = ibc_dir / name
                    if p.exists():
                        try:
                            p.unlink()
                        except Exception as e:
                            logger.warning("[CLEANUP] Failed removing %s: %s", p, e)

        # Emissions window slice (per member)
        if cleanup_cfg.delete_window_emissions:
            for mem in range(self.no_mems):
                emi_dir = self.paths.path_data / f"RUN_{mem}/EMISSION_{mem}"
                p = emi_dir / f"AEMISSIONS.{start_ts}_{end_ts}_ITA7.nc"
                if p.exists():
                    try:
                        p.unlink()
                    except Exception as e:
                        logger.warning("[CLEANUP] Failed removing %s: %s", p, e)

        # Meteo window slice (shared)
        if cleanup_cfg.delete_window_meteo:
            atm_dir = self.paths.path_data / "ATM"
            meteo_slice = atm_dir / f"exdomout.{start_ts}_{end_ts}_ITA7.nc"
            if meteo_slice.exists():
                try:
                    meteo_slice.unlink()
                except Exception as e:
                    logger.warning("[CLEANUP] Failed removing %s: %s", meteo_slice, e)

    def _cleanup_trim_outputs(self) -> None:
        """
        Trim CHIMERE window outputs (out/end) to keep only requested variables.
        """
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        cleanup_cfg = getattr(self.config, "cleanup", None)
        if cleanup_cfg is None or not cleanup_cfg.enabled:
            return

        start_time = self.current_window.start_time
        run_hours = self.current_window.run_hours

        def _default_keep_list(raw: list[str] | None) -> list[str]:
            if raw:
                base = list(raw)
            else:
                base = []
            # Always try to keep basic coords/time if present.
            for v in ("lon", "lat", "Times", "Time"):
                if v not in base:
                    base.insert(0, v)
            # Preserve ensemble member coordinate if present in file.
            if "member" not in base:
                base.append("member")
            return base

        if cleanup_cfg.trim_out:
            keep_out = _default_keep_list(cleanup_cfg.keep_out_vars)
            for mem in range(self.no_mems):
                out_path = self.paths.get_chimere_output_path(
                    self.model_type, mem, start_time, "out", run_hours
                )
                self._subset_netcdf_inplace(out_path, keep_out)

        if cleanup_cfg.trim_end:
            keep_end = _default_keep_list(cleanup_cfg.keep_end_vars)
            for mem in range(self.no_mems):
                end_path = self.paths.get_chimere_output_path(
                    self.model_type, mem, start_time, "end", run_hours
                )
                self._subset_netcdf_inplace(end_path, keep_end)

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
            target = Path(lines[1])
        else:
            target = Path(lines[-1])

        if not target.is_absolute():
            target = (list_path.parent / target).resolve()
        return target

    def _resolve_boun_daily_list(self, daily_list_name: str) -> Path:
        ibc_dir = getattr(self.paths, "chimere_input_ibc_dir", None)
        if ibc_dir is None:
            raise FatalPipelineError(
                "paths.chimere_input_ibc_dir is not configured (required for CHIMERE IBC inputs)"
            )

        candidate = ibc_dir / daily_list_name
        if candidate.exists():
            return candidate

        raise FatalPipelineError(f"Daily BOUN list not found: {candidate}")

    def _get_same_day_window_bounds(self):
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        start_time = self.current_window.start_time
        end_time = self.current_window.end_time
        if end_time <= start_time:
            raise FatalPipelineError(
                f"Invalid assimilation window: start={start_time} end={end_time}"
            )

        daily_start = start_time.replace(hour=0, minute=0, second=0, microsecond=0)
        daily_end = daily_start + timedelta(days=1)
        if end_time > daily_end:
            raise FatalPipelineError(
                "Cross-day assimilation windows are not supported yet: "
                f"start={start_time} end={end_time}"
            )

        start_ts = start_time.strftime("%Y%m%d%H")
        end_ts = end_time.strftime("%Y%m%d%H")

        start_index = start_time.hour
        if self.current_window.run_hours < 1:
            raise FatalPipelineError(
                f"Invalid run_hours for assimilation window: {self.current_window.run_hours}"
            )
        end_index = start_index + self.current_window.run_hours
        return (
            start_time,
            end_time,
            daily_start,
            daily_end,
            start_ts,
            end_ts,
            start_index,
            end_index,
        )

    def _slice_time_window(self, source_path: Path, output_path: Path, start_index: int, end_index: int):
        try:
            subprocess.run(
                [
                    "ncks",
                    "-O",
                    "-d",
                    f"Time,{start_index},{end_index}",
                    str(source_path),
                    str(output_path),
                ],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise FatalPipelineError(
                f"ncks failed extracting Time,{start_index},{end_index} from {source_path}"
            ) from e

    def update_ibc_inputs(self):
        """
        Update INI/BOUN list files and extract same-day BOUN files for CHIMERE runs.
        """
        _, _, daily_start, daily_end, start_ts, end_ts, start_index, end_index = (
            self._get_same_day_window_bounds()
        )

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

            window_boun_nc = (
                mem_ibc_dir / f"BOUN_CONCS.{start_ts}_{end_ts}_ITA7.nc"
            )
            self._slice_time_window(
                daily_boun_path,
                window_boun_nc,
                start_index,
                end_index,
            )

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

    def update_meteo_inputs(self):
        _, _, daily_start, daily_end, start_ts, end_ts, start_index, end_index = (
            self._get_same_day_window_bounds()
        )

        meteo_src_dir = getattr(self.paths, "chimere_input_atm_dir", None)
        if meteo_src_dir is None:
            raise FatalPipelineError(
                "paths.chimere_input_atm_dir is not configured (required for CHIMERE ATM inputs)"
            )
        meteo_out_dir = self.paths.path_data / "ATM"
        meteo_daily_name = (
            f"exdomout.{daily_start:%Y%m%d%H}_{daily_end:%Y%m%d%H}_ITA7.nc"
        )
        meteo_daily_path = meteo_src_dir / meteo_daily_name
        if not meteo_daily_path.exists():
            raise FatalPipelineError(f"Daily meteo netcdf not found: {meteo_daily_path}")

        meteo_out_dir.mkdir(parents=True, exist_ok=True)
        window_meteo_nc = meteo_out_dir / f"exdomout.{start_ts}_{end_ts}_ITA7.nc"
        self._slice_time_window(
            meteo_daily_path,
            window_meteo_nc,
            start_index,
            end_index,
        )

    def update_emission_inputs(self):
        _, _, daily_start, daily_end, start_ts, end_ts, start_index, end_index = (
            self._get_same_day_window_bounds()
        )
        daily_stamp = f"{daily_start:%Y%m%d%H}_{daily_end:%Y%m%d%H}"

        logger.info(
            "[EMISSIONS] Preparing files for %s -> %s from daily window %s (Time=%d:%d)",
            start_ts,
            end_ts,
            daily_stamp,
            start_index,
            end_index,
        )

        for mem in range(self.no_mems):
            perturbed_root = self.paths.path_data / f"RUN_{mem}/EMISSION_{mem}"
            if perturbed_root is None:
                raise FatalPipelineError("paths.path_perturbed_emi is not configured")


            mem_emi_dir = self.paths.path_data / f"RUN_{mem}" / f"EMISSION_{mem}"
            mem_emi_dir.mkdir(parents=True, exist_ok=True)

            emi_daily_path = mem_emi_dir / f"AEMISSIONS.{daily_stamp}_ITA7.nc"
            if not emi_daily_path.exists():
                raise FatalPipelineError(
                    f"Daily emission netcdf not found for member {mem}: {emi_daily_path}"
                )

            window_emi_nc = mem_emi_dir / f"AEMISSIONS.{start_ts}_{end_ts}_ITA7.nc"
            logger.info(
                "[EMISSIONS] Member %s extracting emission file %s from %s",
                mem,
                window_emi_nc.name,
                emi_daily_path,
            )
            self._slice_time_window(
                emi_daily_path,
                window_emi_nc,
                start_index,
                end_index,
            )

    def replace_perturb_into_original_emissions(self):
        current_time = self.time_manager.current_time
        if current_time.hour != 0:
            logger.debug(
                "[EMISSIONS] Skipping daily emission replacement at %s because hour != 00",
                current_time.strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        perturbed_root = self.paths.path_perturbed_emi
        if perturbed_root is None:
            raise FatalPipelineError("paths.path_perturbed_emi is not configured")
        emissions_src_dir = getattr(self.paths, "chimere_input_emissions_dir", None)
        if emissions_src_dir is None:
            raise FatalPipelineError(
                "paths.chimere_input_emissions_dir is not configured (required for CHIMERE emissions inputs)"
            )
        emi_pairs: list[tuple[str, str]] = []
        if self.emi_perturbations is not None:
            if not isinstance(self.emi_perturbations, dict):
                raise FatalPipelineError("assimilation.emi_perturbations must be a mapping")
            if not self.emi_perturbations:
                raise FatalPipelineError("assimilation.emi_perturbations is empty")
            for key, value in self.emi_perturbations.items():
                emi_var = str(key).strip()
                emi_dir = str(value).strip()
                if not emi_var:
                    raise FatalPipelineError("assimilation.emi_perturbations contains an empty variable key")
                if not emi_dir:
                    raise FatalPipelineError(
                        f"assimilation.emi_perturbations contains an empty dir for variable {emi_var}"
                    )
                emi_pairs.append((emi_var, emi_dir))
        else:
            if self.emi_var is None:
                raise FatalPipelineError("assimilation.emi_var is not configured")
            if self.emi_perturbation_dir is None:
                raise FatalPipelineError("assimilation.emi_perturbation_dir is not configured")

            raw_emi_vars = (
                [self.emi_var] if isinstance(self.emi_var, str) else list(self.emi_var)
            )
            raw_emi_dirs = (
                [self.emi_perturbation_dir]
                if isinstance(self.emi_perturbation_dir, str)
                else list(self.emi_perturbation_dir)
            )
            emi_vars = [str(v).strip() for v in raw_emi_vars]
            emi_dirs = [str(d).strip() for d in raw_emi_dirs]

            empty_vars = [idx for idx, v in enumerate(emi_vars) if not v]
            empty_dirs = [idx for idx, d in enumerate(emi_dirs) if not d]
            if empty_vars:
                raise FatalPipelineError(
                    "assimilation.emi_var contains empty entries at indices "
                    + ",".join(map(str, empty_vars))
                )
            if empty_dirs:
                raise FatalPipelineError(
                    "assimilation.emi_perturbation_dir contains empty entries at indices "
                    + ",".join(map(str, empty_dirs))
                )
            if len(emi_vars) != len(emi_dirs):
                raise FatalPipelineError(
                    "assimilation.emi_var and assimilation.emi_perturbation_dir must have the same length "
                    f"(got emi_var={len(emi_vars)} and emi_perturbation_dir={len(emi_dirs)})"
                )
            emi_pairs = list(zip(emi_vars, emi_dirs))

        daily_start = current_time.replace(hour=0, minute=0, second=0, microsecond=0)
        daily_end = daily_start + timedelta(days=1)
        date = daily_start.strftime("%Y%m%d%H")
        datep1 = daily_end.strftime("%Y%m%d%H")
        daily_stamp = f"{date}_{datep1}"
        file_original_name = f"AEMISSIONS.{daily_stamp}_ITA7.nc"

        base_file = emissions_src_dir / file_original_name
        if not base_file.exists():
            raise FatalPipelineError(f"Base emission file not found: {base_file}")

        logger.info(
            "[EMISSIONS] Building daily perturbed emission files for %s using base=%s perturbations=%s variables=%s",
            daily_stamp,
            base_file,
            perturbed_root,
            ",".join([p[0] for p in emi_pairs]),
        )

        generated_files: list[Path] = []
        for mem in range(self.no_mems):
            file_dest_dir = self.paths.path_data / f"RUN_{mem}" / f"EMISSION_{mem}"
            file_dest_dir.mkdir(parents=True, exist_ok=True)
            output_file = file_dest_dir / file_original_name
            perturbed_files: list[tuple[str, Path]] = []
            for emi_var, emi_dir in emi_pairs:
                perturbed_file = (
                    perturbed_root
                    / f"emi_{mem}"
                    / emi_dir
                    / f"AEMISSIONS.{daily_stamp}_ITA7_{mem}.nc"
                )
                if not perturbed_file.exists():
                    raise FatalPipelineError(
                        f"Perturbed emission file not found for member {mem}, variable {emi_var}: {perturbed_file}"
                    )
                perturbed_files.append((emi_var, perturbed_file))

            shutil.copy2(base_file, output_file)

            with xr.open_dataset(base_file) as ds_base:
                for emi_var, perturbed_file in perturbed_files:

                    logger.info(
                        "[EMISSIONS] Member %s replacing %s into %s from %s",
                        mem,
                        emi_var,
                        output_file,
                        perturbed_file,
                    )

                    with xr.open_dataset(perturbed_file) as ds_perturbed:
                        if emi_var not in ds_base:
                            raise FatalPipelineError(
                                f"Variable {emi_var} not found in base emission file {base_file}"
                            )
                        if emi_var not in ds_perturbed:
                            raise FatalPipelineError(
                                f"Variable {emi_var} not found in perturbed emission file {perturbed_file}"
                            )

                        orig = ds_base[emi_var].values
                        pert = ds_perturbed[emi_var].values
                        if np.array_equal(orig, pert):
                            logger.warning(
                                "[EMISSIONS] Member %s has identical %s values in %s",
                                mem,
                                emi_var,
                                perturbed_file,
                            )

                    try:
                        subprocess.run(
                            [
                                "ncks",
                                "-A",
                                "-v",
                                emi_var,
                                str(perturbed_file),
                                str(output_file),
                            ],
                            check=True,
                        )
                    except subprocess.CalledProcessError as e:
                        raise FatalPipelineError(
                            f"ncks failed replacing {emi_var} from {perturbed_file} into {output_file}"
                        ) from e

            generated_files.append(output_file)

        self._generated_daily_emission_files = generated_files
        logger.info(
            "[EMISSIONS] Generated %d daily perturbed emission files for %s",
            len(generated_files),
            daily_stamp,
        )
        self._generated_daily_emission_stamp = daily_stamp

    def finalize_step(self):
        """
        Cleanup + YAML update.
        """
    
        modify_yaml_date(
            self.config._config_path,
            self.time_manager.simulated_time.strftime("%Y-%m-%d %H:00:00"),
        )


        if (
            self._generated_daily_emission_files
            and self.time_manager.current_time.hour == 0
        ):
            for file_path in self._generated_daily_emission_files:
                if file_path.exists():
                    file_path.unlink()
            self._generated_daily_emission_files = []
            self._generated_daily_emission_stamp = None


        for mem in range(self.no_mems):
            run_dir = self.paths.path_data / f"RUN_{mem}"

            ts = self.time_manager.current_time.strftime('%Y%m%d%H')
            tmp_dir = run_dir / f"tmp{ts}-{self.case_dir}"

            if tmp_dir.exists():
                if not tmp_dir.is_dir():
                    raise RuntimeError(f"[CLEANUP] Expected directory, got file: {tmp_dir}")

                logger.info(f"[CLEANUP] Removing tmp directory: {tmp_dir}")
                shutil.rmtree(tmp_dir)

        cleanup_cfg = getattr(self.config, "cleanup", None)
        if cleanup_cfg is not None and cleanup_cfg.enabled:
            if cleanup_cfg.delete_window_ibc or cleanup_cfg.delete_window_emissions or cleanup_cfg.delete_window_meteo:
                self._cleanup_window_inputs()
            self._cleanup_trim_outputs()

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
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

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
                    "@mimesi_nhours_list": str(self.current_window.run_hours),
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
            run_hours = self.current_window.run_hours
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

        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        target_time = self.current_window.end_time
        if not self.current_window.has_assimilation:
            logger.info(
                "[DART] No satellite data for %s, skipping after_model",
                target_time.strftime("%Y%m%d%H"),
            )
            self._clear_pending_orbit()
            return

        self._pending_orbit_filename = self.current_window.orbit_filename
        self._pending_orbit_start = self.current_window.obs_time
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
                self.current_window.run_hours,
            )

            if not out_file.exists():
                logger.error("No CHIMERE out.*.nc files found after run_model().")
                return

            logger.info(
                "Post-processing CHIMERE output (mem %s), out_file %s",
                mem,
                out_file
            )
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
            _run_cdo_suppress_5d(
                ["selname,pres,temp,psfc,lat,lon,Times", temp_out_psfc, tmp_pres]
            )

            t1 = self.current_window.start_time.strftime("%Y%m%d%H")
            tp = self.current_window.end_time.strftime("%Y%m%d%H")
            out_file_with_mem = to_dart_dir / f"out.{t1}_{tp}_{mem}.nc"
            _run_cdo_suppress_5d(
                [f"selname,{self.ass_var}", tmp_ts, out_file_with_mem]
            )

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
        if self.scheduler != Scheduler.SLURM:
            raise FatalPipelineError(
                f"Chimere2017 DART pipeline only supports SLURM, got {self.scheduler}"
            )
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        self.output_sim_folder = self.paths.dart_posteriors_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        Path(self.output_sim_folder).mkdir(parents=True, exist_ok=True)
        filter_cores = (
            self.config.monitoring.cores
            if self.config.monitoring is not None
            else 20
        )
        input_start_time = self.current_window.start_time
        input_end_time = self.current_window.end_time
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
            "SBATCH_PARTITION": self.cineca_queue,
            "DEST_LOG_PATH": self.output_sim_folder,
        }

        replace_nml_template(
            self.paths.base_path / self.paths.path_submit_bsh / "templates/run_mimesi-filter_template.bsh",
            entries_tbr_dict=submit_filter_entries,
            output_nml_path=self.paths.path_submit_bsh / "submit_filter.bsh",
        )
        
        submit_filter_script = self.paths.path_submit_bsh / "submit_filter.bsh"
        job_ids = submit_and_wait_cineca(
            CommandSpec(
                command="sbatch",
                directory=self.paths.path_submit_bsh,
                args=[str(submit_filter_script)],
            )
        )
        time.sleep(10)
        self.monitor_job_dart(job_ids)

    def monitor_job_dart(self, job_id):
        job_ids = job_id if isinstance(job_id, (list, tuple)) else [job_id]
        logger.info(f"Monitoring SLURM jobs {job_ids}")

        while True:
            active_jobs = {}
            for jid in job_ids:
                job_state = self._get_slurm_job_state(str(jid))
                if job_state not in {"COMPLETED", "NOT_FOUND"}:
                    active_jobs[str(jid)] = job_state

            if not active_jobs:
                print("Job completed successfully.")
                self.move_analysis_files()
                self.move_preassim_files()
                replace_priorinflation(
                    self.paths.base_path,
                    self.paths.path_filter,
                    self.time_manager.simulated_time.strftime("%Y%m%d%H"),
                )
                break

            status_summary = ", ".join(
                f"{jid}:{state}" for jid, state in active_jobs.items()
            )
            print(f"Jobs still active: {status_summary}. Waiting...")
            time.sleep(10)

    def _get_slurm_job_state(self, job_id: str) -> str:
        result = subprocess.run(
            ["squeue", "-j", job_id, "-h", "-o", "%T"],
            capture_output=True,
            text=True,
        )

        state = result.stdout.strip()
        if not state:
            return "COMPLETED"
        return state.splitlines()[0].strip().upper()

    def move_analysis_files(self):
        analysis_sim_folder = self.paths.dart_analysis_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        analysis_sim_folder.mkdir(parents=True, exist_ok=True)

        for filename in os.listdir(self.paths.path_filter):
            if not filename.startswith("analysis_"):
                continue

            src = Path(self.paths.path_filter) / filename
            dest = analysis_sim_folder / filename
            try:
                shutil.move(src, dest)
            except shutil.Error:
                print(
                    f"Failed to move '{filename}' to '{analysis_sim_folder}' because it already exists."
                )

    def move_preassim_files(self):
        preassim_sim_folder = self.paths.dart_preassim_dir(
            self.time_manager.simulated_time.strftime("%Y%m%d%H")
        )
        preassim_sim_folder.mkdir(parents=True, exist_ok=True)

        for filename in os.listdir(self.paths.path_filter):
            if not filename.startswith("preassim_"):
                continue

            src = Path(self.paths.path_filter) / filename
            dest = preassim_sim_folder / filename
            try:
                shutil.move(src, dest)
            except shutil.Error:
                print(
                    f"Failed to move '{filename}' to '{preassim_sim_folder}' because it already exists."
                )

    def after_assimilation(self):
        
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled, skipping after_assimilation")
            return
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        assimilation_time = self.time_manager.simulated_time
        if assimilation_time is None:
            logger.info("[DART] simulated_time is not set, skipping after_assimilation")
            return

        if not self.current_window.has_assimilation:
            logger.info(
                "[DART] Current window has no assimilation, skipping after_assimilation"
            )
            return

        prior_start_time = self.current_window.start_time
        posterior_dir = self.paths.dart_posteriors_dir(
            assimilation_time.strftime("%Y%m%d%H")
        )

        if not posterior_dir.exists():
            logger.info(
                "[DART] No posterior directory for %s, skipping after_assimilation",
                assimilation_time.strftime("%Y%m%d%H"),
            )
            return

        for mem in range(self.no_mems):
            posterior_file = (
                posterior_dir
                / f"out.posterior_{assimilation_time.strftime('%Y%m%d%H')}_{mem}.nc"
            )
            if not posterior_file.exists():
                logger.info(
                    "[DART] Posterior file missing for mem %s: %s. Skipping member.",
                    mem,
                    posterior_file,
                )
                continue
            
            restart_chimere_folder = self.paths.chimere_output_runs_dir(mem)
            restart_chimere_folder.mkdir(parents=True, exist_ok=True)
            restart_backup_dir = restart_chimere_folder / "restart_backup"
            restart_backup_dir.mkdir(parents=True, exist_ok=True)

            restart_from_chimere_file = self.paths.get_chimere_output_path(
                self.model_type,
                mem,
                prior_start_time,
                "end",
                self.current_window.run_hours,
            )
            if not restart_from_chimere_file.exists():
                raise FatalPipelineError(
                    "CHIMERE restart file not found for "
                    f"mem {mem}: {restart_from_chimere_file}"
                )
            

            backup_restart_file = restart_backup_dir / restart_from_chimere_file.name
            if backup_restart_file.exists():
                backup_restart_file.unlink()
            shutil.copy2(restart_from_chimere_file, backup_restart_file)
            result_tmp = restart_chimere_folder / f"{restart_from_chimere_file.stem}.tmp.nc"
            result_var_tmp = (
                restart_chimere_folder / f"{restart_from_chimere_file.stem}.{self.ass_var}.tmp.nc"
            )
            try:
                with xr.open_dataset(backup_restart_file) as ds_restart:
                    ds_restart = ds_restart.load()
                with xr.open_dataset(posterior_file) as ds_posterior:
                    ds_posterior = ds_posterior.load()

                if self.ass_var not in ds_restart:
                    raise FatalPipelineError(
                        f"Variable {self.ass_var} not found in CHIMERE restart file {backup_restart_file}"
                    )
                if self.ass_var not in ds_posterior:
                    raise FatalPipelineError(
                        f"Variable {self.ass_var} not found in DART posterior file {posterior_file}"
                    )

                prior_var = ds_restart[self.ass_var]
                posterior_var = ds_posterior[self.ass_var]

                time_dim = prior_var.dims[0]
                if time_dim != "Time":
                    raise FatalPipelineError(
                        f"Unexpected time dim for {self.ass_var}: {time_dim}"
                    )

                prior_spatial = prior_var.shape[1:]
                posterior_spatial = posterior_var.shape[1:]
                if prior_spatial != posterior_spatial:
                    raise FatalPipelineError(
                        f"Shape mismatch for {self.ass_var} (non-time dims): "
                        f"CHIMERE shape={prior_var.shape}, "
                        f"DART shape={posterior_var.shape}"
                    )

                if posterior_var.shape[0] != 1:
                    raise FatalPipelineError(
                        f"Unexpected DART time dimension for {self.ass_var}: "
                        f"{posterior_var.shape[0]} (expected 1)"
                    )

                if prior_var.shape[0] < 1:
                    raise FatalPipelineError(
                        f"Unexpected CHIMERE time dimension for {self.ass_var}: "
                        f"{prior_var.shape[0]}"
                    )

                for dim_name in prior_var.dims:
                    if dim_name in ds_restart.coords and dim_name in ds_posterior.coords:
                        if not ds_restart[dim_name].identical(ds_posterior[dim_name]):
                            raise FatalPipelineError(
                                f"Coordinate mismatch for {self.ass_var} on dim {dim_name}: "
                                f"CHIMERE coord differs from DART coord"
                            )
                #conversion
                # posterior is in ppb 
                # airm is molec/m3
                # from ppb to molec(species)/molec(air) is 1e-9
                # from molec/m3 to molec/cm3 is 1e-6
                # factor is 1e-15
                updated_values = (
                    posterior_var.values[0, :, :, :]
                    * ds_restart["airm"].values[-1, :, :, :]
                    * 1e-15
                )
                ds_update = ds_restart[[self.ass_var]].copy(deep=True)
                ds_update[self.ass_var].values[-1, :, :, :] = updated_values
                ds_update.to_netcdf(
                    result_var_tmp,
                    format="NETCDF3_64BIT",
                    engine="netcdf4",
                )
                shutil.copy2(backup_restart_file, result_tmp)
                subprocess.run(
                    [
                        "ncks",
                        "-A",
                        "-v",
                        self.ass_var,
                        str(result_var_tmp),
                        str(result_tmp),
                    ],
                    check=True,
                )
                os.replace(result_tmp, restart_from_chimere_file)
                logger.info(
                    "[DART] Updated CHIMERE prior for mem %s with posterior %s -> %s",
                    mem,
                    posterior_file.name,
                    restart_from_chimere_file.name,
                )
            finally:
                if result_var_tmp.exists():
                    result_var_tmp.unlink()

    def run_assimilation_if_needed(self):
        """
        Run DART assimilation only if:
        - assimilation is enabled
        - satellite observations are available
        """
        if not self.run_assimilation_flag:
            logger.info("[DART] Assimilation disabled by config")
            return
        if self.current_window is None:
            raise FatalPipelineError("Assimilation window is not initialized")

        orbit_filename = None
        if self.current_window.has_assimilation and self.current_window.orbit_filename:
            orbit_filename = self.current_window.orbit_filename
            self.time_manager.sat_obs = self.current_window.obs_time
        elif (
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
