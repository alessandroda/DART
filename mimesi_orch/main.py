import argparse
import time
import logging
import yaml

from pathlib import Path
from orchestrator_utils import TimeManager
from config_models import AppConfig
from pipelines.factory import build_pipeline


def _resolve_config_path(base_path: Path, value: Path) -> Path:
    if value.is_absolute():
        return value
    return base_path / value


def _resolve_log_dir(config: AppConfig) -> Path:
    path_data = _resolve_config_path(config.paths.base_path, config.paths.path_data)
    runtime = config.runtime
    if runtime is None:
        return path_data / "mimesi_orchestrator_logs"

    log_dir_name = runtime.log_dir_name
    if not log_dir_name:
        if runtime.profile == "irene":
            log_dir_name = "irene_orchestrator_logs"
        else:
            log_dir_name = "mimesi_orchestrator_logs"

    if runtime.log_dir_mode == "cwd":
        return Path.cwd() / log_dir_name
    return path_data / log_dir_name


parser = argparse.ArgumentParser(description="Python orchestrator")
parser.add_argument(
    "-c", "--conf", type=str, required=True, help="Path to the YAML config file"
)
args = parser.parse_args()

CONFIG_PATH = args.conf

with open(CONFIG_PATH, "r") as f:
    cfg = yaml.safe_load(f)

config = AppConfig.model_validate(cfg)
config._config_path = Path(CONFIG_PATH)

LOG_DIR = _resolve_log_dir(config)
LOG_DIR.mkdir(parents=True, exist_ok=True)

logfile = LOG_DIR / f"{config.assimilation.model_type.value}_DART_{time.strftime('%Y%m%d_%H%M%S')}.log"

log_datefmt = "%Y-%m-%d %H:%M:%S"
file_formatter = logging.Formatter(
    "%(asctime)s %(levelname)s %(name)s [%(filename)s:%(lineno)d] %(message)s",
    datefmt=log_datefmt,
)
console_formatter = logging.Formatter(
    "%(asctime)s %(levelname)s %(message)s",
    datefmt=log_datefmt,
)

file_handler = logging.FileHandler(logfile)
file_handler.setFormatter(file_formatter)

console_handler = logging.StreamHandler()
console_handler.setFormatter(console_formatter)

logging.basicConfig(
    level=logging.INFO,
    handlers=[file_handler, console_handler],
)

logger = logging.getLogger(__name__)
logger.info(
    "\n"
    "███    ███ ██ ███    ███ ███████ ███████ ██\n"
    "████  ████ ██ ████  ████ ██      ██      ██\n"
    "██ ████ ██ ██ ██ ████ ██ █████   ███████ ██\n"
    "██  ██  ██ ██ ██  ██  ██ ██           ██ ██\n"
    "██      ██ ██ ██      ██ ███████ ███████ ██"
    "\n"
    "\n"
    "\n"
)
logger.info("[PIPELINE] =======================================")
logger.info(f"[PIPELINE] Starting {config.assimilation.model_type.value}–DART orchestrator")
logger.info(f"[PIPELINE] Config file: {CONFIG_PATH}")
logger.info(f"[PIPELINE] Run dir: {config.paths.path_data}")
logger.info(f"[PIPELINE] Log file: {logfile}")
logger.info("[PIPELINE] =======================================")


time_manager = TimeManager(
    start_time=config.time.start_time,
    end_time=config.time.end_time,
    dt_seconds=config.time.dt_seconds,
)

pipeline = build_pipeline(config, time_manager)
pipeline.run_pipeline()
