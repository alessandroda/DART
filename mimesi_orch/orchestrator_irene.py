import argparse
import time
import logging
import yaml

from pathlib import Path
from orchestrator_utils import TimeManager
from pipelines.chimere_irene_pipeline import ChimereV2023DartPipeline
from config_models import AppConfig
from paths import PathManager


parser = argparse.ArgumentParser(description="Python orchestrator for CHIMERE-DART")
parser.add_argument(
    "-c", "--conf", type=str, required=True, help="Path to the YAML config file"
)
args = parser.parse_args()

CONFIG_PATH = args.conf

with open(CONFIG_PATH, "r") as f:
    cfg = yaml.safe_load(f)

config = AppConfig.model_validate(cfg)
path_manager = PathManager(config.paths)

LOG_DIR = Path(config.paths.run_dir) / "irene_orchestrator_logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

logfile = LOG_DIR / f"{config.assimilation.model_type.value}_DART_{time.strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(processName)s/%(threadName)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(logfile),
        logging.StreamHandler(),
    ],
)

logger = logging.getLogger(__name__)

logger.info(f"Starting {config.assimilation.model_type.value}–DART orchestrator")
logger.info(f"Config file: {CONFIG_PATH}")
logger.info(f"Run dir: {config.paths.run_dir}")
logger.info(f"Log file: {logfile}")

time_manager = TimeManager(
    start_time=config.time.start_time,
    end_time=config.time.end_time,
    dt_seconds=config.time.dt_seconds,
)
path_manager = PathManager(config.paths)
pipeline = ChimereV2023DartPipeline(time_manager, path_manager, config)
pipeline.run_pipeline()
