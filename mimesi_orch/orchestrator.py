import argparse
import time
import logging
import yaml
from mimesi_orch.orchestrator_utils import PathManager, TimeManager
from mimesi_orch.pipelines.farm_pipeline import FarmToDartPipeline

logging.basicConfig(
    filename=f'logs_orchestrator/farm_to_dart_{time.strftime("%Y%m%d_%H%M%S")}.log',
    format="%(asctime)s [%(processName)s/%(threadName)s] %(levelname)s: %(message)s",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Python orchestrator for FARM-DART")
parser.add_argument(
    "-c", "--conf", type=str, required=True, help="Path to the YAML config file"
)
args = parser.parse_args()

CONFIG_PATH = args.conf


def load_config(file_path):
    print(f"Using config file: {file_path}")
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


config = load_config(CONFIG_PATH)

time_manager = TimeManager(
    start_time=config["time"]["start_time"],
    end_time=config["time"]["end_time"],
    dt_seconds=config["time"]["dt_seconds"],
)

path_manager = PathManager(
    base_path=config["paths"]["base_path"],
    env_python=config["paths"]["env_python"],
    listing_file=config["paths"]["listing_file"],
    run_submit_farm_template=config["paths"]["run_submit_farm_template"],
    path_submit_bsh=config["paths"]["path_submit_bsh"],
    path_filter=config["paths"]["path_filter"],
    path_data=config["paths"]["path_data"],
    run_submit_replace_perturbations=config["paths"][
        "run_submit_replace_perturbations"
    ],
    log_paths=True,
)


pipeline = FarmToDartPipeline(time_manager, path_manager, config)
pipeline.run_pipeline()
