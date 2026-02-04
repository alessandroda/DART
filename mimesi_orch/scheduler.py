import logging
import os
import re
from pathlib import Path
from mimesi_types import Scheduler
import subprocess
import time

logger = logging.getLogger(__name__)

def submit_sbatch(command: Path, workdir: Path) -> list[str]:
    original_directory = os.getcwd()
    try:
        os.chdir(workdir)

        subprocess.run(["chmod", "+x", command], check=True)
        result = subprocess.run(
            ["sbatch", str(command)],
            capture_output=True,
            text=True,
            check=True,
        )

        # Typical output:
        # "Submitted batch job 123456"
        job_id = result.stdout.strip().split()[-1]
        return [job_id]

    finally:
        os.chdir(original_directory)


def submit_bsub(command: Path, workdir: Path) -> list[str]:
    original_directory = os.getcwd()
    try:
        os.chdir(workdir)

        subprocess.run(["chmod", "+x", command], check=True)

        result = subprocess.run(
            [str(command)],
            capture_output=True,
            text=True,
            check=True,
        )

        # Example:
        # "Job <987654> is submitted to queue <normal>"
        job_ids = re.findall(r"Job <(\d+)>", result.stdout)
        return job_ids

    finally:
        os.chdir(original_directory)


def submit_job(
    scheduler: Scheduler,
    command: Path,
    workdir: Path,
) -> list[str]:

    if scheduler == Scheduler.SLURM:
        return submit_sbatch(command, workdir)

    elif scheduler == Scheduler.LSF:
        return submit_bsub(command, workdir)

    else:
        raise ValueError(f"Unsupported scheduler {scheduler}")


def check_job_status_slurm(job_id: str) -> bool:
    result = subprocess.run(
        ["squeue", "-j", job_id, "-h"],
        capture_output=True,
        text=True,
    )

    # If squeue prints something → job is still in queue
    return bool(result.stdout.strip())


def check_job_status_lsf(job_id: str) -> str:
    result = subprocess.run(
        ["bjobs", "-noheader", job_id],
        capture_output=True,
        text=True,
    )

    if "DONE" in result.stdout:
        return "COMPLETED"
    if "EXIT" in result.stdout:
        return "FAILED"
    if result.stdout.strip() == "":
        return "COMPLETED"

    return "RUNNING"


def check_job_status(scheduler: Scheduler, job_id: str) -> str:
    if scheduler == Scheduler.SLURM:
        return check_job_status_slurm(job_id)
    elif scheduler == Scheduler.LSF:
        return check_job_status_lsf(job_id)
    else:
        raise ValueError("Unknown scheduler")


def slurm_job_is_running(job_id: str) -> bool:
    result = subprocess.run(
        ["squeue", "-j", job_id, "-h"],
        capture_output=True,
        text=True,
    )

    # If squeue prints something → job is still in queue
    return bool(result.stdout.strip())


def wait_for_slurm_jobs(job_ids: list[str], poll_interval: int = 30):
    while True:
        running = [job_id for job_id in job_ids if slurm_job_is_running(job_id)]

        if not running:
            logger.info("All SLURM jobs have finished")
            return

        logger.info(f"Jobs still running: {running}")
        time.sleep(poll_interval)
