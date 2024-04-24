import subprocess


def submit_slurm_job(case, option):
    job_name = f"DART_{case}"
    error_file = f"err_{case}_{option}.log"
    output_file = f"out_{case}_{option}.log"
    match option:
        case "filter":
            command_execute = f"mpirun -np 16 ./filter"

    slurm_script = f"""#!/bin/sh


#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=16
#SBATCH --error=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/errs/{error_file}
#SBATCH --output=/mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work/outs/{output_file}


conda activate /home/dausilio/miniconda3/envs/dartenv
cd /mnt/mumbai_n4r5/dausilio/projects/DART/models/FARM/work
{command_execute}
"""

    slurm_script_file = f"submit_job_{case}.sh"

    with open(slurm_script_file, "w") as file:
        file.write(slurm_script)

    subprocess.run(["sbatch", slurm_script_file])

    print(f"Job submitted for {case}")


submit_slurm_job("DART_rf_ens", "filter")
