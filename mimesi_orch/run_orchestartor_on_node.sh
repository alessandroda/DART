#!/bin/bash
#MSUB -r orchestrator              # Job name
#MSUB -n 1                      # Number of tasks to use
#MSUB -c 64                      # Number of cores (or threads) per task to use
#MSUB -T 86400                    # Elapsed time limit in seconds of the job (default: 7200)
#MSUB -o err_out/job_mpi_%I.o         # Standard output. %I is the job id
#MSUB -e err_out/job_mpi_%I.e         # Error output. %I is the job id
#MSUB -A gen7232                 # Project ID
#MSUB -q rome                    # Partition name (see ccc_mpinfo)
#MSUB -m work,scratch,store

source /ccc/cont003/home/lisa/demoling/CHIMDART_env/bin/activate
ccc_mprun python -u orchestrator_irene.py -c config_pydantic_irene.yaml

# execute: ccc_msub run_orchestartor_on_node.sh
