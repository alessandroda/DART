#!/bin/bash
#MSUB -r CHIM-DART             # Job name
#MSUB -n 1                      # Number of tasks to use
#MSUB -c 8                       # Number of cores (or threads) per task to use
#MSUB -T 86400                    # Elapsed time limit in seconds of the job (default: 7200)
#MSUB -o job_mpi_%I.o         # Standard output. %I is the job id
#MSUB -e job_mpi_%I.e         # Error output. %I is the job id
#MSUB -A gen7232                 # Project ID
#MSUB -q rome                    # Partition name (see ccc_mpinfo)
#MSUB -m work,scratch,store

module unload intel/20.0.0
module load intel/20.0.4 mpi/openmpi/4.0.5 netcdf-fortran/4.5.3 nco
python -u orchestrator_irene.py -c config_pydantic_irene.yaml #without ccc_mprun
# submit the job with ccc_msub

