#!/bin/bash
#MSUB -r CHIM-DART             # Job name
#MSUB -n 1                      # Number of tasks to use
#MSUB -c 8                       # Number of cores (or threads) per task to use
#MSUB -T 86400                    # Elapsed time limit in seconds of the job (default: 7200)
#MSUB -e $CCCSCRATCHDIR/OUT_orch_chimdart/mimesi_orchestrator_logs/job_mpi_%I.e         # Error output. %I is the job id
#MSUB -A gen7232                 # Project ID
#MSUB -q rome                    # Partition name (see ccc_mpinfo)
#MSUB -m work,scratch,store

module unload intel/20.0.0
module load intel/20.0.4 mpi/openmpi/4.0.5 netcdf-fortran/4.5.3 nco
python -u main.py -c config/config_irene_IM.yaml #config_pydantic_irene.yaml #without ccc_mprun
# submit the job with ccc_msub

