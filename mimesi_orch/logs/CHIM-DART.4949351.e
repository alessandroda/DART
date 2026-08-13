+ SCRIPT_PID=2456946
+ /bin/bash -x /tmp/tmp.gK2MQZFK3S
+ set +x
+ unset _mlshdbg
+ '[' 1 = 1 ']'
+ case "$-" in
+ set +x
+ unset _mlshdbg
+ module unload intel/20.0.0
+ local _mlredir=0
+ '[' -n '' ']'
+ case " $@ " in
+ '[' 0 -eq 0 ']'
+ _module_raw unload intel/20.0.0
+ unset _mlshdbg
+ '[' 1 = 1 ']'
+ case "$-" in
+ set +x
unload module proj/8.0.1
unload module ncview/2.1.7
unload module nco/5.0.1
unload module eccodes/2.23.0
unload module netcdf-c/4.7.4
unload module libaec/1.1.2
unload module jasper/2.0.14
unload module hdf5/1.12.0
unload module szip/2.1
unload module intel/20.0.0
unload module udunits/2.2.28
unload module mkl/20.0.0
unload module feature/mkl/vector/amd
unload module feature/mkl/single_node
unload module feature/mkl/sequential
unload module feature/mkl/lp64
unload module fortran/intel/20.0.0
unload module c/intel/20.0.0
unload module c++/intel/20.0.0
unload module licsrv/intel
unload module flavor/buildcompiler/intel/20
unload module flavor/hdf5/serial
+ unset _mlshdbg
+ return 0
+ module load intel/20.0.4 mpi/openmpi/4.0.5 netcdf-fortran/4.5.3 nco
+ local _mlredir=0
+ '[' -n '' ']'
+ case " $@ " in
+ '[' 0 -eq 0 ']'
+ _module_raw load intel/20.0.4 mpi/openmpi/4.0.5 netcdf-fortran/4.5.3 nco
+ unset _mlshdbg
+ '[' 1 = 1 ']'
+ case "$-" in
+ set +x
load module flavor/buildcompiler/intel/20
load module licsrv/intel
load module c++/intel/20.0.4
load module c/intel/20.0.4
load module fortran/intel/20.0.4
load module feature/mkl/lp64
load module feature/mkl/sequential
load module feature/mkl/single_node
load module feature/mkl/vector/amd
load module mkl/20.0.4
load module intel/20.0.4
load module flavor/buildmpi/openmpi/4
load module feature/openmpi/mpi_compiler/intel
load module feature/openmpi/io/standard
load module feature/openmpi/net/auto
load module flavor/libccc_user/hwloc2
load module hwloc/2.2.0
load module feature/system/pmix_v4
load module pmix/3.1.5
load module flavor/ucx/cuda-11.6
load module flavor/cuda/nvhpc-222
load module cuda/11.6
load module ucx/1.18.1
load module mpi/openmpi/4.0.5
load module flavor/hdf5/serial
load module szip/2.1
load module hdf5/1.12.0
load module netcdf-c/4.7.4
load module netcdf-fortran/4.5.3
load module jasper/2.0.14
load module libaec/1.1.2
load module eccodes/2.23.0
load module nco/5.0.1
+ unset _mlshdbg
+ return 0
+ python -u main.py -c config/config_irene_IM.yaml
2026-06-22 16:48:55 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-22 16:48:55 INFO [PIPELINE] =======================================
2026-06-22 16:48:55 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-22 16:48:55 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-22 16:48:55 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-22 16:48:55 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260622_164855.log
2026-06-22 16:48:55 INFO [PIPELINE] =======================================
2026-06-22 16:48:55 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-22 16:48:55 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-22 16:48:55 INFO [STEP] ---- TIME LOOP START ----
2026-06-22 16:48:55 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-22 16:48:55 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-06-22 16:48:55 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-22 16:48:55 INFO >> Checking links...
2026-06-22 16:49:05 INFO >> All links are good for ENS1  ...
2026-06-22 16:49:05 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-22 16:49:05 INFO >> Checking links...
2026-06-22 16:49:15 INFO >> All links are good for ENS2  ...
2026-06-22 16:49:15 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-22 16:49:15 INFO >> Checking links...
2026-06-22 16:49:25 INFO >> All links are good for ENS3  ...
2026-06-22 16:49:25 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-22 16:49:25 INFO >> Checking links...
2026-06-22 16:49:35 INFO >> All links are good for ENS4  ...
2026-06-22 16:49:35 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-22 16:49:35 INFO >> Checking links...
2026-06-22 16:49:45 INFO >> All links are good for ENS5  ...
2026-06-22 16:49:45 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-22 16:49:45 INFO >> Checking links...
2026-06-22 16:49:55 INFO >> All links are good for ENS6  ...
2026-06-22 16:49:55 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-22 16:49:55 INFO >> Checking links...
2026-06-22 16:50:05 INFO >> All links are good for ENS7  ...
2026-06-22 16:50:05 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-22 16:50:05 INFO >> Checking links...
2026-06-22 16:50:15 INFO >> All links are good for ENS8  ...
2026-06-22 16:50:15 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-22 16:50:15 INFO >> Checking links...
2026-06-22 16:50:26 INFO >> All links are good for ENS9  ...
2026-06-22 16:50:26 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-22 16:50:26 INFO >> Checking links...
2026-06-22 16:50:36 INFO >> All links are good for ENS10  ...
2026-06-22 16:50:36 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-22 16:50:36 INFO >> Checking links...
[2026-06-22T16:50:43.346] error: *** JOB 4949351 ON irene4050 CANCELLED AT 2026-06-22T16:50:43 DUE to SIGNAL Terminated ***
