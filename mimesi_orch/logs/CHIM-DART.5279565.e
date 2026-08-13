+ SCRIPT_PID=1728247
+ /bin/bash -x /tmp/tmp.ORZy2STecd
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
+ python -u main.py -c config/config_irene_IM_cp2.yaml
2026-07-26 12:25:31 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-26 12:25:31 INFO [PIPELINE] =======================================
2026-07-26 12:25:31 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-26 12:25:31 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-26 12:25:31 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-26 12:25:31 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260726_122531.log
2026-07-26 12:25:31 INFO [PIPELINE] =======================================
2026-07-26 12:25:31 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-26 12:25:31 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-26 12:25:31 INFO [STEP] ---- TIME LOOP START ----
2026-07-26 12:25:31 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:25:31 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-26 12:25:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:39 INFO Hourly dataset computed and listing created
2026-07-26 12:25:45 INFO Hourly dataset computed
2026-07-26 12:25:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:46 INFO Hourly dataset computed and listing created
2026-07-26 12:25:48 INFO Hourly dataset computed
2026-07-26 12:25:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:49 INFO Hourly dataset computed and listing created
2026-07-26 12:25:49 INFO Hourly dataset computed
2026-07-26 12:25:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:50 INFO Hourly dataset computed and listing created
2026-07-26 12:25:51 INFO Hourly dataset computed
2026-07-26 12:25:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:52 INFO Hourly dataset computed and listing created
2026-07-26 12:25:53 INFO Hourly dataset computed
2026-07-26 12:25:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:54 INFO Hourly dataset computed and listing created
2026-07-26 12:25:55 INFO Hourly dataset computed
2026-07-26 12:25:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:56 INFO Hourly dataset computed and listing created
2026-07-26 12:25:57 INFO Hourly dataset computed
2026-07-26 12:25:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:25:58 INFO Hourly dataset computed and listing created
2026-07-26 12:25:59 INFO Hourly dataset computed
2026-07-26 12:25:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:00 INFO Hourly dataset computed and listing created
2026-07-26 12:26:01 INFO Hourly dataset computed
2026-07-26 12:26:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:02 INFO Hourly dataset computed and listing created
2026-07-26 12:26:03 INFO Hourly dataset computed
2026-07-26 12:26:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:04 INFO Hourly dataset computed and listing created
2026-07-26 12:26:04 INFO Hourly dataset computed
2026-07-26 12:26:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:05 INFO Hourly dataset computed and listing created
2026-07-26 12:26:06 INFO Hourly dataset computed
2026-07-26 12:26:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:07 INFO Hourly dataset computed and listing created
2026-07-26 12:26:08 INFO Hourly dataset computed
2026-07-26 12:26:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:09 INFO Hourly dataset computed and listing created
2026-07-26 12:26:10 INFO Hourly dataset computed
2026-07-26 12:26:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:26:11 INFO Hourly dataset computed and listing created
2026-07-26 12:26:11 INFO Hourly dataset computed
2026-07-26 12:26:11 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-26 12:26:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:26:11 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-26 12:26:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:26:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:11 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:26:12 INFO Queuing job for member 1...
2026-07-26 12:26:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:12 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:26:12 INFO Found: ['5279579']
2026-07-26 12:26:17 INFO [TGCC-IRENE] Submitted job with ID:['5279579']
2026-07-26 12:26:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:26:17 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-26 12:26:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:26:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:17 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:26:17 INFO Queuing job for member 2...
2026-07-26 12:26:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:17 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:26:18 INFO Found: ['5279582']
2026-07-26 12:26:23 INFO [TGCC-IRENE] Submitted job with ID:['5279582']
2026-07-26 12:26:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:26:23 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-26 12:26:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:26:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:23 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:26:23 INFO Queuing job for member 3...
2026-07-26 12:26:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:23 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:26:24 INFO Found: ['5279583']
2026-07-26 12:26:29 INFO [TGCC-IRENE] Submitted job with ID:['5279583']
2026-07-26 12:26:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:26:29 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-26 12:26:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:26:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:29 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:26:29 INFO Queuing job for member 4...
2026-07-26 12:26:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:29 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:26:30 INFO Found: ['5279586']
2026-07-26 12:26:35 INFO [TGCC-IRENE] Submitted job with ID:['5279586']
2026-07-26 12:26:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:26:35 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-26 12:26:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:26:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:35 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:26:35 INFO Queuing job for member 5...
2026-07-26 12:26:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:35 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:26:35 INFO Found: ['5279587']
2026-07-26 12:26:40 INFO [TGCC-IRENE] Submitted job with ID:['5279587']
2026-07-26 12:26:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:26:40 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-26 12:26:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:26:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:40 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:26:40 INFO Queuing job for member 6...
2026-07-26 12:26:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:40 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:26:41 INFO Found: ['5279588']
2026-07-26 12:26:46 INFO [TGCC-IRENE] Submitted job with ID:['5279588']
2026-07-26 12:26:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:26:46 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-26 12:26:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:26:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:46 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:26:46 INFO Queuing job for member 7...
2026-07-26 12:26:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:46 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:26:47 INFO Found: ['5279589']
2026-07-26 12:26:52 INFO [TGCC-IRENE] Submitted job with ID:['5279589']
2026-07-26 12:26:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 12:26:52 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-26 12:26:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 12:26:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:52 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 12:26:52 INFO Queuing job for member 8...
2026-07-26 12:26:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:52 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 12:26:54 INFO Found: ['5279590']
2026-07-26 12:26:59 INFO [TGCC-IRENE] Submitted job with ID:['5279590']
2026-07-26 12:26:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:26:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 12:26:59 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-26 12:26:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 12:26:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:26:59 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 12:26:59 INFO Queuing job for member 9...
2026-07-26 12:26:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:26:59 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 12:27:02 INFO Found: ['5279591']
2026-07-26 12:27:07 INFO [TGCC-IRENE] Submitted job with ID:['5279591']
2026-07-26 12:27:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 12:27:07 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-26 12:27:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 12:27:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:07 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 12:27:07 INFO Queuing job for member 10...
2026-07-26 12:27:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:07 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 12:27:09 INFO Found: ['5279594']
2026-07-26 12:27:14 INFO [TGCC-IRENE] Submitted job with ID:['5279594']
2026-07-26 12:27:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 12:27:14 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-26 12:27:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 12:27:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:14 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 12:27:14 INFO Queuing job for member 11...
2026-07-26 12:27:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:14 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 12:27:17 INFO Found: ['5279595']
2026-07-26 12:27:22 INFO [TGCC-IRENE] Submitted job with ID:['5279595']
2026-07-26 12:27:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 12:27:22 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-26 12:27:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 12:27:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:22 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 12:27:22 INFO Queuing job for member 12...
2026-07-26 12:27:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:22 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 12:27:24 INFO Found: ['5279596']
2026-07-26 12:27:29 INFO [TGCC-IRENE] Submitted job with ID:['5279596']
2026-07-26 12:27:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 12:27:29 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-26 12:27:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 12:27:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:29 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 12:27:29 INFO Queuing job for member 13...
2026-07-26 12:27:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:29 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 12:27:32 INFO Found: ['5279597']
2026-07-26 12:27:37 INFO [TGCC-IRENE] Submitted job with ID:['5279597']
2026-07-26 12:27:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 12:27:37 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-26 12:27:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 12:27:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:37 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 12:27:37 INFO Queuing job for member 14...
2026-07-26 12:27:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:37 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 12:27:40 INFO Found: ['5279598']
2026-07-26 12:27:45 INFO [TGCC-IRENE] Submitted job with ID:['5279598']
2026-07-26 12:27:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:27:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 12:27:45 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-26 12:27:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 12:27:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:27:45 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 12:27:45 INFO Queuing job for member 15...
2026-07-26 12:27:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:27:45 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 12:27:47 INFO Found: ['5279599']
2026-07-26 12:27:52 INFO [TGCC-IRENE] Submitted job with ID:['5279599']
2026-07-26 12:27:52 INFO Checking job status ...
2026-07-26 12:27:52 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:27:52 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:27:52 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:28:07 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:28:07 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:28:08 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:28:10 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:28:10 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:28:25 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:28:25 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:28:25 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:28:40 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:28:40 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:28:40 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:28:55 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:28:55 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:28:56 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:28:56 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:28:56 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:29:11 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:29:11 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:29:11 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:29:26 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:29:26 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:29:27 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:29:27 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:29:27 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:29:27 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:29:29 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:29:29 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:29:44 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:29:44 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:29:44 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:29:59 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279588: status RUNNING/PENDING
2026-07-26 12:29:59 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:30:02 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:30:02 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:30:17 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279588: status FINISHED
2026-07-26 12:30:17 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:30:17 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:30:17 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:30:32 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279588: status FINISHED
2026-07-26 12:30:32 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:30:32 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:30:32 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:30:47 INFO None 5279579: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279582: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279583: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279588: status FINISHED
2026-07-26 12:30:47 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:30:47 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:30:48 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:30:48 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:30:48 INFO Jobs still running: ['5279579', '5279582', '5279583', '5279586', '5279587', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:31:03 INFO None 5279579: status FINISHED
2026-07-26 12:31:03 INFO None 5279582: status FINISHED
2026-07-26 12:31:03 INFO None 5279583: status FINISHED
2026-07-26 12:31:03 INFO None 5279586: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279587: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279588: status FINISHED
2026-07-26 12:31:03 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:31:03 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:31:03 INFO Jobs still running: ['5279586', '5279587', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:31:18 INFO None 5279579: status FINISHED
2026-07-26 12:31:18 INFO None 5279582: status FINISHED
2026-07-26 12:31:18 INFO None 5279583: status FINISHED
2026-07-26 12:31:18 INFO None 5279586: status FINISHED
2026-07-26 12:31:18 INFO None 5279587: status FINISHED
2026-07-26 12:31:18 INFO None 5279588: status FINISHED
2026-07-26 12:31:18 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:31:18 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:31:18 INFO Jobs still running: ['5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:31:33 INFO None 5279579: status FINISHED
2026-07-26 12:31:33 INFO None 5279582: status FINISHED
2026-07-26 12:31:33 INFO None 5279583: status FINISHED
2026-07-26 12:31:33 INFO None 5279586: status FINISHED
2026-07-26 12:31:33 INFO None 5279587: status FINISHED
2026-07-26 12:31:35 INFO None 5279588: status FINISHED
2026-07-26 12:31:35 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:31:35 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:31:35 INFO Jobs still running: ['5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:31:50 INFO None 5279579: status FINISHED
2026-07-26 12:31:50 INFO None 5279582: status FINISHED
2026-07-26 12:31:50 INFO None 5279583: status FINISHED
2026-07-26 12:31:50 INFO None 5279586: status FINISHED
2026-07-26 12:31:50 INFO None 5279587: status FINISHED
2026-07-26 12:31:51 INFO None 5279588: status FINISHED
2026-07-26 12:31:51 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:31:51 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:31:51 INFO Jobs still running: ['5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:32:06 INFO None 5279579: status FINISHED
2026-07-26 12:32:06 INFO None 5279582: status FINISHED
2026-07-26 12:32:06 INFO None 5279583: status FINISHED
2026-07-26 12:32:06 INFO None 5279586: status FINISHED
2026-07-26 12:32:06 INFO None 5279587: status FINISHED
2026-07-26 12:32:06 INFO None 5279588: status FINISHED
2026-07-26 12:32:08 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:32:08 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:32:08 INFO Jobs still running: ['5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:32:23 INFO None 5279579: status FINISHED
2026-07-26 12:32:23 INFO None 5279582: status FINISHED
2026-07-26 12:32:23 INFO None 5279583: status FINISHED
2026-07-26 12:32:23 INFO None 5279586: status FINISHED
2026-07-26 12:32:23 INFO None 5279587: status FINISHED
2026-07-26 12:32:23 INFO None 5279588: status FINISHED
2026-07-26 12:32:23 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279591: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:32:23 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:32:23 INFO Jobs still running: ['5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:32:38 INFO None 5279579: status FINISHED
2026-07-26 12:32:38 INFO None 5279582: status FINISHED
2026-07-26 12:32:38 INFO None 5279583: status FINISHED
2026-07-26 12:32:38 INFO None 5279586: status FINISHED
2026-07-26 12:32:38 INFO None 5279587: status FINISHED
2026-07-26 12:32:38 INFO None 5279588: status FINISHED
2026-07-26 12:32:38 INFO None 5279589: status RUNNING/PENDING
2026-07-26 12:32:38 INFO None 5279590: status RUNNING/PENDING
2026-07-26 12:32:38 INFO None 5279591: status FINISHED
2026-07-26 12:32:38 INFO None 5279594: status RUNNING/PENDING
2026-07-26 12:32:38 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:32:39 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:32:39 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:32:39 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:32:39 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:32:39 INFO Jobs still running: ['5279589', '5279590', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:32:54 INFO None 5279579: status FINISHED
2026-07-26 12:32:54 INFO None 5279582: status FINISHED
2026-07-26 12:32:54 INFO None 5279583: status FINISHED
2026-07-26 12:32:54 INFO None 5279586: status FINISHED
2026-07-26 12:32:54 INFO None 5279587: status FINISHED
2026-07-26 12:32:54 INFO None 5279588: status FINISHED
2026-07-26 12:32:54 INFO None 5279589: status FINISHED
2026-07-26 12:32:54 INFO None 5279590: status FINISHED
2026-07-26 12:32:54 INFO None 5279591: status FINISHED
2026-07-26 12:32:54 INFO None 5279594: status FINISHED
2026-07-26 12:32:54 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:32:54 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:32:54 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:32:54 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:32:54 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:32:54 INFO Jobs still running: ['5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:33:09 INFO None 5279579: status FINISHED
2026-07-26 12:33:09 INFO None 5279582: status FINISHED
2026-07-26 12:33:09 INFO None 5279583: status FINISHED
2026-07-26 12:33:09 INFO None 5279586: status FINISHED
2026-07-26 12:33:09 INFO None 5279587: status FINISHED
2026-07-26 12:33:09 INFO None 5279588: status FINISHED
2026-07-26 12:33:09 INFO None 5279589: status FINISHED
2026-07-26 12:33:09 INFO None 5279590: status FINISHED
2026-07-26 12:33:09 INFO None 5279591: status FINISHED
2026-07-26 12:33:09 INFO None 5279594: status FINISHED
2026-07-26 12:33:09 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:33:09 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:33:09 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:33:09 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:33:09 INFO None 5279599: status RUNNING/PENDING
2026-07-26 12:33:09 INFO Jobs still running: ['5279595', '5279596', '5279597', '5279598', '5279599']. Waiting...
2026-07-26 12:33:25 INFO None 5279579: status FINISHED
2026-07-26 12:33:25 INFO None 5279582: status FINISHED
2026-07-26 12:33:25 INFO None 5279583: status FINISHED
2026-07-26 12:33:25 INFO None 5279586: status FINISHED
2026-07-26 12:33:25 INFO None 5279587: status FINISHED
2026-07-26 12:33:25 INFO None 5279588: status FINISHED
2026-07-26 12:33:25 INFO None 5279589: status FINISHED
2026-07-26 12:33:25 INFO None 5279590: status FINISHED
2026-07-26 12:33:25 INFO None 5279591: status FINISHED
2026-07-26 12:33:25 INFO None 5279594: status FINISHED
2026-07-26 12:33:25 INFO None 5279595: status RUNNING/PENDING
2026-07-26 12:33:25 INFO None 5279596: status RUNNING/PENDING
2026-07-26 12:33:25 INFO None 5279597: status RUNNING/PENDING
2026-07-26 12:33:25 INFO None 5279598: status RUNNING/PENDING
2026-07-26 12:33:25 INFO None 5279599: status FINISHED
2026-07-26 12:33:25 INFO Jobs still running: ['5279595', '5279596', '5279597', '5279598']. Waiting...
2026-07-26 12:33:40 INFO None 5279579: status FINISHED
2026-07-26 12:33:40 INFO None 5279582: status FINISHED
2026-07-26 12:33:40 INFO None 5279583: status FINISHED
2026-07-26 12:33:40 INFO None 5279586: status FINISHED
2026-07-26 12:33:40 INFO None 5279587: status FINISHED
2026-07-26 12:33:40 INFO None 5279588: status FINISHED
2026-07-26 12:33:40 INFO None 5279589: status FINISHED
2026-07-26 12:33:40 INFO None 5279590: status FINISHED
2026-07-26 12:33:40 INFO None 5279591: status FINISHED
2026-07-26 12:33:40 INFO None 5279594: status FINISHED
2026-07-26 12:33:40 INFO None 5279595: status FINISHED
2026-07-26 12:33:40 INFO None 5279596: status FINISHED
2026-07-26 12:33:40 INFO None 5279597: status FINISHED
2026-07-26 12:33:40 INFO None 5279598: status FINISHED
2026-07-26 12:33:40 INFO None 5279599: status FINISHED
2026-07-26 12:33:40 INFO Jobs ['5279579', '5279582', '5279583', '5279586', '5279587', '5279588', '5279589', '5279590', '5279591', '5279594', '5279595', '5279596', '5279597', '5279598', '5279599'] have finished
2026-07-26 12:33:40 INFO Checking restart files were created ...
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-26 12:33:40 INFO  Run_model() completed successfully.
2026-07-26 12:33:40 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:33:40 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-26 12:33:40 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 12:33:40 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-26 12:33:40 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:33:40 INFO ---------->>> Running process_satellite_data()
2026-07-26 12:33:40 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-26 12:33:40 INFO ---------->>> Running run_obs_converter()
2026-07-26 12:33:40 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-26 12:33:40 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-26 12:33:40 INFO ---------->>> Running DART
2026-07-26 12:33:40 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 12:33:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-26 12:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-26 12:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-26 12:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-26 12:33:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-26 12:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-26 12:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-26 12:33:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-26 12:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-26 12:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-26 12:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-26 12:33:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-26 12:33:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-26 12:33:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-26 12:33:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-26 12:33:45 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 12:33:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 12:33:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 12:33:45 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 12:33:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 12:33:45 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 12:33:59 INFO Found: []
2026-07-26 12:33:59 INFO No job id returned by command ./run_filter.bsh
2026-07-26 12:33:59 INFO No monitoring will be performed
2026-07-26 12:33:59 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-26 12:33:59 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:33:59 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:33:59 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:33:59 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 12:34:00 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 12:34:02 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 12:34:02 INFO run_dart() is DONE.
2026-07-26 12:34:02 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 12:34:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-26 12:34:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-26 12:34:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-26 12:34:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-26 12:34:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-26 12:34:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-26 12:34:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-26 12:34:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-26 12:34:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-26 12:34:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-26 12:34:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:34:57 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-26 12:35:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:35:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-26 12:35:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:35:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-26 12:35:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:35:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-26 12:35:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:35:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-26 12:35:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:35:23 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 12:35:23 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:35:23 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:35:23 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-26 12:35:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:25 INFO Hourly dataset computed and listing created
2026-07-26 12:35:29 INFO Hourly dataset computed
2026-07-26 12:35:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:31 INFO Hourly dataset computed and listing created
2026-07-26 12:35:33 INFO Hourly dataset computed
2026-07-26 12:35:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:34 INFO Hourly dataset computed and listing created
2026-07-26 12:35:36 INFO Hourly dataset computed
2026-07-26 12:35:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:37 INFO Hourly dataset computed and listing created
2026-07-26 12:35:39 INFO Hourly dataset computed
2026-07-26 12:35:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:41 INFO Hourly dataset computed and listing created
2026-07-26 12:35:43 INFO Hourly dataset computed
2026-07-26 12:35:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:44 INFO Hourly dataset computed and listing created
2026-07-26 12:35:46 INFO Hourly dataset computed
2026-07-26 12:35:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:47 INFO Hourly dataset computed and listing created
2026-07-26 12:35:48 INFO Hourly dataset computed
2026-07-26 12:35:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:49 INFO Hourly dataset computed and listing created
2026-07-26 12:35:50 INFO Hourly dataset computed
2026-07-26 12:35:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:51 INFO Hourly dataset computed and listing created
2026-07-26 12:35:52 INFO Hourly dataset computed
2026-07-26 12:35:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:53 INFO Hourly dataset computed and listing created
2026-07-26 12:35:54 INFO Hourly dataset computed
2026-07-26 12:35:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:55 INFO Hourly dataset computed and listing created
2026-07-26 12:35:55 INFO Hourly dataset computed
2026-07-26 12:35:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:56 INFO Hourly dataset computed and listing created
2026-07-26 12:35:57 INFO Hourly dataset computed
2026-07-26 12:35:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:35:58 INFO Hourly dataset computed and listing created
2026-07-26 12:35:59 INFO Hourly dataset computed
2026-07-26 12:35:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:36:00 INFO Hourly dataset computed and listing created
2026-07-26 12:36:00 INFO Hourly dataset computed
2026-07-26 12:36:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:36:01 INFO Hourly dataset computed and listing created
2026-07-26 12:36:02 INFO Hourly dataset computed
2026-07-26 12:36:02 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-26 12:36:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:36:02 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-26 12:36:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:36:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:02 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:36:02 INFO Queuing job for member 1...
2026-07-26 12:36:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:02 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:36:03 INFO Found: ['5279619']
2026-07-26 12:36:08 INFO [TGCC-IRENE] Submitted job with ID:['5279619']
2026-07-26 12:36:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:36:08 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-26 12:36:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:36:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:08 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:36:08 INFO Queuing job for member 2...
2026-07-26 12:36:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:08 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:36:10 INFO Found: ['5279621']
2026-07-26 12:36:15 INFO [TGCC-IRENE] Submitted job with ID:['5279621']
2026-07-26 12:36:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:36:15 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-26 12:36:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:36:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:15 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:36:15 INFO Queuing job for member 3...
2026-07-26 12:36:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:15 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:36:18 INFO Found: ['5279622']
2026-07-26 12:36:23 INFO [TGCC-IRENE] Submitted job with ID:['5279622']
2026-07-26 12:36:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:36:23 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-26 12:36:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:36:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:23 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:36:23 INFO Queuing job for member 4...
2026-07-26 12:36:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:23 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:36:25 INFO Found: ['5279623']
2026-07-26 12:36:30 INFO [TGCC-IRENE] Submitted job with ID:['5279623']
2026-07-26 12:36:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:36:30 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-26 12:36:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:36:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:30 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:36:30 INFO Queuing job for member 5...
2026-07-26 12:36:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:30 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:36:33 INFO Found: ['5279624']
2026-07-26 12:36:38 INFO [TGCC-IRENE] Submitted job with ID:['5279624']
2026-07-26 12:36:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:36:38 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-26 12:36:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:36:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:38 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:36:38 INFO Queuing job for member 6...
2026-07-26 12:36:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:38 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:36:40 INFO Found: ['5279625']
2026-07-26 12:36:45 INFO [TGCC-IRENE] Submitted job with ID:['5279625']
2026-07-26 12:36:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:36:45 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-26 12:36:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:36:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:45 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:36:45 INFO Queuing job for member 7...
2026-07-26 12:36:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:45 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:36:48 INFO Found: ['5279626']
2026-07-26 12:36:53 INFO [TGCC-IRENE] Submitted job with ID:['5279626']
2026-07-26 12:36:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:36:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 12:36:53 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-26 12:36:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 12:36:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:36:53 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 12:36:53 INFO Queuing job for member 8...
2026-07-26 12:36:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:36:53 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 12:36:55 INFO Found: ['5279627']
2026-07-26 12:37:00 INFO [TGCC-IRENE] Submitted job with ID:['5279627']
2026-07-26 12:37:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 12:37:00 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-26 12:37:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 12:37:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:00 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 12:37:00 INFO Queuing job for member 9...
2026-07-26 12:37:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:00 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 12:37:03 INFO Found: ['5279628']
2026-07-26 12:37:08 INFO [TGCC-IRENE] Submitted job with ID:['5279628']
2026-07-26 12:37:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 12:37:08 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-26 12:37:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 12:37:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:08 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 12:37:08 INFO Queuing job for member 10...
2026-07-26 12:37:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:08 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 12:37:09 INFO Found: ['5279630']
2026-07-26 12:37:14 INFO [TGCC-IRENE] Submitted job with ID:['5279630']
2026-07-26 12:37:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 12:37:14 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-26 12:37:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 12:37:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:14 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 12:37:14 INFO Queuing job for member 11...
2026-07-26 12:37:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:14 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 12:37:14 INFO Found: ['5279631']
2026-07-26 12:37:19 INFO [TGCC-IRENE] Submitted job with ID:['5279631']
2026-07-26 12:37:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 12:37:19 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-26 12:37:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 12:37:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:19 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 12:37:19 INFO Queuing job for member 12...
2026-07-26 12:37:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:19 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 12:37:20 INFO Found: ['5279632']
2026-07-26 12:37:25 INFO [TGCC-IRENE] Submitted job with ID:['5279632']
2026-07-26 12:37:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 12:37:25 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-26 12:37:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 12:37:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:25 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 12:37:25 INFO Queuing job for member 13...
2026-07-26 12:37:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:25 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 12:37:26 INFO Found: ['5279633']
2026-07-26 12:37:31 INFO [TGCC-IRENE] Submitted job with ID:['5279633']
2026-07-26 12:37:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 12:37:31 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-26 12:37:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 12:37:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:31 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 12:37:31 INFO Queuing job for member 14...
2026-07-26 12:37:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:31 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 12:37:32 INFO Found: ['5279634']
2026-07-26 12:37:37 INFO [TGCC-IRENE] Submitted job with ID:['5279634']
2026-07-26 12:37:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:37:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 12:37:37 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-26 12:37:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 12:37:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:37:37 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 12:37:37 INFO Queuing job for member 15...
2026-07-26 12:37:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:37:37 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 12:37:37 INFO Found: ['5279635']
2026-07-26 12:37:42 INFO [TGCC-IRENE] Submitted job with ID:['5279635']
2026-07-26 12:37:42 INFO Checking job status ...
2026-07-26 12:37:42 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:37:42 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:37:43 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:37:43 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:37:43 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:37:43 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:37:43 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:37:58 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:37:58 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:37:58 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:38:13 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:38:13 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:38:13 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:38:13 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:38:13 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:38:13 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:38:14 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:38:16 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:38:16 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:38:31 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:38:31 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:38:31 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:38:46 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:38:46 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:38:48 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:38:48 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:39:03 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:39:03 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:39:04 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:39:04 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:39:19 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:39:19 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:39:21 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:39:21 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:39:36 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:39:36 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:39:36 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:39:51 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:39:51 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:39:51 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:40:08 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:40:08 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:40:08 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:40:23 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:40:23 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:40:25 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:40:25 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:40:25 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:40:40 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:40:40 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:40:40 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:40:55 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:40:55 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:40:55 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:40:55 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:40:55 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:40:56 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:40:58 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:40:58 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:41:13 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:41:13 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:41:13 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:41:28 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:41:28 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:41:28 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279627: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:41:30 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:41:31 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:41:31 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:41:31 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:41:46 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279624: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279625: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279626: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279627: status FINISHED
2026-07-26 12:41:46 INFO None 5279628: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:41:46 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:41:46 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:42:01 INFO None 5279619: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279621: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279622: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279624: status FINISHED
2026-07-26 12:42:01 INFO None 5279625: status FINISHED
2026-07-26 12:42:01 INFO None 5279626: status FINISHED
2026-07-26 12:42:01 INFO None 5279627: status FINISHED
2026-07-26 12:42:01 INFO None 5279628: status FINISHED
2026-07-26 12:42:01 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:42:01 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:42:01 INFO Jobs still running: ['5279619', '5279621', '5279622', '5279623', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:42:16 INFO None 5279619: status FINISHED
2026-07-26 12:42:16 INFO None 5279621: status FINISHED
2026-07-26 12:42:16 INFO None 5279622: status FINISHED
2026-07-26 12:42:16 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279624: status FINISHED
2026-07-26 12:42:16 INFO None 5279625: status FINISHED
2026-07-26 12:42:16 INFO None 5279626: status FINISHED
2026-07-26 12:42:16 INFO None 5279627: status FINISHED
2026-07-26 12:42:16 INFO None 5279628: status FINISHED
2026-07-26 12:42:16 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279634: status RUNNING/PENDING
2026-07-26 12:42:16 INFO None 5279635: status RUNNING/PENDING
2026-07-26 12:42:16 INFO Jobs still running: ['5279623', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635']. Waiting...
2026-07-26 12:42:33 INFO None 5279619: status FINISHED
2026-07-26 12:42:33 INFO None 5279621: status FINISHED
2026-07-26 12:42:33 INFO None 5279622: status FINISHED
2026-07-26 12:42:33 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:42:33 INFO None 5279624: status FINISHED
2026-07-26 12:42:33 INFO None 5279625: status FINISHED
2026-07-26 12:42:33 INFO None 5279626: status FINISHED
2026-07-26 12:42:33 INFO None 5279627: status FINISHED
2026-07-26 12:42:33 INFO None 5279628: status FINISHED
2026-07-26 12:42:33 INFO None 5279630: status RUNNING/PENDING
2026-07-26 12:42:33 INFO None 5279631: status RUNNING/PENDING
2026-07-26 12:42:33 INFO None 5279632: status RUNNING/PENDING
2026-07-26 12:42:33 INFO None 5279633: status RUNNING/PENDING
2026-07-26 12:42:33 INFO None 5279634: status FINISHED
2026-07-26 12:42:33 INFO None 5279635: status FINISHED
2026-07-26 12:42:33 INFO Jobs still running: ['5279623', '5279630', '5279631', '5279632', '5279633']. Waiting...
2026-07-26 12:42:48 INFO None 5279619: status FINISHED
2026-07-26 12:42:48 INFO None 5279621: status FINISHED
2026-07-26 12:42:48 INFO None 5279622: status FINISHED
2026-07-26 12:42:48 INFO None 5279623: status RUNNING/PENDING
2026-07-26 12:42:48 INFO None 5279624: status FINISHED
2026-07-26 12:42:48 INFO None 5279625: status FINISHED
2026-07-26 12:42:48 INFO None 5279626: status FINISHED
2026-07-26 12:42:48 INFO None 5279627: status FINISHED
2026-07-26 12:42:48 INFO None 5279628: status FINISHED
2026-07-26 12:42:48 INFO None 5279630: status FINISHED
2026-07-26 12:42:48 INFO None 5279631: status FINISHED
2026-07-26 12:42:48 INFO None 5279632: status FINISHED
2026-07-26 12:42:48 INFO None 5279633: status FINISHED
2026-07-26 12:42:48 INFO None 5279634: status FINISHED
2026-07-26 12:42:48 INFO None 5279635: status FINISHED
2026-07-26 12:42:48 INFO Jobs still running: ['5279623']. Waiting...
2026-07-26 12:43:03 INFO None 5279619: status FINISHED
2026-07-26 12:43:03 INFO None 5279621: status FINISHED
2026-07-26 12:43:03 INFO None 5279622: status FINISHED
2026-07-26 12:43:05 INFO None 5279623: status FINISHED
2026-07-26 12:43:06 INFO None 5279624: status FINISHED
2026-07-26 12:43:06 INFO None 5279625: status FINISHED
2026-07-26 12:43:06 INFO None 5279626: status FINISHED
2026-07-26 12:43:06 INFO None 5279627: status FINISHED
2026-07-26 12:43:06 INFO None 5279628: status FINISHED
2026-07-26 12:43:06 INFO None 5279630: status FINISHED
2026-07-26 12:43:06 INFO None 5279631: status FINISHED
2026-07-26 12:43:06 INFO None 5279632: status FINISHED
2026-07-26 12:43:06 INFO None 5279633: status FINISHED
2026-07-26 12:43:06 INFO None 5279634: status FINISHED
2026-07-26 12:43:06 INFO None 5279635: status FINISHED
2026-07-26 12:43:06 INFO Jobs ['5279619', '5279621', '5279622', '5279623', '5279624', '5279625', '5279626', '5279627', '5279628', '5279630', '5279631', '5279632', '5279633', '5279634', '5279635'] have finished
2026-07-26 12:43:06 INFO Checking restart files were created ...
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-26 12:43:06 INFO  Run_model() completed successfully.
2026-07-26 12:43:06 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:43:06 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-26 12:43:06 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 12:43:06 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-26 12:43:06 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:43:06 INFO ---------->>> Running process_satellite_data()
2026-07-26 12:43:06 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-26 12:43:06 INFO ---------->>> Running run_obs_converter()
2026-07-26 12:43:06 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-26 12:43:06 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-26 12:43:06 INFO ---------->>> Running DART
2026-07-26 12:43:06 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 12:43:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-26 12:43:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-26 12:43:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-26 12:43:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-26 12:43:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-26 12:43:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-26 12:43:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-26 12:43:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-26 12:43:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-26 12:43:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-26 12:43:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-26 12:43:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-26 12:43:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-26 12:43:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-26 12:43:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-26 12:43:11 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 12:43:11 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 12:43:11 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 12:43:11 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 12:43:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 12:43:11 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 12:43:27 INFO Found: []
2026-07-26 12:43:27 INFO No job id returned by command ./run_filter.bsh
2026-07-26 12:43:27 INFO No monitoring will be performed
2026-07-26 12:43:27 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-26 12:43:27 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:43:27 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 12:43:27 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 12:43:27 INFO run_dart() is DONE.
2026-07-26 12:43:27 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 12:43:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-26 12:43:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-26 12:43:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-26 12:43:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-26 12:43:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-26 12:43:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-26 12:43:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:43:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-26 12:44:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-26 12:44:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-26 12:44:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-26 12:44:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-26 12:44:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-26 12:44:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-26 12:44:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-26 12:44:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-26 12:44:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:44:49 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 12:44:49 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:44:49 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:44:49 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-26 12:44:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:44:50 INFO Hourly dataset computed and listing created
2026-07-26 12:44:55 INFO Hourly dataset computed
2026-07-26 12:44:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:44:56 INFO Hourly dataset computed and listing created
2026-07-26 12:44:57 INFO Hourly dataset computed
2026-07-26 12:44:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:44:58 INFO Hourly dataset computed and listing created
2026-07-26 12:44:59 INFO Hourly dataset computed
2026-07-26 12:44:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:00 INFO Hourly dataset computed and listing created
2026-07-26 12:45:00 INFO Hourly dataset computed
2026-07-26 12:45:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:01 INFO Hourly dataset computed and listing created
2026-07-26 12:45:02 INFO Hourly dataset computed
2026-07-26 12:45:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:03 INFO Hourly dataset computed and listing created
2026-07-26 12:45:03 INFO Hourly dataset computed
2026-07-26 12:45:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:04 INFO Hourly dataset computed and listing created
2026-07-26 12:45:05 INFO Hourly dataset computed
2026-07-26 12:45:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:06 INFO Hourly dataset computed and listing created
2026-07-26 12:45:06 INFO Hourly dataset computed
2026-07-26 12:45:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:07 INFO Hourly dataset computed and listing created
2026-07-26 12:45:08 INFO Hourly dataset computed
2026-07-26 12:45:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:09 INFO Hourly dataset computed and listing created
2026-07-26 12:45:10 INFO Hourly dataset computed
2026-07-26 12:45:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:11 INFO Hourly dataset computed and listing created
2026-07-26 12:45:11 INFO Hourly dataset computed
2026-07-26 12:45:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:12 INFO Hourly dataset computed and listing created
2026-07-26 12:45:13 INFO Hourly dataset computed
2026-07-26 12:45:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:13 INFO Hourly dataset computed and listing created
2026-07-26 12:45:14 INFO Hourly dataset computed
2026-07-26 12:45:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:15 INFO Hourly dataset computed and listing created
2026-07-26 12:45:16 INFO Hourly dataset computed
2026-07-26 12:45:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:45:17 INFO Hourly dataset computed and listing created
2026-07-26 12:45:17 INFO Hourly dataset computed
2026-07-26 12:45:17 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-26 12:45:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:45:17 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-26 12:45:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:45:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:17 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:45:17 INFO Queuing job for member 1...
2026-07-26 12:45:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:17 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:45:18 INFO Found: ['5279648']
2026-07-26 12:45:23 INFO [TGCC-IRENE] Submitted job with ID:['5279648']
2026-07-26 12:45:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:45:23 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-26 12:45:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:45:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:23 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:45:23 INFO Queuing job for member 2...
2026-07-26 12:45:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:23 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:45:25 INFO Found: ['5279649']
2026-07-26 12:45:30 INFO [TGCC-IRENE] Submitted job with ID:['5279649']
2026-07-26 12:45:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:45:30 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-26 12:45:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:45:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:30 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:45:30 INFO Queuing job for member 3...
2026-07-26 12:45:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:30 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:45:33 INFO Found: ['5279650']
2026-07-26 12:45:38 INFO [TGCC-IRENE] Submitted job with ID:['5279650']
2026-07-26 12:45:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:45:38 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-26 12:45:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:45:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:38 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:45:38 INFO Queuing job for member 4...
2026-07-26 12:45:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:38 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:45:40 INFO Found: ['5279651']
2026-07-26 12:45:45 INFO [TGCC-IRENE] Submitted job with ID:['5279651']
2026-07-26 12:45:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:45:45 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-26 12:45:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:45:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:45 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:45:45 INFO Queuing job for member 5...
2026-07-26 12:45:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:45 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:45:48 INFO Found: ['5279652']
2026-07-26 12:45:53 INFO [TGCC-IRENE] Submitted job with ID:['5279652']
2026-07-26 12:45:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:45:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:45:53 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-26 12:45:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:45:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:45:53 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:45:53 INFO Queuing job for member 6...
2026-07-26 12:45:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:45:53 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:45:55 INFO Found: ['5279653']
2026-07-26 12:46:00 INFO [TGCC-IRENE] Submitted job with ID:['5279653']
2026-07-26 12:46:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:46:00 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-26 12:46:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:46:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:00 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:46:00 INFO Queuing job for member 7...
2026-07-26 12:46:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:00 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:46:03 INFO Found: ['5279654']
2026-07-26 12:46:08 INFO [TGCC-IRENE] Submitted job with ID:['5279654']
2026-07-26 12:46:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 12:46:08 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-26 12:46:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 12:46:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:08 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 12:46:08 INFO Queuing job for member 8...
2026-07-26 12:46:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:08 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 12:46:09 INFO Found: ['5279656']
2026-07-26 12:46:14 INFO [TGCC-IRENE] Submitted job with ID:['5279656']
2026-07-26 12:46:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 12:46:14 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-26 12:46:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 12:46:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:14 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 12:46:14 INFO Queuing job for member 9...
2026-07-26 12:46:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:14 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 12:46:14 INFO Found: ['5279657']
2026-07-26 12:46:19 INFO [TGCC-IRENE] Submitted job with ID:['5279657']
2026-07-26 12:46:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 12:46:19 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-26 12:46:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 12:46:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:19 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 12:46:19 INFO Queuing job for member 10...
2026-07-26 12:46:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:19 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 12:46:20 INFO Found: ['5279658']
2026-07-26 12:46:25 INFO [TGCC-IRENE] Submitted job with ID:['5279658']
2026-07-26 12:46:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 12:46:25 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-26 12:46:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 12:46:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:25 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 12:46:25 INFO Queuing job for member 11...
2026-07-26 12:46:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:25 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 12:46:26 INFO Found: ['5279659']
2026-07-26 12:46:31 INFO [TGCC-IRENE] Submitted job with ID:['5279659']
2026-07-26 12:46:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 12:46:31 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-26 12:46:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 12:46:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:31 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 12:46:31 INFO Queuing job for member 12...
2026-07-26 12:46:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:31 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 12:46:32 INFO Found: ['5279660']
2026-07-26 12:46:37 INFO [TGCC-IRENE] Submitted job with ID:['5279660']
2026-07-26 12:46:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 12:46:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-26 12:46:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 12:46:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 12:46:37 INFO Queuing job for member 13...
2026-07-26 12:46:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 12:46:37 INFO Found: ['5279661']
2026-07-26 12:46:42 INFO [TGCC-IRENE] Submitted job with ID:['5279661']
2026-07-26 12:46:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 12:46:42 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-26 12:46:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 12:46:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:42 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 12:46:42 INFO Queuing job for member 14...
2026-07-26 12:46:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:42 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 12:46:43 INFO Found: ['5279662']
2026-07-26 12:46:48 INFO [TGCC-IRENE] Submitted job with ID:['5279662']
2026-07-26 12:46:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:46:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 12:46:48 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-26 12:46:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 12:46:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:46:48 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 12:46:48 INFO Queuing job for member 15...
2026-07-26 12:46:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:46:48 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 12:46:49 INFO Found: ['5279663']
2026-07-26 12:46:54 INFO [TGCC-IRENE] Submitted job with ID:['5279663']
2026-07-26 12:46:54 INFO Checking job status ...
2026-07-26 12:46:55 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:46:55 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:46:56 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:46:56 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:47:11 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:47:11 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:47:11 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:47:26 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:47:28 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:47:28 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:47:43 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:47:43 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:47:44 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:47:44 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:47:59 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:47:59 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:47:59 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:48:01 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:48:01 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:48:16 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:48:16 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:48:16 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:48:31 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:48:31 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:48:31 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:48:46 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:48:46 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:48:47 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:48:47 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:49:02 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:49:02 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:49:02 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:49:18 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:49:18 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:49:18 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:49:33 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279652: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279656: status RUNNING/PENDING
2026-07-26 12:49:33 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:49:35 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:49:35 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:49:50 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279652: status FINISHED
2026-07-26 12:49:50 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279656: status FINISHED
2026-07-26 12:49:50 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279659: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279660: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279661: status RUNNING/PENDING
2026-07-26 12:49:50 INFO None 5279662: status RUNNING/PENDING
2026-07-26 12:49:51 INFO None 5279663: status RUNNING/PENDING
2026-07-26 12:49:51 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279653', '5279654', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663']. Waiting...
2026-07-26 12:50:06 INFO None 5279648: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279649: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279650: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279652: status FINISHED
2026-07-26 12:50:06 INFO None 5279653: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279654: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279656: status FINISHED
2026-07-26 12:50:06 INFO None 5279657: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279658: status RUNNING/PENDING
2026-07-26 12:50:06 INFO None 5279659: status FINISHED
2026-07-26 12:50:08 INFO None 5279660: status FINISHED
2026-07-26 12:50:08 INFO None 5279661: status FINISHED
2026-07-26 12:50:08 INFO None 5279662: status FINISHED
2026-07-26 12:50:08 INFO None 5279663: status FINISHED
2026-07-26 12:50:08 INFO Jobs still running: ['5279648', '5279649', '5279650', '5279651', '5279653', '5279654', '5279657', '5279658']. Waiting...
2026-07-26 12:50:23 INFO None 5279648: status FINISHED
2026-07-26 12:50:23 INFO None 5279649: status FINISHED
2026-07-26 12:50:23 INFO None 5279650: status FINISHED
2026-07-26 12:50:23 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:50:23 INFO None 5279652: status FINISHED
2026-07-26 12:50:23 INFO None 5279653: status FINISHED
2026-07-26 12:50:23 INFO None 5279654: status FINISHED
2026-07-26 12:50:23 INFO None 5279656: status FINISHED
2026-07-26 12:50:23 INFO None 5279657: status FINISHED
2026-07-26 12:50:23 INFO None 5279658: status FINISHED
2026-07-26 12:50:23 INFO None 5279659: status FINISHED
2026-07-26 12:50:23 INFO None 5279660: status FINISHED
2026-07-26 12:50:23 INFO None 5279661: status FINISHED
2026-07-26 12:50:23 INFO None 5279662: status FINISHED
2026-07-26 12:50:23 INFO None 5279663: status FINISHED
2026-07-26 12:50:23 INFO Jobs still running: ['5279651']. Waiting...
2026-07-26 12:50:38 INFO None 5279648: status FINISHED
2026-07-26 12:50:38 INFO None 5279649: status FINISHED
2026-07-26 12:50:38 INFO None 5279650: status FINISHED
2026-07-26 12:50:38 INFO None 5279651: status RUNNING/PENDING
2026-07-26 12:50:38 INFO None 5279652: status FINISHED
2026-07-26 12:50:38 INFO None 5279653: status FINISHED
2026-07-26 12:50:38 INFO None 5279654: status FINISHED
2026-07-26 12:50:38 INFO None 5279656: status FINISHED
2026-07-26 12:50:38 INFO None 5279657: status FINISHED
2026-07-26 12:50:38 INFO None 5279658: status FINISHED
2026-07-26 12:50:38 INFO None 5279659: status FINISHED
2026-07-26 12:50:38 INFO None 5279660: status FINISHED
2026-07-26 12:50:38 INFO None 5279661: status FINISHED
2026-07-26 12:50:38 INFO None 5279662: status FINISHED
2026-07-26 12:50:38 INFO None 5279663: status FINISHED
2026-07-26 12:50:38 INFO Jobs still running: ['5279651']. Waiting...
2026-07-26 12:50:53 INFO None 5279648: status FINISHED
2026-07-26 12:50:53 INFO None 5279649: status FINISHED
2026-07-26 12:50:53 INFO None 5279650: status FINISHED
2026-07-26 12:50:53 INFO None 5279651: status FINISHED
2026-07-26 12:50:53 INFO None 5279652: status FINISHED
2026-07-26 12:50:54 INFO None 5279653: status FINISHED
2026-07-26 12:50:54 INFO None 5279654: status FINISHED
2026-07-26 12:50:54 INFO None 5279656: status FINISHED
2026-07-26 12:50:54 INFO None 5279657: status FINISHED
2026-07-26 12:50:54 INFO None 5279658: status FINISHED
2026-07-26 12:50:54 INFO None 5279659: status FINISHED
2026-07-26 12:50:54 INFO None 5279660: status FINISHED
2026-07-26 12:50:54 INFO None 5279661: status FINISHED
2026-07-26 12:50:54 INFO None 5279662: status FINISHED
2026-07-26 12:50:54 INFO None 5279663: status FINISHED
2026-07-26 12:50:54 INFO Jobs ['5279648', '5279649', '5279650', '5279651', '5279652', '5279653', '5279654', '5279656', '5279657', '5279658', '5279659', '5279660', '5279661', '5279662', '5279663'] have finished
2026-07-26 12:50:54 INFO Checking restart files were created ...
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-26 12:50:54 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-26 12:50:54 INFO  Run_model() completed successfully.
2026-07-26 12:50:54 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:50:54 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-26 12:50:54 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 12:50:54 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-26 12:50:54 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:50:54 INFO ---------->>> Running process_satellite_data()
2026-07-26 12:50:54 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-26 12:50:54 INFO ---------->>> Running run_obs_converter()
2026-07-26 12:50:54 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-26 12:50:54 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-26 12:50:54 INFO ---------->>> Running DART
2026-07-26 12:50:54 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 12:50:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-26 12:50:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-26 12:50:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-26 12:50:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-26 12:50:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-26 12:50:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-26 12:50:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-26 12:50:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-26 12:50:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-26 12:50:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-26 12:50:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-26 12:50:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-26 12:50:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-26 12:50:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-26 12:50:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-26 12:50:59 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 12:50:59 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 12:50:59 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 12:50:59 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 12:50:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 12:50:59 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 12:51:07 INFO Found: []
2026-07-26 12:51:07 INFO No job id returned by command ./run_filter.bsh
2026-07-26 12:51:07 INFO No monitoring will be performed
2026-07-26 12:51:07 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-26 12:51:07 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:07 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:07 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020614'
2026-07-26 12:51:08 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 12:51:08 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 12:51:08 INFO run_dart() is DONE.
2026-07-26 12:51:08 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 12:51:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-26 12:51:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-26 12:51:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-26 12:51:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-26 12:51:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-26 12:51:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-26 12:51:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-26 12:51:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-26 12:51:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-26 12:51:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-26 12:51:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-26 12:51:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-26 12:51:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-26 12:51:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:51:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-26 12:52:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:52:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-26 12:52:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:52:04 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 12:52:04 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:52:04 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:52:04 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-26 12:52:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:52:05 INFO Hourly dataset computed and listing created
2026-07-26 12:52:23 INFO Hourly dataset computed
2026-07-26 12:52:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:52:24 INFO Hourly dataset computed and listing created
2026-07-26 12:52:41 INFO Hourly dataset computed
2026-07-26 12:52:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:52:42 INFO Hourly dataset computed and listing created
2026-07-26 12:52:58 INFO Hourly dataset computed
2026-07-26 12:52:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:52:59 INFO Hourly dataset computed and listing created
2026-07-26 12:53:15 INFO Hourly dataset computed
2026-07-26 12:53:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:53:17 INFO Hourly dataset computed and listing created
2026-07-26 12:53:32 INFO Hourly dataset computed
2026-07-26 12:53:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:53:34 INFO Hourly dataset computed and listing created
2026-07-26 12:53:49 INFO Hourly dataset computed
2026-07-26 12:53:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:53:50 INFO Hourly dataset computed and listing created
2026-07-26 12:54:06 INFO Hourly dataset computed
2026-07-26 12:54:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:54:08 INFO Hourly dataset computed and listing created
2026-07-26 12:54:24 INFO Hourly dataset computed
2026-07-26 12:54:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:54:26 INFO Hourly dataset computed and listing created
2026-07-26 12:54:42 INFO Hourly dataset computed
2026-07-26 12:54:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:54:43 INFO Hourly dataset computed and listing created
2026-07-26 12:54:59 INFO Hourly dataset computed
2026-07-26 12:54:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:55:01 INFO Hourly dataset computed and listing created
2026-07-26 12:55:04 INFO Hourly dataset computed
2026-07-26 12:55:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:55:05 INFO Hourly dataset computed and listing created
2026-07-26 12:55:07 INFO Hourly dataset computed
2026-07-26 12:55:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:55:09 INFO Hourly dataset computed and listing created
2026-07-26 12:55:11 INFO Hourly dataset computed
2026-07-26 12:55:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:55:12 INFO Hourly dataset computed and listing created
2026-07-26 12:55:15 INFO Hourly dataset computed
2026-07-26 12:55:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:55:16 INFO Hourly dataset computed and listing created
2026-07-26 12:55:19 INFO Hourly dataset computed
2026-07-26 12:55:19 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-26 12:55:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:55:19 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-26 12:55:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:55:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:19 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:55:19 INFO Queuing job for member 1...
2026-07-26 12:55:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:19 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:55:20 INFO Found: ['5279686']
2026-07-26 12:55:25 INFO [TGCC-IRENE] Submitted job with ID:['5279686']
2026-07-26 12:55:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:55:25 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-26 12:55:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:55:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:25 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:55:25 INFO Queuing job for member 2...
2026-07-26 12:55:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:25 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:55:26 INFO Found: ['5279688']
2026-07-26 12:55:31 INFO [TGCC-IRENE] Submitted job with ID:['5279688']
2026-07-26 12:55:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:55:31 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-26 12:55:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:55:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:31 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:55:31 INFO Queuing job for member 3...
2026-07-26 12:55:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:31 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:55:31 INFO Found: ['5279691']
2026-07-26 12:55:36 INFO [TGCC-IRENE] Submitted job with ID:['5279691']
2026-07-26 12:55:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:55:36 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-26 12:55:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:55:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:36 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:55:36 INFO Queuing job for member 4...
2026-07-26 12:55:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:36 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:55:37 INFO Found: ['5279693']
2026-07-26 12:55:42 INFO [TGCC-IRENE] Submitted job with ID:['5279693']
2026-07-26 12:55:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:55:42 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-26 12:55:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:55:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:42 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:55:42 INFO Queuing job for member 5...
2026-07-26 12:55:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:42 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:55:45 INFO Found: ['5279694']
2026-07-26 12:55:50 INFO [TGCC-IRENE] Submitted job with ID:['5279694']
2026-07-26 12:55:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:55:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-26 12:55:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:55:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:55:50 INFO Queuing job for member 6...
2026-07-26 12:55:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:55:52 INFO Found: ['5279695']
2026-07-26 12:55:57 INFO [TGCC-IRENE] Submitted job with ID:['5279695']
2026-07-26 12:55:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:55:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:55:57 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-26 12:55:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:55:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:55:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:55:57 INFO Queuing job for member 7...
2026-07-26 12:55:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:55:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:56:00 INFO Found: ['5279697']
2026-07-26 12:56:05 INFO [TGCC-IRENE] Submitted job with ID:['5279697']
2026-07-26 12:56:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 12:56:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-26 12:56:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 12:56:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 12:56:05 INFO Queuing job for member 8...
2026-07-26 12:56:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 12:56:07 INFO Found: ['5279699']
2026-07-26 12:56:12 INFO [TGCC-IRENE] Submitted job with ID:['5279699']
2026-07-26 12:56:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 12:56:12 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-26 12:56:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 12:56:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:12 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 12:56:12 INFO Queuing job for member 9...
2026-07-26 12:56:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 12:56:15 INFO Found: ['5279700']
2026-07-26 12:56:20 INFO [TGCC-IRENE] Submitted job with ID:['5279700']
2026-07-26 12:56:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 12:56:20 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-26 12:56:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 12:56:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:20 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 12:56:20 INFO Queuing job for member 10...
2026-07-26 12:56:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:20 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 12:56:22 INFO Found: ['5279701']
2026-07-26 12:56:27 INFO [TGCC-IRENE] Submitted job with ID:['5279701']
2026-07-26 12:56:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 12:56:27 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-26 12:56:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 12:56:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:27 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 12:56:27 INFO Queuing job for member 11...
2026-07-26 12:56:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:27 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 12:56:30 INFO Found: ['5279702']
2026-07-26 12:56:35 INFO [TGCC-IRENE] Submitted job with ID:['5279702']
2026-07-26 12:56:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 12:56:35 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-26 12:56:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 12:56:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:35 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 12:56:35 INFO Queuing job for member 12...
2026-07-26 12:56:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:35 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 12:56:37 INFO Found: ['5279703']
2026-07-26 12:56:42 INFO [TGCC-IRENE] Submitted job with ID:['5279703']
2026-07-26 12:56:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 12:56:42 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-26 12:56:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 12:56:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:42 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 12:56:42 INFO Queuing job for member 13...
2026-07-26 12:56:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:42 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 12:56:45 INFO Found: ['5279704']
2026-07-26 12:56:50 INFO [TGCC-IRENE] Submitted job with ID:['5279704']
2026-07-26 12:56:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 12:56:50 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-26 12:56:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 12:56:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:50 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 12:56:50 INFO Queuing job for member 14...
2026-07-26 12:56:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:50 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 12:56:52 INFO Found: ['5279706']
2026-07-26 12:56:57 INFO [TGCC-IRENE] Submitted job with ID:['5279706']
2026-07-26 12:56:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:56:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 12:56:57 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-26 12:56:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 12:56:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:56:57 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 12:56:57 INFO Queuing job for member 15...
2026-07-26 12:56:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:56:57 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 12:57:00 INFO Found: ['5279708']
2026-07-26 12:57:05 INFO [TGCC-IRENE] Submitted job with ID:['5279708']
2026-07-26 12:57:05 INFO Checking job status ...
2026-07-26 12:57:05 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:57:05 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:57:05 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:57:20 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:57:20 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:57:20 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:57:35 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:57:35 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:57:36 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:57:36 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:57:36 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:57:36 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:57:36 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:57:51 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:57:51 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:57:53 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:57:53 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:57:53 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:57:53 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:57:53 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:58:08 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:58:08 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:58:08 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:58:23 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:58:23 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:58:25 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:58:25 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:58:25 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:58:40 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:58:41 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:58:41 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:58:56 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:58:56 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:58:58 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:58:58 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:59:13 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:59:13 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:59:13 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:59:28 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:59:28 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:59:29 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:59:29 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:59:44 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:59:44 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:59:44 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 12:59:59 INFO None 5279686: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279688: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279691: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279693: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279694: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279695: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279697: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279699: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279700: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279701: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279702: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279703: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279704: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279706: status RUNNING/PENDING
2026-07-26 12:59:59 INFO None 5279708: status RUNNING/PENDING
2026-07-26 12:59:59 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:00:15 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:00:15 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:00:15 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:00:30 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:00:30 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:00:32 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:00:32 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:00:32 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:00:32 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:00:32 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:00:32 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:00:47 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:00:47 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:00:47 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:01:05 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:01:05 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:01:05 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:01:20 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:01:20 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:01:22 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:01:22 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:01:23 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:01:23 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:01:38 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:01:38 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:01:38 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:01:53 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:01:53 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:01:53 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:02:08 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:02:08 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:02:08 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:02:25 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:02:25 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:02:25 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:02:40 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:02:40 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:02:41 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:02:41 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:02:41 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:02:41 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:02:41 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:02:43 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:02:43 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:02:58 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:02:58 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:02:58 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:03:13 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:03:13 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:03:15 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:03:15 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:03:30 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:03:30 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:03:31 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:03:31 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:03:31 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:03:31 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:03:31 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:03:46 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:03:46 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:03:47 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:03:47 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:03:47 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:04:02 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:04:02 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:04:02 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:04:18 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:04:18 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:04:18 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:04:34 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:04:34 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:04:34 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:04:49 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:04:50 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:04:52 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:04:52 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:05:07 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279693: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:05:07 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:05:07 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:05:22 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279693: status FINISHED
2026-07-26 13:05:22 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:05:22 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:05:24 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:05:24 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:05:24 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:05:24 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:05:24 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:05:39 INFO None 5279686: status RUNNING/PENDING
2026-07-26 13:05:39 INFO None 5279688: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279691: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279693: status FINISHED
2026-07-26 13:05:40 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:05:40 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:05:40 INFO Jobs still running: ['5279686', '5279688', '5279691', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:05:55 INFO None 5279686: status FINISHED
2026-07-26 13:05:55 INFO None 5279688: status FINISHED
2026-07-26 13:05:55 INFO None 5279691: status FINISHED
2026-07-26 13:05:55 INFO None 5279693: status FINISHED
2026-07-26 13:05:55 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:05:55 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:05:57 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:05:57 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:05:57 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:05:57 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:06:12 INFO None 5279686: status FINISHED
2026-07-26 13:06:12 INFO None 5279688: status FINISHED
2026-07-26 13:06:12 INFO None 5279691: status FINISHED
2026-07-26 13:06:12 INFO None 5279693: status FINISHED
2026-07-26 13:06:12 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:06:12 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:06:12 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:06:27 INFO None 5279686: status FINISHED
2026-07-26 13:06:27 INFO None 5279688: status FINISHED
2026-07-26 13:06:27 INFO None 5279691: status FINISHED
2026-07-26 13:06:27 INFO None 5279693: status FINISHED
2026-07-26 13:06:27 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:06:27 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:06:27 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:06:27 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:06:27 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:06:28 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:06:28 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:06:43 INFO None 5279686: status FINISHED
2026-07-26 13:06:43 INFO None 5279688: status FINISHED
2026-07-26 13:06:43 INFO None 5279691: status FINISHED
2026-07-26 13:06:43 INFO None 5279693: status FINISHED
2026-07-26 13:06:43 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:06:43 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:06:43 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:06:59 INFO None 5279686: status FINISHED
2026-07-26 13:06:59 INFO None 5279688: status FINISHED
2026-07-26 13:06:59 INFO None 5279691: status FINISHED
2026-07-26 13:06:59 INFO None 5279693: status FINISHED
2026-07-26 13:06:59 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:06:59 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:06:59 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:07:14 INFO None 5279686: status FINISHED
2026-07-26 13:07:14 INFO None 5279688: status FINISHED
2026-07-26 13:07:14 INFO None 5279691: status FINISHED
2026-07-26 13:07:14 INFO None 5279693: status FINISHED
2026-07-26 13:07:14 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:07:14 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:07:16 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:07:16 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:07:16 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:07:31 INFO None 5279686: status FINISHED
2026-07-26 13:07:31 INFO None 5279688: status FINISHED
2026-07-26 13:07:31 INFO None 5279691: status FINISHED
2026-07-26 13:07:31 INFO None 5279693: status FINISHED
2026-07-26 13:07:31 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279697: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279701: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:07:31 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:07:31 INFO Jobs still running: ['5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:07:46 INFO None 5279686: status FINISHED
2026-07-26 13:07:46 INFO None 5279688: status FINISHED
2026-07-26 13:07:46 INFO None 5279691: status FINISHED
2026-07-26 13:07:47 INFO None 5279693: status FINISHED
2026-07-26 13:07:47 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279697: status FINISHED
2026-07-26 13:07:47 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279701: status FINISHED
2026-07-26 13:07:47 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:07:47 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:07:47 INFO Jobs still running: ['5279694', '5279695', '5279699', '5279700', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:08:02 INFO None 5279686: status FINISHED
2026-07-26 13:08:02 INFO None 5279688: status FINISHED
2026-07-26 13:08:02 INFO None 5279691: status FINISHED
2026-07-26 13:08:04 INFO None 5279693: status FINISHED
2026-07-26 13:08:04 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279695: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279697: status FINISHED
2026-07-26 13:08:04 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279701: status FINISHED
2026-07-26 13:08:04 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:08:04 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:08:04 INFO Jobs still running: ['5279694', '5279695', '5279699', '5279700', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:08:19 INFO None 5279686: status FINISHED
2026-07-26 13:08:19 INFO None 5279688: status FINISHED
2026-07-26 13:08:19 INFO None 5279691: status FINISHED
2026-07-26 13:08:19 INFO None 5279693: status FINISHED
2026-07-26 13:08:19 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279695: status FINISHED
2026-07-26 13:08:19 INFO None 5279697: status FINISHED
2026-07-26 13:08:19 INFO None 5279699: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279700: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279701: status FINISHED
2026-07-26 13:08:19 INFO None 5279702: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:08:19 INFO None 5279708: status RUNNING/PENDING
2026-07-26 13:08:19 INFO Jobs still running: ['5279694', '5279699', '5279700', '5279702', '5279703', '5279704', '5279706', '5279708']. Waiting...
2026-07-26 13:08:34 INFO None 5279686: status FINISHED
2026-07-26 13:08:34 INFO None 5279688: status FINISHED
2026-07-26 13:08:34 INFO None 5279691: status FINISHED
2026-07-26 13:08:34 INFO None 5279693: status FINISHED
2026-07-26 13:08:34 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:08:34 INFO None 5279695: status FINISHED
2026-07-26 13:08:34 INFO None 5279697: status FINISHED
2026-07-26 13:08:34 INFO None 5279699: status FINISHED
2026-07-26 13:08:34 INFO None 5279700: status FINISHED
2026-07-26 13:08:35 INFO None 5279701: status FINISHED
2026-07-26 13:08:35 INFO None 5279702: status FINISHED
2026-07-26 13:08:35 INFO None 5279703: status RUNNING/PENDING
2026-07-26 13:08:35 INFO None 5279704: status RUNNING/PENDING
2026-07-26 13:08:35 INFO None 5279706: status RUNNING/PENDING
2026-07-26 13:08:35 INFO None 5279708: status FINISHED
2026-07-26 13:08:35 INFO Jobs still running: ['5279694', '5279703', '5279704', '5279706']. Waiting...
2026-07-26 13:08:50 INFO None 5279686: status FINISHED
2026-07-26 13:08:50 INFO None 5279688: status FINISHED
2026-07-26 13:08:50 INFO None 5279691: status FINISHED
2026-07-26 13:08:50 INFO None 5279693: status FINISHED
2026-07-26 13:08:50 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:08:50 INFO None 5279695: status FINISHED
2026-07-26 13:08:50 INFO None 5279697: status FINISHED
2026-07-26 13:08:50 INFO None 5279699: status FINISHED
2026-07-26 13:08:50 INFO None 5279700: status FINISHED
2026-07-26 13:08:50 INFO None 5279701: status FINISHED
2026-07-26 13:08:50 INFO None 5279702: status FINISHED
2026-07-26 13:08:50 INFO None 5279703: status FINISHED
2026-07-26 13:08:50 INFO None 5279704: status FINISHED
2026-07-26 13:08:50 INFO None 5279706: status FINISHED
2026-07-26 13:08:50 INFO None 5279708: status FINISHED
2026-07-26 13:08:50 INFO Jobs still running: ['5279694']. Waiting...
2026-07-26 13:09:06 INFO None 5279686: status FINISHED
2026-07-26 13:09:06 INFO None 5279688: status FINISHED
2026-07-26 13:09:06 INFO None 5279691: status FINISHED
2026-07-26 13:09:06 INFO None 5279693: status FINISHED
2026-07-26 13:09:06 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:09:06 INFO None 5279695: status FINISHED
2026-07-26 13:09:07 INFO None 5279697: status FINISHED
2026-07-26 13:09:07 INFO None 5279699: status FINISHED
2026-07-26 13:09:07 INFO None 5279700: status FINISHED
2026-07-26 13:09:07 INFO None 5279701: status FINISHED
2026-07-26 13:09:07 INFO None 5279702: status FINISHED
2026-07-26 13:09:07 INFO None 5279703: status FINISHED
2026-07-26 13:09:07 INFO None 5279704: status FINISHED
2026-07-26 13:09:07 INFO None 5279706: status FINISHED
2026-07-26 13:09:07 INFO None 5279708: status FINISHED
2026-07-26 13:09:07 INFO Jobs still running: ['5279694']. Waiting...
2026-07-26 13:09:22 INFO None 5279686: status FINISHED
2026-07-26 13:09:22 INFO None 5279688: status FINISHED
2026-07-26 13:09:22 INFO None 5279691: status FINISHED
2026-07-26 13:09:22 INFO None 5279693: status FINISHED
2026-07-26 13:09:22 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:09:22 INFO None 5279695: status FINISHED
2026-07-26 13:09:22 INFO None 5279697: status FINISHED
2026-07-26 13:09:22 INFO None 5279699: status FINISHED
2026-07-26 13:09:22 INFO None 5279700: status FINISHED
2026-07-26 13:09:22 INFO None 5279701: status FINISHED
2026-07-26 13:09:22 INFO None 5279702: status FINISHED
2026-07-26 13:09:22 INFO None 5279703: status FINISHED
2026-07-26 13:09:22 INFO None 5279704: status FINISHED
2026-07-26 13:09:24 INFO None 5279706: status FINISHED
2026-07-26 13:09:24 INFO None 5279708: status FINISHED
2026-07-26 13:09:24 INFO Jobs still running: ['5279694']. Waiting...
2026-07-26 13:09:39 INFO None 5279686: status FINISHED
2026-07-26 13:09:39 INFO None 5279688: status FINISHED
2026-07-26 13:09:39 INFO None 5279691: status FINISHED
2026-07-26 13:09:39 INFO None 5279693: status FINISHED
2026-07-26 13:09:39 INFO None 5279694: status RUNNING/PENDING
2026-07-26 13:09:39 INFO None 5279695: status FINISHED
2026-07-26 13:09:39 INFO None 5279697: status FINISHED
2026-07-26 13:09:39 INFO None 5279699: status FINISHED
2026-07-26 13:09:39 INFO None 5279700: status FINISHED
2026-07-26 13:09:39 INFO None 5279701: status FINISHED
2026-07-26 13:09:39 INFO None 5279702: status FINISHED
2026-07-26 13:09:39 INFO None 5279703: status FINISHED
2026-07-26 13:09:39 INFO None 5279704: status FINISHED
2026-07-26 13:09:39 INFO None 5279706: status FINISHED
2026-07-26 13:09:39 INFO None 5279708: status FINISHED
2026-07-26 13:09:39 INFO Jobs still running: ['5279694']. Waiting...
2026-07-26 13:09:54 INFO None 5279686: status FINISHED
2026-07-26 13:09:54 INFO None 5279688: status FINISHED
2026-07-26 13:09:54 INFO None 5279691: status FINISHED
2026-07-26 13:09:54 INFO None 5279693: status FINISHED
2026-07-26 13:09:54 INFO None 5279694: status FINISHED
2026-07-26 13:09:54 INFO None 5279695: status FINISHED
2026-07-26 13:09:54 INFO None 5279697: status FINISHED
2026-07-26 13:09:54 INFO None 5279699: status FINISHED
2026-07-26 13:09:54 INFO None 5279700: status FINISHED
2026-07-26 13:09:54 INFO None 5279701: status FINISHED
2026-07-26 13:09:54 INFO None 5279702: status FINISHED
2026-07-26 13:09:54 INFO None 5279703: status FINISHED
2026-07-26 13:09:55 INFO None 5279704: status FINISHED
2026-07-26 13:09:55 INFO None 5279706: status FINISHED
2026-07-26 13:09:57 INFO None 5279708: status FINISHED
2026-07-26 13:09:57 INFO Jobs ['5279686', '5279688', '5279691', '5279693', '5279694', '5279695', '5279697', '5279699', '5279700', '5279701', '5279702', '5279703', '5279704', '5279706', '5279708'] have finished
2026-07-26 13:09:57 INFO Checking restart files were created ...
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-26 13:09:57 INFO  Run_model() completed successfully.
2026-07-26 13:09:57 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:09:57 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-26 13:09:57 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 13:09:57 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-26 13:09:57 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:09:57 INFO ---------->>> Running process_satellite_data()
2026-07-26 13:09:57 INFO [DART] No satellite data found, skipping assimilation
2026-07-26 13:09:57 INFO after_assimilation() skipped
2026-07-26 13:09:57 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 13:09:57 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:09:57 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:09:57 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-26 13:09:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:09:58 INFO Hourly dataset computed and listing created
2026-07-26 13:10:03 INFO Hourly dataset computed
2026-07-26 13:10:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:04 INFO Hourly dataset computed and listing created
2026-07-26 13:10:05 INFO Hourly dataset computed
2026-07-26 13:10:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:06 INFO Hourly dataset computed and listing created
2026-07-26 13:10:06 INFO Hourly dataset computed
2026-07-26 13:10:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:07 INFO Hourly dataset computed and listing created
2026-07-26 13:10:08 INFO Hourly dataset computed
2026-07-26 13:10:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:09 INFO Hourly dataset computed and listing created
2026-07-26 13:10:09 INFO Hourly dataset computed
2026-07-26 13:10:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:10 INFO Hourly dataset computed and listing created
2026-07-26 13:10:11 INFO Hourly dataset computed
2026-07-26 13:10:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:12 INFO Hourly dataset computed and listing created
2026-07-26 13:10:12 INFO Hourly dataset computed
2026-07-26 13:10:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:13 INFO Hourly dataset computed and listing created
2026-07-26 13:10:14 INFO Hourly dataset computed
2026-07-26 13:10:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:14 INFO Hourly dataset computed and listing created
2026-07-26 13:10:15 INFO Hourly dataset computed
2026-07-26 13:10:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:16 INFO Hourly dataset computed and listing created
2026-07-26 13:10:17 INFO Hourly dataset computed
2026-07-26 13:10:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:17 INFO Hourly dataset computed and listing created
2026-07-26 13:10:18 INFO Hourly dataset computed
2026-07-26 13:10:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:19 INFO Hourly dataset computed and listing created
2026-07-26 13:10:19 INFO Hourly dataset computed
2026-07-26 13:10:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:20 INFO Hourly dataset computed and listing created
2026-07-26 13:10:21 INFO Hourly dataset computed
2026-07-26 13:10:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:22 INFO Hourly dataset computed and listing created
2026-07-26 13:10:22 INFO Hourly dataset computed
2026-07-26 13:10:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:10:23 INFO Hourly dataset computed and listing created
2026-07-26 13:10:24 INFO Hourly dataset computed
2026-07-26 13:10:24 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-26 13:10:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 13:10:24 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-26 13:10:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 13:10:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:24 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 13:10:24 INFO Queuing job for member 1...
2026-07-26 13:10:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:24 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 13:10:25 INFO Found: ['5279744']
2026-07-26 13:10:30 INFO [TGCC-IRENE] Submitted job with ID:['5279744']
2026-07-26 13:10:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 13:10:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-26 13:10:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 13:10:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 13:10:30 INFO Queuing job for member 2...
2026-07-26 13:10:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 13:10:30 INFO Found: ['5279745']
2026-07-26 13:10:35 INFO [TGCC-IRENE] Submitted job with ID:['5279745']
2026-07-26 13:10:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 13:10:35 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-26 13:10:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 13:10:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:35 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 13:10:35 INFO Queuing job for member 3...
2026-07-26 13:10:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:35 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 13:10:36 INFO Found: ['5279747']
2026-07-26 13:10:41 INFO [TGCC-IRENE] Submitted job with ID:['5279747']
2026-07-26 13:10:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 13:10:41 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-26 13:10:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 13:10:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:41 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 13:10:41 INFO Queuing job for member 4...
2026-07-26 13:10:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:41 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 13:10:42 INFO Found: ['5279748']
2026-07-26 13:10:47 INFO [TGCC-IRENE] Submitted job with ID:['5279748']
2026-07-26 13:10:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 13:10:47 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-26 13:10:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 13:10:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:47 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 13:10:47 INFO Queuing job for member 5...
2026-07-26 13:10:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:47 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 13:10:48 INFO Found: ['5279749']
2026-07-26 13:10:53 INFO [TGCC-IRENE] Submitted job with ID:['5279749']
2026-07-26 13:10:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 13:10:53 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-26 13:10:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 13:10:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:53 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 13:10:53 INFO Queuing job for member 6...
2026-07-26 13:10:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:53 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 13:10:53 INFO Found: ['5279750']
2026-07-26 13:10:58 INFO [TGCC-IRENE] Submitted job with ID:['5279750']
2026-07-26 13:10:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:10:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 13:10:58 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-26 13:10:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 13:10:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:10:58 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 13:10:58 INFO Queuing job for member 7...
2026-07-26 13:10:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:10:58 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 13:10:59 INFO Found: ['5279751']
2026-07-26 13:11:04 INFO [TGCC-IRENE] Submitted job with ID:['5279751']
2026-07-26 13:11:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 13:11:04 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-26 13:11:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 13:11:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:04 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 13:11:04 INFO Queuing job for member 8...
2026-07-26 13:11:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:04 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 13:11:05 INFO Found: ['5279753']
2026-07-26 13:11:10 INFO [TGCC-IRENE] Submitted job with ID:['5279753']
2026-07-26 13:11:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 13:11:10 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-26 13:11:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 13:11:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:10 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 13:11:10 INFO Queuing job for member 9...
2026-07-26 13:11:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:10 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 13:11:11 INFO Found: ['5279754']
2026-07-26 13:11:16 INFO [TGCC-IRENE] Submitted job with ID:['5279754']
2026-07-26 13:11:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 13:11:16 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-26 13:11:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 13:11:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:16 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 13:11:16 INFO Queuing job for member 10...
2026-07-26 13:11:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:16 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 13:11:17 INFO Found: ['5279755']
2026-07-26 13:11:22 INFO [TGCC-IRENE] Submitted job with ID:['5279755']
2026-07-26 13:11:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 13:11:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-26 13:11:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 13:11:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 13:11:22 INFO Queuing job for member 11...
2026-07-26 13:11:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:22 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 13:11:25 INFO Found: ['5279756']
2026-07-26 13:11:30 INFO [TGCC-IRENE] Submitted job with ID:['5279756']
2026-07-26 13:11:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 13:11:30 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-26 13:11:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 13:11:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:30 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 13:11:30 INFO Queuing job for member 12...
2026-07-26 13:11:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:30 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 13:11:32 INFO Found: ['5279757']
2026-07-26 13:11:37 INFO [TGCC-IRENE] Submitted job with ID:['5279757']
2026-07-26 13:11:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 13:11:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-26 13:11:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 13:11:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 13:11:37 INFO Queuing job for member 13...
2026-07-26 13:11:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 13:11:40 INFO Found: ['5279758']
2026-07-26 13:11:45 INFO [TGCC-IRENE] Submitted job with ID:['5279758']
2026-07-26 13:11:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 13:11:45 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-26 13:11:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 13:11:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:45 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 13:11:45 INFO Queuing job for member 14...
2026-07-26 13:11:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:45 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 13:11:47 INFO Found: ['5279759']
2026-07-26 13:11:52 INFO [TGCC-IRENE] Submitted job with ID:['5279759']
2026-07-26 13:11:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:11:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 13:11:52 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-26 13:11:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 13:11:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:11:52 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 13:11:52 INFO Queuing job for member 15...
2026-07-26 13:11:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:11:52 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 13:11:55 INFO Found: ['5279760']
2026-07-26 13:12:00 INFO [TGCC-IRENE] Submitted job with ID:['5279760']
2026-07-26 13:12:00 INFO Checking job status ...
2026-07-26 13:12:00 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:12:00 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:12:00 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:12:15 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:12:15 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:12:17 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:12:17 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:12:32 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:12:32 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:12:33 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:12:33 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:12:48 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:12:48 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:12:48 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:13:03 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:13:03 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:13:03 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:13:18 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:13:18 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:13:19 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:13:19 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:13:19 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:13:19 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:13:19 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:13:19 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:13:34 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:13:34 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:13:35 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:13:35 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:13:50 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279749: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279750: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:13:50 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:13:52 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:13:52 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:13:52 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:14:07 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279749: status FINISHED
2026-07-26 13:14:07 INFO None 5279750: status FINISHED
2026-07-26 13:14:07 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:14:07 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:14:07 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:14:22 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279749: status FINISHED
2026-07-26 13:14:22 INFO None 5279750: status FINISHED
2026-07-26 13:14:22 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279753: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279754: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:14:22 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:14:24 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:14:24 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:14:25 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:14:25 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:14:40 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279749: status FINISHED
2026-07-26 13:14:40 INFO None 5279750: status FINISHED
2026-07-26 13:14:40 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279753: status FINISHED
2026-07-26 13:14:40 INFO None 5279754: status FINISHED
2026-07-26 13:14:40 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:14:40 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:14:40 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279751', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:14:55 INFO None 5279744: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279745: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279747: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279748: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279749: status FINISHED
2026-07-26 13:14:55 INFO None 5279750: status FINISHED
2026-07-26 13:14:55 INFO None 5279751: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279753: status FINISHED
2026-07-26 13:14:55 INFO None 5279754: status FINISHED
2026-07-26 13:14:55 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:14:55 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:14:55 INFO Jobs still running: ['5279744', '5279745', '5279747', '5279748', '5279751', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:15:10 INFO None 5279744: status FINISHED
2026-07-26 13:15:10 INFO None 5279745: status FINISHED
2026-07-26 13:15:10 INFO None 5279747: status FINISHED
2026-07-26 13:15:10 INFO None 5279748: status FINISHED
2026-07-26 13:15:10 INFO None 5279749: status FINISHED
2026-07-26 13:15:10 INFO None 5279750: status FINISHED
2026-07-26 13:15:10 INFO None 5279751: status FINISHED
2026-07-26 13:15:10 INFO None 5279753: status FINISHED
2026-07-26 13:15:10 INFO None 5279754: status FINISHED
2026-07-26 13:15:10 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:15:10 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:15:10 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:15:10 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:15:10 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:15:10 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:15:10 INFO Jobs still running: ['5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:15:25 INFO None 5279744: status FINISHED
2026-07-26 13:15:25 INFO None 5279745: status FINISHED
2026-07-26 13:15:25 INFO None 5279747: status FINISHED
2026-07-26 13:15:25 INFO None 5279748: status FINISHED
2026-07-26 13:15:25 INFO None 5279749: status FINISHED
2026-07-26 13:15:25 INFO None 5279750: status FINISHED
2026-07-26 13:15:26 INFO None 5279751: status FINISHED
2026-07-26 13:15:26 INFO None 5279753: status FINISHED
2026-07-26 13:15:26 INFO None 5279754: status FINISHED
2026-07-26 13:15:26 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:15:26 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:15:26 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:15:26 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:15:26 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:15:26 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:15:26 INFO Jobs still running: ['5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:15:41 INFO None 5279744: status FINISHED
2026-07-26 13:15:41 INFO None 5279745: status FINISHED
2026-07-26 13:15:41 INFO None 5279747: status FINISHED
2026-07-26 13:15:42 INFO None 5279748: status FINISHED
2026-07-26 13:15:42 INFO None 5279749: status FINISHED
2026-07-26 13:15:42 INFO None 5279750: status FINISHED
2026-07-26 13:15:42 INFO None 5279751: status FINISHED
2026-07-26 13:15:42 INFO None 5279753: status FINISHED
2026-07-26 13:15:42 INFO None 5279754: status FINISHED
2026-07-26 13:15:42 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:15:42 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:15:42 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:15:42 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:15:42 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:15:42 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:15:42 INFO Jobs still running: ['5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:15:57 INFO None 5279744: status FINISHED
2026-07-26 13:15:57 INFO None 5279745: status FINISHED
2026-07-26 13:15:57 INFO None 5279747: status FINISHED
2026-07-26 13:15:57 INFO None 5279748: status FINISHED
2026-07-26 13:15:57 INFO None 5279749: status FINISHED
2026-07-26 13:15:57 INFO None 5279750: status FINISHED
2026-07-26 13:15:57 INFO None 5279751: status FINISHED
2026-07-26 13:15:57 INFO None 5279753: status FINISHED
2026-07-26 13:15:57 INFO None 5279754: status FINISHED
2026-07-26 13:15:57 INFO None 5279755: status RUNNING/PENDING
2026-07-26 13:15:57 INFO None 5279756: status RUNNING/PENDING
2026-07-26 13:15:57 INFO None 5279757: status RUNNING/PENDING
2026-07-26 13:15:57 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:15:59 INFO None 5279759: status RUNNING/PENDING
2026-07-26 13:15:59 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:15:59 INFO Jobs still running: ['5279755', '5279756', '5279757', '5279758', '5279759', '5279760']. Waiting...
2026-07-26 13:16:14 INFO None 5279744: status FINISHED
2026-07-26 13:16:14 INFO None 5279745: status FINISHED
2026-07-26 13:16:14 INFO None 5279747: status FINISHED
2026-07-26 13:16:14 INFO None 5279748: status FINISHED
2026-07-26 13:16:14 INFO None 5279749: status FINISHED
2026-07-26 13:16:14 INFO None 5279750: status FINISHED
2026-07-26 13:16:14 INFO None 5279751: status FINISHED
2026-07-26 13:16:14 INFO None 5279753: status FINISHED
2026-07-26 13:16:14 INFO None 5279754: status FINISHED
2026-07-26 13:16:14 INFO None 5279755: status FINISHED
2026-07-26 13:16:14 INFO None 5279756: status FINISHED
2026-07-26 13:16:14 INFO None 5279757: status FINISHED
2026-07-26 13:16:14 INFO None 5279758: status RUNNING/PENDING
2026-07-26 13:16:14 INFO None 5279759: status FINISHED
2026-07-26 13:16:14 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:16:14 INFO Jobs still running: ['5279758', '5279760']. Waiting...
2026-07-26 13:16:29 INFO None 5279744: status FINISHED
2026-07-26 13:16:29 INFO None 5279745: status FINISHED
2026-07-26 13:16:29 INFO None 5279747: status FINISHED
2026-07-26 13:16:29 INFO None 5279748: status FINISHED
2026-07-26 13:16:29 INFO None 5279749: status FINISHED
2026-07-26 13:16:29 INFO None 5279750: status FINISHED
2026-07-26 13:16:29 INFO None 5279751: status FINISHED
2026-07-26 13:16:30 INFO None 5279753: status FINISHED
2026-07-26 13:16:30 INFO None 5279754: status FINISHED
2026-07-26 13:16:30 INFO None 5279755: status FINISHED
2026-07-26 13:16:30 INFO None 5279756: status FINISHED
2026-07-26 13:16:30 INFO None 5279757: status FINISHED
2026-07-26 13:16:32 INFO None 5279758: status FINISHED
2026-07-26 13:16:32 INFO None 5279759: status FINISHED
2026-07-26 13:16:32 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:16:32 INFO Jobs still running: ['5279760']. Waiting...
2026-07-26 13:16:47 INFO None 5279744: status FINISHED
2026-07-26 13:16:47 INFO None 5279745: status FINISHED
2026-07-26 13:16:47 INFO None 5279747: status FINISHED
2026-07-26 13:16:47 INFO None 5279748: status FINISHED
2026-07-26 13:16:47 INFO None 5279749: status FINISHED
2026-07-26 13:16:47 INFO None 5279750: status FINISHED
2026-07-26 13:16:47 INFO None 5279751: status FINISHED
2026-07-26 13:16:47 INFO None 5279753: status FINISHED
2026-07-26 13:16:47 INFO None 5279754: status FINISHED
2026-07-26 13:16:47 INFO None 5279755: status FINISHED
2026-07-26 13:16:47 INFO None 5279756: status FINISHED
2026-07-26 13:16:47 INFO None 5279757: status FINISHED
2026-07-26 13:16:47 INFO None 5279758: status FINISHED
2026-07-26 13:16:47 INFO None 5279759: status FINISHED
2026-07-26 13:16:47 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:16:47 INFO Jobs still running: ['5279760']. Waiting...
2026-07-26 13:17:02 INFO None 5279744: status FINISHED
2026-07-26 13:17:02 INFO None 5279745: status FINISHED
2026-07-26 13:17:02 INFO None 5279747: status FINISHED
2026-07-26 13:17:02 INFO None 5279748: status FINISHED
2026-07-26 13:17:02 INFO None 5279749: status FINISHED
2026-07-26 13:17:02 INFO None 5279750: status FINISHED
2026-07-26 13:17:02 INFO None 5279751: status FINISHED
2026-07-26 13:17:02 INFO None 5279753: status FINISHED
2026-07-26 13:17:02 INFO None 5279754: status FINISHED
2026-07-26 13:17:02 INFO None 5279755: status FINISHED
2026-07-26 13:17:02 INFO None 5279756: status FINISHED
2026-07-26 13:17:02 INFO None 5279757: status FINISHED
2026-07-26 13:17:02 INFO None 5279758: status FINISHED
2026-07-26 13:17:02 INFO None 5279759: status FINISHED
2026-07-26 13:17:02 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:17:02 INFO Jobs still running: ['5279760']. Waiting...
2026-07-26 13:17:17 INFO None 5279744: status FINISHED
2026-07-26 13:17:17 INFO None 5279745: status FINISHED
2026-07-26 13:17:17 INFO None 5279747: status FINISHED
2026-07-26 13:17:17 INFO None 5279748: status FINISHED
2026-07-26 13:17:17 INFO None 5279749: status FINISHED
2026-07-26 13:17:17 INFO None 5279750: status FINISHED
2026-07-26 13:17:17 INFO None 5279751: status FINISHED
2026-07-26 13:17:17 INFO None 5279753: status FINISHED
2026-07-26 13:17:17 INFO None 5279754: status FINISHED
2026-07-26 13:17:17 INFO None 5279755: status FINISHED
2026-07-26 13:17:17 INFO None 5279756: status FINISHED
2026-07-26 13:17:17 INFO None 5279757: status FINISHED
2026-07-26 13:17:17 INFO None 5279758: status FINISHED
2026-07-26 13:17:17 INFO None 5279759: status FINISHED
2026-07-26 13:17:17 INFO None 5279760: status RUNNING/PENDING
2026-07-26 13:17:17 INFO Jobs still running: ['5279760']. Waiting...
2026-07-26 13:17:33 INFO None 5279744: status FINISHED
2026-07-26 13:17:33 INFO None 5279745: status FINISHED
2026-07-26 13:17:33 INFO None 5279747: status FINISHED
2026-07-26 13:17:33 INFO None 5279748: status FINISHED
2026-07-26 13:17:33 INFO None 5279749: status FINISHED
2026-07-26 13:17:33 INFO None 5279750: status FINISHED
2026-07-26 13:17:33 INFO None 5279751: status FINISHED
2026-07-26 13:17:33 INFO None 5279753: status FINISHED
2026-07-26 13:17:33 INFO None 5279754: status FINISHED
2026-07-26 13:17:33 INFO None 5279755: status FINISHED
2026-07-26 13:17:33 INFO None 5279756: status FINISHED
2026-07-26 13:17:33 INFO None 5279757: status FINISHED
2026-07-26 13:17:33 INFO None 5279758: status FINISHED
2026-07-26 13:17:33 INFO None 5279759: status FINISHED
2026-07-26 13:17:33 INFO None 5279760: status FINISHED
2026-07-26 13:17:33 INFO Jobs ['5279744', '5279745', '5279747', '5279748', '5279749', '5279750', '5279751', '5279753', '5279754', '5279755', '5279756', '5279757', '5279758', '5279759', '5279760'] have finished
2026-07-26 13:17:33 INFO Checking restart files were created ...
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-26 13:17:33 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-26 13:17:33 INFO  Run_model() completed successfully.
2026-07-26 13:17:33 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:17:33 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-26 13:17:33 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 13:17:33 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-26 13:17:33 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:17:33 INFO ---------->>> Running process_satellite_data()
2026-07-26 13:17:33 INFO [DART] No satellite data found, skipping assimilation
2026-07-26 13:17:33 INFO after_assimilation() skipped
2026-07-26 13:17:33 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 13:17:33 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:17:33 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:17:33 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-26 13:17:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:17:35 INFO Hourly dataset computed and listing created
2026-07-26 13:17:51 INFO Hourly dataset computed
2026-07-26 13:17:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:17:52 INFO Hourly dataset computed and listing created
2026-07-26 13:18:03 INFO Hourly dataset computed
2026-07-26 13:18:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:04 INFO Hourly dataset computed and listing created
2026-07-26 13:18:08 INFO Hourly dataset computed
2026-07-26 13:18:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:09 INFO Hourly dataset computed and listing created
2026-07-26 13:18:11 INFO Hourly dataset computed
2026-07-26 13:18:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:12 INFO Hourly dataset computed and listing created
2026-07-26 13:18:14 INFO Hourly dataset computed
2026-07-26 13:18:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:15 INFO Hourly dataset computed and listing created
2026-07-26 13:18:17 INFO Hourly dataset computed
2026-07-26 13:18:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:18 INFO Hourly dataset computed and listing created
2026-07-26 13:18:20 INFO Hourly dataset computed
2026-07-26 13:18:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:22 INFO Hourly dataset computed and listing created
2026-07-26 13:18:24 INFO Hourly dataset computed
2026-07-26 13:18:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:25 INFO Hourly dataset computed and listing created
2026-07-26 13:18:27 INFO Hourly dataset computed
2026-07-26 13:18:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:28 INFO Hourly dataset computed and listing created
2026-07-26 13:18:30 INFO Hourly dataset computed
2026-07-26 13:18:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:31 INFO Hourly dataset computed and listing created
2026-07-26 13:18:33 INFO Hourly dataset computed
2026-07-26 13:18:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:34 INFO Hourly dataset computed and listing created
2026-07-26 13:18:36 INFO Hourly dataset computed
2026-07-26 13:18:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:37 INFO Hourly dataset computed and listing created
2026-07-26 13:18:39 INFO Hourly dataset computed
2026-07-26 13:18:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:40 INFO Hourly dataset computed and listing created
2026-07-26 13:18:42 INFO Hourly dataset computed
2026-07-26 13:18:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:18:43 INFO Hourly dataset computed and listing created
2026-07-26 13:18:46 INFO Hourly dataset computed
2026-07-26 13:18:46 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-26 13:18:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:18:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 13:18:46 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-26 13:18:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 13:18:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:18:46 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 13:18:46 INFO Queuing job for member 1...
2026-07-26 13:18:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:18:46 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 13:18:48 INFO Found: ['5279772']
2026-07-26 13:18:53 INFO [TGCC-IRENE] Submitted job with ID:['5279772']
2026-07-26 13:18:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:18:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 13:18:53 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-26 13:18:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 13:18:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:18:53 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 13:18:53 INFO Queuing job for member 2...
2026-07-26 13:18:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:18:53 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 13:18:54 INFO Found: ['5279773']
2026-07-26 13:18:59 INFO [TGCC-IRENE] Submitted job with ID:['5279773']
2026-07-26 13:18:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:18:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 13:18:59 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-26 13:18:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 13:18:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:18:59 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 13:18:59 INFO Queuing job for member 3...
2026-07-26 13:18:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:18:59 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 13:19:00 INFO Found: ['5279775']
2026-07-26 13:19:05 INFO [TGCC-IRENE] Submitted job with ID:['5279775']
2026-07-26 13:19:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 13:19:05 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-26 13:19:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 13:19:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:05 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 13:19:05 INFO Queuing job for member 4...
2026-07-26 13:19:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:05 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 13:19:06 INFO Found: ['5279777']
2026-07-26 13:19:11 INFO [TGCC-IRENE] Submitted job with ID:['5279777']
2026-07-26 13:19:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 13:19:11 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-26 13:19:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 13:19:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:11 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 13:19:11 INFO Queuing job for member 5...
2026-07-26 13:19:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:11 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 13:19:11 INFO Found: ['5279778']
2026-07-26 13:19:16 INFO [TGCC-IRENE] Submitted job with ID:['5279778']
2026-07-26 13:19:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 13:19:16 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-26 13:19:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 13:19:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:16 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 13:19:16 INFO Queuing job for member 6...
2026-07-26 13:19:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:16 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 13:19:17 INFO Found: ['5279779']
2026-07-26 13:19:22 INFO [TGCC-IRENE] Submitted job with ID:['5279779']
2026-07-26 13:19:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 13:19:22 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-26 13:19:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 13:19:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:22 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 13:19:22 INFO Queuing job for member 7...
2026-07-26 13:19:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:22 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 13:19:23 INFO Found: ['5279780']
2026-07-26 13:19:28 INFO [TGCC-IRENE] Submitted job with ID:['5279780']
2026-07-26 13:19:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 13:19:28 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-26 13:19:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 13:19:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:28 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 13:19:28 INFO Queuing job for member 8...
2026-07-26 13:19:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:28 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 13:19:29 INFO Found: ['5279782']
2026-07-26 13:19:34 INFO [TGCC-IRENE] Submitted job with ID:['5279782']
2026-07-26 13:19:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 13:19:34 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-26 13:19:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 13:19:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:34 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 13:19:34 INFO Queuing job for member 9...
2026-07-26 13:19:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:34 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 13:19:34 INFO Found: ['5279783']
2026-07-26 13:19:39 INFO [TGCC-IRENE] Submitted job with ID:['5279783']
2026-07-26 13:19:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 13:19:39 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-26 13:19:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 13:19:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:39 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 13:19:39 INFO Queuing job for member 10...
2026-07-26 13:19:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:39 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 13:19:41 INFO Found: ['5279784']
2026-07-26 13:19:46 INFO [TGCC-IRENE] Submitted job with ID:['5279784']
2026-07-26 13:19:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 13:19:46 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-26 13:19:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 13:19:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:46 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 13:19:46 INFO Queuing job for member 11...
2026-07-26 13:19:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:46 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 13:19:49 INFO Found: ['5279785']
2026-07-26 13:19:54 INFO [TGCC-IRENE] Submitted job with ID:['5279785']
2026-07-26 13:19:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:19:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 13:19:54 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-26 13:19:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 13:19:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:19:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 13:19:54 INFO Queuing job for member 12...
2026-07-26 13:19:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:19:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 13:19:56 INFO Found: ['5279786']
2026-07-26 13:20:01 INFO [TGCC-IRENE] Submitted job with ID:['5279786']
2026-07-26 13:20:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:20:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 13:20:01 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-26 13:20:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 13:20:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:20:01 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 13:20:01 INFO Queuing job for member 13...
2026-07-26 13:20:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:20:01 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 13:20:04 INFO Found: ['5279789']
2026-07-26 13:20:09 INFO [TGCC-IRENE] Submitted job with ID:['5279789']
2026-07-26 13:20:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:20:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 13:20:09 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-26 13:20:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 13:20:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:20:09 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 13:20:09 INFO Queuing job for member 14...
2026-07-26 13:20:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:20:09 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 13:20:11 INFO Found: ['5279790']
2026-07-26 13:20:16 INFO [TGCC-IRENE] Submitted job with ID:['5279790']
2026-07-26 13:20:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:20:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 13:20:16 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-26 13:20:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 13:20:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:20:16 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 13:20:16 INFO Queuing job for member 15...
2026-07-26 13:20:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:20:16 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 13:20:19 INFO Found: ['5279791']
2026-07-26 13:20:24 INFO [TGCC-IRENE] Submitted job with ID:['5279791']
2026-07-26 13:20:24 INFO Checking job status ...
2026-07-26 13:20:24 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:20:24 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:20:24 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:20:39 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:20:39 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:20:39 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:20:39 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:20:39 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:20:39 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:20:41 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:20:41 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:20:56 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:20:56 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:20:57 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:20:57 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:21:12 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:21:12 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:21:12 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:21:27 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:21:27 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:21:27 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:21:43 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:21:43 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:21:43 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:21:58 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:21:58 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:22:00 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:22:00 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:22:00 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:22:15 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:22:15 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:22:15 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:22:15 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:22:16 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:22:16 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:22:31 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:22:31 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:22:33 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:22:33 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:22:48 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:22:48 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:22:48 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:23:03 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:23:03 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:23:04 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:23:04 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:23:04 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:23:04 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:23:04 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:23:04 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:23:19 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:23:19 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:23:19 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:23:34 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:23:34 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:23:34 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:23:50 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:23:50 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:23:50 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:24:05 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:24:05 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:24:07 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:24:07 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:24:22 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:24:22 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:24:22 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:24:37 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:24:37 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:24:37 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:24:37 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:24:38 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:24:38 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:24:53 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:24:53 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:24:55 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:24:55 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:25:10 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:25:10 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:25:10 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:25:25 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:25:25 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:25:26 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:25:26 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:25:26 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:25:26 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:25:26 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:25:41 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:25:41 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:25:41 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:25:56 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:25:56 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:25:56 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:26:11 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:26:11 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:26:11 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:26:11 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:26:11 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:26:11 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:26:12 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:26:12 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:26:27 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:26:27 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:26:27 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:26:42 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:26:42 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:26:44 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:26:44 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:26:59 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:26:59 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:27:00 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:27:00 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:27:00 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:27:00 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:27:15 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:27:15 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:27:17 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:27:17 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:27:32 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:27:32 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:27:32 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:27:47 INFO None 5279772: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279775: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:27:47 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:27:47 INFO Jobs still running: ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:28:02 INFO None 5279772: status FINISHED
2026-07-26 13:28:03 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279775: status FINISHED
2026-07-26 13:28:03 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:28:03 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:28:03 INFO Jobs still running: ['5279773', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:28:18 INFO None 5279772: status FINISHED
2026-07-26 13:28:18 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279775: status FINISHED
2026-07-26 13:28:18 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279778: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279783: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:28:18 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:28:18 INFO Jobs still running: ['5279773', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:28:33 INFO None 5279772: status FINISHED
2026-07-26 13:28:33 INFO None 5279773: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279775: status FINISHED
2026-07-26 13:28:34 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279778: status FINISHED
2026-07-26 13:28:34 INFO None 5279779: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279780: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279782: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279783: status FINISHED
2026-07-26 13:28:34 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:28:34 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:28:34 INFO Jobs still running: ['5279773', '5279777', '5279779', '5279780', '5279782', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:28:49 INFO None 5279772: status FINISHED
2026-07-26 13:28:49 INFO None 5279773: status FINISHED
2026-07-26 13:28:49 INFO None 5279775: status FINISHED
2026-07-26 13:28:49 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:28:49 INFO None 5279778: status FINISHED
2026-07-26 13:28:51 INFO None 5279779: status FINISHED
2026-07-26 13:28:51 INFO None 5279780: status FINISHED
2026-07-26 13:28:51 INFO None 5279782: status FINISHED
2026-07-26 13:28:51 INFO None 5279783: status FINISHED
2026-07-26 13:28:51 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:28:51 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:28:51 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:28:51 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:28:51 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:28:51 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:28:51 INFO Jobs still running: ['5279777', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:29:06 INFO None 5279772: status FINISHED
2026-07-26 13:29:06 INFO None 5279773: status FINISHED
2026-07-26 13:29:06 INFO None 5279775: status FINISHED
2026-07-26 13:29:06 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279778: status FINISHED
2026-07-26 13:29:06 INFO None 5279779: status FINISHED
2026-07-26 13:29:06 INFO None 5279780: status FINISHED
2026-07-26 13:29:06 INFO None 5279782: status FINISHED
2026-07-26 13:29:06 INFO None 5279783: status FINISHED
2026-07-26 13:29:06 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:29:06 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:29:06 INFO Jobs still running: ['5279777', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:29:21 INFO None 5279772: status FINISHED
2026-07-26 13:29:21 INFO None 5279773: status FINISHED
2026-07-26 13:29:21 INFO None 5279775: status FINISHED
2026-07-26 13:29:23 INFO None 5279777: status RUNNING/PENDING
2026-07-26 13:29:23 INFO None 5279778: status FINISHED
2026-07-26 13:29:24 INFO None 5279779: status FINISHED
2026-07-26 13:29:24 INFO None 5279780: status FINISHED
2026-07-26 13:29:24 INFO None 5279782: status FINISHED
2026-07-26 13:29:24 INFO None 5279783: status FINISHED
2026-07-26 13:29:24 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:29:24 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:29:24 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:29:24 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:29:24 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:29:24 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:29:24 INFO Jobs still running: ['5279777', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:29:39 INFO None 5279772: status FINISHED
2026-07-26 13:29:39 INFO None 5279773: status FINISHED
2026-07-26 13:29:39 INFO None 5279775: status FINISHED
2026-07-26 13:29:39 INFO None 5279777: status FINISHED
2026-07-26 13:29:39 INFO None 5279778: status FINISHED
2026-07-26 13:29:39 INFO None 5279779: status FINISHED
2026-07-26 13:29:39 INFO None 5279780: status FINISHED
2026-07-26 13:29:39 INFO None 5279782: status FINISHED
2026-07-26 13:29:39 INFO None 5279783: status FINISHED
2026-07-26 13:29:39 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:29:39 INFO None 5279785: status RUNNING/PENDING
2026-07-26 13:29:39 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:29:39 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:29:39 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:29:39 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:29:39 INFO Jobs still running: ['5279784', '5279785', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:29:54 INFO None 5279772: status FINISHED
2026-07-26 13:29:54 INFO None 5279773: status FINISHED
2026-07-26 13:29:54 INFO None 5279775: status FINISHED
2026-07-26 13:29:54 INFO None 5279777: status FINISHED
2026-07-26 13:29:54 INFO None 5279778: status FINISHED
2026-07-26 13:29:54 INFO None 5279779: status FINISHED
2026-07-26 13:29:54 INFO None 5279780: status FINISHED
2026-07-26 13:29:54 INFO None 5279782: status FINISHED
2026-07-26 13:29:54 INFO None 5279783: status FINISHED
2026-07-26 13:29:54 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:29:54 INFO None 5279785: status FINISHED
2026-07-26 13:29:54 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:29:54 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:29:54 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:29:54 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:29:54 INFO Jobs still running: ['5279784', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:30:11 INFO None 5279772: status FINISHED
2026-07-26 13:30:11 INFO None 5279773: status FINISHED
2026-07-26 13:30:11 INFO None 5279775: status FINISHED
2026-07-26 13:30:11 INFO None 5279777: status FINISHED
2026-07-26 13:30:11 INFO None 5279778: status FINISHED
2026-07-26 13:30:11 INFO None 5279779: status FINISHED
2026-07-26 13:30:11 INFO None 5279780: status FINISHED
2026-07-26 13:30:11 INFO None 5279782: status FINISHED
2026-07-26 13:30:11 INFO None 5279783: status FINISHED
2026-07-26 13:30:11 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:30:11 INFO None 5279785: status FINISHED
2026-07-26 13:30:11 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:30:11 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:30:11 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:30:11 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:30:11 INFO Jobs still running: ['5279784', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:30:26 INFO None 5279772: status FINISHED
2026-07-26 13:30:26 INFO None 5279773: status FINISHED
2026-07-26 13:30:26 INFO None 5279775: status FINISHED
2026-07-26 13:30:26 INFO None 5279777: status FINISHED
2026-07-26 13:30:26 INFO None 5279778: status FINISHED
2026-07-26 13:30:26 INFO None 5279779: status FINISHED
2026-07-26 13:30:26 INFO None 5279780: status FINISHED
2026-07-26 13:30:27 INFO None 5279782: status FINISHED
2026-07-26 13:30:27 INFO None 5279783: status FINISHED
2026-07-26 13:30:27 INFO None 5279784: status RUNNING/PENDING
2026-07-26 13:30:27 INFO None 5279785: status FINISHED
2026-07-26 13:30:27 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:30:27 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:30:27 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:30:27 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:30:27 INFO Jobs still running: ['5279784', '5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:30:42 INFO None 5279772: status FINISHED
2026-07-26 13:30:42 INFO None 5279773: status FINISHED
2026-07-26 13:30:44 INFO None 5279775: status FINISHED
2026-07-26 13:30:44 INFO None 5279777: status FINISHED
2026-07-26 13:30:44 INFO None 5279778: status FINISHED
2026-07-26 13:30:44 INFO None 5279779: status FINISHED
2026-07-26 13:30:44 INFO None 5279780: status FINISHED
2026-07-26 13:30:44 INFO None 5279782: status FINISHED
2026-07-26 13:30:44 INFO None 5279783: status FINISHED
2026-07-26 13:30:44 INFO None 5279784: status FINISHED
2026-07-26 13:30:44 INFO None 5279785: status FINISHED
2026-07-26 13:30:44 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:30:44 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:30:44 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:30:44 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:30:44 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:30:59 INFO None 5279772: status FINISHED
2026-07-26 13:30:59 INFO None 5279773: status FINISHED
2026-07-26 13:30:59 INFO None 5279775: status FINISHED
2026-07-26 13:30:59 INFO None 5279777: status FINISHED
2026-07-26 13:30:59 INFO None 5279778: status FINISHED
2026-07-26 13:30:59 INFO None 5279779: status FINISHED
2026-07-26 13:30:59 INFO None 5279780: status FINISHED
2026-07-26 13:30:59 INFO None 5279782: status FINISHED
2026-07-26 13:30:59 INFO None 5279783: status FINISHED
2026-07-26 13:30:59 INFO None 5279784: status FINISHED
2026-07-26 13:30:59 INFO None 5279785: status FINISHED
2026-07-26 13:30:59 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:30:59 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:30:59 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:30:59 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:30:59 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:31:14 INFO None 5279772: status FINISHED
2026-07-26 13:31:14 INFO None 5279773: status FINISHED
2026-07-26 13:31:14 INFO None 5279775: status FINISHED
2026-07-26 13:31:14 INFO None 5279777: status FINISHED
2026-07-26 13:31:16 INFO None 5279778: status FINISHED
2026-07-26 13:31:16 INFO None 5279779: status FINISHED
2026-07-26 13:31:16 INFO None 5279780: status FINISHED
2026-07-26 13:31:16 INFO None 5279782: status FINISHED
2026-07-26 13:31:16 INFO None 5279783: status FINISHED
2026-07-26 13:31:16 INFO None 5279784: status FINISHED
2026-07-26 13:31:16 INFO None 5279785: status FINISHED
2026-07-26 13:31:16 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:31:16 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:31:16 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:31:17 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:31:17 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:31:32 INFO None 5279772: status FINISHED
2026-07-26 13:31:32 INFO None 5279773: status FINISHED
2026-07-26 13:31:32 INFO None 5279775: status FINISHED
2026-07-26 13:31:32 INFO None 5279777: status FINISHED
2026-07-26 13:31:32 INFO None 5279778: status FINISHED
2026-07-26 13:31:32 INFO None 5279779: status FINISHED
2026-07-26 13:31:32 INFO None 5279780: status FINISHED
2026-07-26 13:31:32 INFO None 5279782: status FINISHED
2026-07-26 13:31:32 INFO None 5279783: status FINISHED
2026-07-26 13:31:32 INFO None 5279784: status FINISHED
2026-07-26 13:31:32 INFO None 5279785: status FINISHED
2026-07-26 13:31:32 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:31:32 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:31:32 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:31:32 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:31:32 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:31:47 INFO None 5279772: status FINISHED
2026-07-26 13:31:47 INFO None 5279773: status FINISHED
2026-07-26 13:31:47 INFO None 5279775: status FINISHED
2026-07-26 13:31:47 INFO None 5279777: status FINISHED
2026-07-26 13:31:47 INFO None 5279778: status FINISHED
2026-07-26 13:31:47 INFO None 5279779: status FINISHED
2026-07-26 13:31:47 INFO None 5279780: status FINISHED
2026-07-26 13:31:47 INFO None 5279782: status FINISHED
2026-07-26 13:31:47 INFO None 5279783: status FINISHED
2026-07-26 13:31:47 INFO None 5279784: status FINISHED
2026-07-26 13:31:47 INFO None 5279785: status FINISHED
2026-07-26 13:31:47 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:31:47 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:31:47 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:31:47 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:31:47 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:32:02 INFO None 5279772: status FINISHED
2026-07-26 13:32:02 INFO None 5279773: status FINISHED
2026-07-26 13:32:02 INFO None 5279775: status FINISHED
2026-07-26 13:32:02 INFO None 5279777: status FINISHED
2026-07-26 13:32:02 INFO None 5279778: status FINISHED
2026-07-26 13:32:02 INFO None 5279779: status FINISHED
2026-07-26 13:32:02 INFO None 5279780: status FINISHED
2026-07-26 13:32:02 INFO None 5279782: status FINISHED
2026-07-26 13:32:02 INFO None 5279783: status FINISHED
2026-07-26 13:32:02 INFO None 5279784: status FINISHED
2026-07-26 13:32:02 INFO None 5279785: status FINISHED
2026-07-26 13:32:02 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:32:02 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:32:02 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:32:02 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:32:02 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:32:17 INFO None 5279772: status FINISHED
2026-07-26 13:32:17 INFO None 5279773: status FINISHED
2026-07-26 13:32:17 INFO None 5279775: status FINISHED
2026-07-26 13:32:17 INFO None 5279777: status FINISHED
2026-07-26 13:32:17 INFO None 5279778: status FINISHED
2026-07-26 13:32:17 INFO None 5279779: status FINISHED
2026-07-26 13:32:18 INFO None 5279780: status FINISHED
2026-07-26 13:32:18 INFO None 5279782: status FINISHED
2026-07-26 13:32:18 INFO None 5279783: status FINISHED
2026-07-26 13:32:18 INFO None 5279784: status FINISHED
2026-07-26 13:32:18 INFO None 5279785: status FINISHED
2026-07-26 13:32:18 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:32:18 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:32:18 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:32:18 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:32:18 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:32:34 INFO None 5279772: status FINISHED
2026-07-26 13:32:34 INFO None 5279773: status FINISHED
2026-07-26 13:32:34 INFO None 5279775: status FINISHED
2026-07-26 13:32:34 INFO None 5279777: status FINISHED
2026-07-26 13:32:34 INFO None 5279778: status FINISHED
2026-07-26 13:32:34 INFO None 5279779: status FINISHED
2026-07-26 13:32:34 INFO None 5279780: status FINISHED
2026-07-26 13:32:34 INFO None 5279782: status FINISHED
2026-07-26 13:32:34 INFO None 5279783: status FINISHED
2026-07-26 13:32:34 INFO None 5279784: status FINISHED
2026-07-26 13:32:34 INFO None 5279785: status FINISHED
2026-07-26 13:32:34 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:32:34 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:32:34 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:32:34 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:32:34 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:32:49 INFO None 5279772: status FINISHED
2026-07-26 13:32:49 INFO None 5279773: status FINISHED
2026-07-26 13:32:49 INFO None 5279775: status FINISHED
2026-07-26 13:32:49 INFO None 5279777: status FINISHED
2026-07-26 13:32:49 INFO None 5279778: status FINISHED
2026-07-26 13:32:49 INFO None 5279779: status FINISHED
2026-07-26 13:32:49 INFO None 5279780: status FINISHED
2026-07-26 13:32:49 INFO None 5279782: status FINISHED
2026-07-26 13:32:49 INFO None 5279783: status FINISHED
2026-07-26 13:32:50 INFO None 5279784: status FINISHED
2026-07-26 13:32:50 INFO None 5279785: status FINISHED
2026-07-26 13:32:50 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:32:50 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:32:50 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:32:52 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:32:52 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:33:07 INFO None 5279772: status FINISHED
2026-07-26 13:33:07 INFO None 5279773: status FINISHED
2026-07-26 13:33:07 INFO None 5279775: status FINISHED
2026-07-26 13:33:07 INFO None 5279777: status FINISHED
2026-07-26 13:33:07 INFO None 5279778: status FINISHED
2026-07-26 13:33:07 INFO None 5279779: status FINISHED
2026-07-26 13:33:07 INFO None 5279780: status FINISHED
2026-07-26 13:33:07 INFO None 5279782: status FINISHED
2026-07-26 13:33:07 INFO None 5279783: status FINISHED
2026-07-26 13:33:07 INFO None 5279784: status FINISHED
2026-07-26 13:33:07 INFO None 5279785: status FINISHED
2026-07-26 13:33:07 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:33:07 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:33:07 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:33:07 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:33:07 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:33:22 INFO None 5279772: status FINISHED
2026-07-26 13:33:22 INFO None 5279773: status FINISHED
2026-07-26 13:33:22 INFO None 5279775: status FINISHED
2026-07-26 13:33:22 INFO None 5279777: status FINISHED
2026-07-26 13:33:22 INFO None 5279778: status FINISHED
2026-07-26 13:33:22 INFO None 5279779: status FINISHED
2026-07-26 13:33:22 INFO None 5279780: status FINISHED
2026-07-26 13:33:22 INFO None 5279782: status FINISHED
2026-07-26 13:33:22 INFO None 5279783: status FINISHED
2026-07-26 13:33:22 INFO None 5279784: status FINISHED
2026-07-26 13:33:22 INFO None 5279785: status FINISHED
2026-07-26 13:33:22 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:33:22 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:33:22 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:33:22 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:33:22 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:33:37 INFO None 5279772: status FINISHED
2026-07-26 13:33:37 INFO None 5279773: status FINISHED
2026-07-26 13:33:37 INFO None 5279775: status FINISHED
2026-07-26 13:33:39 INFO None 5279777: status FINISHED
2026-07-26 13:33:39 INFO None 5279778: status FINISHED
2026-07-26 13:33:39 INFO None 5279779: status FINISHED
2026-07-26 13:33:39 INFO None 5279780: status FINISHED
2026-07-26 13:33:39 INFO None 5279782: status FINISHED
2026-07-26 13:33:39 INFO None 5279783: status FINISHED
2026-07-26 13:33:39 INFO None 5279784: status FINISHED
2026-07-26 13:33:39 INFO None 5279785: status FINISHED
2026-07-26 13:33:39 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:33:39 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:33:39 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:33:39 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:33:39 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:33:54 INFO None 5279772: status FINISHED
2026-07-26 13:33:54 INFO None 5279773: status FINISHED
2026-07-26 13:33:55 INFO None 5279775: status FINISHED
2026-07-26 13:33:55 INFO None 5279777: status FINISHED
2026-07-26 13:33:55 INFO None 5279778: status FINISHED
2026-07-26 13:33:55 INFO None 5279779: status FINISHED
2026-07-26 13:33:55 INFO None 5279780: status FINISHED
2026-07-26 13:33:55 INFO None 5279782: status FINISHED
2026-07-26 13:33:55 INFO None 5279783: status FINISHED
2026-07-26 13:33:55 INFO None 5279784: status FINISHED
2026-07-26 13:33:55 INFO None 5279785: status FINISHED
2026-07-26 13:33:55 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:33:55 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:33:55 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:33:55 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:33:55 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:34:10 INFO None 5279772: status FINISHED
2026-07-26 13:34:10 INFO None 5279773: status FINISHED
2026-07-26 13:34:10 INFO None 5279775: status FINISHED
2026-07-26 13:34:10 INFO None 5279777: status FINISHED
2026-07-26 13:34:10 INFO None 5279778: status FINISHED
2026-07-26 13:34:10 INFO None 5279779: status FINISHED
2026-07-26 13:34:10 INFO None 5279780: status FINISHED
2026-07-26 13:34:10 INFO None 5279782: status FINISHED
2026-07-26 13:34:10 INFO None 5279783: status FINISHED
2026-07-26 13:34:10 INFO None 5279784: status FINISHED
2026-07-26 13:34:10 INFO None 5279785: status FINISHED
2026-07-26 13:34:10 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:34:10 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:34:10 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:34:10 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:34:10 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:34:25 INFO None 5279772: status FINISHED
2026-07-26 13:34:25 INFO None 5279773: status FINISHED
2026-07-26 13:34:25 INFO None 5279775: status FINISHED
2026-07-26 13:34:25 INFO None 5279777: status FINISHED
2026-07-26 13:34:25 INFO None 5279778: status FINISHED
2026-07-26 13:34:25 INFO None 5279779: status FINISHED
2026-07-26 13:34:25 INFO None 5279780: status FINISHED
2026-07-26 13:34:25 INFO None 5279782: status FINISHED
2026-07-26 13:34:25 INFO None 5279783: status FINISHED
2026-07-26 13:34:25 INFO None 5279784: status FINISHED
2026-07-26 13:34:25 INFO None 5279785: status FINISHED
2026-07-26 13:34:25 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:34:25 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:34:25 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:34:25 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:34:25 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:34:40 INFO None 5279772: status FINISHED
2026-07-26 13:34:40 INFO None 5279773: status FINISHED
2026-07-26 13:34:40 INFO None 5279775: status FINISHED
2026-07-26 13:34:40 INFO None 5279777: status FINISHED
2026-07-26 13:34:40 INFO None 5279778: status FINISHED
2026-07-26 13:34:40 INFO None 5279779: status FINISHED
2026-07-26 13:34:40 INFO None 5279780: status FINISHED
2026-07-26 13:34:40 INFO None 5279782: status FINISHED
2026-07-26 13:34:40 INFO None 5279783: status FINISHED
2026-07-26 13:34:40 INFO None 5279784: status FINISHED
2026-07-26 13:34:41 INFO None 5279785: status FINISHED
2026-07-26 13:34:41 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:34:41 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:34:41 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:34:41 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:34:41 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:34:56 INFO None 5279772: status FINISHED
2026-07-26 13:34:56 INFO None 5279773: status FINISHED
2026-07-26 13:34:56 INFO None 5279775: status FINISHED
2026-07-26 13:34:56 INFO None 5279777: status FINISHED
2026-07-26 13:34:56 INFO None 5279778: status FINISHED
2026-07-26 13:34:56 INFO None 5279779: status FINISHED
2026-07-26 13:34:56 INFO None 5279780: status FINISHED
2026-07-26 13:34:56 INFO None 5279782: status FINISHED
2026-07-26 13:34:56 INFO None 5279783: status FINISHED
2026-07-26 13:34:56 INFO None 5279784: status FINISHED
2026-07-26 13:34:56 INFO None 5279785: status FINISHED
2026-07-26 13:34:56 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:34:58 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:34:58 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:34:58 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:34:58 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:35:13 INFO None 5279772: status FINISHED
2026-07-26 13:35:13 INFO None 5279773: status FINISHED
2026-07-26 13:35:13 INFO None 5279775: status FINISHED
2026-07-26 13:35:13 INFO None 5279777: status FINISHED
2026-07-26 13:35:13 INFO None 5279778: status FINISHED
2026-07-26 13:35:13 INFO None 5279779: status FINISHED
2026-07-26 13:35:13 INFO None 5279780: status FINISHED
2026-07-26 13:35:13 INFO None 5279782: status FINISHED
2026-07-26 13:35:13 INFO None 5279783: status FINISHED
2026-07-26 13:35:13 INFO None 5279784: status FINISHED
2026-07-26 13:35:13 INFO None 5279785: status FINISHED
2026-07-26 13:35:13 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:35:13 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:35:13 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:35:13 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:35:13 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:35:28 INFO None 5279772: status FINISHED
2026-07-26 13:35:28 INFO None 5279773: status FINISHED
2026-07-26 13:35:28 INFO None 5279775: status FINISHED
2026-07-26 13:35:28 INFO None 5279777: status FINISHED
2026-07-26 13:35:28 INFO None 5279778: status FINISHED
2026-07-26 13:35:28 INFO None 5279779: status FINISHED
2026-07-26 13:35:28 INFO None 5279780: status FINISHED
2026-07-26 13:35:28 INFO None 5279782: status FINISHED
2026-07-26 13:35:28 INFO None 5279783: status FINISHED
2026-07-26 13:35:28 INFO None 5279784: status FINISHED
2026-07-26 13:35:28 INFO None 5279785: status FINISHED
2026-07-26 13:35:28 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:35:28 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:35:28 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:35:30 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:35:30 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:35:46 INFO None 5279772: status FINISHED
2026-07-26 13:35:46 INFO None 5279773: status FINISHED
2026-07-26 13:35:46 INFO None 5279775: status FINISHED
2026-07-26 13:35:46 INFO None 5279777: status FINISHED
2026-07-26 13:35:46 INFO None 5279778: status FINISHED
2026-07-26 13:35:46 INFO None 5279779: status FINISHED
2026-07-26 13:35:46 INFO None 5279780: status FINISHED
2026-07-26 13:35:46 INFO None 5279782: status FINISHED
2026-07-26 13:35:46 INFO None 5279783: status FINISHED
2026-07-26 13:35:46 INFO None 5279784: status FINISHED
2026-07-26 13:35:46 INFO None 5279785: status FINISHED
2026-07-26 13:35:46 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:35:46 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:35:46 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:35:46 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:35:46 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:36:01 INFO None 5279772: status FINISHED
2026-07-26 13:36:01 INFO None 5279773: status FINISHED
2026-07-26 13:36:01 INFO None 5279775: status FINISHED
2026-07-26 13:36:01 INFO None 5279777: status FINISHED
2026-07-26 13:36:01 INFO None 5279778: status FINISHED
2026-07-26 13:36:01 INFO None 5279779: status FINISHED
2026-07-26 13:36:01 INFO None 5279780: status FINISHED
2026-07-26 13:36:01 INFO None 5279782: status FINISHED
2026-07-26 13:36:01 INFO None 5279783: status FINISHED
2026-07-26 13:36:01 INFO None 5279784: status FINISHED
2026-07-26 13:36:01 INFO None 5279785: status FINISHED
2026-07-26 13:36:01 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:36:01 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:36:01 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:36:01 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:36:01 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:36:16 INFO None 5279772: status FINISHED
2026-07-26 13:36:16 INFO None 5279773: status FINISHED
2026-07-26 13:36:16 INFO None 5279775: status FINISHED
2026-07-26 13:36:16 INFO None 5279777: status FINISHED
2026-07-26 13:36:16 INFO None 5279778: status FINISHED
2026-07-26 13:36:16 INFO None 5279779: status FINISHED
2026-07-26 13:36:16 INFO None 5279780: status FINISHED
2026-07-26 13:36:16 INFO None 5279782: status FINISHED
2026-07-26 13:36:16 INFO None 5279783: status FINISHED
2026-07-26 13:36:16 INFO None 5279784: status FINISHED
2026-07-26 13:36:16 INFO None 5279785: status FINISHED
2026-07-26 13:36:16 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:36:16 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:36:16 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:36:16 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:36:16 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:36:31 INFO None 5279772: status FINISHED
2026-07-26 13:36:31 INFO None 5279773: status FINISHED
2026-07-26 13:36:31 INFO None 5279775: status FINISHED
2026-07-26 13:36:31 INFO None 5279777: status FINISHED
2026-07-26 13:36:32 INFO None 5279778: status FINISHED
2026-07-26 13:36:32 INFO None 5279779: status FINISHED
2026-07-26 13:36:32 INFO None 5279780: status FINISHED
2026-07-26 13:36:32 INFO None 5279782: status FINISHED
2026-07-26 13:36:32 INFO None 5279783: status FINISHED
2026-07-26 13:36:32 INFO None 5279784: status FINISHED
2026-07-26 13:36:32 INFO None 5279785: status FINISHED
2026-07-26 13:36:32 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:36:32 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:36:32 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:36:32 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:36:32 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:36:47 INFO None 5279772: status FINISHED
2026-07-26 13:36:47 INFO None 5279773: status FINISHED
2026-07-26 13:36:47 INFO None 5279775: status FINISHED
2026-07-26 13:36:47 INFO None 5279777: status FINISHED
2026-07-26 13:36:47 INFO None 5279778: status FINISHED
2026-07-26 13:36:47 INFO None 5279779: status FINISHED
2026-07-26 13:36:47 INFO None 5279780: status FINISHED
2026-07-26 13:36:47 INFO None 5279782: status FINISHED
2026-07-26 13:36:49 INFO None 5279783: status FINISHED
2026-07-26 13:36:49 INFO None 5279784: status FINISHED
2026-07-26 13:36:49 INFO None 5279785: status FINISHED
2026-07-26 13:36:49 INFO None 5279786: status RUNNING/PENDING
2026-07-26 13:36:49 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:36:49 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:36:49 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:36:49 INFO Jobs still running: ['5279786', '5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:37:04 INFO None 5279772: status FINISHED
2026-07-26 13:37:04 INFO None 5279773: status FINISHED
2026-07-26 13:37:04 INFO None 5279775: status FINISHED
2026-07-26 13:37:04 INFO None 5279777: status FINISHED
2026-07-26 13:37:04 INFO None 5279778: status FINISHED
2026-07-26 13:37:04 INFO None 5279779: status FINISHED
2026-07-26 13:37:04 INFO None 5279780: status FINISHED
2026-07-26 13:37:04 INFO None 5279782: status FINISHED
2026-07-26 13:37:04 INFO None 5279783: status FINISHED
2026-07-26 13:37:04 INFO None 5279784: status FINISHED
2026-07-26 13:37:04 INFO None 5279785: status FINISHED
2026-07-26 13:37:04 INFO None 5279786: status FINISHED
2026-07-26 13:37:04 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:37:04 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:37:04 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:37:04 INFO Jobs still running: ['5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:37:19 INFO None 5279772: status FINISHED
2026-07-26 13:37:19 INFO None 5279773: status FINISHED
2026-07-26 13:37:19 INFO None 5279775: status FINISHED
2026-07-26 13:37:19 INFO None 5279777: status FINISHED
2026-07-26 13:37:19 INFO None 5279778: status FINISHED
2026-07-26 13:37:19 INFO None 5279779: status FINISHED
2026-07-26 13:37:19 INFO None 5279780: status FINISHED
2026-07-26 13:37:19 INFO None 5279782: status FINISHED
2026-07-26 13:37:19 INFO None 5279783: status FINISHED
2026-07-26 13:37:21 INFO None 5279784: status FINISHED
2026-07-26 13:37:22 INFO None 5279785: status FINISHED
2026-07-26 13:37:22 INFO None 5279786: status FINISHED
2026-07-26 13:37:22 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:37:22 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:37:22 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:37:22 INFO Jobs still running: ['5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:37:37 INFO None 5279772: status FINISHED
2026-07-26 13:37:37 INFO None 5279773: status FINISHED
2026-07-26 13:37:37 INFO None 5279775: status FINISHED
2026-07-26 13:37:37 INFO None 5279777: status FINISHED
2026-07-26 13:37:37 INFO None 5279778: status FINISHED
2026-07-26 13:37:37 INFO None 5279779: status FINISHED
2026-07-26 13:37:37 INFO None 5279780: status FINISHED
2026-07-26 13:37:37 INFO None 5279782: status FINISHED
2026-07-26 13:37:37 INFO None 5279783: status FINISHED
2026-07-26 13:37:37 INFO None 5279784: status FINISHED
2026-07-26 13:37:37 INFO None 5279785: status FINISHED
2026-07-26 13:37:37 INFO None 5279786: status FINISHED
2026-07-26 13:37:37 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:37:37 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:37:37 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:37:37 INFO Jobs still running: ['5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:37:52 INFO None 5279772: status FINISHED
2026-07-26 13:37:52 INFO None 5279773: status FINISHED
2026-07-26 13:37:52 INFO None 5279775: status FINISHED
2026-07-26 13:37:52 INFO None 5279777: status FINISHED
2026-07-26 13:37:52 INFO None 5279778: status FINISHED
2026-07-26 13:37:52 INFO None 5279779: status FINISHED
2026-07-26 13:37:52 INFO None 5279780: status FINISHED
2026-07-26 13:37:52 INFO None 5279782: status FINISHED
2026-07-26 13:37:52 INFO None 5279783: status FINISHED
2026-07-26 13:37:52 INFO None 5279784: status FINISHED
2026-07-26 13:37:52 INFO None 5279785: status FINISHED
2026-07-26 13:37:52 INFO None 5279786: status FINISHED
2026-07-26 13:37:52 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:37:54 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:37:54 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:37:54 INFO Jobs still running: ['5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:38:09 INFO None 5279772: status FINISHED (not in squeue)
2026-07-26 13:38:09 INFO None 5279773: status FINISHED
2026-07-26 13:38:09 INFO None 5279775: status FINISHED (not in squeue)
2026-07-26 13:38:09 INFO None 5279777: status FINISHED (not in squeue)
2026-07-26 13:38:09 INFO None 5279778: status FINISHED (not in squeue)
2026-07-26 13:38:09 INFO None 5279779: status FINISHED (not in squeue)
2026-07-26 13:38:09 INFO None 5279780: status FINISHED
2026-07-26 13:38:09 INFO None 5279782: status FINISHED
2026-07-26 13:38:09 INFO None 5279783: status FINISHED
2026-07-26 13:38:09 INFO None 5279784: status FINISHED
2026-07-26 13:38:09 INFO None 5279785: status FINISHED
2026-07-26 13:38:09 INFO None 5279786: status FINISHED
2026-07-26 13:38:09 INFO None 5279789: status RUNNING/PENDING
2026-07-26 13:38:09 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:38:09 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:38:09 INFO Jobs still running: ['5279789', '5279790', '5279791']. Waiting...
2026-07-26 13:38:24 INFO None 5279772: status FINISHED (not in squeue)
2026-07-26 13:38:24 INFO None 5279773: status FINISHED
2026-07-26 13:38:24 INFO None 5279775: status FINISHED (not in squeue)
2026-07-26 13:38:24 INFO None 5279777: status FINISHED (not in squeue)
2026-07-26 13:38:24 INFO None 5279778: status FINISHED (not in squeue)
2026-07-26 13:38:24 INFO None 5279779: status FINISHED (not in squeue)
2026-07-26 13:38:25 INFO None 5279780: status FINISHED
2026-07-26 13:38:25 INFO None 5279782: status FINISHED
2026-07-26 13:38:25 INFO None 5279783: status FINISHED
2026-07-26 13:38:25 INFO None 5279784: status FINISHED
2026-07-26 13:38:25 INFO None 5279785: status FINISHED
2026-07-26 13:38:25 INFO None 5279786: status FINISHED
2026-07-26 13:38:25 INFO None 5279789: status FINISHED
2026-07-26 13:38:25 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:38:25 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:38:25 INFO Jobs still running: ['5279790', '5279791']. Waiting...
2026-07-26 13:38:40 INFO None 5279772: status FINISHED (not in squeue)
2026-07-26 13:38:40 INFO None 5279773: status FINISHED
2026-07-26 13:38:40 INFO None 5279775: status FINISHED (not in squeue)
2026-07-26 13:38:40 INFO None 5279777: status FINISHED (not in squeue)
2026-07-26 13:38:40 INFO None 5279778: status FINISHED (not in squeue)
2026-07-26 13:38:40 INFO None 5279779: status FINISHED (not in squeue)
2026-07-26 13:38:40 INFO None 5279780: status FINISHED
2026-07-26 13:38:40 INFO None 5279782: status FINISHED
2026-07-26 13:38:40 INFO None 5279783: status FINISHED
2026-07-26 13:38:40 INFO None 5279784: status FINISHED
2026-07-26 13:38:40 INFO None 5279785: status FINISHED
2026-07-26 13:38:40 INFO None 5279786: status FINISHED
2026-07-26 13:38:40 INFO None 5279789: status FINISHED
2026-07-26 13:38:40 INFO None 5279790: status RUNNING/PENDING
2026-07-26 13:38:40 INFO None 5279791: status RUNNING/PENDING
2026-07-26 13:38:40 INFO Jobs still running: ['5279790', '5279791']. Waiting...
2026-07-26 13:38:56 INFO None 5279772: status FINISHED (not in squeue)
2026-07-26 13:38:56 INFO None 5279773: status FINISHED
2026-07-26 13:38:56 INFO None 5279775: status FINISHED (not in squeue)
2026-07-26 13:38:56 INFO None 5279777: status FINISHED (not in squeue)
2026-07-26 13:38:56 INFO None 5279778: status FINISHED (not in squeue)
2026-07-26 13:38:56 INFO None 5279779: status FINISHED (not in squeue)
2026-07-26 13:38:56 INFO None 5279780: status FINISHED
2026-07-26 13:38:56 INFO None 5279782: status FINISHED
2026-07-26 13:38:56 INFO None 5279783: status FINISHED
2026-07-26 13:38:56 INFO None 5279784: status FINISHED
2026-07-26 13:38:56 INFO None 5279785: status FINISHED
2026-07-26 13:38:56 INFO None 5279786: status FINISHED
2026-07-26 13:38:56 INFO None 5279789: status FINISHED
2026-07-26 13:38:56 INFO None 5279790: status FINISHED
2026-07-26 13:38:56 INFO None 5279791: status FINISHED
2026-07-26 13:38:56 INFO Jobs ['5279772', '5279773', '5279775', '5279777', '5279778', '5279779', '5279780', '5279782', '5279783', '5279784', '5279785', '5279786', '5279789', '5279790', '5279791'] have finished
2026-07-26 13:38:56 INFO Checking restart files were created ...
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-26 13:38:56 INFO  Run_model() completed successfully.
2026-07-26 13:38:56 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:38:56 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-26 13:38:56 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 13:38:56 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-26 13:38:56 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:38:56 INFO ---------->>> Running process_satellite_data()
2026-07-26 13:38:56 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-26 13:38:56 INFO ---------->>> Running run_obs_converter()
2026-07-26 13:38:56 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-26 13:38:56 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-26 13:38:56 INFO ---------->>> Running DART
2026-07-26 13:38:56 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 13:38:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-26 13:38:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-26 13:38:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-26 13:38:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-26 13:38:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-26 13:38:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-26 13:38:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-26 13:38:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-26 13:38:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-26 13:38:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-26 13:39:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-26 13:39:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-26 13:39:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-26 13:39:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-26 13:39:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-26 13:39:02 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 13:39:02 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 13:39:02 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 13:39:02 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 13:39:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 13:39:02 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 13:39:08 INFO Found: []
2026-07-26 13:39:08 INFO No job id returned by command ./run_filter.bsh
2026-07-26 13:39:08 INFO No monitoring will be performed
2026-07-26 13:39:08 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-26 13:39:08 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:08 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020709'
2026-07-26 13:39:09 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 13:39:09 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 13:39:09 INFO run_dart() is DONE.
2026-07-26 13:39:09 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 13:39:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-26 13:39:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-26 13:39:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:39:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-26 13:39:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-26 13:39:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:39:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-26 13:39:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-26 13:39:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:39:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-26 13:39:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-26 13:40:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:40:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-26 13:40:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-26 13:40:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:40:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-26 13:40:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-26 13:40:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:40:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-26 13:40:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-26 13:40:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:40:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-26 13:40:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-26 13:41:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:41:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-26 13:41:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-26 13:41:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:41:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-26 13:41:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-26 13:41:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:41:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-26 13:41:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-26 13:41:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:41:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-26 13:41:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-26 13:41:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:41:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-26 13:41:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-26 13:42:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:42:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-26 13:42:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-26 13:42:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:42:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-26 13:42:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-26 13:42:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:42:40 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 13:42:40 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:42:40 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:42:41 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-26 13:42:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:42 INFO Hourly dataset computed and listing created
2026-07-26 13:42:48 INFO Hourly dataset computed
2026-07-26 13:42:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:49 INFO Hourly dataset computed and listing created
2026-07-26 13:42:50 INFO Hourly dataset computed
2026-07-26 13:42:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:51 INFO Hourly dataset computed and listing created
2026-07-26 13:42:51 INFO Hourly dataset computed
2026-07-26 13:42:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:52 INFO Hourly dataset computed and listing created
2026-07-26 13:42:53 INFO Hourly dataset computed
2026-07-26 13:42:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:54 INFO Hourly dataset computed and listing created
2026-07-26 13:42:55 INFO Hourly dataset computed
2026-07-26 13:42:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:56 INFO Hourly dataset computed and listing created
2026-07-26 13:42:57 INFO Hourly dataset computed
2026-07-26 13:42:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:57 INFO Hourly dataset computed and listing created
2026-07-26 13:42:58 INFO Hourly dataset computed
2026-07-26 13:42:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:42:59 INFO Hourly dataset computed and listing created
2026-07-26 13:43:00 INFO Hourly dataset computed
2026-07-26 13:43:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:01 INFO Hourly dataset computed and listing created
2026-07-26 13:43:02 INFO Hourly dataset computed
2026-07-26 13:43:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:03 INFO Hourly dataset computed and listing created
2026-07-26 13:43:03 INFO Hourly dataset computed
2026-07-26 13:43:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:04 INFO Hourly dataset computed and listing created
2026-07-26 13:43:05 INFO Hourly dataset computed
2026-07-26 13:43:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:06 INFO Hourly dataset computed and listing created
2026-07-26 13:43:07 INFO Hourly dataset computed
2026-07-26 13:43:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:08 INFO Hourly dataset computed and listing created
2026-07-26 13:43:09 INFO Hourly dataset computed
2026-07-26 13:43:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:09 INFO Hourly dataset computed and listing created
2026-07-26 13:43:10 INFO Hourly dataset computed
2026-07-26 13:43:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:43:11 INFO Hourly dataset computed and listing created
2026-07-26 13:43:12 INFO Hourly dataset computed
2026-07-26 13:43:12 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-26 13:43:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 13:43:12 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-26 13:43:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 13:43:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:12 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 13:43:12 INFO Queuing job for member 1...
2026-07-26 13:43:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:12 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 13:43:14 INFO Found: ['5279844']
2026-07-26 13:43:19 INFO [TGCC-IRENE] Submitted job with ID:['5279844']
2026-07-26 13:43:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 13:43:19 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-26 13:43:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 13:43:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:19 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 13:43:19 INFO Queuing job for member 2...
2026-07-26 13:43:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:19 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 13:43:21 INFO Found: ['5279845']
2026-07-26 13:43:26 INFO [TGCC-IRENE] Submitted job with ID:['5279845']
2026-07-26 13:43:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 13:43:26 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-26 13:43:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 13:43:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:26 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 13:43:26 INFO Queuing job for member 3...
2026-07-26 13:43:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:26 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 13:43:29 INFO Found: ['5279846']
2026-07-26 13:43:34 INFO [TGCC-IRENE] Submitted job with ID:['5279846']
2026-07-26 13:43:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 13:43:34 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-26 13:43:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 13:43:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:34 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 13:43:34 INFO Queuing job for member 4...
2026-07-26 13:43:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:34 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 13:43:36 INFO Found: ['5279847']
2026-07-26 13:43:41 INFO [TGCC-IRENE] Submitted job with ID:['5279847']
2026-07-26 13:43:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 13:43:41 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-26 13:43:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 13:43:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:41 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 13:43:41 INFO Queuing job for member 5...
2026-07-26 13:43:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:41 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 13:43:44 INFO Found: ['5279849']
2026-07-26 13:43:49 INFO [TGCC-IRENE] Submitted job with ID:['5279849']
2026-07-26 13:43:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 13:43:49 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-26 13:43:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 13:43:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:49 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 13:43:49 INFO Queuing job for member 6...
2026-07-26 13:43:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:49 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 13:43:51 INFO Found: ['5279850']
2026-07-26 13:43:56 INFO [TGCC-IRENE] Submitted job with ID:['5279850']
2026-07-26 13:43:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:43:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 13:43:56 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-26 13:43:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 13:43:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:43:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 13:43:57 INFO Queuing job for member 7...
2026-07-26 13:43:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:43:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 13:43:59 INFO Found: ['5279851']
2026-07-26 13:44:04 INFO [TGCC-IRENE] Submitted job with ID:['5279851']
2026-07-26 13:44:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 13:44:04 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-26 13:44:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 13:44:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:04 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 13:44:04 INFO Queuing job for member 8...
2026-07-26 13:44:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:04 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 13:44:07 INFO Found: ['5279853']
2026-07-26 13:44:12 INFO [TGCC-IRENE] Submitted job with ID:['5279853']
2026-07-26 13:44:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 13:44:12 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-26 13:44:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 13:44:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:12 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 13:44:12 INFO Queuing job for member 9...
2026-07-26 13:44:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 13:44:14 INFO Found: ['5279854']
2026-07-26 13:44:19 INFO [TGCC-IRENE] Submitted job with ID:['5279854']
2026-07-26 13:44:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 13:44:19 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-26 13:44:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 13:44:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:19 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 13:44:19 INFO Queuing job for member 10...
2026-07-26 13:44:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:19 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 13:44:21 INFO Found: ['5279855']
2026-07-26 13:44:26 INFO [TGCC-IRENE] Submitted job with ID:['5279855']
2026-07-26 13:44:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 13:44:26 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-26 13:44:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 13:44:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:26 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 13:44:26 INFO Queuing job for member 11...
2026-07-26 13:44:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:26 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 13:44:27 INFO Found: ['5279857']
2026-07-26 13:44:32 INFO [TGCC-IRENE] Submitted job with ID:['5279857']
2026-07-26 13:44:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 13:44:32 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-26 13:44:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 13:44:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:32 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 13:44:32 INFO Queuing job for member 12...
2026-07-26 13:44:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:32 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 13:44:33 INFO Found: ['5279858']
2026-07-26 13:44:38 INFO [TGCC-IRENE] Submitted job with ID:['5279858']
2026-07-26 13:44:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 13:44:38 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-26 13:44:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 13:44:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:38 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 13:44:38 INFO Queuing job for member 13...
2026-07-26 13:44:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:38 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 13:44:39 INFO Found: ['5279859']
2026-07-26 13:44:44 INFO [TGCC-IRENE] Submitted job with ID:['5279859']
2026-07-26 13:44:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 13:44:44 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-26 13:44:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 13:44:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:44 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 13:44:44 INFO Queuing job for member 14...
2026-07-26 13:44:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:44 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 13:44:44 INFO Found: ['5279860']
2026-07-26 13:44:49 INFO [TGCC-IRENE] Submitted job with ID:['5279860']
2026-07-26 13:44:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:44:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 13:44:49 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-26 13:44:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 13:44:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:44:49 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 13:44:49 INFO Queuing job for member 15...
2026-07-26 13:44:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:44:49 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 13:44:50 INFO Found: ['5279861']
2026-07-26 13:44:55 INFO [TGCC-IRENE] Submitted job with ID:['5279861']
2026-07-26 13:44:55 INFO Checking job status ...
2026-07-26 13:44:55 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:44:55 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:44:55 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:45:11 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:45:11 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:45:11 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:45:26 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:45:26 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:45:26 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:45:41 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:45:41 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:45:43 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:45:44 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:45:44 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:45:59 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:45:59 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:46:01 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:46:01 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:46:16 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:46:16 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:46:16 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:46:31 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:46:31 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:46:32 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:46:32 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:46:47 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:46:47 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:46:47 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:47:02 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:47:02 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:47:02 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:47:17 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:47:17 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:47:18 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:47:18 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:47:18 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:47:18 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:47:18 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:47:33 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:47:33 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:47:35 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:47:35 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:47:35 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:47:35 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:47:35 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:47:35 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:47:50 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:47:50 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:47:50 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:48:05 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:48:05 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:48:07 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:48:07 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:48:07 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:48:07 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:48:07 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:48:23 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279853: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279855: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:48:23 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:48:23 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:48:38 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279851: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279853: status FINISHED
2026-07-26 13:48:38 INFO None 5279854: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279855: status FINISHED
2026-07-26 13:48:38 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:48:38 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:48:38 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279854', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:48:53 INFO None 5279844: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279845: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279846: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279847: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279849: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279850: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279851: status FINISHED
2026-07-26 13:48:53 INFO None 5279853: status FINISHED
2026-07-26 13:48:53 INFO None 5279854: status FINISHED
2026-07-26 13:48:53 INFO None 5279855: status FINISHED
2026-07-26 13:48:53 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:48:53 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:48:53 INFO Jobs still running: ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:49:08 INFO None 5279844: status FINISHED
2026-07-26 13:49:08 INFO None 5279845: status FINISHED
2026-07-26 13:49:08 INFO None 5279846: status FINISHED
2026-07-26 13:49:08 INFO None 5279847: status FINISHED
2026-07-26 13:49:08 INFO None 5279849: status FINISHED
2026-07-26 13:49:08 INFO None 5279850: status FINISHED
2026-07-26 13:49:08 INFO None 5279851: status FINISHED
2026-07-26 13:49:09 INFO None 5279853: status FINISHED
2026-07-26 13:49:09 INFO None 5279854: status FINISHED
2026-07-26 13:49:09 INFO None 5279855: status FINISHED
2026-07-26 13:49:09 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:49:09 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:49:09 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:49:09 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:49:09 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:49:09 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:49:25 INFO None 5279844: status FINISHED
2026-07-26 13:49:25 INFO None 5279845: status FINISHED
2026-07-26 13:49:25 INFO None 5279846: status FINISHED
2026-07-26 13:49:25 INFO None 5279847: status FINISHED
2026-07-26 13:49:25 INFO None 5279849: status FINISHED
2026-07-26 13:49:25 INFO None 5279850: status FINISHED
2026-07-26 13:49:25 INFO None 5279851: status FINISHED
2026-07-26 13:49:25 INFO None 5279853: status FINISHED
2026-07-26 13:49:25 INFO None 5279854: status FINISHED
2026-07-26 13:49:25 INFO None 5279855: status FINISHED
2026-07-26 13:49:25 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:49:25 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:49:25 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:49:25 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:49:25 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:49:25 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:49:40 INFO None 5279844: status FINISHED
2026-07-26 13:49:40 INFO None 5279845: status FINISHED
2026-07-26 13:49:40 INFO None 5279846: status FINISHED
2026-07-26 13:49:41 INFO None 5279847: status FINISHED
2026-07-26 13:49:41 INFO None 5279849: status FINISHED
2026-07-26 13:49:41 INFO None 5279850: status FINISHED
2026-07-26 13:49:41 INFO None 5279851: status FINISHED
2026-07-26 13:49:41 INFO None 5279853: status FINISHED
2026-07-26 13:49:41 INFO None 5279854: status FINISHED
2026-07-26 13:49:41 INFO None 5279855: status FINISHED
2026-07-26 13:49:41 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:49:41 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:49:41 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:49:41 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:49:41 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:49:41 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:49:56 INFO None 5279844: status FINISHED
2026-07-26 13:49:56 INFO None 5279845: status FINISHED
2026-07-26 13:49:56 INFO None 5279846: status FINISHED
2026-07-26 13:49:58 INFO None 5279847: status FINISHED
2026-07-26 13:49:58 INFO None 5279849: status FINISHED
2026-07-26 13:49:58 INFO None 5279850: status FINISHED
2026-07-26 13:49:58 INFO None 5279851: status FINISHED
2026-07-26 13:49:58 INFO None 5279853: status FINISHED
2026-07-26 13:49:58 INFO None 5279854: status FINISHED
2026-07-26 13:49:58 INFO None 5279855: status FINISHED
2026-07-26 13:49:58 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:49:58 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:49:58 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:49:58 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:49:58 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:49:58 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:50:13 INFO None 5279844: status FINISHED
2026-07-26 13:50:13 INFO None 5279845: status FINISHED
2026-07-26 13:50:13 INFO None 5279846: status FINISHED
2026-07-26 13:50:13 INFO None 5279847: status FINISHED
2026-07-26 13:50:13 INFO None 5279849: status FINISHED
2026-07-26 13:50:13 INFO None 5279850: status FINISHED
2026-07-26 13:50:13 INFO None 5279851: status FINISHED
2026-07-26 13:50:13 INFO None 5279853: status FINISHED
2026-07-26 13:50:13 INFO None 5279854: status FINISHED
2026-07-26 13:50:13 INFO None 5279855: status FINISHED
2026-07-26 13:50:13 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:50:13 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:50:13 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:50:13 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:50:13 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:50:13 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:50:28 INFO None 5279844: status FINISHED
2026-07-26 13:50:28 INFO None 5279845: status FINISHED
2026-07-26 13:50:28 INFO None 5279846: status FINISHED
2026-07-26 13:50:28 INFO None 5279847: status FINISHED
2026-07-26 13:50:28 INFO None 5279849: status FINISHED
2026-07-26 13:50:28 INFO None 5279850: status FINISHED
2026-07-26 13:50:30 INFO None 5279851: status FINISHED
2026-07-26 13:50:30 INFO None 5279853: status FINISHED
2026-07-26 13:50:31 INFO None 5279854: status FINISHED
2026-07-26 13:50:31 INFO None 5279855: status FINISHED
2026-07-26 13:50:31 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:50:31 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:50:31 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:50:31 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:50:31 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:50:31 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:50:46 INFO None 5279844: status FINISHED
2026-07-26 13:50:46 INFO None 5279845: status FINISHED
2026-07-26 13:50:46 INFO None 5279846: status FINISHED
2026-07-26 13:50:46 INFO None 5279847: status FINISHED
2026-07-26 13:50:46 INFO None 5279849: status FINISHED
2026-07-26 13:50:46 INFO None 5279850: status FINISHED
2026-07-26 13:50:46 INFO None 5279851: status FINISHED
2026-07-26 13:50:46 INFO None 5279853: status FINISHED
2026-07-26 13:50:46 INFO None 5279854: status FINISHED
2026-07-26 13:50:46 INFO None 5279855: status FINISHED
2026-07-26 13:50:46 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:50:46 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:50:46 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:50:46 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:50:46 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:50:46 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:51:01 INFO None 5279844: status FINISHED
2026-07-26 13:51:01 INFO None 5279845: status FINISHED
2026-07-26 13:51:01 INFO None 5279846: status FINISHED
2026-07-26 13:51:01 INFO None 5279847: status FINISHED
2026-07-26 13:51:01 INFO None 5279849: status FINISHED
2026-07-26 13:51:01 INFO None 5279850: status FINISHED
2026-07-26 13:51:01 INFO None 5279851: status FINISHED
2026-07-26 13:51:01 INFO None 5279853: status FINISHED
2026-07-26 13:51:01 INFO None 5279854: status FINISHED
2026-07-26 13:51:01 INFO None 5279855: status FINISHED
2026-07-26 13:51:01 INFO None 5279857: status RUNNING/PENDING
2026-07-26 13:51:01 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:51:01 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:51:01 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:51:01 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:51:01 INFO Jobs still running: ['5279857', '5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:51:16 INFO None 5279844: status FINISHED
2026-07-26 13:51:16 INFO None 5279845: status FINISHED
2026-07-26 13:51:16 INFO None 5279846: status FINISHED
2026-07-26 13:51:16 INFO None 5279847: status FINISHED
2026-07-26 13:51:16 INFO None 5279849: status FINISHED
2026-07-26 13:51:16 INFO None 5279850: status FINISHED
2026-07-26 13:51:16 INFO None 5279851: status FINISHED
2026-07-26 13:51:16 INFO None 5279853: status FINISHED
2026-07-26 13:51:16 INFO None 5279854: status FINISHED
2026-07-26 13:51:17 INFO None 5279855: status FINISHED
2026-07-26 13:51:17 INFO None 5279857: status FINISHED
2026-07-26 13:51:17 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:51:17 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:51:17 INFO None 5279860: status RUNNING/PENDING
2026-07-26 13:51:17 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:51:17 INFO Jobs still running: ['5279858', '5279859', '5279860', '5279861']. Waiting...
2026-07-26 13:51:32 INFO None 5279844: status FINISHED
2026-07-26 13:51:32 INFO None 5279845: status FINISHED
2026-07-26 13:51:32 INFO None 5279846: status FINISHED
2026-07-26 13:51:32 INFO None 5279847: status FINISHED
2026-07-26 13:51:32 INFO None 5279849: status FINISHED
2026-07-26 13:51:32 INFO None 5279850: status FINISHED
2026-07-26 13:51:32 INFO None 5279851: status FINISHED
2026-07-26 13:51:32 INFO None 5279853: status FINISHED
2026-07-26 13:51:32 INFO None 5279854: status FINISHED
2026-07-26 13:51:32 INFO None 5279855: status FINISHED
2026-07-26 13:51:32 INFO None 5279857: status FINISHED
2026-07-26 13:51:32 INFO None 5279858: status RUNNING/PENDING
2026-07-26 13:51:32 INFO None 5279859: status RUNNING/PENDING
2026-07-26 13:51:32 INFO None 5279860: status FINISHED
2026-07-26 13:51:32 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:51:32 INFO Jobs still running: ['5279858', '5279859', '5279861']. Waiting...
2026-07-26 13:51:47 INFO None 5279844: status FINISHED
2026-07-26 13:51:47 INFO None 5279845: status FINISHED
2026-07-26 13:51:47 INFO None 5279846: status FINISHED
2026-07-26 13:51:47 INFO None 5279847: status FINISHED
2026-07-26 13:51:47 INFO None 5279849: status FINISHED
2026-07-26 13:51:49 INFO None 5279850: status FINISHED
2026-07-26 13:51:49 INFO None 5279851: status FINISHED
2026-07-26 13:51:49 INFO None 5279853: status FINISHED
2026-07-26 13:51:49 INFO None 5279854: status FINISHED
2026-07-26 13:51:49 INFO None 5279855: status FINISHED
2026-07-26 13:51:49 INFO None 5279857: status FINISHED
2026-07-26 13:51:49 INFO None 5279858: status FINISHED
2026-07-26 13:51:49 INFO None 5279859: status FINISHED
2026-07-26 13:51:49 INFO None 5279860: status FINISHED
2026-07-26 13:51:49 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:51:49 INFO Jobs still running: ['5279861']. Waiting...
2026-07-26 13:52:04 INFO None 5279844: status FINISHED
2026-07-26 13:52:04 INFO None 5279845: status FINISHED
2026-07-26 13:52:05 INFO None 5279846: status FINISHED
2026-07-26 13:52:05 INFO None 5279847: status FINISHED
2026-07-26 13:52:05 INFO None 5279849: status FINISHED
2026-07-26 13:52:05 INFO None 5279850: status FINISHED
2026-07-26 13:52:05 INFO None 5279851: status FINISHED
2026-07-26 13:52:05 INFO None 5279853: status FINISHED
2026-07-26 13:52:05 INFO None 5279854: status FINISHED
2026-07-26 13:52:05 INFO None 5279855: status FINISHED
2026-07-26 13:52:05 INFO None 5279857: status FINISHED
2026-07-26 13:52:05 INFO None 5279858: status FINISHED
2026-07-26 13:52:05 INFO None 5279859: status FINISHED
2026-07-26 13:52:05 INFO None 5279860: status FINISHED
2026-07-26 13:52:05 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:52:05 INFO Jobs still running: ['5279861']. Waiting...
2026-07-26 13:52:22 INFO None 5279844: status FINISHED
2026-07-26 13:52:22 INFO None 5279845: status FINISHED
2026-07-26 13:52:22 INFO None 5279846: status FINISHED
2026-07-26 13:52:22 INFO None 5279847: status FINISHED
2026-07-26 13:52:22 INFO None 5279849: status FINISHED
2026-07-26 13:52:22 INFO None 5279850: status FINISHED
2026-07-26 13:52:22 INFO None 5279851: status FINISHED
2026-07-26 13:52:22 INFO None 5279853: status FINISHED
2026-07-26 13:52:22 INFO None 5279854: status FINISHED
2026-07-26 13:52:22 INFO None 5279855: status FINISHED
2026-07-26 13:52:22 INFO None 5279857: status FINISHED
2026-07-26 13:52:22 INFO None 5279858: status FINISHED
2026-07-26 13:52:22 INFO None 5279859: status FINISHED
2026-07-26 13:52:22 INFO None 5279860: status FINISHED
2026-07-26 13:52:22 INFO None 5279861: status RUNNING/PENDING
2026-07-26 13:52:22 INFO Jobs still running: ['5279861']. Waiting...
2026-07-26 13:52:37 INFO None 5279844: status FINISHED
2026-07-26 13:52:37 INFO None 5279845: status FINISHED
2026-07-26 13:52:37 INFO None 5279846: status FINISHED
2026-07-26 13:52:37 INFO None 5279847: status FINISHED
2026-07-26 13:52:37 INFO None 5279849: status FINISHED
2026-07-26 13:52:37 INFO None 5279850: status FINISHED
2026-07-26 13:52:37 INFO None 5279851: status FINISHED
2026-07-26 13:52:37 INFO None 5279853: status FINISHED
2026-07-26 13:52:37 INFO None 5279854: status FINISHED
2026-07-26 13:52:37 INFO None 5279855: status FINISHED
2026-07-26 13:52:37 INFO None 5279857: status FINISHED
2026-07-26 13:52:37 INFO None 5279858: status FINISHED
2026-07-26 13:52:37 INFO None 5279859: status FINISHED
2026-07-26 13:52:37 INFO None 5279860: status FINISHED
2026-07-26 13:52:37 INFO None 5279861: status FINISHED
2026-07-26 13:52:37 INFO Jobs ['5279844', '5279845', '5279846', '5279847', '5279849', '5279850', '5279851', '5279853', '5279854', '5279855', '5279857', '5279858', '5279859', '5279860', '5279861'] have finished
2026-07-26 13:52:37 INFO Checking restart files were created ...
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-26 13:52:37 INFO  Run_model() completed successfully.
2026-07-26 13:52:37 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:52:37 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-26 13:52:37 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 13:52:37 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-26 13:52:37 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:52:37 INFO ---------->>> Running process_satellite_data()
2026-07-26 13:52:38 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-26 13:52:38 INFO ---------->>> Running run_obs_converter()
2026-07-26 13:52:38 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-26 13:52:38 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-26 13:52:38 INFO ---------->>> Running DART
2026-07-26 13:52:38 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 13:52:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-26 13:52:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-26 13:52:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-26 13:52:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-26 13:52:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-26 13:52:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-26 13:52:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-26 13:52:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-26 13:52:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-26 13:52:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-26 13:52:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-26 13:52:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-26 13:52:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-26 13:52:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-26 13:52:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-26 13:52:42 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 13:52:42 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 13:52:42 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 13:52:42 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 13:52:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 13:52:42 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 13:52:55 INFO Found: []
2026-07-26 13:52:55 INFO No job id returned by command ./run_filter.bsh
2026-07-26 13:52:55 INFO No monitoring will be performed
2026-07-26 13:52:55 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-26 13:52:55 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020711'
2026-07-26 13:52:55 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 13:52:55 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 13:52:55 INFO run_dart() is DONE.
2026-07-26 13:52:55 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 13:52:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-26 13:53:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-26 13:53:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-26 13:53:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-26 13:53:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-26 13:53:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-26 13:53:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-26 13:53:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-26 13:53:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-26 13:53:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-26 13:53:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-26 13:53:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:53:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-26 13:54:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:54:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-26 13:54:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:54:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-26 13:54:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:54:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-26 13:54:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 13:54:17 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 13:54:17 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:54:17 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 13:54:17 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-26 13:54:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:18 INFO Hourly dataset computed and listing created
2026-07-26 13:54:22 INFO Hourly dataset computed
2026-07-26 13:54:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:23 INFO Hourly dataset computed and listing created
2026-07-26 13:54:24 INFO Hourly dataset computed
2026-07-26 13:54:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:25 INFO Hourly dataset computed and listing created
2026-07-26 13:54:26 INFO Hourly dataset computed
2026-07-26 13:54:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:27 INFO Hourly dataset computed and listing created
2026-07-26 13:54:27 INFO Hourly dataset computed
2026-07-26 13:54:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:28 INFO Hourly dataset computed and listing created
2026-07-26 13:54:29 INFO Hourly dataset computed
2026-07-26 13:54:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:30 INFO Hourly dataset computed and listing created
2026-07-26 13:54:30 INFO Hourly dataset computed
2026-07-26 13:54:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:31 INFO Hourly dataset computed and listing created
2026-07-26 13:54:32 INFO Hourly dataset computed
2026-07-26 13:54:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:32 INFO Hourly dataset computed and listing created
2026-07-26 13:54:33 INFO Hourly dataset computed
2026-07-26 13:54:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:34 INFO Hourly dataset computed and listing created
2026-07-26 13:54:35 INFO Hourly dataset computed
2026-07-26 13:54:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:36 INFO Hourly dataset computed and listing created
2026-07-26 13:54:36 INFO Hourly dataset computed
2026-07-26 13:54:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:37 INFO Hourly dataset computed and listing created
2026-07-26 13:54:38 INFO Hourly dataset computed
2026-07-26 13:54:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:38 INFO Hourly dataset computed and listing created
2026-07-26 13:54:39 INFO Hourly dataset computed
2026-07-26 13:54:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:40 INFO Hourly dataset computed and listing created
2026-07-26 13:54:41 INFO Hourly dataset computed
2026-07-26 13:54:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:42 INFO Hourly dataset computed and listing created
2026-07-26 13:54:42 INFO Hourly dataset computed
2026-07-26 13:54:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 13:54:43 INFO Hourly dataset computed and listing created
2026-07-26 13:54:44 INFO Hourly dataset computed
2026-07-26 13:54:44 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-26 13:54:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:54:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 13:54:44 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-26 13:54:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 13:54:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:54:44 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 13:54:44 INFO Queuing job for member 1...
2026-07-26 13:54:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:54:44 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 13:54:45 INFO Found: ['5279898']
2026-07-26 13:54:50 INFO [TGCC-IRENE] Submitted job with ID:['5279898']
2026-07-26 13:54:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:54:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 13:54:50 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-26 13:54:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 13:54:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:54:50 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 13:54:50 INFO Queuing job for member 2...
2026-07-26 13:54:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:54:50 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 13:54:50 INFO Found: ['5279899']
2026-07-26 13:54:55 INFO [TGCC-IRENE] Submitted job with ID:['5279899']
2026-07-26 13:54:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:54:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 13:54:55 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-26 13:54:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 13:54:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:54:55 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 13:54:55 INFO Queuing job for member 3...
2026-07-26 13:54:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:54:55 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 13:54:56 INFO Found: ['5279900']
2026-07-26 13:55:01 INFO [TGCC-IRENE] Submitted job with ID:['5279900']
2026-07-26 13:55:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 13:55:01 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-26 13:55:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 13:55:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:01 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 13:55:01 INFO Queuing job for member 4...
2026-07-26 13:55:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:01 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 13:55:02 INFO Found: ['5279902']
2026-07-26 13:55:07 INFO [TGCC-IRENE] Submitted job with ID:['5279902']
2026-07-26 13:55:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 13:55:07 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-26 13:55:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 13:55:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:07 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 13:55:07 INFO Queuing job for member 5...
2026-07-26 13:55:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:07 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 13:55:07 INFO Found: ['5279903']
2026-07-26 13:55:12 INFO [TGCC-IRENE] Submitted job with ID:['5279903']
2026-07-26 13:55:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 13:55:12 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-26 13:55:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 13:55:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:12 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 13:55:12 INFO Queuing job for member 6...
2026-07-26 13:55:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:12 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 13:55:13 INFO Found: ['5279904']
2026-07-26 13:55:18 INFO [TGCC-IRENE] Submitted job with ID:['5279904']
2026-07-26 13:55:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 13:55:18 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-26 13:55:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 13:55:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:18 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 13:55:18 INFO Queuing job for member 7...
2026-07-26 13:55:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:18 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 13:55:21 INFO Found: ['5279905']
2026-07-26 13:55:26 INFO [TGCC-IRENE] Submitted job with ID:['5279905']
2026-07-26 13:55:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 13:55:26 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-26 13:55:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 13:55:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:26 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 13:55:26 INFO Queuing job for member 8...
2026-07-26 13:55:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:26 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 13:55:28 INFO Found: ['5279906']
2026-07-26 13:55:33 INFO [TGCC-IRENE] Submitted job with ID:['5279906']
2026-07-26 13:55:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 13:55:33 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-26 13:55:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 13:55:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:33 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 13:55:33 INFO Queuing job for member 9...
2026-07-26 13:55:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:33 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 13:55:36 INFO Found: ['5279907']
2026-07-26 13:55:41 INFO [TGCC-IRENE] Submitted job with ID:['5279907']
2026-07-26 13:55:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 13:55:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-26 13:55:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 13:55:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 13:55:41 INFO Queuing job for member 10...
2026-07-26 13:55:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 13:55:43 INFO Found: ['5279908']
2026-07-26 13:55:48 INFO [TGCC-IRENE] Submitted job with ID:['5279908']
2026-07-26 13:55:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 13:55:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-26 13:55:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 13:55:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 13:55:48 INFO Queuing job for member 11...
2026-07-26 13:55:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 13:55:51 INFO Found: ['5279909']
2026-07-26 13:55:56 INFO [TGCC-IRENE] Submitted job with ID:['5279909']
2026-07-26 13:55:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:55:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 13:55:56 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-26 13:55:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 13:55:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:55:56 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 13:55:56 INFO Queuing job for member 12...
2026-07-26 13:55:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:55:56 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 13:55:58 INFO Found: ['5279910']
2026-07-26 13:56:03 INFO [TGCC-IRENE] Submitted job with ID:['5279910']
2026-07-26 13:56:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:56:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 13:56:03 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-26 13:56:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 13:56:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:56:03 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 13:56:03 INFO Queuing job for member 13...
2026-07-26 13:56:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:56:03 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 13:56:06 INFO Found: ['5279912']
2026-07-26 13:56:11 INFO [TGCC-IRENE] Submitted job with ID:['5279912']
2026-07-26 13:56:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:56:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 13:56:11 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-26 13:56:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 13:56:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:56:11 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 13:56:11 INFO Queuing job for member 14...
2026-07-26 13:56:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:56:11 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 13:56:13 INFO Found: ['5279913']
2026-07-26 13:56:18 INFO [TGCC-IRENE] Submitted job with ID:['5279913']
2026-07-26 13:56:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 13:56:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 13:56:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-26 13:56:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 13:56:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 13:56:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 13:56:18 INFO Queuing job for member 15...
2026-07-26 13:56:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 13:56:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 13:56:21 INFO Found: ['5279914']
2026-07-26 13:56:26 INFO [TGCC-IRENE] Submitted job with ID:['5279914']
2026-07-26 13:56:26 INFO Checking job status ...
2026-07-26 13:56:26 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:56:26 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:56:26 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:56:41 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:56:41 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:56:41 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:56:56 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:56:56 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:56:57 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:56:57 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:56:57 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:57:12 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:57:12 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:57:12 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:57:27 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:57:27 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:57:27 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:57:42 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:57:42 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:57:45 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:57:45 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:58:00 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279900: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:58:00 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:58:00 INFO Jobs still running: ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:58:15 INFO None 5279898: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279899: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279900: status FINISHED
2026-07-26 13:58:15 INFO None 5279902: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279903: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:58:15 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:58:17 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:58:17 INFO Jobs still running: ['5279898', '5279899', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:58:32 INFO None 5279898: status FINISHED
2026-07-26 13:58:32 INFO None 5279899: status FINISHED
2026-07-26 13:58:32 INFO None 5279900: status FINISHED
2026-07-26 13:58:32 INFO None 5279902: status FINISHED
2026-07-26 13:58:32 INFO None 5279903: status FINISHED
2026-07-26 13:58:32 INFO None 5279904: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:58:32 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:58:32 INFO Jobs still running: ['5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:58:48 INFO None 5279898: status FINISHED
2026-07-26 13:58:48 INFO None 5279899: status FINISHED
2026-07-26 13:58:48 INFO None 5279900: status FINISHED
2026-07-26 13:58:48 INFO None 5279902: status FINISHED
2026-07-26 13:58:48 INFO None 5279903: status FINISHED
2026-07-26 13:58:48 INFO None 5279904: status FINISHED
2026-07-26 13:58:48 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:58:48 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:58:48 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:59:03 INFO None 5279898: status FINISHED
2026-07-26 13:59:03 INFO None 5279899: status FINISHED
2026-07-26 13:59:03 INFO None 5279900: status FINISHED
2026-07-26 13:59:03 INFO None 5279902: status FINISHED
2026-07-26 13:59:03 INFO None 5279903: status FINISHED
2026-07-26 13:59:03 INFO None 5279904: status FINISHED
2026-07-26 13:59:03 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:59:03 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:59:03 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:59:18 INFO None 5279898: status FINISHED
2026-07-26 13:59:18 INFO None 5279899: status FINISHED
2026-07-26 13:59:18 INFO None 5279900: status FINISHED
2026-07-26 13:59:18 INFO None 5279902: status FINISHED
2026-07-26 13:59:18 INFO None 5279903: status FINISHED
2026-07-26 13:59:18 INFO None 5279904: status FINISHED
2026-07-26 13:59:18 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:59:18 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:59:18 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:59:34 INFO None 5279898: status FINISHED
2026-07-26 13:59:34 INFO None 5279899: status FINISHED
2026-07-26 13:59:34 INFO None 5279900: status FINISHED
2026-07-26 13:59:34 INFO None 5279902: status FINISHED
2026-07-26 13:59:34 INFO None 5279903: status FINISHED
2026-07-26 13:59:34 INFO None 5279904: status FINISHED
2026-07-26 13:59:34 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:59:34 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:59:34 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 13:59:49 INFO None 5279898: status FINISHED
2026-07-26 13:59:49 INFO None 5279899: status FINISHED
2026-07-26 13:59:49 INFO None 5279900: status FINISHED
2026-07-26 13:59:49 INFO None 5279902: status FINISHED
2026-07-26 13:59:51 INFO None 5279903: status FINISHED
2026-07-26 13:59:51 INFO None 5279904: status FINISHED
2026-07-26 13:59:51 INFO None 5279905: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279906: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279907: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279908: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279909: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279910: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279912: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279913: status RUNNING/PENDING
2026-07-26 13:59:51 INFO None 5279914: status RUNNING/PENDING
2026-07-26 13:59:51 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:00:06 INFO None 5279898: status FINISHED
2026-07-26 14:00:06 INFO None 5279899: status FINISHED
2026-07-26 14:00:06 INFO None 5279900: status FINISHED
2026-07-26 14:00:06 INFO None 5279902: status FINISHED
2026-07-26 14:00:06 INFO None 5279903: status FINISHED
2026-07-26 14:00:06 INFO None 5279904: status FINISHED
2026-07-26 14:00:06 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:00:06 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:00:06 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:00:06 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:00:06 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:00:09 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:00:09 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:00:09 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:00:09 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:00:09 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:00:24 INFO None 5279898: status FINISHED
2026-07-26 14:00:24 INFO None 5279899: status FINISHED
2026-07-26 14:00:24 INFO None 5279900: status FINISHED
2026-07-26 14:00:24 INFO None 5279902: status FINISHED
2026-07-26 14:00:24 INFO None 5279903: status FINISHED
2026-07-26 14:00:24 INFO None 5279904: status FINISHED
2026-07-26 14:00:24 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:00:24 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:00:24 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:00:39 INFO None 5279898: status FINISHED
2026-07-26 14:00:39 INFO None 5279899: status FINISHED
2026-07-26 14:00:39 INFO None 5279900: status FINISHED
2026-07-26 14:00:39 INFO None 5279902: status FINISHED
2026-07-26 14:00:39 INFO None 5279903: status FINISHED
2026-07-26 14:00:39 INFO None 5279904: status FINISHED
2026-07-26 14:00:39 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:00:39 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:00:39 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:00:39 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:00:39 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:00:39 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:00:41 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:00:41 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:00:41 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:00:41 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:00:56 INFO None 5279898: status FINISHED
2026-07-26 14:00:56 INFO None 5279899: status FINISHED
2026-07-26 14:00:56 INFO None 5279900: status FINISHED
2026-07-26 14:00:56 INFO None 5279902: status FINISHED
2026-07-26 14:00:56 INFO None 5279903: status FINISHED
2026-07-26 14:00:56 INFO None 5279904: status FINISHED
2026-07-26 14:00:56 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:00:56 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:00:56 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:00:57 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:00:57 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:01:12 INFO None 5279898: status FINISHED
2026-07-26 14:01:12 INFO None 5279899: status FINISHED
2026-07-26 14:01:12 INFO None 5279900: status FINISHED
2026-07-26 14:01:12 INFO None 5279902: status FINISHED
2026-07-26 14:01:12 INFO None 5279903: status FINISHED
2026-07-26 14:01:12 INFO None 5279904: status FINISHED
2026-07-26 14:01:12 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:01:12 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:01:12 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:01:27 INFO None 5279898: status FINISHED
2026-07-26 14:01:27 INFO None 5279899: status FINISHED
2026-07-26 14:01:27 INFO None 5279900: status FINISHED
2026-07-26 14:01:27 INFO None 5279902: status FINISHED
2026-07-26 14:01:27 INFO None 5279903: status FINISHED
2026-07-26 14:01:27 INFO None 5279904: status FINISHED
2026-07-26 14:01:27 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:01:27 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:01:27 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:01:43 INFO None 5279898: status FINISHED
2026-07-26 14:01:43 INFO None 5279899: status FINISHED
2026-07-26 14:01:43 INFO None 5279900: status FINISHED
2026-07-26 14:01:43 INFO None 5279902: status FINISHED
2026-07-26 14:01:43 INFO None 5279903: status FINISHED
2026-07-26 14:01:43 INFO None 5279904: status FINISHED
2026-07-26 14:01:43 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:01:43 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:01:43 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:01:58 INFO None 5279898: status FINISHED
2026-07-26 14:01:58 INFO None 5279899: status FINISHED
2026-07-26 14:01:58 INFO None 5279900: status FINISHED
2026-07-26 14:01:58 INFO None 5279902: status FINISHED
2026-07-26 14:01:58 INFO None 5279903: status FINISHED
2026-07-26 14:01:58 INFO None 5279904: status FINISHED
2026-07-26 14:01:58 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:01:58 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:02:00 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:02:00 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:02:15 INFO None 5279898: status FINISHED
2026-07-26 14:02:15 INFO None 5279899: status FINISHED
2026-07-26 14:02:15 INFO None 5279900: status FINISHED
2026-07-26 14:02:15 INFO None 5279902: status FINISHED
2026-07-26 14:02:16 INFO None 5279903: status FINISHED
2026-07-26 14:02:16 INFO None 5279904: status FINISHED
2026-07-26 14:02:16 INFO None 5279905: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279906: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279907: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279908: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279909: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:02:16 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:02:16 INFO Jobs still running: ['5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:02:31 INFO None 5279898: status FINISHED
2026-07-26 14:02:31 INFO None 5279899: status FINISHED
2026-07-26 14:02:31 INFO None 5279900: status FINISHED
2026-07-26 14:02:31 INFO None 5279902: status FINISHED
2026-07-26 14:02:31 INFO None 5279903: status FINISHED
2026-07-26 14:02:31 INFO None 5279904: status FINISHED
2026-07-26 14:02:31 INFO None 5279905: status FINISHED
2026-07-26 14:02:31 INFO None 5279906: status FINISHED
2026-07-26 14:02:31 INFO None 5279907: status FINISHED
2026-07-26 14:02:31 INFO None 5279908: status FINISHED
2026-07-26 14:02:33 INFO None 5279909: status FINISHED
2026-07-26 14:02:33 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:02:33 INFO None 5279912: status RUNNING/PENDING
2026-07-26 14:02:33 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:02:33 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:02:33 INFO Jobs still running: ['5279910', '5279912', '5279913', '5279914']. Waiting...
2026-07-26 14:02:48 INFO None 5279898: status FINISHED
2026-07-26 14:02:48 INFO None 5279899: status FINISHED
2026-07-26 14:02:48 INFO None 5279900: status FINISHED
2026-07-26 14:02:48 INFO None 5279902: status FINISHED
2026-07-26 14:02:48 INFO None 5279903: status FINISHED
2026-07-26 14:02:48 INFO None 5279904: status FINISHED
2026-07-26 14:02:48 INFO None 5279905: status FINISHED
2026-07-26 14:02:48 INFO None 5279906: status FINISHED
2026-07-26 14:02:48 INFO None 5279907: status FINISHED
2026-07-26 14:02:48 INFO None 5279908: status FINISHED
2026-07-26 14:02:48 INFO None 5279909: status FINISHED
2026-07-26 14:02:48 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:02:48 INFO None 5279912: status FINISHED
2026-07-26 14:02:48 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:02:48 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:02:48 INFO Jobs still running: ['5279910', '5279913', '5279914']. Waiting...
2026-07-26 14:03:03 INFO None 5279898: status FINISHED
2026-07-26 14:03:03 INFO None 5279899: status FINISHED
2026-07-26 14:03:03 INFO None 5279900: status FINISHED
2026-07-26 14:03:03 INFO None 5279902: status FINISHED
2026-07-26 14:03:03 INFO None 5279903: status FINISHED
2026-07-26 14:03:03 INFO None 5279904: status FINISHED
2026-07-26 14:03:03 INFO None 5279905: status FINISHED
2026-07-26 14:03:03 INFO None 5279906: status FINISHED
2026-07-26 14:03:03 INFO None 5279907: status FINISHED
2026-07-26 14:03:03 INFO None 5279908: status FINISHED
2026-07-26 14:03:03 INFO None 5279909: status FINISHED
2026-07-26 14:03:04 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:03:04 INFO None 5279912: status FINISHED
2026-07-26 14:03:04 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:03:04 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:03:04 INFO Jobs still running: ['5279910', '5279913', '5279914']. Waiting...
2026-07-26 14:03:19 INFO None 5279898: status FINISHED
2026-07-26 14:03:19 INFO None 5279899: status FINISHED
2026-07-26 14:03:19 INFO None 5279900: status FINISHED
2026-07-26 14:03:19 INFO None 5279902: status FINISHED
2026-07-26 14:03:19 INFO None 5279903: status FINISHED
2026-07-26 14:03:19 INFO None 5279904: status FINISHED
2026-07-26 14:03:19 INFO None 5279905: status FINISHED
2026-07-26 14:03:19 INFO None 5279906: status FINISHED
2026-07-26 14:03:19 INFO None 5279907: status FINISHED
2026-07-26 14:03:19 INFO None 5279908: status FINISHED
2026-07-26 14:03:19 INFO None 5279909: status FINISHED
2026-07-26 14:03:19 INFO None 5279910: status RUNNING/PENDING
2026-07-26 14:03:19 INFO None 5279912: status FINISHED
2026-07-26 14:03:19 INFO None 5279913: status RUNNING/PENDING
2026-07-26 14:03:19 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:03:19 INFO Jobs still running: ['5279910', '5279913', '5279914']. Waiting...
2026-07-26 14:03:34 INFO None 5279898: status FINISHED
2026-07-26 14:03:34 INFO None 5279899: status FINISHED
2026-07-26 14:03:34 INFO None 5279900: status FINISHED
2026-07-26 14:03:34 INFO None 5279902: status FINISHED
2026-07-26 14:03:34 INFO None 5279903: status FINISHED
2026-07-26 14:03:34 INFO None 5279904: status FINISHED
2026-07-26 14:03:34 INFO None 5279905: status FINISHED
2026-07-26 14:03:34 INFO None 5279906: status FINISHED
2026-07-26 14:03:34 INFO None 5279907: status FINISHED
2026-07-26 14:03:34 INFO None 5279908: status FINISHED
2026-07-26 14:03:34 INFO None 5279909: status FINISHED
2026-07-26 14:03:34 INFO None 5279910: status FINISHED
2026-07-26 14:03:34 INFO None 5279912: status FINISHED
2026-07-26 14:03:34 INFO None 5279913: status FINISHED
2026-07-26 14:03:34 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:03:34 INFO Jobs still running: ['5279914']. Waiting...
2026-07-26 14:03:49 INFO None 5279898: status FINISHED
2026-07-26 14:03:49 INFO None 5279899: status FINISHED
2026-07-26 14:03:49 INFO None 5279900: status FINISHED
2026-07-26 14:03:49 INFO None 5279902: status FINISHED
2026-07-26 14:03:49 INFO None 5279903: status FINISHED
2026-07-26 14:03:49 INFO None 5279904: status FINISHED
2026-07-26 14:03:49 INFO None 5279905: status FINISHED
2026-07-26 14:03:49 INFO None 5279906: status FINISHED
2026-07-26 14:03:51 INFO None 5279907: status FINISHED
2026-07-26 14:03:51 INFO None 5279908: status FINISHED
2026-07-26 14:03:51 INFO None 5279909: status FINISHED
2026-07-26 14:03:51 INFO None 5279910: status FINISHED
2026-07-26 14:03:51 INFO None 5279912: status FINISHED
2026-07-26 14:03:52 INFO None 5279913: status FINISHED
2026-07-26 14:03:52 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:03:52 INFO Jobs still running: ['5279914']. Waiting...
2026-07-26 14:04:07 INFO None 5279898: status FINISHED
2026-07-26 14:04:07 INFO None 5279899: status FINISHED
2026-07-26 14:04:07 INFO None 5279900: status FINISHED
2026-07-26 14:04:07 INFO None 5279902: status FINISHED
2026-07-26 14:04:07 INFO None 5279903: status FINISHED
2026-07-26 14:04:07 INFO None 5279904: status FINISHED
2026-07-26 14:04:07 INFO None 5279905: status FINISHED
2026-07-26 14:04:07 INFO None 5279906: status FINISHED
2026-07-26 14:04:07 INFO None 5279907: status FINISHED
2026-07-26 14:04:07 INFO None 5279908: status FINISHED
2026-07-26 14:04:07 INFO None 5279909: status FINISHED
2026-07-26 14:04:07 INFO None 5279910: status FINISHED
2026-07-26 14:04:07 INFO None 5279912: status FINISHED
2026-07-26 14:04:07 INFO None 5279913: status FINISHED
2026-07-26 14:04:07 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:04:07 INFO Jobs still running: ['5279914']. Waiting...
2026-07-26 14:04:22 INFO None 5279898: status FINISHED
2026-07-26 14:04:22 INFO None 5279899: status FINISHED
2026-07-26 14:04:22 INFO None 5279900: status FINISHED
2026-07-26 14:04:22 INFO None 5279902: status FINISHED
2026-07-26 14:04:22 INFO None 5279903: status FINISHED
2026-07-26 14:04:22 INFO None 5279904: status FINISHED
2026-07-26 14:04:22 INFO None 5279905: status FINISHED
2026-07-26 14:04:22 INFO None 5279906: status FINISHED
2026-07-26 14:04:22 INFO None 5279907: status FINISHED
2026-07-26 14:04:22 INFO None 5279908: status FINISHED
2026-07-26 14:04:22 INFO None 5279909: status FINISHED
2026-07-26 14:04:22 INFO None 5279910: status FINISHED
2026-07-26 14:04:24 INFO None 5279912: status FINISHED
2026-07-26 14:04:24 INFO None 5279913: status FINISHED
2026-07-26 14:04:24 INFO None 5279914: status RUNNING/PENDING
2026-07-26 14:04:24 INFO Jobs still running: ['5279914']. Waiting...
2026-07-26 14:04:39 INFO None 5279898: status FINISHED
2026-07-26 14:04:39 INFO None 5279899: status FINISHED
2026-07-26 14:04:39 INFO None 5279900: status FINISHED
2026-07-26 14:04:39 INFO None 5279902: status FINISHED
2026-07-26 14:04:39 INFO None 5279903: status FINISHED
2026-07-26 14:04:39 INFO None 5279904: status FINISHED
2026-07-26 14:04:39 INFO None 5279905: status FINISHED
2026-07-26 14:04:39 INFO None 5279906: status FINISHED
2026-07-26 14:04:39 INFO None 5279907: status FINISHED
2026-07-26 14:04:39 INFO None 5279908: status FINISHED
2026-07-26 14:04:39 INFO None 5279909: status FINISHED
2026-07-26 14:04:39 INFO None 5279910: status FINISHED
2026-07-26 14:04:39 INFO None 5279912: status FINISHED
2026-07-26 14:04:39 INFO None 5279913: status FINISHED
2026-07-26 14:04:39 INFO None 5279914: status FINISHED
2026-07-26 14:04:39 INFO Jobs ['5279898', '5279899', '5279900', '5279902', '5279903', '5279904', '5279905', '5279906', '5279907', '5279908', '5279909', '5279910', '5279912', '5279913', '5279914'] have finished
2026-07-26 14:04:39 INFO Checking restart files were created ...
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-26 14:04:39 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-26 14:04:39 INFO  Run_model() completed successfully.
2026-07-26 14:04:39 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:04:39 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-26 14:04:39 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 14:04:39 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-26 14:04:39 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:04:39 INFO ---------->>> Running process_satellite_data()
2026-07-26 14:04:40 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-26 14:04:40 INFO ---------->>> Running run_obs_converter()
2026-07-26 14:04:40 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-26 14:04:40 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-26 14:04:40 INFO ---------->>> Running DART
2026-07-26 14:04:40 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 14:04:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-26 14:04:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-26 14:04:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-26 14:04:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-26 14:04:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-26 14:04:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-26 14:04:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-26 14:04:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-26 14:04:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-26 14:04:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-26 14:04:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-26 14:04:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-26 14:04:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-26 14:04:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-26 14:04:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-26 14:04:44 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 14:04:44 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 14:04:44 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 14:04:44 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 14:04:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 14:04:44 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 14:05:01 INFO Found: []
2026-07-26 14:05:01 INFO No job id returned by command ./run_filter.bsh
2026-07-26 14:05:01 INFO No monitoring will be performed
2026-07-26 14:05:01 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-26 14:05:01 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:01 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020712'
2026-07-26 14:05:02 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 14:05:02 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 14:05:02 INFO run_dart() is DONE.
2026-07-26 14:05:02 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 14:05:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-26 14:05:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-26 14:05:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-26 14:05:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:14 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-26 14:05:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-26 14:05:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-26 14:05:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-26 14:05:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-26 14:05:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-26 14:05:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-26 14:05:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-26 14:05:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-26 14:05:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-26 14:05:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-26 14:05:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-26 14:05:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:05:59 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 14:05:59 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:05:59 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:06:00 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-26 14:06:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:01 INFO Hourly dataset computed and listing created
2026-07-26 14:06:06 INFO Hourly dataset computed
2026-07-26 14:06:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:07 INFO Hourly dataset computed and listing created
2026-07-26 14:06:09 INFO Hourly dataset computed
2026-07-26 14:06:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:10 INFO Hourly dataset computed and listing created
2026-07-26 14:06:12 INFO Hourly dataset computed
2026-07-26 14:06:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:13 INFO Hourly dataset computed and listing created
2026-07-26 14:06:13 INFO Hourly dataset computed
2026-07-26 14:06:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:14 INFO Hourly dataset computed and listing created
2026-07-26 14:06:15 INFO Hourly dataset computed
2026-07-26 14:06:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:16 INFO Hourly dataset computed and listing created
2026-07-26 14:06:17 INFO Hourly dataset computed
2026-07-26 14:06:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:18 INFO Hourly dataset computed and listing created
2026-07-26 14:06:18 INFO Hourly dataset computed
2026-07-26 14:06:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:19 INFO Hourly dataset computed and listing created
2026-07-26 14:06:20 INFO Hourly dataset computed
2026-07-26 14:06:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:21 INFO Hourly dataset computed and listing created
2026-07-26 14:06:22 INFO Hourly dataset computed
2026-07-26 14:06:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:23 INFO Hourly dataset computed and listing created
2026-07-26 14:06:24 INFO Hourly dataset computed
2026-07-26 14:06:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:25 INFO Hourly dataset computed and listing created
2026-07-26 14:06:25 INFO Hourly dataset computed
2026-07-26 14:06:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:26 INFO Hourly dataset computed and listing created
2026-07-26 14:06:27 INFO Hourly dataset computed
2026-07-26 14:06:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:28 INFO Hourly dataset computed and listing created
2026-07-26 14:06:29 INFO Hourly dataset computed
2026-07-26 14:06:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:30 INFO Hourly dataset computed and listing created
2026-07-26 14:06:31 INFO Hourly dataset computed
2026-07-26 14:06:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:06:32 INFO Hourly dataset computed and listing created
2026-07-26 14:06:32 INFO Hourly dataset computed
2026-07-26 14:06:32 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-26 14:06:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:06:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 14:06:32 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-26 14:06:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 14:06:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:06:32 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 14:06:32 INFO Queuing job for member 1...
2026-07-26 14:06:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:06:32 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 14:06:35 INFO Found: ['5279938']
2026-07-26 14:06:40 INFO [TGCC-IRENE] Submitted job with ID:['5279938']
2026-07-26 14:06:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:06:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 14:06:40 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-26 14:06:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 14:06:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:06:40 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 14:06:40 INFO Queuing job for member 2...
2026-07-26 14:06:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:06:40 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 14:06:42 INFO Found: ['5279939']
2026-07-26 14:06:47 INFO [TGCC-IRENE] Submitted job with ID:['5279939']
2026-07-26 14:06:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:06:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 14:06:47 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-26 14:06:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 14:06:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:06:47 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 14:06:47 INFO Queuing job for member 3...
2026-07-26 14:06:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:06:47 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 14:06:50 INFO Found: ['5279940']
2026-07-26 14:06:55 INFO [TGCC-IRENE] Submitted job with ID:['5279940']
2026-07-26 14:06:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:06:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 14:06:55 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-26 14:06:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 14:06:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:06:55 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 14:06:55 INFO Queuing job for member 4...
2026-07-26 14:06:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:06:55 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 14:06:55 INFO Found: ['5279941']
2026-07-26 14:07:00 INFO [TGCC-IRENE] Submitted job with ID:['5279941']
2026-07-26 14:07:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 14:07:00 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-26 14:07:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 14:07:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:00 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 14:07:00 INFO Queuing job for member 5...
2026-07-26 14:07:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:00 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 14:07:01 INFO Found: ['5279942']
2026-07-26 14:07:06 INFO [TGCC-IRENE] Submitted job with ID:['5279942']
2026-07-26 14:07:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 14:07:06 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-26 14:07:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 14:07:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:06 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 14:07:06 INFO Queuing job for member 6...
2026-07-26 14:07:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:06 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 14:07:07 INFO Found: ['5279944']
2026-07-26 14:07:12 INFO [TGCC-IRENE] Submitted job with ID:['5279944']
2026-07-26 14:07:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 14:07:12 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-26 14:07:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 14:07:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:12 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 14:07:12 INFO Queuing job for member 7...
2026-07-26 14:07:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:12 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 14:07:13 INFO Found: ['5279945']
2026-07-26 14:07:18 INFO [TGCC-IRENE] Submitted job with ID:['5279945']
2026-07-26 14:07:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 14:07:18 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-26 14:07:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 14:07:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:18 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 14:07:18 INFO Queuing job for member 8...
2026-07-26 14:07:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:18 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 14:07:18 INFO Found: ['5279946']
2026-07-26 14:07:23 INFO [TGCC-IRENE] Submitted job with ID:['5279946']
2026-07-26 14:07:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 14:07:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-26 14:07:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 14:07:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 14:07:23 INFO Queuing job for member 9...
2026-07-26 14:07:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 14:07:24 INFO Found: ['5279947']
2026-07-26 14:07:29 INFO [TGCC-IRENE] Submitted job with ID:['5279947']
2026-07-26 14:07:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 14:07:29 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-26 14:07:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 14:07:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:29 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 14:07:29 INFO Queuing job for member 10...
2026-07-26 14:07:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:29 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 14:07:30 INFO Found: ['5279948']
2026-07-26 14:07:35 INFO [TGCC-IRENE] Submitted job with ID:['5279948']
2026-07-26 14:07:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 14:07:35 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-26 14:07:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 14:07:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:35 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 14:07:35 INFO Queuing job for member 11...
2026-07-26 14:07:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:35 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 14:07:36 INFO Found: ['5279949']
2026-07-26 14:07:41 INFO [TGCC-IRENE] Submitted job with ID:['5279949']
2026-07-26 14:07:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 14:07:41 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-26 14:07:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 14:07:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:41 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 14:07:41 INFO Queuing job for member 12...
2026-07-26 14:07:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:41 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 14:07:42 INFO Found: ['5279950']
2026-07-26 14:07:47 INFO [TGCC-IRENE] Submitted job with ID:['5279950']
2026-07-26 14:07:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 14:07:47 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-26 14:07:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 14:07:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:47 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 14:07:47 INFO Queuing job for member 13...
2026-07-26 14:07:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:47 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 14:07:50 INFO Found: ['5279951']
2026-07-26 14:07:55 INFO [TGCC-IRENE] Submitted job with ID:['5279951']
2026-07-26 14:07:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:07:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 14:07:55 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-26 14:07:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 14:07:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:07:55 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 14:07:55 INFO Queuing job for member 14...
2026-07-26 14:07:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:07:55 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 14:07:57 INFO Found: ['5279952']
2026-07-26 14:08:02 INFO [TGCC-IRENE] Submitted job with ID:['5279952']
2026-07-26 14:08:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:08:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 14:08:02 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-26 14:08:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 14:08:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:08:02 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 14:08:02 INFO Queuing job for member 15...
2026-07-26 14:08:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:08:02 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 14:08:05 INFO Found: ['5279954']
2026-07-26 14:08:10 INFO [TGCC-IRENE] Submitted job with ID:['5279954']
2026-07-26 14:08:10 INFO Checking job status ...
2026-07-26 14:08:10 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:08:10 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:08:10 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:08:25 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:08:25 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:08:27 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:08:27 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:08:42 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:08:42 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:08:43 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:08:43 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:08:58 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:08:58 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:08:59 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:08:59 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:09:14 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:09:14 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:09:14 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:09:14 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:09:14 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:09:15 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:09:15 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:09:30 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:09:30 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:09:30 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:09:46 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:09:46 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:09:47 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:09:47 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:09:47 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:09:47 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:09:47 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:09:47 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:10:02 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:10:02 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:10:04 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:10:04 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:10:19 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:10:19 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:10:19 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:10:34 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:10:34 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:10:36 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:10:37 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:10:37 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:10:37 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:10:52 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:10:52 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:10:52 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:11:07 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:11:07 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279940: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279942: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279946: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279947: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:11:09 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:11:09 INFO Jobs still running: ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:11:24 INFO None 5279938: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279939: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279940: status FINISHED
2026-07-26 14:11:24 INFO None 5279941: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279942: status FINISHED
2026-07-26 14:11:24 INFO None 5279944: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279945: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279946: status FINISHED
2026-07-26 14:11:24 INFO None 5279947: status FINISHED
2026-07-26 14:11:24 INFO None 5279948: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:11:24 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:11:24 INFO Jobs still running: ['5279938', '5279939', '5279941', '5279944', '5279945', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:11:39 INFO None 5279938: status FINISHED
2026-07-26 14:11:39 INFO None 5279939: status FINISHED
2026-07-26 14:11:40 INFO None 5279940: status FINISHED
2026-07-26 14:11:40 INFO None 5279941: status FINISHED
2026-07-26 14:11:40 INFO None 5279942: status FINISHED
2026-07-26 14:11:40 INFO None 5279944: status FINISHED
2026-07-26 14:11:40 INFO None 5279945: status FINISHED
2026-07-26 14:11:40 INFO None 5279946: status FINISHED
2026-07-26 14:11:40 INFO None 5279947: status FINISHED
2026-07-26 14:11:40 INFO None 5279948: status FINISHED
2026-07-26 14:11:40 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:11:40 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:11:40 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:11:40 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:11:40 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:11:40 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:11:56 INFO None 5279938: status FINISHED
2026-07-26 14:11:56 INFO None 5279939: status FINISHED
2026-07-26 14:11:56 INFO None 5279940: status FINISHED
2026-07-26 14:11:56 INFO None 5279941: status FINISHED
2026-07-26 14:11:56 INFO None 5279942: status FINISHED
2026-07-26 14:11:56 INFO None 5279944: status FINISHED
2026-07-26 14:11:56 INFO None 5279945: status FINISHED
2026-07-26 14:11:56 INFO None 5279946: status FINISHED
2026-07-26 14:11:56 INFO None 5279947: status FINISHED
2026-07-26 14:11:56 INFO None 5279948: status FINISHED
2026-07-26 14:11:56 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:11:56 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:11:56 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:11:56 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:11:56 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:11:56 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:12:11 INFO None 5279938: status FINISHED
2026-07-26 14:12:11 INFO None 5279939: status FINISHED
2026-07-26 14:12:11 INFO None 5279940: status FINISHED
2026-07-26 14:12:11 INFO None 5279941: status FINISHED
2026-07-26 14:12:11 INFO None 5279942: status FINISHED
2026-07-26 14:12:11 INFO None 5279944: status FINISHED
2026-07-26 14:12:13 INFO None 5279945: status FINISHED
2026-07-26 14:12:14 INFO None 5279946: status FINISHED
2026-07-26 14:12:14 INFO None 5279947: status FINISHED
2026-07-26 14:12:14 INFO None 5279948: status FINISHED
2026-07-26 14:12:14 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:12:14 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:12:14 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:12:14 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:12:14 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:12:14 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:12:29 INFO None 5279938: status FINISHED
2026-07-26 14:12:29 INFO None 5279939: status FINISHED
2026-07-26 14:12:29 INFO None 5279940: status FINISHED
2026-07-26 14:12:29 INFO None 5279941: status FINISHED
2026-07-26 14:12:29 INFO None 5279942: status FINISHED
2026-07-26 14:12:29 INFO None 5279944: status FINISHED
2026-07-26 14:12:29 INFO None 5279945: status FINISHED
2026-07-26 14:12:29 INFO None 5279946: status FINISHED
2026-07-26 14:12:29 INFO None 5279947: status FINISHED
2026-07-26 14:12:29 INFO None 5279948: status FINISHED
2026-07-26 14:12:29 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:12:31 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:12:31 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:12:31 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:12:31 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:12:31 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:12:46 INFO None 5279938: status FINISHED
2026-07-26 14:12:46 INFO None 5279939: status FINISHED
2026-07-26 14:12:46 INFO None 5279940: status FINISHED
2026-07-26 14:12:46 INFO None 5279941: status FINISHED
2026-07-26 14:12:46 INFO None 5279942: status FINISHED
2026-07-26 14:12:46 INFO None 5279944: status FINISHED
2026-07-26 14:12:46 INFO None 5279945: status FINISHED
2026-07-26 14:12:46 INFO None 5279946: status FINISHED
2026-07-26 14:12:46 INFO None 5279947: status FINISHED
2026-07-26 14:12:46 INFO None 5279948: status FINISHED
2026-07-26 14:12:46 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:12:46 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:12:46 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:12:46 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:12:46 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:12:46 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:13:01 INFO None 5279938: status FINISHED
2026-07-26 14:13:01 INFO None 5279939: status FINISHED
2026-07-26 14:13:01 INFO None 5279940: status FINISHED
2026-07-26 14:13:01 INFO None 5279941: status FINISHED
2026-07-26 14:13:01 INFO None 5279942: status FINISHED
2026-07-26 14:13:01 INFO None 5279944: status FINISHED
2026-07-26 14:13:02 INFO None 5279945: status FINISHED
2026-07-26 14:13:02 INFO None 5279946: status FINISHED
2026-07-26 14:13:02 INFO None 5279947: status FINISHED
2026-07-26 14:13:02 INFO None 5279948: status FINISHED
2026-07-26 14:13:02 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:13:02 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:13:02 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:13:04 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:13:04 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:13:04 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:13:19 INFO None 5279938: status FINISHED
2026-07-26 14:13:19 INFO None 5279939: status FINISHED
2026-07-26 14:13:19 INFO None 5279940: status FINISHED
2026-07-26 14:13:19 INFO None 5279941: status FINISHED
2026-07-26 14:13:19 INFO None 5279942: status FINISHED
2026-07-26 14:13:19 INFO None 5279944: status FINISHED
2026-07-26 14:13:19 INFO None 5279945: status FINISHED
2026-07-26 14:13:19 INFO None 5279946: status FINISHED
2026-07-26 14:13:19 INFO None 5279947: status FINISHED
2026-07-26 14:13:19 INFO None 5279948: status FINISHED
2026-07-26 14:13:19 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:13:19 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:13:19 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:13:19 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:13:19 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:13:19 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:13:34 INFO None 5279938: status FINISHED
2026-07-26 14:13:34 INFO None 5279939: status FINISHED
2026-07-26 14:13:34 INFO None 5279940: status FINISHED
2026-07-26 14:13:34 INFO None 5279941: status FINISHED
2026-07-26 14:13:34 INFO None 5279942: status FINISHED
2026-07-26 14:13:34 INFO None 5279944: status FINISHED
2026-07-26 14:13:34 INFO None 5279945: status FINISHED
2026-07-26 14:13:34 INFO None 5279946: status FINISHED
2026-07-26 14:13:34 INFO None 5279947: status FINISHED
2026-07-26 14:13:34 INFO None 5279948: status FINISHED
2026-07-26 14:13:34 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:13:34 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:13:34 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:13:34 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:13:34 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:13:34 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:13:49 INFO None 5279938: status FINISHED
2026-07-26 14:13:49 INFO None 5279939: status FINISHED
2026-07-26 14:13:49 INFO None 5279940: status FINISHED
2026-07-26 14:13:49 INFO None 5279941: status FINISHED
2026-07-26 14:13:49 INFO None 5279942: status FINISHED
2026-07-26 14:13:49 INFO None 5279944: status FINISHED
2026-07-26 14:13:49 INFO None 5279945: status FINISHED
2026-07-26 14:13:49 INFO None 5279946: status FINISHED
2026-07-26 14:13:50 INFO None 5279947: status FINISHED
2026-07-26 14:13:50 INFO None 5279948: status FINISHED
2026-07-26 14:13:50 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:13:50 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:13:50 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:13:50 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:13:50 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:13:50 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:14:05 INFO None 5279938: status FINISHED
2026-07-26 14:14:05 INFO None 5279939: status FINISHED
2026-07-26 14:14:05 INFO None 5279940: status FINISHED
2026-07-26 14:14:05 INFO None 5279941: status FINISHED
2026-07-26 14:14:05 INFO None 5279942: status FINISHED
2026-07-26 14:14:05 INFO None 5279944: status FINISHED
2026-07-26 14:14:05 INFO None 5279945: status FINISHED
2026-07-26 14:14:05 INFO None 5279946: status FINISHED
2026-07-26 14:14:05 INFO None 5279947: status FINISHED
2026-07-26 14:14:05 INFO None 5279948: status FINISHED
2026-07-26 14:14:05 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:14:05 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:14:05 INFO None 5279951: status RUNNING/PENDING
2026-07-26 14:14:05 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:14:05 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:14:05 INFO Jobs still running: ['5279949', '5279950', '5279951', '5279952', '5279954']. Waiting...
2026-07-26 14:14:20 INFO None 5279938: status FINISHED
2026-07-26 14:14:20 INFO None 5279939: status FINISHED
2026-07-26 14:14:20 INFO None 5279940: status FINISHED
2026-07-26 14:14:20 INFO None 5279941: status FINISHED
2026-07-26 14:14:20 INFO None 5279942: status FINISHED
2026-07-26 14:14:20 INFO None 5279944: status FINISHED
2026-07-26 14:14:20 INFO None 5279945: status FINISHED
2026-07-26 14:14:20 INFO None 5279946: status FINISHED
2026-07-26 14:14:20 INFO None 5279947: status FINISHED
2026-07-26 14:14:20 INFO None 5279948: status FINISHED
2026-07-26 14:14:20 INFO None 5279949: status RUNNING/PENDING
2026-07-26 14:14:20 INFO None 5279950: status RUNNING/PENDING
2026-07-26 14:14:20 INFO None 5279951: status FINISHED
2026-07-26 14:14:20 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:14:20 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:14:20 INFO Jobs still running: ['5279949', '5279950', '5279952', '5279954']. Waiting...
2026-07-26 14:14:35 INFO None 5279938: status FINISHED
2026-07-26 14:14:35 INFO None 5279939: status FINISHED
2026-07-26 14:14:35 INFO None 5279940: status FINISHED
2026-07-26 14:14:37 INFO None 5279941: status FINISHED
2026-07-26 14:14:37 INFO None 5279942: status FINISHED
2026-07-26 14:14:37 INFO None 5279944: status FINISHED
2026-07-26 14:14:37 INFO None 5279945: status FINISHED
2026-07-26 14:14:37 INFO None 5279946: status FINISHED
2026-07-26 14:14:37 INFO None 5279947: status FINISHED
2026-07-26 14:14:37 INFO None 5279948: status FINISHED
2026-07-26 14:14:37 INFO None 5279949: status FINISHED
2026-07-26 14:14:37 INFO None 5279950: status FINISHED
2026-07-26 14:14:37 INFO None 5279951: status FINISHED
2026-07-26 14:14:37 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:14:37 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:14:37 INFO Jobs still running: ['5279952', '5279954']. Waiting...
2026-07-26 14:14:52 INFO None 5279938: status FINISHED
2026-07-26 14:14:53 INFO None 5279939: status FINISHED
2026-07-26 14:14:53 INFO None 5279940: status FINISHED
2026-07-26 14:14:53 INFO None 5279941: status FINISHED
2026-07-26 14:14:53 INFO None 5279942: status FINISHED
2026-07-26 14:14:53 INFO None 5279944: status FINISHED
2026-07-26 14:14:53 INFO None 5279945: status FINISHED
2026-07-26 14:14:53 INFO None 5279946: status FINISHED
2026-07-26 14:14:53 INFO None 5279947: status FINISHED
2026-07-26 14:14:53 INFO None 5279948: status FINISHED
2026-07-26 14:14:53 INFO None 5279949: status FINISHED
2026-07-26 14:14:53 INFO None 5279950: status FINISHED
2026-07-26 14:14:53 INFO None 5279951: status FINISHED
2026-07-26 14:14:53 INFO None 5279952: status RUNNING/PENDING
2026-07-26 14:14:53 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:14:53 INFO Jobs still running: ['5279952', '5279954']. Waiting...
2026-07-26 14:15:08 INFO None 5279938: status FINISHED
2026-07-26 14:15:08 INFO None 5279939: status FINISHED
2026-07-26 14:15:08 INFO None 5279940: status FINISHED
2026-07-26 14:15:08 INFO None 5279941: status FINISHED
2026-07-26 14:15:08 INFO None 5279942: status FINISHED
2026-07-26 14:15:10 INFO None 5279944: status FINISHED
2026-07-26 14:15:10 INFO None 5279945: status FINISHED
2026-07-26 14:15:10 INFO None 5279946: status FINISHED
2026-07-26 14:15:10 INFO None 5279947: status FINISHED
2026-07-26 14:15:10 INFO None 5279948: status FINISHED
2026-07-26 14:15:10 INFO None 5279949: status FINISHED
2026-07-26 14:15:10 INFO None 5279950: status FINISHED
2026-07-26 14:15:10 INFO None 5279951: status FINISHED
2026-07-26 14:15:10 INFO None 5279952: status FINISHED
2026-07-26 14:15:10 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:15:10 INFO Jobs still running: ['5279954']. Waiting...
2026-07-26 14:15:25 INFO None 5279938: status FINISHED
2026-07-26 14:15:25 INFO None 5279939: status FINISHED
2026-07-26 14:15:25 INFO None 5279940: status FINISHED
2026-07-26 14:15:25 INFO None 5279941: status FINISHED
2026-07-26 14:15:25 INFO None 5279942: status FINISHED
2026-07-26 14:15:25 INFO None 5279944: status FINISHED
2026-07-26 14:15:25 INFO None 5279945: status FINISHED
2026-07-26 14:15:25 INFO None 5279946: status FINISHED
2026-07-26 14:15:25 INFO None 5279947: status FINISHED
2026-07-26 14:15:25 INFO None 5279948: status FINISHED
2026-07-26 14:15:25 INFO None 5279949: status FINISHED
2026-07-26 14:15:25 INFO None 5279950: status FINISHED
2026-07-26 14:15:25 INFO None 5279951: status FINISHED
2026-07-26 14:15:25 INFO None 5279952: status FINISHED
2026-07-26 14:15:25 INFO None 5279954: status RUNNING/PENDING
2026-07-26 14:15:25 INFO Jobs still running: ['5279954']. Waiting...
2026-07-26 14:15:40 INFO None 5279938: status FINISHED
2026-07-26 14:15:40 INFO None 5279939: status FINISHED
2026-07-26 14:15:40 INFO None 5279940: status FINISHED
2026-07-26 14:15:40 INFO None 5279941: status FINISHED
2026-07-26 14:15:40 INFO None 5279942: status FINISHED
2026-07-26 14:15:40 INFO None 5279944: status FINISHED
2026-07-26 14:15:40 INFO None 5279945: status FINISHED
2026-07-26 14:15:40 INFO None 5279946: status FINISHED
2026-07-26 14:15:40 INFO None 5279947: status FINISHED
2026-07-26 14:15:40 INFO None 5279948: status FINISHED
2026-07-26 14:15:41 INFO None 5279949: status FINISHED
2026-07-26 14:15:41 INFO None 5279950: status FINISHED
2026-07-26 14:15:41 INFO None 5279951: status FINISHED
2026-07-26 14:15:41 INFO None 5279952: status FINISHED
2026-07-26 14:15:41 INFO None 5279954: status FINISHED
2026-07-26 14:15:41 INFO Jobs ['5279938', '5279939', '5279940', '5279941', '5279942', '5279944', '5279945', '5279946', '5279947', '5279948', '5279949', '5279950', '5279951', '5279952', '5279954'] have finished
2026-07-26 14:15:41 INFO Checking restart files were created ...
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc(1002685915 bytes)
2026-07-26 14:15:41 INFO  Run_model() completed successfully.
2026-07-26 14:15:41 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:15:41 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 14:00:00 days=153073 seconds=50400
2026-07-26 14:15:41 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 14:15:41 INFO [TIME] increment current_time 2020-02-07 12:00:00 -> 2020-02-07 14:00:00
2026-07-26 14:15:41 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:15:41 INFO ---------->>> Running process_satellite_data()
2026-07-26 14:15:41 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12019.nc
2026-07-26 14:15:41 INFO ---------->>> Running run_obs_converter()
2026-07-26 14:15:41 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-26 14:15:41 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-26 14:15:41 INFO ---------->>> Running DART
2026-07-26 14:15:41 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 14:15:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020714_1_out_toDART.nc
2026-07-26 14:15:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020714_1_out_toDART.nc
2026-07-26 14:15:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020714_1_out_toDART.nc
2026-07-26 14:15:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020714_1_out_toDART.nc
2026-07-26 14:15:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020714_1_out_toDART.nc
2026-07-26 14:15:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020714_1_out_toDART.nc
2026-07-26 14:15:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020714_1_out_toDART.nc
2026-07-26 14:15:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020714_1_out_toDART.nc
2026-07-26 14:15:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020714_1_out_toDART.nc
2026-07-26 14:15:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020714_1_out_toDART.nc
2026-07-26 14:15:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020714_1_out_toDART.nc
2026-07-26 14:15:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020714_1_out_toDART.nc
2026-07-26 14:15:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020714_1_out_toDART.nc
2026-07-26 14:15:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020714_1_out_toDART.nc
2026-07-26 14:15:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020714_1_out_toDART.nc
2026-07-26 14:15:46 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 14:15:46 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 14:15:46 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 14:15:46 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 14:15:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 14:15:46 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 14:15:56 INFO Found: []
2026-07-26 14:15:56 INFO No job id returned by command ./run_filter.bsh
2026-07-26 14:15:56 INFO No monitoring will be performed
2026-07-26 14:15:56 INFO Moving DART output files to analysis and preassim directories for date 2020020714 if present ...
2026-07-26 14:15:56 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020714'
2026-07-26 14:15:56 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 14:15:57 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 14:15:57 INFO run_dart() is DONE.
2026-07-26 14:15:57 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 14:15:57 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-26 14:16:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-26 14:16:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-26 14:16:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-26 14:16:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-26 14:16:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-26 14:16:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-26 14:16:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-26 14:16:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-26 14:16:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-26 14:16:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-26 14:16:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:16:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-26 14:17:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:17:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-26 14:17:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:17:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-26 14:17:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:17:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-26 14:17:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 14:17:17 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 14:17:17 INFO [TIME] step_end current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:17:17 INFO [TIME] step_start current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:17:17 INFO [TIME] window start=2020-02-07 14:00:00 end=2020-02-08 00:00:00 run_hours=10 has_assimilation=False
2026-07-26 14:17:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:17:19 INFO Hourly dataset computed and listing created
2026-07-26 14:17:37 INFO Hourly dataset computed
2026-07-26 14:17:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:17:38 INFO Hourly dataset computed and listing created
2026-07-26 14:17:51 INFO Hourly dataset computed
2026-07-26 14:17:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:17:52 INFO Hourly dataset computed and listing created
2026-07-26 14:18:05 INFO Hourly dataset computed
2026-07-26 14:18:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:18:07 INFO Hourly dataset computed and listing created
2026-07-26 14:18:19 INFO Hourly dataset computed
2026-07-26 14:18:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:18:20 INFO Hourly dataset computed and listing created
2026-07-26 14:18:32 INFO Hourly dataset computed
2026-07-26 14:18:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:18:34 INFO Hourly dataset computed and listing created
2026-07-26 14:18:46 INFO Hourly dataset computed
2026-07-26 14:18:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:18:47 INFO Hourly dataset computed and listing created
2026-07-26 14:18:57 INFO Hourly dataset computed
2026-07-26 14:18:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:18:58 INFO Hourly dataset computed and listing created
2026-07-26 14:19:08 INFO Hourly dataset computed
2026-07-26 14:19:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:19:09 INFO Hourly dataset computed and listing created
2026-07-26 14:19:19 INFO Hourly dataset computed
2026-07-26 14:19:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:19:20 INFO Hourly dataset computed and listing created
2026-07-26 14:19:29 INFO Hourly dataset computed
2026-07-26 14:19:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:19:30 INFO Hourly dataset computed and listing created
2026-07-26 14:19:40 INFO Hourly dataset computed
2026-07-26 14:19:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:19:41 INFO Hourly dataset computed and listing created
2026-07-26 14:19:52 INFO Hourly dataset computed
2026-07-26 14:19:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:19:53 INFO Hourly dataset computed and listing created
2026-07-26 14:20:02 INFO Hourly dataset computed
2026-07-26 14:20:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:20:04 INFO Hourly dataset computed and listing created
2026-07-26 14:20:13 INFO Hourly dataset computed
2026-07-26 14:20:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 14:20:15 INFO Hourly dataset computed and listing created
2026-07-26 14:20:24 INFO Hourly dataset computed
2026-07-26 14:20:24 INFO ---------->>> Running CHIMERE model from 2020-02-07 14:00:00 to 2020-02-08 00:00:00
2026-07-26 14:20:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 14:20:24 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-26 14:20:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 14:20:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:24 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 14:20:24 INFO Queuing job for member 1...
2026-07-26 14:20:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:24 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 14:20:25 INFO Found: ['5279980']
2026-07-26 14:20:30 INFO [TGCC-IRENE] Submitted job with ID:['5279980']
2026-07-26 14:20:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 14:20:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-26 14:20:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 14:20:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 14:20:30 INFO Queuing job for member 2...
2026-07-26 14:20:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 14:20:31 INFO Found: ['5279981']
2026-07-26 14:20:36 INFO [TGCC-IRENE] Submitted job with ID:['5279981']
2026-07-26 14:20:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 14:20:36 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-26 14:20:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 14:20:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:36 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 14:20:36 INFO Queuing job for member 3...
2026-07-26 14:20:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:36 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 14:20:38 INFO Found: ['5279982']
2026-07-26 14:20:43 INFO [TGCC-IRENE] Submitted job with ID:['5279982']
2026-07-26 14:20:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 14:20:43 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-26 14:20:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 14:20:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:43 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 14:20:43 INFO Queuing job for member 4...
2026-07-26 14:20:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:43 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 14:20:46 INFO Found: ['5279983']
2026-07-26 14:20:51 INFO [TGCC-IRENE] Submitted job with ID:['5279983']
2026-07-26 14:20:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 14:20:51 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-26 14:20:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 14:20:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:51 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 14:20:51 INFO Queuing job for member 5...
2026-07-26 14:20:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:51 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 14:20:53 INFO Found: ['5279984']
2026-07-26 14:20:58 INFO [TGCC-IRENE] Submitted job with ID:['5279984']
2026-07-26 14:20:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:20:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 14:20:58 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-26 14:20:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 14:20:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:20:58 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 14:20:58 INFO Queuing job for member 6...
2026-07-26 14:20:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:20:58 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 14:21:01 INFO Found: ['5279985']
2026-07-26 14:21:06 INFO [TGCC-IRENE] Submitted job with ID:['5279985']
2026-07-26 14:21:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 14:21:06 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-26 14:21:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 14:21:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:06 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 14:21:06 INFO Queuing job for member 7...
2026-07-26 14:21:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:06 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 14:21:08 INFO Found: ['5279989']
2026-07-26 14:21:13 INFO [TGCC-IRENE] Submitted job with ID:['5279989']
2026-07-26 14:21:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 14:21:13 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-26 14:21:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 14:21:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:13 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 14:21:13 INFO Queuing job for member 8...
2026-07-26 14:21:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:13 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 14:21:16 INFO Found: ['5279991']
2026-07-26 14:21:21 INFO [TGCC-IRENE] Submitted job with ID:['5279991']
2026-07-26 14:21:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 14:21:21 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-26 14:21:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 14:21:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:21 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 14:21:21 INFO Queuing job for member 9...
2026-07-26 14:21:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:21 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 14:21:23 INFO Found: ['5279993']
2026-07-26 14:21:28 INFO [TGCC-IRENE] Submitted job with ID:['5279993']
2026-07-26 14:21:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 14:21:28 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-26 14:21:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 14:21:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:28 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 14:21:28 INFO Queuing job for member 10...
2026-07-26 14:21:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:28 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 14:21:31 INFO Found: ['5279994']
2026-07-26 14:21:36 INFO [TGCC-IRENE] Submitted job with ID:['5279994']
2026-07-26 14:21:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 14:21:36 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-26 14:21:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 14:21:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:36 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 14:21:36 INFO Queuing job for member 11...
2026-07-26 14:21:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:36 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 14:21:38 INFO Found: ['5279995']
2026-07-26 14:21:43 INFO [TGCC-IRENE] Submitted job with ID:['5279995']
2026-07-26 14:21:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 14:21:43 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-26 14:21:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 14:21:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:43 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 14:21:43 INFO Queuing job for member 12...
2026-07-26 14:21:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:43 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 14:21:46 INFO Found: ['5279996']
2026-07-26 14:21:51 INFO [TGCC-IRENE] Submitted job with ID:['5279996']
2026-07-26 14:21:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 14:21:51 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-26 14:21:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 14:21:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:51 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 14:21:51 INFO Queuing job for member 13...
2026-07-26 14:21:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:51 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 14:21:53 INFO Found: ['5279997']
2026-07-26 14:21:58 INFO [TGCC-IRENE] Submitted job with ID:['5279997']
2026-07-26 14:21:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:21:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 14:21:58 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-26 14:21:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 14:21:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:21:58 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 14:21:58 INFO Queuing job for member 14...
2026-07-26 14:21:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:21:58 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 14:22:01 INFO Found: ['5279998']
2026-07-26 14:22:06 INFO [TGCC-IRENE] Submitted job with ID:['5279998']
2026-07-26 14:22:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 14:22:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 14:22:06 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-26 14:22:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 14:22:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 14:22:06 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 14:22:06 INFO Queuing job for member 15...
2026-07-26 14:22:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 14:22:06 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 14:22:07 INFO Found: ['5280000']
2026-07-26 14:22:12 INFO [TGCC-IRENE] Submitted job with ID:['5280000']
2026-07-26 14:22:12 INFO Checking job status ...
2026-07-26 14:22:12 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:22:12 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:22:12 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:22:27 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:22:27 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:22:27 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:22:42 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:22:42 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:22:42 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:22:59 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:22:59 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:22:59 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:23:14 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:23:14 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:23:15 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:23:17 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:23:17 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:23:32 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:23:32 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:23:32 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:23:47 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:23:47 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:23:47 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:23:47 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:23:47 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:23:47 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:23:49 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:23:49 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:23:50 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:23:50 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:24:05 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:24:05 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:24:07 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:24:07 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:24:07 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:24:22 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:24:22 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:24:22 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:24:37 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:24:37 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:24:38 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:24:38 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:24:53 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:24:53 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:24:53 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:25:08 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:25:08 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:25:08 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:25:24 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:25:24 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:25:26 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:25:26 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:25:41 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:25:41 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:25:41 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:25:56 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:25:56 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:25:58 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:25:58 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:25:58 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:26:13 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:26:13 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:26:13 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:26:14 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:26:14 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:26:29 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:26:29 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:26:31 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:26:31 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:26:31 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:26:31 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:26:31 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:26:46 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:26:46 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:26:46 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:27:01 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:27:01 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:27:01 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:27:18 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:27:18 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:27:18 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:27:33 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:27:33 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:27:35 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:27:35 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:27:50 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:27:50 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:27:50 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:28:05 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:28:05 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:28:05 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:28:06 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:28:08 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:28:08 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:28:23 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:28:23 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:28:23 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:28:38 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:28:38 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:28:40 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:28:40 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:28:55 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:28:55 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:28:56 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:28:56 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:29:11 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:29:11 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:29:11 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:29:26 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:29:26 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:29:26 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:29:42 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:29:42 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:29:42 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:29:57 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279981: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:29:57 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:29:58 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:29:58 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:30:00 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:30:00 INFO Jobs still running: ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:30:15 INFO None 5279980: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279981: status FINISHED
2026-07-26 14:30:15 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:30:15 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:30:15 INFO Jobs still running: ['5279980', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:30:32 INFO None 5279980: status FINISHED
2026-07-26 14:30:32 INFO None 5279981: status FINISHED
2026-07-26 14:30:32 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:30:32 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:30:32 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:30:48 INFO None 5279980: status FINISHED
2026-07-26 14:30:48 INFO None 5279981: status FINISHED
2026-07-26 14:30:48 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:30:48 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:30:50 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:30:50 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:30:50 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:31:05 INFO None 5279980: status FINISHED
2026-07-26 14:31:05 INFO None 5279981: status FINISHED
2026-07-26 14:31:05 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:31:05 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:31:05 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:31:20 INFO None 5279980: status FINISHED
2026-07-26 14:31:20 INFO None 5279981: status FINISHED
2026-07-26 14:31:20 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:31:20 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:31:20 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:31:37 INFO None 5279980: status FINISHED
2026-07-26 14:31:37 INFO None 5279981: status FINISHED
2026-07-26 14:31:37 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:31:37 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:31:37 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:31:52 INFO None 5279980: status FINISHED
2026-07-26 14:31:52 INFO None 5279981: status FINISHED
2026-07-26 14:31:52 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:31:52 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:31:54 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:31:54 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:31:54 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:32:09 INFO None 5279980: status FINISHED
2026-07-26 14:32:09 INFO None 5279981: status FINISHED
2026-07-26 14:32:09 INFO None 5279982: status RUNNING/PENDING
2026-07-26 14:32:09 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:32:09 INFO None 5279984: status RUNNING/PENDING
2026-07-26 14:32:09 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:32:10 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:32:10 INFO Jobs still running: ['5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:32:25 INFO None 5279980: status FINISHED
2026-07-26 14:32:25 INFO None 5279981: status FINISHED
2026-07-26 14:32:25 INFO None 5279982: status FINISHED
2026-07-26 14:32:25 INFO None 5279983: status RUNNING/PENDING
2026-07-26 14:32:25 INFO None 5279984: status FINISHED
2026-07-26 14:32:27 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:32:27 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:32:27 INFO Jobs still running: ['5279983', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:32:42 INFO None 5279980: status FINISHED
2026-07-26 14:32:42 INFO None 5279981: status FINISHED
2026-07-26 14:32:42 INFO None 5279982: status FINISHED
2026-07-26 14:32:42 INFO None 5279983: status FINISHED
2026-07-26 14:32:42 INFO None 5279984: status FINISHED
2026-07-26 14:32:42 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:32:42 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:32:42 INFO Jobs still running: ['5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:32:57 INFO None 5279980: status FINISHED
2026-07-26 14:32:57 INFO None 5279981: status FINISHED
2026-07-26 14:32:57 INFO None 5279982: status FINISHED
2026-07-26 14:32:57 INFO None 5279983: status FINISHED
2026-07-26 14:32:57 INFO None 5279984: status FINISHED
2026-07-26 14:32:57 INFO None 5279985: status RUNNING/PENDING
2026-07-26 14:32:57 INFO None 5279989: status RUNNING/PENDING
2026-07-26 14:32:57 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:32:57 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:32:57 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:32:57 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:32:58 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:32:58 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:32:58 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:32:58 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:32:58 INFO Jobs still running: ['5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:33:13 INFO None 5279980: status FINISHED
2026-07-26 14:33:13 INFO None 5279981: status FINISHED
2026-07-26 14:33:13 INFO None 5279982: status FINISHED
2026-07-26 14:33:13 INFO None 5279983: status FINISHED
2026-07-26 14:33:13 INFO None 5279984: status FINISHED
2026-07-26 14:33:13 INFO None 5279985: status FINISHED
2026-07-26 14:33:13 INFO None 5279989: status FINISHED
2026-07-26 14:33:13 INFO None 5279991: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279993: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:33:13 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:33:13 INFO Jobs still running: ['5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:33:28 INFO None 5279980: status FINISHED
2026-07-26 14:33:28 INFO None 5279981: status FINISHED
2026-07-26 14:33:28 INFO None 5279982: status FINISHED
2026-07-26 14:33:28 INFO None 5279983: status FINISHED
2026-07-26 14:33:28 INFO None 5279984: status FINISHED
2026-07-26 14:33:28 INFO None 5279985: status FINISHED
2026-07-26 14:33:28 INFO None 5279989: status FINISHED
2026-07-26 14:33:28 INFO None 5279991: status FINISHED
2026-07-26 14:33:28 INFO None 5279993: status FINISHED
2026-07-26 14:33:28 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:33:28 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:33:28 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:33:28 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:33:28 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:33:28 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:33:28 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:33:44 INFO None 5279980: status FINISHED
2026-07-26 14:33:44 INFO None 5279981: status FINISHED
2026-07-26 14:33:44 INFO None 5279982: status FINISHED
2026-07-26 14:33:44 INFO None 5279983: status FINISHED
2026-07-26 14:33:44 INFO None 5279984: status FINISHED
2026-07-26 14:33:44 INFO None 5279985: status FINISHED
2026-07-26 14:33:44 INFO None 5279989: status FINISHED
2026-07-26 14:33:44 INFO None 5279991: status FINISHED
2026-07-26 14:33:44 INFO None 5279993: status FINISHED
2026-07-26 14:33:44 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:33:44 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:33:44 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:33:44 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:33:44 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:33:44 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:33:44 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:33:59 INFO None 5279980: status FINISHED
2026-07-26 14:33:59 INFO None 5279981: status FINISHED
2026-07-26 14:33:59 INFO None 5279982: status FINISHED
2026-07-26 14:33:59 INFO None 5279983: status FINISHED
2026-07-26 14:33:59 INFO None 5279984: status FINISHED
2026-07-26 14:33:59 INFO None 5279985: status FINISHED
2026-07-26 14:33:59 INFO None 5279989: status FINISHED
2026-07-26 14:33:59 INFO None 5279991: status FINISHED
2026-07-26 14:33:59 INFO None 5279993: status FINISHED
2026-07-26 14:33:59 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:33:59 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:33:59 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:34:01 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:34:01 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:34:01 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:34:01 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:34:16 INFO None 5279980: status FINISHED
2026-07-26 14:34:16 INFO None 5279981: status FINISHED
2026-07-26 14:34:16 INFO None 5279982: status FINISHED
2026-07-26 14:34:16 INFO None 5279983: status FINISHED
2026-07-26 14:34:16 INFO None 5279984: status FINISHED
2026-07-26 14:34:17 INFO None 5279985: status FINISHED
2026-07-26 14:34:17 INFO None 5279989: status FINISHED
2026-07-26 14:34:17 INFO None 5279991: status FINISHED
2026-07-26 14:34:17 INFO None 5279993: status FINISHED
2026-07-26 14:34:17 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:34:17 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:34:17 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:34:17 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:34:17 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:34:17 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:34:17 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:34:32 INFO None 5279980: status FINISHED
2026-07-26 14:34:32 INFO None 5279981: status FINISHED
2026-07-26 14:34:32 INFO None 5279982: status FINISHED
2026-07-26 14:34:32 INFO None 5279983: status FINISHED
2026-07-26 14:34:32 INFO None 5279984: status FINISHED
2026-07-26 14:34:32 INFO None 5279985: status FINISHED
2026-07-26 14:34:32 INFO None 5279989: status FINISHED
2026-07-26 14:34:32 INFO None 5279991: status FINISHED
2026-07-26 14:34:32 INFO None 5279993: status FINISHED
2026-07-26 14:34:32 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:34:32 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:34:32 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:34:32 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:34:32 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:34:32 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:34:32 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:34:47 INFO None 5279980: status FINISHED
2026-07-26 14:34:47 INFO None 5279981: status FINISHED
2026-07-26 14:34:47 INFO None 5279982: status FINISHED
2026-07-26 14:34:49 INFO None 5279983: status FINISHED
2026-07-26 14:34:49 INFO None 5279984: status FINISHED
2026-07-26 14:34:49 INFO None 5279985: status FINISHED
2026-07-26 14:34:49 INFO None 5279989: status FINISHED
2026-07-26 14:34:49 INFO None 5279991: status FINISHED
2026-07-26 14:34:49 INFO None 5279993: status FINISHED
2026-07-26 14:34:49 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:34:49 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:34:49 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:34:49 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:34:49 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:34:49 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:34:49 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:35:04 INFO None 5279980: status FINISHED
2026-07-26 14:35:04 INFO None 5279981: status FINISHED
2026-07-26 14:35:04 INFO None 5279982: status FINISHED
2026-07-26 14:35:04 INFO None 5279983: status FINISHED
2026-07-26 14:35:04 INFO None 5279984: status FINISHED
2026-07-26 14:35:04 INFO None 5279985: status FINISHED
2026-07-26 14:35:04 INFO None 5279989: status FINISHED
2026-07-26 14:35:04 INFO None 5279991: status FINISHED
2026-07-26 14:35:04 INFO None 5279993: status FINISHED
2026-07-26 14:35:04 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:35:04 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:35:04 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:35:05 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:35:05 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:35:05 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:35:05 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:35:20 INFO None 5279980: status FINISHED
2026-07-26 14:35:20 INFO None 5279981: status FINISHED
2026-07-26 14:35:20 INFO None 5279982: status FINISHED
2026-07-26 14:35:20 INFO None 5279983: status FINISHED
2026-07-26 14:35:20 INFO None 5279984: status FINISHED
2026-07-26 14:35:20 INFO None 5279985: status FINISHED
2026-07-26 14:35:20 INFO None 5279989: status FINISHED
2026-07-26 14:35:20 INFO None 5279991: status FINISHED
2026-07-26 14:35:20 INFO None 5279993: status FINISHED
2026-07-26 14:35:20 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:35:20 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:35:20 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:35:20 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:35:20 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:35:20 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:35:20 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:35:35 INFO None 5279980: status FINISHED
2026-07-26 14:35:35 INFO None 5279981: status FINISHED
2026-07-26 14:35:35 INFO None 5279982: status FINISHED
2026-07-26 14:35:35 INFO None 5279983: status FINISHED
2026-07-26 14:35:35 INFO None 5279984: status FINISHED
2026-07-26 14:35:35 INFO None 5279985: status FINISHED
2026-07-26 14:35:35 INFO None 5279989: status FINISHED
2026-07-26 14:35:35 INFO None 5279991: status FINISHED
2026-07-26 14:35:35 INFO None 5279993: status FINISHED
2026-07-26 14:35:35 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:35:35 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:35:35 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:35:35 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:35:35 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:35:35 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:35:35 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:35:52 INFO None 5279980: status FINISHED
2026-07-26 14:35:52 INFO None 5279981: status FINISHED
2026-07-26 14:35:52 INFO None 5279982: status FINISHED
2026-07-26 14:35:52 INFO None 5279983: status FINISHED
2026-07-26 14:35:52 INFO None 5279984: status FINISHED
2026-07-26 14:35:52 INFO None 5279985: status FINISHED
2026-07-26 14:35:52 INFO None 5279989: status FINISHED
2026-07-26 14:35:52 INFO None 5279991: status FINISHED
2026-07-26 14:35:52 INFO None 5279993: status FINISHED
2026-07-26 14:35:52 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:35:52 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:35:52 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:35:52 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:35:52 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:35:52 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:35:52 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:36:07 INFO None 5279980: status FINISHED
2026-07-26 14:36:07 INFO None 5279981: status FINISHED
2026-07-26 14:36:07 INFO None 5279982: status FINISHED
2026-07-26 14:36:07 INFO None 5279983: status FINISHED
2026-07-26 14:36:07 INFO None 5279984: status FINISHED
2026-07-26 14:36:07 INFO None 5279985: status FINISHED
2026-07-26 14:36:07 INFO None 5279989: status FINISHED
2026-07-26 14:36:07 INFO None 5279991: status FINISHED
2026-07-26 14:36:07 INFO None 5279993: status FINISHED
2026-07-26 14:36:07 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:36:07 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:36:07 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:36:07 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:36:07 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:36:07 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:36:07 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:36:22 INFO None 5279980: status FINISHED
2026-07-26 14:36:22 INFO None 5279981: status FINISHED
2026-07-26 14:36:22 INFO None 5279982: status FINISHED
2026-07-26 14:36:22 INFO None 5279983: status FINISHED
2026-07-26 14:36:24 INFO None 5279984: status FINISHED
2026-07-26 14:36:25 INFO None 5279985: status FINISHED
2026-07-26 14:36:25 INFO None 5279989: status FINISHED
2026-07-26 14:36:25 INFO None 5279991: status FINISHED
2026-07-26 14:36:25 INFO None 5279993: status FINISHED
2026-07-26 14:36:25 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:36:25 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:36:25 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:36:25 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:36:25 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:36:25 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:36:25 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:36:40 INFO None 5279980: status FINISHED
2026-07-26 14:36:40 INFO None 5279981: status FINISHED
2026-07-26 14:36:40 INFO None 5279982: status FINISHED
2026-07-26 14:36:40 INFO None 5279983: status FINISHED
2026-07-26 14:36:40 INFO None 5279984: status FINISHED
2026-07-26 14:36:40 INFO None 5279985: status FINISHED
2026-07-26 14:36:40 INFO None 5279989: status FINISHED
2026-07-26 14:36:40 INFO None 5279991: status FINISHED
2026-07-26 14:36:40 INFO None 5279993: status FINISHED
2026-07-26 14:36:40 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:36:40 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:36:40 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:36:40 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:36:40 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:36:40 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:36:40 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:36:55 INFO None 5279980: status FINISHED
2026-07-26 14:36:55 INFO None 5279981: status FINISHED
2026-07-26 14:36:55 INFO None 5279982: status FINISHED
2026-07-26 14:36:55 INFO None 5279983: status FINISHED
2026-07-26 14:36:55 INFO None 5279984: status FINISHED
2026-07-26 14:36:57 INFO None 5279985: status FINISHED
2026-07-26 14:36:57 INFO None 5279989: status FINISHED
2026-07-26 14:36:57 INFO None 5279991: status FINISHED
2026-07-26 14:36:57 INFO None 5279993: status FINISHED
2026-07-26 14:36:57 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:36:57 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:36:57 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:36:57 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:36:57 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:36:57 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:36:57 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:37:12 INFO None 5279980: status FINISHED
2026-07-26 14:37:12 INFO None 5279981: status FINISHED
2026-07-26 14:37:12 INFO None 5279982: status FINISHED
2026-07-26 14:37:12 INFO None 5279983: status FINISHED
2026-07-26 14:37:12 INFO None 5279984: status FINISHED
2026-07-26 14:37:12 INFO None 5279985: status FINISHED
2026-07-26 14:37:12 INFO None 5279989: status FINISHED
2026-07-26 14:37:12 INFO None 5279991: status FINISHED
2026-07-26 14:37:12 INFO None 5279993: status FINISHED
2026-07-26 14:37:12 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:37:12 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:37:12 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:37:12 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:37:12 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:37:13 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:37:13 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:37:28 INFO None 5279980: status FINISHED
2026-07-26 14:37:28 INFO None 5279981: status FINISHED
2026-07-26 14:37:28 INFO None 5279982: status FINISHED
2026-07-26 14:37:28 INFO None 5279983: status FINISHED
2026-07-26 14:37:28 INFO None 5279984: status FINISHED
2026-07-26 14:37:28 INFO None 5279985: status FINISHED
2026-07-26 14:37:28 INFO None 5279989: status FINISHED
2026-07-26 14:37:28 INFO None 5279991: status FINISHED
2026-07-26 14:37:28 INFO None 5279993: status FINISHED
2026-07-26 14:37:28 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:37:28 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:37:28 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:37:28 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:37:28 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:37:28 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:37:28 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:37:43 INFO None 5279980: status FINISHED
2026-07-26 14:37:43 INFO None 5279981: status FINISHED
2026-07-26 14:37:43 INFO None 5279982: status FINISHED
2026-07-26 14:37:43 INFO None 5279983: status FINISHED
2026-07-26 14:37:43 INFO None 5279984: status FINISHED
2026-07-26 14:37:43 INFO None 5279985: status FINISHED
2026-07-26 14:37:43 INFO None 5279989: status FINISHED
2026-07-26 14:37:43 INFO None 5279991: status FINISHED
2026-07-26 14:37:43 INFO None 5279993: status FINISHED
2026-07-26 14:37:43 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:37:43 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:37:43 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:37:43 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:37:43 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:37:43 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:37:43 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:37:59 INFO None 5279980: status FINISHED
2026-07-26 14:37:59 INFO None 5279981: status FINISHED
2026-07-26 14:37:59 INFO None 5279982: status FINISHED
2026-07-26 14:37:59 INFO None 5279983: status FINISHED
2026-07-26 14:37:59 INFO None 5279984: status FINISHED
2026-07-26 14:37:59 INFO None 5279985: status FINISHED
2026-07-26 14:37:59 INFO None 5279989: status FINISHED
2026-07-26 14:37:59 INFO None 5279991: status FINISHED
2026-07-26 14:37:59 INFO None 5279993: status FINISHED
2026-07-26 14:37:59 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:37:59 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:37:59 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:37:59 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:37:59 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:37:59 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:37:59 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:38:14 INFO None 5279980: status FINISHED
2026-07-26 14:38:14 INFO None 5279981: status FINISHED
2026-07-26 14:38:14 INFO None 5279982: status FINISHED
2026-07-26 14:38:14 INFO None 5279983: status FINISHED
2026-07-26 14:38:14 INFO None 5279984: status FINISHED
2026-07-26 14:38:14 INFO None 5279985: status FINISHED
2026-07-26 14:38:14 INFO None 5279989: status FINISHED
2026-07-26 14:38:14 INFO None 5279991: status FINISHED
2026-07-26 14:38:14 INFO None 5279993: status FINISHED
2026-07-26 14:38:14 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:38:14 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:38:14 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:38:14 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:38:15 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:38:15 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:38:15 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:38:30 INFO None 5279980: status FINISHED
2026-07-26 14:38:30 INFO None 5279981: status FINISHED
2026-07-26 14:38:30 INFO None 5279982: status FINISHED
2026-07-26 14:38:32 INFO None 5279983: status FINISHED
2026-07-26 14:38:32 INFO None 5279984: status FINISHED
2026-07-26 14:38:32 INFO None 5279985: status FINISHED
2026-07-26 14:38:32 INFO None 5279989: status FINISHED
2026-07-26 14:38:32 INFO None 5279991: status FINISHED
2026-07-26 14:38:32 INFO None 5279993: status FINISHED
2026-07-26 14:38:32 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:38:32 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:38:32 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:38:32 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:38:32 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:38:32 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:38:32 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:38:47 INFO None 5279980: status FINISHED
2026-07-26 14:38:47 INFO None 5279981: status FINISHED
2026-07-26 14:38:47 INFO None 5279982: status FINISHED
2026-07-26 14:38:47 INFO None 5279983: status FINISHED
2026-07-26 14:38:47 INFO None 5279984: status FINISHED
2026-07-26 14:38:47 INFO None 5279985: status FINISHED
2026-07-26 14:38:47 INFO None 5279989: status FINISHED
2026-07-26 14:38:47 INFO None 5279991: status FINISHED
2026-07-26 14:38:47 INFO None 5279993: status FINISHED
2026-07-26 14:38:47 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:38:47 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:38:47 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:38:47 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:38:47 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:38:47 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:38:47 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:39:02 INFO None 5279980: status FINISHED
2026-07-26 14:39:02 INFO None 5279981: status FINISHED
2026-07-26 14:39:02 INFO None 5279982: status FINISHED
2026-07-26 14:39:04 INFO None 5279983: status FINISHED
2026-07-26 14:39:04 INFO None 5279984: status FINISHED
2026-07-26 14:39:04 INFO None 5279985: status FINISHED
2026-07-26 14:39:04 INFO None 5279989: status FINISHED
2026-07-26 14:39:04 INFO None 5279991: status FINISHED
2026-07-26 14:39:04 INFO None 5279993: status FINISHED
2026-07-26 14:39:04 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:39:04 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:39:04 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:39:04 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:39:04 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:39:04 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:39:04 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:39:19 INFO None 5279980: status FINISHED
2026-07-26 14:39:20 INFO None 5279981: status FINISHED
2026-07-26 14:39:20 INFO None 5279982: status FINISHED
2026-07-26 14:39:20 INFO None 5279983: status FINISHED
2026-07-26 14:39:20 INFO None 5279984: status FINISHED
2026-07-26 14:39:20 INFO None 5279985: status FINISHED
2026-07-26 14:39:20 INFO None 5279989: status FINISHED
2026-07-26 14:39:20 INFO None 5279991: status FINISHED
2026-07-26 14:39:20 INFO None 5279993: status FINISHED
2026-07-26 14:39:20 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:39:20 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:39:20 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:39:20 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:39:20 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:39:20 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:39:20 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:39:35 INFO None 5279980: status FINISHED
2026-07-26 14:39:35 INFO None 5279981: status FINISHED
2026-07-26 14:39:35 INFO None 5279982: status FINISHED
2026-07-26 14:39:35 INFO None 5279983: status FINISHED
2026-07-26 14:39:35 INFO None 5279984: status FINISHED
2026-07-26 14:39:35 INFO None 5279985: status FINISHED
2026-07-26 14:39:35 INFO None 5279989: status FINISHED
2026-07-26 14:39:35 INFO None 5279991: status FINISHED
2026-07-26 14:39:35 INFO None 5279993: status FINISHED
2026-07-26 14:39:35 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:39:35 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:39:35 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:39:35 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:39:35 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:39:35 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:39:35 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:39:50 INFO None 5279980: status FINISHED
2026-07-26 14:39:50 INFO None 5279981: status FINISHED
2026-07-26 14:39:50 INFO None 5279982: status FINISHED
2026-07-26 14:39:50 INFO None 5279983: status FINISHED
2026-07-26 14:39:50 INFO None 5279984: status FINISHED
2026-07-26 14:39:50 INFO None 5279985: status FINISHED
2026-07-26 14:39:50 INFO None 5279989: status FINISHED
2026-07-26 14:39:50 INFO None 5279991: status FINISHED
2026-07-26 14:39:50 INFO None 5279993: status FINISHED
2026-07-26 14:39:50 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:39:50 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:39:50 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:39:50 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:39:50 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:39:50 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:39:50 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:40:06 INFO None 5279980: status FINISHED
2026-07-26 14:40:06 INFO None 5279981: status FINISHED
2026-07-26 14:40:06 INFO None 5279982: status FINISHED
2026-07-26 14:40:06 INFO None 5279983: status FINISHED
2026-07-26 14:40:06 INFO None 5279984: status FINISHED
2026-07-26 14:40:06 INFO None 5279985: status FINISHED
2026-07-26 14:40:06 INFO None 5279989: status FINISHED
2026-07-26 14:40:06 INFO None 5279991: status FINISHED
2026-07-26 14:40:06 INFO None 5279993: status FINISHED
2026-07-26 14:40:06 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:40:06 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:40:06 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:40:06 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:40:06 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:40:06 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:40:06 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:40:21 INFO None 5279980: status FINISHED
2026-07-26 14:40:21 INFO None 5279981: status FINISHED
2026-07-26 14:40:21 INFO None 5279982: status FINISHED
2026-07-26 14:40:21 INFO None 5279983: status FINISHED
2026-07-26 14:40:21 INFO None 5279984: status FINISHED
2026-07-26 14:40:21 INFO None 5279985: status FINISHED
2026-07-26 14:40:21 INFO None 5279989: status FINISHED
2026-07-26 14:40:21 INFO None 5279991: status FINISHED
2026-07-26 14:40:21 INFO None 5279993: status FINISHED
2026-07-26 14:40:21 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:40:21 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:40:23 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:40:23 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:40:23 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:40:23 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:40:23 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:40:38 INFO None 5279980: status FINISHED
2026-07-26 14:40:38 INFO None 5279981: status FINISHED
2026-07-26 14:40:38 INFO None 5279982: status FINISHED
2026-07-26 14:40:38 INFO None 5279983: status FINISHED
2026-07-26 14:40:38 INFO None 5279984: status FINISHED
2026-07-26 14:40:38 INFO None 5279985: status FINISHED
2026-07-26 14:40:38 INFO None 5279989: status FINISHED
2026-07-26 14:40:38 INFO None 5279991: status FINISHED
2026-07-26 14:40:38 INFO None 5279993: status FINISHED
2026-07-26 14:40:39 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:40:39 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:40:39 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:40:39 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:40:39 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:40:39 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:40:39 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:40:54 INFO None 5279980: status FINISHED
2026-07-26 14:40:54 INFO None 5279981: status FINISHED
2026-07-26 14:40:54 INFO None 5279982: status FINISHED
2026-07-26 14:40:54 INFO None 5279983: status FINISHED
2026-07-26 14:40:54 INFO None 5279984: status FINISHED
2026-07-26 14:40:54 INFO None 5279985: status FINISHED
2026-07-26 14:40:56 INFO None 5279989: status FINISHED
2026-07-26 14:40:56 INFO None 5279991: status FINISHED
2026-07-26 14:40:56 INFO None 5279993: status FINISHED
2026-07-26 14:40:56 INFO None 5279994: status RUNNING/PENDING
2026-07-26 14:40:56 INFO None 5279995: status RUNNING/PENDING
2026-07-26 14:40:56 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:40:56 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:40:56 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:40:56 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:40:56 INFO Jobs still running: ['5279994', '5279995', '5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:41:11 INFO None 5279980: status FINISHED
2026-07-26 14:41:11 INFO None 5279981: status FINISHED
2026-07-26 14:41:11 INFO None 5279982: status FINISHED
2026-07-26 14:41:11 INFO None 5279983: status FINISHED
2026-07-26 14:41:11 INFO None 5279984: status FINISHED
2026-07-26 14:41:11 INFO None 5279985: status FINISHED
2026-07-26 14:41:11 INFO None 5279989: status FINISHED
2026-07-26 14:41:11 INFO None 5279991: status FINISHED
2026-07-26 14:41:11 INFO None 5279993: status FINISHED
2026-07-26 14:41:11 INFO None 5279994: status FINISHED
2026-07-26 14:41:11 INFO None 5279995: status FINISHED
2026-07-26 14:41:11 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:41:11 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:41:11 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:41:11 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:41:11 INFO Jobs still running: ['5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:41:26 INFO None 5279980: status FINISHED
2026-07-26 14:41:26 INFO None 5279981: status FINISHED
2026-07-26 14:41:26 INFO None 5279982: status FINISHED
2026-07-26 14:41:26 INFO None 5279983: status FINISHED
2026-07-26 14:41:26 INFO None 5279984: status FINISHED
2026-07-26 14:41:26 INFO None 5279985: status FINISHED
2026-07-26 14:41:26 INFO None 5279989: status FINISHED
2026-07-26 14:41:26 INFO None 5279991: status FINISHED
2026-07-26 14:41:26 INFO None 5279993: status FINISHED
2026-07-26 14:41:26 INFO None 5279994: status FINISHED
2026-07-26 14:41:28 INFO None 5279995: status FINISHED
2026-07-26 14:41:28 INFO None 5279996: status RUNNING/PENDING
2026-07-26 14:41:28 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:41:28 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:41:28 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:41:28 INFO Jobs still running: ['5279996', '5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:41:43 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:41:43 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:41:43 INFO None 5279982: status FINISHED
2026-07-26 14:41:43 INFO None 5279983: status FINISHED
2026-07-26 14:41:43 INFO None 5279984: status FINISHED
2026-07-26 14:41:43 INFO None 5279985: status FINISHED
2026-07-26 14:41:43 INFO None 5279989: status FINISHED
2026-07-26 14:41:43 INFO None 5279991: status FINISHED
2026-07-26 14:41:43 INFO None 5279993: status FINISHED
2026-07-26 14:41:43 INFO None 5279994: status FINISHED
2026-07-26 14:41:43 INFO None 5279995: status FINISHED
2026-07-26 14:41:43 INFO None 5279996: status FINISHED
2026-07-26 14:41:43 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:41:43 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:41:43 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:41:43 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:41:58 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:41:58 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:41:58 INFO None 5279982: status FINISHED
2026-07-26 14:41:58 INFO None 5279983: status FINISHED
2026-07-26 14:41:58 INFO None 5279984: status FINISHED
2026-07-26 14:41:58 INFO None 5279985: status FINISHED
2026-07-26 14:41:58 INFO None 5279989: status FINISHED
2026-07-26 14:41:58 INFO None 5279991: status FINISHED
2026-07-26 14:41:58 INFO None 5279993: status FINISHED
2026-07-26 14:41:58 INFO None 5279994: status FINISHED
2026-07-26 14:41:58 INFO None 5279995: status FINISHED
2026-07-26 14:41:58 INFO None 5279996: status FINISHED
2026-07-26 14:41:58 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:41:58 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:41:58 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:41:58 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:42:15 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:42:15 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:42:15 INFO None 5279982: status FINISHED
2026-07-26 14:42:15 INFO None 5279983: status FINISHED
2026-07-26 14:42:15 INFO None 5279984: status FINISHED
2026-07-26 14:42:15 INFO None 5279985: status FINISHED
2026-07-26 14:42:15 INFO None 5279989: status FINISHED
2026-07-26 14:42:15 INFO None 5279991: status FINISHED
2026-07-26 14:42:15 INFO None 5279993: status FINISHED
2026-07-26 14:42:15 INFO None 5279994: status FINISHED
2026-07-26 14:42:15 INFO None 5279995: status FINISHED
2026-07-26 14:42:15 INFO None 5279996: status FINISHED
2026-07-26 14:42:15 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:42:15 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:42:15 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:42:15 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:42:30 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:42:30 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:42:30 INFO None 5279982: status FINISHED
2026-07-26 14:42:30 INFO None 5279983: status FINISHED
2026-07-26 14:42:30 INFO None 5279984: status FINISHED
2026-07-26 14:42:30 INFO None 5279985: status FINISHED
2026-07-26 14:42:30 INFO None 5279989: status FINISHED
2026-07-26 14:42:30 INFO None 5279991: status FINISHED
2026-07-26 14:42:30 INFO None 5279993: status FINISHED
2026-07-26 14:42:30 INFO None 5279994: status FINISHED
2026-07-26 14:42:30 INFO None 5279995: status FINISHED
2026-07-26 14:42:30 INFO None 5279996: status FINISHED
2026-07-26 14:42:30 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:42:32 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:42:32 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:42:32 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:42:47 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:42:47 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:42:47 INFO None 5279982: status FINISHED
2026-07-26 14:42:47 INFO None 5279983: status FINISHED
2026-07-26 14:42:48 INFO None 5279984: status FINISHED
2026-07-26 14:42:48 INFO None 5279985: status FINISHED
2026-07-26 14:42:48 INFO None 5279989: status FINISHED
2026-07-26 14:42:48 INFO None 5279991: status FINISHED
2026-07-26 14:42:48 INFO None 5279993: status FINISHED
2026-07-26 14:42:48 INFO None 5279994: status FINISHED
2026-07-26 14:42:48 INFO None 5279995: status FINISHED
2026-07-26 14:42:48 INFO None 5279996: status FINISHED
2026-07-26 14:42:48 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:42:48 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:42:48 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:42:48 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:43:03 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:43:03 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:43:03 INFO None 5279982: status FINISHED
2026-07-26 14:43:03 INFO None 5279983: status FINISHED
2026-07-26 14:43:03 INFO None 5279984: status FINISHED
2026-07-26 14:43:03 INFO None 5279985: status FINISHED
2026-07-26 14:43:03 INFO None 5279989: status FINISHED
2026-07-26 14:43:03 INFO None 5279991: status FINISHED
2026-07-26 14:43:03 INFO None 5279993: status FINISHED
2026-07-26 14:43:03 INFO None 5279994: status FINISHED
2026-07-26 14:43:03 INFO None 5279995: status FINISHED
2026-07-26 14:43:03 INFO None 5279996: status FINISHED
2026-07-26 14:43:03 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:43:03 INFO None 5279998: status RUNNING/PENDING
2026-07-26 14:43:03 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:43:03 INFO Jobs still running: ['5279997', '5279998', '5280000']. Waiting...
2026-07-26 14:43:18 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:43:18 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:43:18 INFO None 5279982: status FINISHED
2026-07-26 14:43:18 INFO None 5279983: status FINISHED
2026-07-26 14:43:18 INFO None 5279984: status FINISHED
2026-07-26 14:43:18 INFO None 5279985: status FINISHED
2026-07-26 14:43:20 INFO None 5279989: status FINISHED
2026-07-26 14:43:20 INFO None 5279991: status FINISHED
2026-07-26 14:43:20 INFO None 5279993: status FINISHED
2026-07-26 14:43:20 INFO None 5279994: status FINISHED
2026-07-26 14:43:20 INFO None 5279995: status FINISHED
2026-07-26 14:43:20 INFO None 5279996: status FINISHED
2026-07-26 14:43:20 INFO None 5279997: status RUNNING/PENDING
2026-07-26 14:43:20 INFO None 5279998: status FINISHED
2026-07-26 14:43:20 INFO None 5280000: status RUNNING/PENDING
2026-07-26 14:43:20 INFO Jobs still running: ['5279997', '5280000']. Waiting...
2026-07-26 14:43:35 INFO None 5279980: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279981: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279982: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279983: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279984: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279985: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279989: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279991: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279993: status FINISHED (not in squeue)
2026-07-26 14:43:35 INFO None 5279994: status FINISHED
2026-07-26 14:43:35 INFO None 5279995: status FINISHED
2026-07-26 14:43:35 INFO None 5279996: status FINISHED
2026-07-26 14:43:35 INFO None 5279997: status FINISHED
2026-07-26 14:43:35 INFO None 5279998: status FINISHED
2026-07-26 14:43:36 INFO None 5280000: status FINISHED
2026-07-26 14:43:36 INFO Jobs ['5279980', '5279981', '5279982', '5279983', '5279984', '5279985', '5279989', '5279991', '5279993', '5279994', '5279995', '5279996', '5279997', '5279998', '5280000'] have finished
2026-07-26 14:43:36 INFO Checking restart files were created ...
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020714_10_ENS1.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020714_10_ENS2.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020714_10_ENS3.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020714_10_ENS4.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020714_10_ENS5.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020714_10_ENS6.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020714_10_ENS7.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020714_10_ENS8.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020714_10_ENS9.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020714_10_ENS10.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020714_10_ENS11.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020714_10_ENS12.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020714_10_ENS13.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020714_10_ENS14.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020714_10_ENS15.nc(3673513755 bytes)
2026-07-26 14:43:36 INFO  Run_model() completed successfully.
2026-07-26 14:43:36 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 14:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:43:36 INFO [TIME] gregorian_conversion simulated_time=2020-02-08 00:00:00 days=153074 seconds=0
2026-07-26 14:43:36 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 14:43:36 INFO [TIME] increment current_time 2020-02-07 14:00:00 -> 2020-02-08 00:00:00
2026-07-26 14:43:36 INFO [TIME] after_increment_before_assimilation current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:43:36 INFO ---------->>> Running process_satellite_data()
2026-07-26 14:43:36 INFO [DART] No satellite data found, skipping assimilation
2026-07-26 14:43:36 INFO after_assimilation() skipped
2026-07-26 14:43:36 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 14:43:36 INFO [TIME] step_end current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 14:43:36 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
