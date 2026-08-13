+ SCRIPT_PID=2564498
+ /bin/bash -x /tmp/tmp.ZhqIrBtrBe
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
2026-07-26 12:14:25 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-26 12:14:25 INFO [PIPELINE] =======================================
2026-07-26 12:14:25 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-26 12:14:25 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-26 12:14:25 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-26 12:14:25 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260726_121425.log
2026-07-26 12:14:25 INFO [PIPELINE] =======================================
2026-07-26 12:14:25 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-26 12:14:25 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-26 12:14:25 INFO [STEP] ---- TIME LOOP START ----
2026-07-26 12:14:25 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:14:25 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-26 12:14:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:33 INFO Hourly dataset computed and listing created
2026-07-26 12:14:38 INFO Hourly dataset computed
2026-07-26 12:14:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:39 INFO Hourly dataset computed and listing created
2026-07-26 12:14:40 INFO Hourly dataset computed
2026-07-26 12:14:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:41 INFO Hourly dataset computed and listing created
2026-07-26 12:14:42 INFO Hourly dataset computed
2026-07-26 12:14:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:43 INFO Hourly dataset computed and listing created
2026-07-26 12:14:43 INFO Hourly dataset computed
2026-07-26 12:14:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:44 INFO Hourly dataset computed and listing created
2026-07-26 12:14:45 INFO Hourly dataset computed
2026-07-26 12:14:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:46 INFO Hourly dataset computed and listing created
2026-07-26 12:14:47 INFO Hourly dataset computed
2026-07-26 12:14:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:48 INFO Hourly dataset computed and listing created
2026-07-26 12:14:49 INFO Hourly dataset computed
2026-07-26 12:14:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:50 INFO Hourly dataset computed and listing created
2026-07-26 12:14:51 INFO Hourly dataset computed
2026-07-26 12:14:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:52 INFO Hourly dataset computed and listing created
2026-07-26 12:14:53 INFO Hourly dataset computed
2026-07-26 12:14:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:54 INFO Hourly dataset computed and listing created
2026-07-26 12:14:55 INFO Hourly dataset computed
2026-07-26 12:14:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:56 INFO Hourly dataset computed and listing created
2026-07-26 12:14:57 INFO Hourly dataset computed
2026-07-26 12:14:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:14:58 INFO Hourly dataset computed and listing created
2026-07-26 12:14:59 INFO Hourly dataset computed
2026-07-26 12:14:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:15:00 INFO Hourly dataset computed and listing created
2026-07-26 12:15:01 INFO Hourly dataset computed
2026-07-26 12:15:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:15:02 INFO Hourly dataset computed and listing created
2026-07-26 12:15:03 INFO Hourly dataset computed
2026-07-26 12:15:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:15:04 INFO Hourly dataset computed and listing created
2026-07-26 12:15:05 INFO Hourly dataset computed
2026-07-26 12:15:05 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-26 12:15:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:15:05 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-26 12:15:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:15:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:05 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:15:05 INFO Queuing job for member 1...
2026-07-26 12:15:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:05 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:15:06 INFO Found: ['5279531']
2026-07-26 12:15:11 INFO [TGCC-IRENE] Submitted job with ID:['5279531']
2026-07-26 12:15:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:15:11 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-26 12:15:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:15:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:11 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:15:11 INFO Queuing job for member 2...
2026-07-26 12:15:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:11 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:15:11 INFO Found: ['5279532']
2026-07-26 12:15:16 INFO [TGCC-IRENE] Submitted job with ID:['5279532']
2026-07-26 12:15:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:15:16 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-26 12:15:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:15:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:16 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:15:16 INFO Queuing job for member 3...
2026-07-26 12:15:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:16 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:15:17 INFO Found: ['5279533']
2026-07-26 12:15:22 INFO [TGCC-IRENE] Submitted job with ID:['5279533']
2026-07-26 12:15:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:15:22 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-26 12:15:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:15:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:22 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:15:22 INFO Queuing job for member 4...
2026-07-26 12:15:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:22 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:15:23 INFO Found: ['5279534']
2026-07-26 12:15:28 INFO [TGCC-IRENE] Submitted job with ID:['5279534']
2026-07-26 12:15:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:15:28 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-26 12:15:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:15:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:28 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:15:28 INFO Queuing job for member 5...
2026-07-26 12:15:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:28 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:15:29 INFO Found: ['5279535']
2026-07-26 12:15:34 INFO [TGCC-IRENE] Submitted job with ID:['5279535']
2026-07-26 12:15:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:15:34 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-26 12:15:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:15:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:34 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:15:34 INFO Queuing job for member 6...
2026-07-26 12:15:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:34 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:15:34 INFO Found: ['5279536']
2026-07-26 12:15:39 INFO [TGCC-IRENE] Submitted job with ID:['5279536']
2026-07-26 12:15:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:15:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:15:39 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-26 12:15:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:15:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:15:39 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:15:39 INFO Queuing job for member 7...
2026-07-26 12:15:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:15:39 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:15:42 INFO Found: ['5279537']
[2026-07-26T12:15:44.787] error: *** JOB 5279524 ON irene4333 CANCELLED AT 2026-07-26T12:15:44 DUE to SIGNAL Terminated ***
