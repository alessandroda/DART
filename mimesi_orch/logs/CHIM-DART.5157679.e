+ /bin/bash -x /tmp/tmp.SGEC2cW1Ul
+ SCRIPT_PID=1592493
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
2026-07-14 20:24:06 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-14 20:24:06 INFO [PIPELINE] =======================================
2026-07-14 20:24:06 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-14 20:24:06 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-14 20:24:06 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-14 20:24:06 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260714_202406.log
2026-07-14 20:24:06 INFO [PIPELINE] =======================================
2026-07-14 20:24:06 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-14 20:24:06 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-14 20:24:06 INFO [STEP] ---- TIME LOOP START ----
2026-07-14 20:24:06 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:24:06 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-14 20:24:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:15 INFO Hourly dataset computed and listing created
2026-07-14 20:24:19 INFO Hourly dataset computed
2026-07-14 20:24:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:20 INFO Hourly dataset computed and listing created
2026-07-14 20:24:21 INFO Hourly dataset computed
2026-07-14 20:24:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:22 INFO Hourly dataset computed and listing created
2026-07-14 20:24:23 INFO Hourly dataset computed
2026-07-14 20:24:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:24 INFO Hourly dataset computed and listing created
2026-07-14 20:24:25 INFO Hourly dataset computed
2026-07-14 20:24:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:26 INFO Hourly dataset computed and listing created
2026-07-14 20:24:27 INFO Hourly dataset computed
2026-07-14 20:24:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:28 INFO Hourly dataset computed and listing created
2026-07-14 20:24:28 INFO Hourly dataset computed
2026-07-14 20:24:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:29 INFO Hourly dataset computed and listing created
2026-07-14 20:24:30 INFO Hourly dataset computed
2026-07-14 20:24:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:31 INFO Hourly dataset computed and listing created
2026-07-14 20:24:32 INFO Hourly dataset computed
2026-07-14 20:24:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:33 INFO Hourly dataset computed and listing created
2026-07-14 20:24:34 INFO Hourly dataset computed
2026-07-14 20:24:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:35 INFO Hourly dataset computed and listing created
2026-07-14 20:24:35 INFO Hourly dataset computed
2026-07-14 20:24:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:36 INFO Hourly dataset computed and listing created
2026-07-14 20:24:37 INFO Hourly dataset computed
2026-07-14 20:24:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:38 INFO Hourly dataset computed and listing created
2026-07-14 20:24:39 INFO Hourly dataset computed
2026-07-14 20:24:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:40 INFO Hourly dataset computed and listing created
2026-07-14 20:24:41 INFO Hourly dataset computed
2026-07-14 20:24:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:42 INFO Hourly dataset computed and listing created
2026-07-14 20:24:42 INFO Hourly dataset computed
2026-07-14 20:24:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:24:44 INFO Hourly dataset computed and listing created
2026-07-14 20:24:48 INFO Hourly dataset computed
2026-07-14 20:24:48 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-14 20:24:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:24:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 20:24:48 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-14 20:24:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 20:24:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:24:48 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 20:24:48 INFO Queuing job for member 1...
2026-07-14 20:24:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:24:48 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 20:24:49 INFO Found: ['5157687']
2026-07-14 20:24:54 INFO [TGCC-IRENE] Submitted job with ID:['5157687']
2026-07-14 20:24:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:24:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 20:24:54 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-14 20:24:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 20:24:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:24:54 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 20:24:54 INFO Queuing job for member 2...
2026-07-14 20:24:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:24:54 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 20:24:55 INFO Found: ['5157689']
2026-07-14 20:25:00 INFO [TGCC-IRENE] Submitted job with ID:['5157689']
2026-07-14 20:25:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 20:25:00 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-14 20:25:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 20:25:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:00 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 20:25:00 INFO Queuing job for member 3...
2026-07-14 20:25:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:00 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 20:25:02 INFO Found: ['5157690']
2026-07-14 20:25:07 INFO [TGCC-IRENE] Submitted job with ID:['5157690']
2026-07-14 20:25:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 20:25:07 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-14 20:25:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 20:25:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:07 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 20:25:07 INFO Queuing job for member 4...
2026-07-14 20:25:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:07 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 20:25:09 INFO Found: ['5157692']
2026-07-14 20:25:14 INFO [TGCC-IRENE] Submitted job with ID:['5157692']
2026-07-14 20:25:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 20:25:14 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-14 20:25:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 20:25:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:15 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 20:25:15 INFO Queuing job for member 5...
2026-07-14 20:25:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:15 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 20:25:15 INFO Found: ['5157694']
2026-07-14 20:25:20 INFO [TGCC-IRENE] Submitted job with ID:['5157694']
2026-07-14 20:25:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 20:25:20 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-14 20:25:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 20:25:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:20 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 20:25:20 INFO Queuing job for member 6...
2026-07-14 20:25:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:20 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 20:25:21 INFO Found: ['5157695']
2026-07-14 20:25:26 INFO [TGCC-IRENE] Submitted job with ID:['5157695']
2026-07-14 20:25:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 20:25:26 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-14 20:25:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 20:25:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:26 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 20:25:26 INFO Queuing job for member 7...
2026-07-14 20:25:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:26 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 20:25:27 INFO Found: ['5157696']
2026-07-14 20:25:32 INFO [TGCC-IRENE] Submitted job with ID:['5157696']
2026-07-14 20:25:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 20:25:32 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-14 20:25:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 20:25:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:32 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 20:25:32 INFO Queuing job for member 8...
2026-07-14 20:25:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:32 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 20:25:33 INFO Found: ['5157697']
2026-07-14 20:25:38 INFO [TGCC-IRENE] Submitted job with ID:['5157697']
2026-07-14 20:25:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 20:25:38 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-14 20:25:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 20:25:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:38 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 20:25:38 INFO Queuing job for member 9...
2026-07-14 20:25:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:38 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 20:25:38 INFO Found: ['5157698']
2026-07-14 20:25:43 INFO [TGCC-IRENE] Submitted job with ID:['5157698']
2026-07-14 20:25:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 20:25:43 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-14 20:25:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 20:25:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:43 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 20:25:44 INFO Queuing job for member 10...
2026-07-14 20:25:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:44 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 20:25:44 INFO Found: ['5157699']
2026-07-14 20:25:49 INFO [TGCC-IRENE] Submitted job with ID:['5157699']
2026-07-14 20:25:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 20:25:49 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-14 20:25:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 20:25:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:49 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 20:25:49 INFO Queuing job for member 11...
2026-07-14 20:25:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:49 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 20:25:50 INFO Found: ['5157700']
2026-07-14 20:25:55 INFO [TGCC-IRENE] Submitted job with ID:['5157700']
2026-07-14 20:25:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:25:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 20:25:55 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-14 20:25:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 20:25:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:25:55 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 20:25:55 INFO Queuing job for member 12...
2026-07-14 20:25:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:25:55 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 20:25:57 INFO Found: ['5157701']
2026-07-14 20:26:02 INFO [TGCC-IRENE] Submitted job with ID:['5157701']
2026-07-14 20:28:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:28:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 20:28:08 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-14 20:28:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 20:28:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:28:08 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 20:28:08 INFO Queuing job for member 13...
2026-07-14 20:28:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:28:08 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 20:28:09 INFO Found: ['5157710']
2026-07-14 20:28:14 INFO [TGCC-IRENE] Submitted job with ID:['5157710']
2026-07-14 20:28:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:28:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 20:28:14 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-14 20:28:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 20:28:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:28:14 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 20:28:14 INFO Queuing job for member 14...
2026-07-14 20:28:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:28:14 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 20:28:15 INFO Found: ['5157711']
2026-07-14 20:28:20 INFO [TGCC-IRENE] Submitted job with ID:['5157711']
2026-07-14 20:28:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:28:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 20:28:20 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-14 20:28:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 20:28:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:28:20 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 20:28:20 INFO Queuing job for member 15...
2026-07-14 20:28:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:28:20 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 20:28:20 INFO Found: ['5157712']
2026-07-14 20:28:25 INFO [TGCC-IRENE] Submitted job with ID:['5157712']
2026-07-14 20:28:25 INFO Checking job status ...
2026-07-14 20:28:25 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:28:25 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:28:26 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:28:26 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:28:26 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:28:26 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:28:26 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:28:41 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:28:41 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:28:41 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:28:57 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:28:57 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:28:57 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:29:12 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:29:12 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:29:12 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:29:12 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:29:12 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:29:13 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:29:13 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:29:28 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:29:28 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:29:28 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:29:43 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:29:43 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:29:43 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:29:58 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:29:58 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:29:58 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:29:58 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:29:58 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:29:59 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:29:59 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:30:14 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:30:14 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:30:14 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:30:29 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157689: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:30:29 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:30:29 INFO Jobs still running: ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:30:46 INFO None 5157687: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157689: status FINISHED
2026-07-14 20:30:46 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:30:46 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:30:46 INFO Jobs still running: ['5157687', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:31:01 INFO None 5157687: status FINISHED
2026-07-14 20:31:01 INFO None 5157689: status FINISHED
2026-07-14 20:31:01 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157695: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157700: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:31:01 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:31:01 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:31:16 INFO None 5157687: status FINISHED
2026-07-14 20:33:14 INFO None 5157689: status FINISHED
2026-07-14 20:33:14 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157695: status FINISHED
2026-07-14 20:33:14 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157700: status FINISHED
2026-07-14 20:33:14 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:33:14 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:33:15 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:33:15 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:33:15 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157696', '5157697', '5157698', '5157699', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:33:30 INFO None 5157687: status FINISHED
2026-07-14 20:33:30 INFO None 5157689: status FINISHED
2026-07-14 20:33:30 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157695: status FINISHED
2026-07-14 20:33:30 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157700: status FINISHED
2026-07-14 20:33:30 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157710: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:33:30 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:33:30 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157696', '5157697', '5157698', '5157699', '5157701', '5157710', '5157711', '5157712']. Waiting...
2026-07-14 20:33:45 INFO None 5157687: status FINISHED
2026-07-14 20:33:45 INFO None 5157689: status FINISHED
2026-07-14 20:33:45 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157695: status FINISHED
2026-07-14 20:33:45 INFO None 5157696: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157697: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157698: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157699: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157700: status FINISHED
2026-07-14 20:33:45 INFO None 5157701: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157710: status FINISHED
2026-07-14 20:33:45 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:33:45 INFO None 5157712: status RUNNING/PENDING
2026-07-14 20:33:45 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157696', '5157697', '5157698', '5157699', '5157701', '5157711', '5157712']. Waiting...
2026-07-14 20:34:00 INFO None 5157687: status FINISHED
2026-07-14 20:34:00 INFO None 5157689: status FINISHED
2026-07-14 20:34:00 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:34:00 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:34:00 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:34:00 INFO None 5157695: status FINISHED
2026-07-14 20:34:00 INFO None 5157696: status FINISHED
2026-07-14 20:34:00 INFO None 5157697: status FINISHED
2026-07-14 20:34:00 INFO None 5157698: status FINISHED
2026-07-14 20:34:00 INFO None 5157699: status FINISHED
2026-07-14 20:34:00 INFO None 5157700: status FINISHED
2026-07-14 20:34:00 INFO None 5157701: status FINISHED
2026-07-14 20:34:00 INFO None 5157710: status FINISHED
2026-07-14 20:34:01 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:34:01 INFO None 5157712: status FINISHED
2026-07-14 20:34:01 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157711']. Waiting...
2026-07-14 20:34:16 INFO None 5157687: status FINISHED
2026-07-14 20:34:16 INFO None 5157689: status FINISHED
2026-07-14 20:34:16 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:34:16 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:34:16 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:34:16 INFO None 5157695: status FINISHED
2026-07-14 20:34:16 INFO None 5157696: status FINISHED
2026-07-14 20:34:16 INFO None 5157697: status FINISHED
2026-07-14 20:34:16 INFO None 5157698: status FINISHED
2026-07-14 20:34:16 INFO None 5157699: status FINISHED
2026-07-14 20:34:16 INFO None 5157700: status FINISHED
2026-07-14 20:34:16 INFO None 5157701: status FINISHED
2026-07-14 20:34:16 INFO None 5157710: status FINISHED
2026-07-14 20:34:16 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:34:16 INFO None 5157712: status FINISHED
2026-07-14 20:34:16 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157711']. Waiting...
2026-07-14 20:34:31 INFO None 5157687: status FINISHED
2026-07-14 20:34:31 INFO None 5157689: status FINISHED
2026-07-14 20:34:31 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:34:31 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:34:31 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:34:31 INFO None 5157695: status FINISHED
2026-07-14 20:34:31 INFO None 5157696: status FINISHED
2026-07-14 20:34:31 INFO None 5157697: status FINISHED
2026-07-14 20:34:31 INFO None 5157698: status FINISHED
2026-07-14 20:34:31 INFO None 5157699: status FINISHED
2026-07-14 20:34:32 INFO None 5157700: status FINISHED
2026-07-14 20:34:32 INFO None 5157701: status FINISHED
2026-07-14 20:34:32 INFO None 5157710: status FINISHED
2026-07-14 20:34:32 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:34:32 INFO None 5157712: status FINISHED
2026-07-14 20:34:32 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157711']. Waiting...
2026-07-14 20:34:47 INFO None 5157687: status FINISHED
2026-07-14 20:34:47 INFO None 5157689: status FINISHED
2026-07-14 20:34:47 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:34:47 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:34:47 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:34:47 INFO None 5157695: status FINISHED
2026-07-14 20:34:47 INFO None 5157696: status FINISHED
2026-07-14 20:34:47 INFO None 5157697: status FINISHED
2026-07-14 20:34:47 INFO None 5157698: status FINISHED
2026-07-14 20:34:47 INFO None 5157699: status FINISHED
2026-07-14 20:34:47 INFO None 5157700: status FINISHED
2026-07-14 20:34:47 INFO None 5157701: status FINISHED
2026-07-14 20:34:47 INFO None 5157710: status FINISHED
2026-07-14 20:34:47 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:34:47 INFO None 5157712: status FINISHED
2026-07-14 20:34:47 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157711']. Waiting...
2026-07-14 20:35:02 INFO None 5157687: status FINISHED
2026-07-14 20:35:02 INFO None 5157689: status FINISHED
2026-07-14 20:35:02 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:35:02 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:35:02 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:35:02 INFO None 5157695: status FINISHED
2026-07-14 20:35:02 INFO None 5157696: status FINISHED
2026-07-14 20:35:02 INFO None 5157697: status FINISHED
2026-07-14 20:35:02 INFO None 5157698: status FINISHED
2026-07-14 20:35:02 INFO None 5157699: status FINISHED
2026-07-14 20:35:02 INFO None 5157700: status FINISHED
2026-07-14 20:35:02 INFO None 5157701: status FINISHED
2026-07-14 20:35:02 INFO None 5157710: status FINISHED
2026-07-14 20:35:02 INFO None 5157711: status RUNNING/PENDING
2026-07-14 20:35:02 INFO None 5157712: status FINISHED
2026-07-14 20:35:02 INFO Jobs still running: ['5157690', '5157692', '5157694', '5157711']. Waiting...
2026-07-14 20:35:17 INFO None 5157687: status FINISHED
2026-07-14 20:35:17 INFO None 5157689: status FINISHED
2026-07-14 20:35:17 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:35:17 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:35:17 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:35:17 INFO None 5157695: status FINISHED
2026-07-14 20:35:17 INFO None 5157696: status FINISHED
2026-07-14 20:35:17 INFO None 5157697: status FINISHED
2026-07-14 20:35:17 INFO None 5157698: status FINISHED
2026-07-14 20:35:17 INFO None 5157699: status FINISHED
2026-07-14 20:35:17 INFO None 5157700: status FINISHED
2026-07-14 20:35:17 INFO None 5157701: status FINISHED
2026-07-14 20:35:17 INFO None 5157710: status FINISHED
2026-07-14 20:35:17 INFO None 5157711: status FINISHED
2026-07-14 20:35:17 INFO None 5157712: status FINISHED
2026-07-14 20:35:17 INFO Jobs still running: ['5157690', '5157692', '5157694']. Waiting...
2026-07-14 20:35:33 INFO None 5157687: status FINISHED
2026-07-14 20:35:33 INFO None 5157689: status FINISHED
2026-07-14 20:35:33 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:35:33 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:35:33 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:35:33 INFO None 5157695: status FINISHED
2026-07-14 20:35:33 INFO None 5157696: status FINISHED
2026-07-14 20:35:33 INFO None 5157697: status FINISHED
2026-07-14 20:35:33 INFO None 5157698: status FINISHED
2026-07-14 20:35:33 INFO None 5157699: status FINISHED
2026-07-14 20:35:33 INFO None 5157700: status FINISHED
2026-07-14 20:35:33 INFO None 5157701: status FINISHED
2026-07-14 20:35:33 INFO None 5157710: status FINISHED
2026-07-14 20:35:33 INFO None 5157711: status FINISHED
2026-07-14 20:35:33 INFO None 5157712: status FINISHED
2026-07-14 20:35:33 INFO Jobs still running: ['5157690', '5157692', '5157694']. Waiting...
2026-07-14 20:35:48 INFO None 5157687: status FINISHED
2026-07-14 20:35:48 INFO None 5157689: status FINISHED
2026-07-14 20:35:48 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:35:48 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:35:48 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:35:48 INFO None 5157695: status FINISHED
2026-07-14 20:35:48 INFO None 5157696: status FINISHED
2026-07-14 20:35:48 INFO None 5157697: status FINISHED
2026-07-14 20:35:48 INFO None 5157698: status FINISHED
2026-07-14 20:35:48 INFO None 5157699: status FINISHED
2026-07-14 20:35:48 INFO None 5157700: status FINISHED
2026-07-14 20:35:48 INFO None 5157701: status FINISHED
2026-07-14 20:35:48 INFO None 5157710: status FINISHED
2026-07-14 20:35:48 INFO None 5157711: status FINISHED
2026-07-14 20:35:48 INFO None 5157712: status FINISHED
2026-07-14 20:35:48 INFO Jobs still running: ['5157690', '5157692', '5157694']. Waiting...
2026-07-14 20:36:03 INFO None 5157687: status FINISHED
2026-07-14 20:36:03 INFO None 5157689: status FINISHED
2026-07-14 20:36:03 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:36:03 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:36:03 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:36:03 INFO None 5157695: status FINISHED
2026-07-14 20:36:03 INFO None 5157696: status FINISHED
2026-07-14 20:36:03 INFO None 5157697: status FINISHED
2026-07-14 20:36:03 INFO None 5157698: status FINISHED
2026-07-14 20:36:03 INFO None 5157699: status FINISHED
2026-07-14 20:36:04 INFO None 5157700: status FINISHED
2026-07-14 20:36:04 INFO None 5157701: status FINISHED
2026-07-14 20:36:04 INFO None 5157710: status FINISHED
2026-07-14 20:36:04 INFO None 5157711: status FINISHED
2026-07-14 20:36:04 INFO None 5157712: status FINISHED
2026-07-14 20:36:04 INFO Jobs still running: ['5157690', '5157692', '5157694']. Waiting...
2026-07-14 20:36:19 INFO None 5157687: status FINISHED
2026-07-14 20:36:19 INFO None 5157689: status FINISHED
2026-07-14 20:36:19 INFO None 5157690: status RUNNING/PENDING
2026-07-14 20:36:19 INFO None 5157692: status RUNNING/PENDING
2026-07-14 20:36:19 INFO None 5157694: status RUNNING/PENDING
2026-07-14 20:36:19 INFO None 5157695: status FINISHED
2026-07-14 20:36:19 INFO None 5157696: status FINISHED
2026-07-14 20:36:19 INFO None 5157697: status FINISHED
2026-07-14 20:36:19 INFO None 5157698: status FINISHED
2026-07-14 20:36:19 INFO None 5157699: status FINISHED
2026-07-14 20:36:19 INFO None 5157700: status FINISHED
2026-07-14 20:36:19 INFO None 5157701: status FINISHED
2026-07-14 20:36:19 INFO None 5157710: status FINISHED
2026-07-14 20:36:19 INFO None 5157711: status FINISHED
2026-07-14 20:36:19 INFO None 5157712: status FINISHED
2026-07-14 20:36:19 INFO Jobs still running: ['5157690', '5157692', '5157694']. Waiting...
2026-07-14 20:36:34 INFO None 5157687: status FINISHED
2026-07-14 20:36:34 INFO None 5157689: status FINISHED
2026-07-14 20:36:34 INFO None 5157690: status FINISHED
2026-07-14 20:36:34 INFO None 5157692: status FINISHED
2026-07-14 20:36:34 INFO None 5157694: status FINISHED
2026-07-14 20:36:34 INFO None 5157695: status FINISHED
2026-07-14 20:36:34 INFO None 5157696: status FINISHED
2026-07-14 20:36:34 INFO None 5157697: status FINISHED
2026-07-14 20:36:34 INFO None 5157698: status FINISHED
2026-07-14 20:36:34 INFO None 5157699: status FINISHED
2026-07-14 20:36:34 INFO None 5157700: status FINISHED
2026-07-14 20:36:34 INFO None 5157701: status FINISHED
2026-07-14 20:36:34 INFO None 5157710: status FINISHED
2026-07-14 20:36:34 INFO None 5157711: status FINISHED
2026-07-14 20:36:34 INFO None 5157712: status FINISHED
2026-07-14 20:36:34 INFO Jobs ['5157687', '5157689', '5157690', '5157692', '5157694', '5157695', '5157696', '5157697', '5157698', '5157699', '5157700', '5157701', '5157710', '5157711', '5157712'] have finished
2026-07-14 20:36:34 INFO Checking restart files were created ...
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-14 20:38:16 INFO  Run_model() completed successfully.
2026-07-14 20:38:16 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:38:16 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-14 20:38:16 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 20:38:16 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-14 20:38:16 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:38:16 INFO ---------->>> Running process_satellite_data()
2026-07-14 20:38:16 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-14 20:38:16 INFO ---------->>> Running run_obs_converter()
2026-07-14 20:38:16 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-14 20:38:16 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-14 20:38:16 INFO ---------->>> Running DART
2026-07-14 20:38:16 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 20:38:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-14 20:38:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-14 20:38:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-14 20:38:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-14 20:38:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-14 20:38:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-14 20:38:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-14 20:38:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-14 20:38:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-14 20:38:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-14 20:38:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-14 20:38:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-14 20:38:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-14 20:38:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-14 20:38:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-14 20:38:21 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 20:38:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 20:38:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 20:38:21 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 20:38:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 20:38:21 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 20:38:35 INFO Found: []
2026-07-14 20:38:35 INFO No job id returned by command ./run_filter.bsh
2026-07-14 20:38:35 INFO No monitoring will be performed
2026-07-14 20:38:35 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-14 20:38:35 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-14 20:38:35 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 20:38:36 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 20:38:36 INFO run_dart() is DONE.
2026-07-14 20:38:36 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 20:38:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-14 20:38:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:38:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-14 20:38:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:38:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-14 20:38:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:38:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-14 20:38:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:38:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-14 20:39:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-14 20:39:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-14 20:39:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-14 20:39:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-14 20:39:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-14 20:39:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-14 20:39:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-14 20:39:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-14 20:39:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-14 20:39:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-14 20:39:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:39:56 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 20:39:56 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:39:56 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:39:57 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-14 20:39:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:39:58 INFO Hourly dataset computed and listing created
2026-07-14 20:40:00 INFO Hourly dataset computed
2026-07-14 20:40:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:01 INFO Hourly dataset computed and listing created
2026-07-14 20:40:01 INFO Hourly dataset computed
2026-07-14 20:40:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:02 INFO Hourly dataset computed and listing created
2026-07-14 20:40:03 INFO Hourly dataset computed
2026-07-14 20:40:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:04 INFO Hourly dataset computed and listing created
2026-07-14 20:40:04 INFO Hourly dataset computed
2026-07-14 20:40:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:05 INFO Hourly dataset computed and listing created
2026-07-14 20:40:06 INFO Hourly dataset computed
2026-07-14 20:40:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:07 INFO Hourly dataset computed and listing created
2026-07-14 20:40:07 INFO Hourly dataset computed
2026-07-14 20:40:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:08 INFO Hourly dataset computed and listing created
2026-07-14 20:40:09 INFO Hourly dataset computed
2026-07-14 20:40:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:10 INFO Hourly dataset computed and listing created
2026-07-14 20:40:10 INFO Hourly dataset computed
2026-07-14 20:40:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:11 INFO Hourly dataset computed and listing created
2026-07-14 20:40:12 INFO Hourly dataset computed
2026-07-14 20:40:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:13 INFO Hourly dataset computed and listing created
2026-07-14 20:40:13 INFO Hourly dataset computed
2026-07-14 20:40:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:14 INFO Hourly dataset computed and listing created
2026-07-14 20:40:15 INFO Hourly dataset computed
2026-07-14 20:40:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:16 INFO Hourly dataset computed and listing created
2026-07-14 20:40:16 INFO Hourly dataset computed
2026-07-14 20:40:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:17 INFO Hourly dataset computed and listing created
2026-07-14 20:40:18 INFO Hourly dataset computed
2026-07-14 20:40:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:19 INFO Hourly dataset computed and listing created
2026-07-14 20:40:19 INFO Hourly dataset computed
2026-07-14 20:40:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:40:20 INFO Hourly dataset computed and listing created
2026-07-14 20:40:21 INFO Hourly dataset computed
2026-07-14 20:40:21 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-14 20:40:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 20:40:21 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-14 20:40:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 20:40:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:21 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 20:40:21 INFO Queuing job for member 1...
2026-07-14 20:40:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:21 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 20:40:22 INFO Found: ['5157740']
2026-07-14 20:40:27 INFO [TGCC-IRENE] Submitted job with ID:['5157740']
2026-07-14 20:40:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 20:40:27 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-14 20:40:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 20:40:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:27 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 20:40:27 INFO Queuing job for member 2...
2026-07-14 20:40:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:27 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 20:40:28 INFO Found: ['5157741']
2026-07-14 20:40:33 INFO [TGCC-IRENE] Submitted job with ID:['5157741']
2026-07-14 20:40:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 20:40:33 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-14 20:40:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 20:40:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:33 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 20:40:33 INFO Queuing job for member 3...
2026-07-14 20:40:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:33 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 20:40:33 INFO Found: ['5157742']
2026-07-14 20:40:38 INFO [TGCC-IRENE] Submitted job with ID:['5157742']
2026-07-14 20:40:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 20:40:38 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-14 20:40:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 20:40:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:38 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 20:40:38 INFO Queuing job for member 4...
2026-07-14 20:40:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:38 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 20:40:39 INFO Found: ['5157743']
2026-07-14 20:40:44 INFO [TGCC-IRENE] Submitted job with ID:['5157743']
2026-07-14 20:40:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 20:40:44 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-14 20:40:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 20:40:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:44 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 20:40:44 INFO Queuing job for member 5...
2026-07-14 20:40:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:44 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 20:40:45 INFO Found: ['5157744']
2026-07-14 20:40:50 INFO [TGCC-IRENE] Submitted job with ID:['5157744']
2026-07-14 20:40:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 20:40:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-14 20:40:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 20:40:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 20:40:50 INFO Queuing job for member 6...
2026-07-14 20:40:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 20:40:50 INFO Found: ['5157745']
2026-07-14 20:40:55 INFO [TGCC-IRENE] Submitted job with ID:['5157745']
2026-07-14 20:40:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:40:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 20:40:55 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-14 20:40:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 20:40:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:40:55 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 20:40:55 INFO Queuing job for member 7...
2026-07-14 20:40:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:40:55 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 20:40:56 INFO Found: ['5157746']
2026-07-14 20:41:01 INFO [TGCC-IRENE] Submitted job with ID:['5157746']
2026-07-14 20:41:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:41:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 20:41:01 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-14 20:41:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 20:41:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:41:01 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 20:41:01 INFO Queuing job for member 8...
2026-07-14 20:41:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:41:01 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 20:43:10 INFO Found: ['5157751']
2026-07-14 20:43:15 INFO [TGCC-IRENE] Submitted job with ID:['5157751']
2026-07-14 20:43:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 20:43:15 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-14 20:43:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 20:43:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:15 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 20:43:15 INFO Queuing job for member 9...
2026-07-14 20:43:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:15 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 20:43:16 INFO Found: ['5157752']
2026-07-14 20:43:21 INFO [TGCC-IRENE] Submitted job with ID:['5157752']
2026-07-14 20:43:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 20:43:21 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-14 20:43:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 20:43:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:21 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 20:43:21 INFO Queuing job for member 10...
2026-07-14 20:43:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:21 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 20:43:22 INFO Found: ['5157753']
2026-07-14 20:43:27 INFO [TGCC-IRENE] Submitted job with ID:['5157753']
2026-07-14 20:43:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 20:43:27 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-14 20:43:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 20:43:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:27 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 20:43:27 INFO Queuing job for member 11...
2026-07-14 20:43:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:27 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 20:43:28 INFO Found: ['5157754']
2026-07-14 20:43:33 INFO [TGCC-IRENE] Submitted job with ID:['5157754']
2026-07-14 20:43:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 20:43:33 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-14 20:43:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 20:43:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:33 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 20:43:33 INFO Queuing job for member 12...
2026-07-14 20:43:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:33 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 20:43:33 INFO Found: ['5157755']
2026-07-14 20:43:38 INFO [TGCC-IRENE] Submitted job with ID:['5157755']
2026-07-14 20:43:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 20:43:38 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-14 20:43:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 20:43:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:38 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 20:43:38 INFO Queuing job for member 13...
2026-07-14 20:43:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:38 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 20:43:39 INFO Found: ['5157756']
2026-07-14 20:43:44 INFO [TGCC-IRENE] Submitted job with ID:['5157756']
2026-07-14 20:43:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 20:43:44 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-14 20:43:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 20:43:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:44 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 20:43:44 INFO Queuing job for member 14...
2026-07-14 20:43:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:44 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 20:43:45 INFO Found: ['5157757']
2026-07-14 20:43:50 INFO [TGCC-IRENE] Submitted job with ID:['5157757']
2026-07-14 20:43:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:43:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 20:43:50 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-14 20:43:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 20:43:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:43:50 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 20:43:50 INFO Queuing job for member 15...
2026-07-14 20:43:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:43:50 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 20:43:52 INFO Found: ['5157758']
2026-07-14 20:43:57 INFO [TGCC-IRENE] Submitted job with ID:['5157758']
2026-07-14 20:43:57 INFO Checking job status ...
2026-07-14 20:43:57 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157745: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157746: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:43:57 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:43:57 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:44:12 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157745: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157746: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:44:12 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:44:12 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:44:28 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157745: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157746: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:44:28 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:44:28 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:44:43 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157745: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157746: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:44:43 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:44:43 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:44:58 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157745: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157746: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:44:58 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:44:58 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:45:13 INFO None 5157740: status RUNNING/PENDING
2026-07-14 20:45:13 INFO None 5157741: status RUNNING/PENDING
2026-07-14 20:45:13 INFO None 5157742: status RUNNING/PENDING
2026-07-14 20:45:13 INFO None 5157743: status RUNNING/PENDING
2026-07-14 20:45:13 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157745: status FINISHED
2026-07-14 20:45:14 INFO None 5157746: status FINISHED
2026-07-14 20:45:14 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:45:14 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:45:14 INFO Jobs still running: ['5157740', '5157741', '5157742', '5157743', '5157744', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:45:29 INFO None 5157740: status FINISHED
2026-07-14 20:45:29 INFO None 5157741: status FINISHED
2026-07-14 20:45:29 INFO None 5157742: status FINISHED
2026-07-14 20:45:29 INFO None 5157743: status FINISHED
2026-07-14 20:45:29 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157745: status FINISHED
2026-07-14 20:45:29 INFO None 5157746: status FINISHED
2026-07-14 20:45:29 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157756: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157757: status RUNNING/PENDING
2026-07-14 20:45:29 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:45:29 INFO Jobs still running: ['5157744', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758']. Waiting...
2026-07-14 20:45:44 INFO None 5157740: status FINISHED
2026-07-14 20:45:44 INFO None 5157741: status FINISHED
2026-07-14 20:45:44 INFO None 5157742: status FINISHED
2026-07-14 20:45:44 INFO None 5157743: status FINISHED
2026-07-14 20:45:44 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:45:44 INFO None 5157745: status FINISHED
2026-07-14 20:45:44 INFO None 5157746: status FINISHED
2026-07-14 20:45:44 INFO None 5157751: status RUNNING/PENDING
2026-07-14 20:45:44 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:45:44 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:45:44 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:45:45 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:45:45 INFO None 5157756: status FINISHED
2026-07-14 20:45:45 INFO None 5157757: status FINISHED
2026-07-14 20:45:45 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:45:45 INFO Jobs still running: ['5157744', '5157751', '5157752', '5157753', '5157754', '5157755', '5157758']. Waiting...
2026-07-14 20:46:00 INFO None 5157740: status FINISHED
2026-07-14 20:46:00 INFO None 5157741: status FINISHED
2026-07-14 20:46:00 INFO None 5157742: status FINISHED
2026-07-14 20:46:00 INFO None 5157743: status FINISHED
2026-07-14 20:46:00 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:46:00 INFO None 5157745: status FINISHED
2026-07-14 20:46:00 INFO None 5157746: status FINISHED
2026-07-14 20:46:00 INFO None 5157751: status FINISHED
2026-07-14 20:46:00 INFO None 5157752: status RUNNING/PENDING
2026-07-14 20:46:00 INFO None 5157753: status RUNNING/PENDING
2026-07-14 20:46:00 INFO None 5157754: status RUNNING/PENDING
2026-07-14 20:46:00 INFO None 5157755: status RUNNING/PENDING
2026-07-14 20:46:00 INFO None 5157756: status FINISHED
2026-07-14 20:46:00 INFO None 5157757: status FINISHED
2026-07-14 20:46:00 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:46:00 INFO Jobs still running: ['5157744', '5157752', '5157753', '5157754', '5157755', '5157758']. Waiting...
2026-07-14 20:46:15 INFO None 5157740: status FINISHED
2026-07-14 20:48:09 INFO None 5157741: status FINISHED
2026-07-14 20:48:09 INFO None 5157742: status FINISHED
2026-07-14 20:48:09 INFO None 5157743: status FINISHED
2026-07-14 20:48:09 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:48:09 INFO None 5157745: status FINISHED
2026-07-14 20:48:09 INFO None 5157746: status FINISHED
2026-07-14 20:48:09 INFO None 5157751: status FINISHED
2026-07-14 20:48:09 INFO None 5157752: status FINISHED
2026-07-14 20:48:09 INFO None 5157753: status FINISHED
2026-07-14 20:48:09 INFO None 5157754: status FINISHED
2026-07-14 20:48:09 INFO None 5157755: status FINISHED
2026-07-14 20:48:09 INFO None 5157756: status FINISHED
2026-07-14 20:48:09 INFO None 5157757: status FINISHED
2026-07-14 20:48:09 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:48:09 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:48:24 INFO None 5157740: status FINISHED
2026-07-14 20:48:25 INFO None 5157741: status FINISHED
2026-07-14 20:48:25 INFO None 5157742: status FINISHED
2026-07-14 20:48:25 INFO None 5157743: status FINISHED
2026-07-14 20:48:25 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:48:25 INFO None 5157745: status FINISHED
2026-07-14 20:48:25 INFO None 5157746: status FINISHED
2026-07-14 20:48:25 INFO None 5157751: status FINISHED
2026-07-14 20:48:25 INFO None 5157752: status FINISHED
2026-07-14 20:48:25 INFO None 5157753: status FINISHED
2026-07-14 20:48:25 INFO None 5157754: status FINISHED
2026-07-14 20:48:25 INFO None 5157755: status FINISHED
2026-07-14 20:48:25 INFO None 5157756: status FINISHED
2026-07-14 20:48:25 INFO None 5157757: status FINISHED
2026-07-14 20:48:25 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:48:25 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:48:41 INFO None 5157740: status FINISHED
2026-07-14 20:48:41 INFO None 5157741: status FINISHED
2026-07-14 20:48:41 INFO None 5157742: status FINISHED
2026-07-14 20:48:41 INFO None 5157743: status FINISHED
2026-07-14 20:48:41 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:48:41 INFO None 5157745: status FINISHED
2026-07-14 20:48:41 INFO None 5157746: status FINISHED
2026-07-14 20:48:41 INFO None 5157751: status FINISHED
2026-07-14 20:48:41 INFO None 5157752: status FINISHED
2026-07-14 20:48:41 INFO None 5157753: status FINISHED
2026-07-14 20:48:41 INFO None 5157754: status FINISHED
2026-07-14 20:48:41 INFO None 5157755: status FINISHED
2026-07-14 20:48:41 INFO None 5157756: status FINISHED
2026-07-14 20:48:41 INFO None 5157757: status FINISHED
2026-07-14 20:48:42 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:48:42 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:48:57 INFO None 5157740: status FINISHED
2026-07-14 20:48:57 INFO None 5157741: status FINISHED
2026-07-14 20:48:57 INFO None 5157742: status FINISHED
2026-07-14 20:48:57 INFO None 5157743: status FINISHED
2026-07-14 20:48:57 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:48:57 INFO None 5157745: status FINISHED
2026-07-14 20:48:57 INFO None 5157746: status FINISHED
2026-07-14 20:48:57 INFO None 5157751: status FINISHED
2026-07-14 20:48:57 INFO None 5157752: status FINISHED
2026-07-14 20:48:57 INFO None 5157753: status FINISHED
2026-07-14 20:48:57 INFO None 5157754: status FINISHED
2026-07-14 20:48:57 INFO None 5157755: status FINISHED
2026-07-14 20:48:57 INFO None 5157756: status FINISHED
2026-07-14 20:48:57 INFO None 5157757: status FINISHED
2026-07-14 20:48:57 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:48:57 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:49:12 INFO None 5157740: status FINISHED
2026-07-14 20:49:12 INFO None 5157741: status FINISHED
2026-07-14 20:49:12 INFO None 5157742: status FINISHED
2026-07-14 20:49:12 INFO None 5157743: status FINISHED
2026-07-14 20:49:12 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:49:12 INFO None 5157745: status FINISHED
2026-07-14 20:49:12 INFO None 5157746: status FINISHED
2026-07-14 20:49:12 INFO None 5157751: status FINISHED
2026-07-14 20:49:12 INFO None 5157752: status FINISHED
2026-07-14 20:49:12 INFO None 5157753: status FINISHED
2026-07-14 20:49:12 INFO None 5157754: status FINISHED
2026-07-14 20:49:12 INFO None 5157755: status FINISHED
2026-07-14 20:49:12 INFO None 5157756: status FINISHED
2026-07-14 20:49:12 INFO None 5157757: status FINISHED
2026-07-14 20:49:12 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:49:12 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:49:27 INFO None 5157740: status FINISHED
2026-07-14 20:49:27 INFO None 5157741: status FINISHED
2026-07-14 20:49:27 INFO None 5157742: status FINISHED
2026-07-14 20:49:27 INFO None 5157743: status FINISHED
2026-07-14 20:49:27 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:49:27 INFO None 5157745: status FINISHED
2026-07-14 20:49:27 INFO None 5157746: status FINISHED
2026-07-14 20:49:27 INFO None 5157751: status FINISHED
2026-07-14 20:49:27 INFO None 5157752: status FINISHED
2026-07-14 20:49:27 INFO None 5157753: status FINISHED
2026-07-14 20:49:27 INFO None 5157754: status FINISHED
2026-07-14 20:49:27 INFO None 5157755: status FINISHED
2026-07-14 20:49:27 INFO None 5157756: status FINISHED
2026-07-14 20:49:27 INFO None 5157757: status FINISHED
2026-07-14 20:49:27 INFO None 5157758: status RUNNING/PENDING
2026-07-14 20:49:27 INFO Jobs still running: ['5157744', '5157758']. Waiting...
2026-07-14 20:49:42 INFO None 5157740: status FINISHED
2026-07-14 20:49:42 INFO None 5157741: status FINISHED
2026-07-14 20:49:42 INFO None 5157742: status FINISHED
2026-07-14 20:49:42 INFO None 5157743: status FINISHED
2026-07-14 20:49:42 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:49:42 INFO None 5157745: status FINISHED
2026-07-14 20:49:43 INFO None 5157746: status FINISHED
2026-07-14 20:49:43 INFO None 5157751: status FINISHED
2026-07-14 20:49:43 INFO None 5157752: status FINISHED
2026-07-14 20:49:43 INFO None 5157753: status FINISHED
2026-07-14 20:49:43 INFO None 5157754: status FINISHED
2026-07-14 20:49:43 INFO None 5157755: status FINISHED
2026-07-14 20:49:43 INFO None 5157756: status FINISHED
2026-07-14 20:49:43 INFO None 5157757: status FINISHED
2026-07-14 20:49:43 INFO None 5157758: status FINISHED
2026-07-14 20:49:43 INFO Jobs still running: ['5157744']. Waiting...
2026-07-14 20:49:58 INFO None 5157740: status FINISHED
2026-07-14 20:49:58 INFO None 5157741: status FINISHED
2026-07-14 20:49:58 INFO None 5157742: status FINISHED
2026-07-14 20:49:58 INFO None 5157743: status FINISHED
2026-07-14 20:49:58 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:49:58 INFO None 5157745: status FINISHED
2026-07-14 20:49:58 INFO None 5157746: status FINISHED
2026-07-14 20:49:58 INFO None 5157751: status FINISHED
2026-07-14 20:49:58 INFO None 5157752: status FINISHED
2026-07-14 20:49:58 INFO None 5157753: status FINISHED
2026-07-14 20:49:58 INFO None 5157754: status FINISHED
2026-07-14 20:49:58 INFO None 5157755: status FINISHED
2026-07-14 20:49:58 INFO None 5157756: status FINISHED
2026-07-14 20:49:58 INFO None 5157757: status FINISHED
2026-07-14 20:49:58 INFO None 5157758: status FINISHED
2026-07-14 20:49:58 INFO Jobs still running: ['5157744']. Waiting...
2026-07-14 20:50:13 INFO None 5157740: status FINISHED
2026-07-14 20:50:13 INFO None 5157741: status FINISHED
2026-07-14 20:50:13 INFO None 5157742: status FINISHED
2026-07-14 20:50:13 INFO None 5157743: status FINISHED
2026-07-14 20:50:13 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:50:13 INFO None 5157745: status FINISHED
2026-07-14 20:50:13 INFO None 5157746: status FINISHED
2026-07-14 20:50:13 INFO None 5157751: status FINISHED
2026-07-14 20:50:13 INFO None 5157752: status FINISHED
2026-07-14 20:50:13 INFO None 5157753: status FINISHED
2026-07-14 20:50:13 INFO None 5157754: status FINISHED
2026-07-14 20:50:13 INFO None 5157755: status FINISHED
2026-07-14 20:50:13 INFO None 5157756: status FINISHED
2026-07-14 20:50:13 INFO None 5157757: status FINISHED
2026-07-14 20:50:13 INFO None 5157758: status FINISHED
2026-07-14 20:50:13 INFO Jobs still running: ['5157744']. Waiting...
2026-07-14 20:50:29 INFO None 5157740: status FINISHED
2026-07-14 20:50:29 INFO None 5157741: status FINISHED
2026-07-14 20:50:29 INFO None 5157742: status FINISHED
2026-07-14 20:50:29 INFO None 5157743: status FINISHED
2026-07-14 20:50:29 INFO None 5157744: status RUNNING/PENDING
2026-07-14 20:50:29 INFO None 5157745: status FINISHED
2026-07-14 20:50:29 INFO None 5157746: status FINISHED
2026-07-14 20:50:29 INFO None 5157751: status FINISHED
2026-07-14 20:50:29 INFO None 5157752: status FINISHED
2026-07-14 20:50:29 INFO None 5157753: status FINISHED
2026-07-14 20:50:29 INFO None 5157754: status FINISHED
2026-07-14 20:50:29 INFO None 5157755: status FINISHED
2026-07-14 20:50:29 INFO None 5157756: status FINISHED
2026-07-14 20:50:29 INFO None 5157757: status FINISHED
2026-07-14 20:50:30 INFO None 5157758: status FINISHED
2026-07-14 20:50:30 INFO Jobs still running: ['5157744']. Waiting...
2026-07-14 20:50:45 INFO None 5157740: status FINISHED
2026-07-14 20:50:45 INFO None 5157741: status FINISHED
2026-07-14 20:50:45 INFO None 5157742: status FINISHED
2026-07-14 20:50:45 INFO None 5157743: status FINISHED
2026-07-14 20:50:45 INFO None 5157744: status FINISHED
2026-07-14 20:50:45 INFO None 5157745: status FINISHED
2026-07-14 20:50:45 INFO None 5157746: status FINISHED
2026-07-14 20:50:45 INFO None 5157751: status FINISHED
2026-07-14 20:50:45 INFO None 5157752: status FINISHED
2026-07-14 20:50:45 INFO None 5157753: status FINISHED
2026-07-14 20:50:45 INFO None 5157754: status FINISHED
2026-07-14 20:50:45 INFO None 5157755: status FINISHED
2026-07-14 20:50:45 INFO None 5157756: status FINISHED
2026-07-14 20:50:45 INFO None 5157757: status FINISHED
2026-07-14 20:50:45 INFO None 5157758: status FINISHED
2026-07-14 20:50:45 INFO Jobs ['5157740', '5157741', '5157742', '5157743', '5157744', '5157745', '5157746', '5157751', '5157752', '5157753', '5157754', '5157755', '5157756', '5157757', '5157758'] have finished
2026-07-14 20:50:45 INFO Checking restart files were created ...
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-14 20:50:45 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-14 20:50:45 INFO  Run_model() completed successfully.
2026-07-14 20:50:45 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:50:45 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-14 20:50:45 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 20:50:45 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-14 20:50:45 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:50:45 INFO ---------->>> Running process_satellite_data()
2026-07-14 20:50:45 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-14 20:50:45 INFO ---------->>> Running run_obs_converter()
2026-07-14 20:50:45 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-14 20:50:45 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-14 20:50:45 INFO ---------->>> Running DART
2026-07-14 20:50:45 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 20:50:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-14 20:50:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-14 20:50:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-14 20:50:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-14 20:50:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-14 20:50:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-14 20:50:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-14 20:50:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-14 20:50:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-14 20:50:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-14 20:50:48 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-14 20:50:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-14 20:50:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-14 20:50:49 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-14 20:50:50 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-14 20:50:50 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 20:50:50 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 20:50:50 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 20:50:50 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 20:50:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 20:50:50 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 20:50:58 INFO Found: []
2026-07-14 20:50:58 INFO No job id returned by command ./run_filter.bsh
2026-07-14 20:50:58 INFO No monitoring will be performed
2026-07-14 20:50:58 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-14 20:50:58 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:58 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-14 20:50:59 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 20:50:59 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 20:50:59 INFO run_dart() is DONE.
2026-07-14 20:50:59 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 20:50:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-14 20:53:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-14 20:53:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-14 20:53:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-14 20:53:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-14 20:53:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-14 20:53:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-14 20:53:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-14 20:53:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-14 20:53:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-14 20:53:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:53:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-14 20:54:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:54:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-14 20:54:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:54:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-14 20:54:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:54:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-14 20:54:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:54:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-14 20:54:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 20:54:18 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 20:54:18 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:54:18 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 20:54:18 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-14 20:54:18 INFO Copying EMIS of next day ...
2026-07-14 20:54:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-14 20:54:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:54:20 INFO Hourly dataset computed and listing created
2026-07-14 20:54:39 INFO Hourly dataset computed
2026-07-14 20:54:39 INFO Copying EMIS of next day ...
2026-07-14 20:54:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-14 20:54:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:54:41 INFO Hourly dataset computed and listing created
2026-07-14 20:54:56 INFO Hourly dataset computed
2026-07-14 20:54:56 INFO Copying EMIS of next day ...
2026-07-14 20:54:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-14 20:54:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:54:58 INFO Hourly dataset computed and listing created
2026-07-14 20:55:08 INFO Hourly dataset computed
2026-07-14 20:55:08 INFO Copying EMIS of next day ...
2026-07-14 20:55:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-14 20:55:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:10 INFO Hourly dataset computed and listing created
2026-07-14 20:55:12 INFO Hourly dataset computed
2026-07-14 20:55:12 INFO Copying EMIS of next day ...
2026-07-14 20:55:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-14 20:55:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:14 INFO Hourly dataset computed and listing created
2026-07-14 20:55:42 INFO Hourly dataset computed
2026-07-14 20:55:42 INFO Copying EMIS of next day ...
2026-07-14 20:55:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-14 20:55:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:44 INFO Hourly dataset computed and listing created
2026-07-14 20:55:47 INFO Hourly dataset computed
2026-07-14 20:55:47 INFO Copying EMIS of next day ...
2026-07-14 20:55:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-14 20:55:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:48 INFO Hourly dataset computed and listing created
2026-07-14 20:55:50 INFO Hourly dataset computed
2026-07-14 20:55:50 INFO Copying EMIS of next day ...
2026-07-14 20:55:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-14 20:55:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:52 INFO Hourly dataset computed and listing created
2026-07-14 20:55:54 INFO Hourly dataset computed
2026-07-14 20:55:54 INFO Copying EMIS of next day ...
2026-07-14 20:55:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-14 20:55:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:55:56 INFO Hourly dataset computed and listing created
2026-07-14 20:55:59 INFO Hourly dataset computed
2026-07-14 20:55:59 INFO Copying EMIS of next day ...
2026-07-14 20:55:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-14 20:55:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:09 INFO Hourly dataset computed and listing created
2026-07-14 20:58:11 INFO Hourly dataset computed
2026-07-14 20:58:11 INFO Copying EMIS of next day ...
2026-07-14 20:58:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-14 20:58:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:12 INFO Hourly dataset computed and listing created
2026-07-14 20:58:15 INFO Hourly dataset computed
2026-07-14 20:58:15 INFO Copying EMIS of next day ...
2026-07-14 20:58:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-14 20:58:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:16 INFO Hourly dataset computed and listing created
2026-07-14 20:58:19 INFO Hourly dataset computed
2026-07-14 20:58:19 INFO Copying EMIS of next day ...
2026-07-14 20:58:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-14 20:58:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:21 INFO Hourly dataset computed and listing created
2026-07-14 20:58:23 INFO Hourly dataset computed
2026-07-14 20:58:23 INFO Copying EMIS of next day ...
2026-07-14 20:58:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-14 20:58:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:24 INFO Hourly dataset computed and listing created
2026-07-14 20:58:27 INFO Hourly dataset computed
2026-07-14 20:58:27 INFO Copying EMIS of next day ...
2026-07-14 20:58:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-14 20:58:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 20:58:29 INFO Hourly dataset computed and listing created
2026-07-14 20:58:48 INFO Hourly dataset computed
2026-07-14 20:58:48 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-14 20:58:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:58:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 20:58:48 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-14 20:58:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 20:58:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:58:48 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 20:58:48 INFO Queuing job for member 1...
2026-07-14 20:58:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:58:48 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 20:58:49 INFO Found: ['5157807']
2026-07-14 20:58:54 INFO [TGCC-IRENE] Submitted job with ID:['5157807']
2026-07-14 20:58:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:58:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 20:58:54 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-14 20:58:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 20:58:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:58:54 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 20:58:54 INFO Queuing job for member 2...
2026-07-14 20:58:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:58:54 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 20:58:56 INFO Found: ['5157808']
2026-07-14 20:59:01 INFO [TGCC-IRENE] Submitted job with ID:['5157808']
2026-07-14 20:59:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 20:59:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-14 20:59:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 20:59:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 20:59:01 INFO Queuing job for member 3...
2026-07-14 20:59:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 20:59:01 INFO Found: ['5157810']
2026-07-14 20:59:06 INFO [TGCC-IRENE] Submitted job with ID:['5157810']
2026-07-14 20:59:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 20:59:06 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-14 20:59:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 20:59:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:06 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 20:59:06 INFO Queuing job for member 4...
2026-07-14 20:59:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:06 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 20:59:07 INFO Found: ['5157811']
2026-07-14 20:59:12 INFO [TGCC-IRENE] Submitted job with ID:['5157811']
2026-07-14 20:59:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 20:59:12 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-14 20:59:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 20:59:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:12 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 20:59:12 INFO Queuing job for member 5...
2026-07-14 20:59:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:12 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 20:59:13 INFO Found: ['5157812']
2026-07-14 20:59:18 INFO [TGCC-IRENE] Submitted job with ID:['5157812']
2026-07-14 20:59:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 20:59:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-14 20:59:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 20:59:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 20:59:18 INFO Queuing job for member 6...
2026-07-14 20:59:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 20:59:19 INFO Found: ['5157813']
2026-07-14 20:59:24 INFO [TGCC-IRENE] Submitted job with ID:['5157813']
2026-07-14 20:59:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 20:59:24 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-14 20:59:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 20:59:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:24 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 20:59:24 INFO Queuing job for member 7...
2026-07-14 20:59:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:24 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 20:59:24 INFO Found: ['5157814']
2026-07-14 20:59:29 INFO [TGCC-IRENE] Submitted job with ID:['5157814']
2026-07-14 20:59:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 20:59:29 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-14 20:59:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 20:59:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:29 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 20:59:29 INFO Queuing job for member 8...
2026-07-14 20:59:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:29 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 20:59:30 INFO Found: ['5157815']
2026-07-14 20:59:35 INFO [TGCC-IRENE] Submitted job with ID:['5157815']
2026-07-14 20:59:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 20:59:35 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-14 20:59:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 20:59:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:35 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 20:59:35 INFO Queuing job for member 9...
2026-07-14 20:59:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:35 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 20:59:36 INFO Found: ['5157817']
2026-07-14 20:59:41 INFO [TGCC-IRENE] Submitted job with ID:['5157817']
2026-07-14 20:59:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 20:59:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-14 20:59:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 20:59:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 20:59:41 INFO Queuing job for member 10...
2026-07-14 20:59:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 20:59:43 INFO Found: ['5157819']
2026-07-14 20:59:48 INFO [TGCC-IRENE] Submitted job with ID:['5157819']
2026-07-14 20:59:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 20:59:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-14 20:59:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 20:59:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 20:59:48 INFO Queuing job for member 11...
2026-07-14 20:59:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 20:59:50 INFO Found: ['5157820']
2026-07-14 20:59:55 INFO [TGCC-IRENE] Submitted job with ID:['5157820']
2026-07-14 20:59:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 20:59:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 20:59:55 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-14 20:59:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 20:59:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 20:59:55 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 20:59:55 INFO Queuing job for member 12...
2026-07-14 20:59:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 20:59:55 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 20:59:56 INFO Found: ['5157821']
2026-07-14 21:00:01 INFO [TGCC-IRENE] Submitted job with ID:['5157821']
2026-07-14 21:00:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 21:00:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 21:00:01 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-14 21:00:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 21:00:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 21:00:01 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 21:00:01 INFO Queuing job for member 13...
2026-07-14 21:00:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 21:00:01 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 21:00:02 INFO Found: ['5157823']
2026-07-14 21:00:07 INFO [TGCC-IRENE] Submitted job with ID:['5157823']
2026-07-14 21:00:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 21:00:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 21:00:07 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-14 21:00:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 21:00:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 21:00:07 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 21:00:07 INFO Queuing job for member 14...
2026-07-14 21:00:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 21:00:07 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 21:00:07 INFO Found: ['5157824']
2026-07-14 21:00:12 INFO [TGCC-IRENE] Submitted job with ID:['5157824']
2026-07-14 21:00:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 21:00:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 21:00:12 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-14 21:00:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 21:00:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 21:00:12 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 21:00:12 INFO Queuing job for member 15...
2026-07-14 21:00:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 21:00:12 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 21:00:13 INFO Found: ['5157825']
2026-07-14 21:00:18 INFO [TGCC-IRENE] Submitted job with ID:['5157825']
2026-07-14 21:00:18 INFO Checking job status ...
2026-07-14 21:00:18 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:00:18 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:00:18 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:00:33 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:00:34 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:00:34 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:00:49 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:00:49 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:00:49 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:01:04 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:01:04 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:01:04 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:01:19 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:01:19 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:01:20 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:01:20 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:01:36 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157820: status RUNNING/PENDING
2026-07-14 21:03:33 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:03:35 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:03:35 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:03:35 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:03:35 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:03:50 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157820: status FINISHED
2026-07-14 21:03:50 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:03:50 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:03:51 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:03:51 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:04:06 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157820: status FINISHED
2026-07-14 21:04:06 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:04:06 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:04:06 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:04:21 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157820: status FINISHED
2026-07-14 21:04:21 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:04:21 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:04:21 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:04:36 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:04:36 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:04:36 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:04:36 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:04:36 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157820: status FINISHED
2026-07-14 21:04:38 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:04:38 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:04:38 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:04:53 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:04:53 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:04:53 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:04:53 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157820: status FINISHED
2026-07-14 21:04:54 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:04:54 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:04:54 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:05:09 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157820: status FINISHED
2026-07-14 21:05:09 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:05:09 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:05:09 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:05:24 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157820: status FINISHED
2026-07-14 21:05:24 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:05:24 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:05:24 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:05:39 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:05:39 INFO None 5157820: status FINISHED
2026-07-14 21:05:40 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:05:40 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:05:40 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:05:40 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:05:40 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:05:55 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157820: status FINISHED
2026-07-14 21:05:55 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:05:55 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:05:55 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:06:10 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157820: status FINISHED
2026-07-14 21:06:10 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:06:10 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:06:10 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:06:26 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157820: status FINISHED
2026-07-14 21:06:26 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:06:26 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:06:27 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:06:27 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:06:27 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:06:42 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157820: status FINISHED
2026-07-14 21:08:39 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:08:39 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:08:39 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:08:54 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:08:54 INFO None 5157820: status FINISHED
2026-07-14 21:08:55 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:08:55 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:08:55 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:08:55 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:08:55 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:09:10 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157820: status FINISHED
2026-07-14 21:09:10 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:09:10 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:09:10 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:09:26 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157820: status FINISHED
2026-07-14 21:09:26 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:09:26 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:09:26 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:09:41 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:09:41 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:09:42 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:09:42 INFO None 5157820: status FINISHED
2026-07-14 21:09:42 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:09:42 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:09:42 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:09:42 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:09:42 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:09:57 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157820: status FINISHED
2026-07-14 21:09:57 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:09:57 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:09:57 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:10:12 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157820: status FINISHED
2026-07-14 21:10:12 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:10:12 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:10:12 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:10:27 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157820: status FINISHED
2026-07-14 21:10:27 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:10:27 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:10:27 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:10:42 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157820: status FINISHED
2026-07-14 21:10:43 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:10:43 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:10:43 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:10:58 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157820: status FINISHED
2026-07-14 21:10:58 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:10:58 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:10:58 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:11:13 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157820: status FINISHED
2026-07-14 21:13:12 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:13:12 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:13:12 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:13:27 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:13:27 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:13:27 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157820: status FINISHED
2026-07-14 21:13:28 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:13:28 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:13:28 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:13:43 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157820: status FINISHED
2026-07-14 21:13:43 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:13:43 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:13:43 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:13:58 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:13:58 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:13:58 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:13:58 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:14:15 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:14:15 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:14:15 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:14:15 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:14:16 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:14:16 INFO None 5157825: status RUNNING/PENDING
2026-07-14 21:14:16 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824', '5157825']. Waiting...
2026-07-14 21:14:31 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:14:31 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:14:31 INFO None 5157825: status FINISHED
2026-07-14 21:14:31 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:14:46 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:14:46 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:14:46 INFO None 5157825: status FINISHED
2026-07-14 21:14:46 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:15:01 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:15:01 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:15:01 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:15:02 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:15:02 INFO None 5157825: status FINISHED
2026-07-14 21:15:02 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:15:17 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157808: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:15:17 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:15:17 INFO None 5157825: status FINISHED
2026-07-14 21:15:17 INFO Jobs still running: ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:15:32 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157808: status FINISHED
2026-07-14 21:15:32 INFO None 5157810: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157811: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157812: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:15:32 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:15:33 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:15:33 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:15:33 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:15:33 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:15:33 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:15:33 INFO None 5157825: status FINISHED
2026-07-14 21:15:33 INFO Jobs still running: ['5157807', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:15:48 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157808: status FINISHED
2026-07-14 21:15:48 INFO None 5157810: status FINISHED
2026-07-14 21:15:48 INFO None 5157811: status FINISHED
2026-07-14 21:15:48 INFO None 5157812: status FINISHED
2026-07-14 21:15:48 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:15:48 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:15:48 INFO None 5157825: status FINISHED
2026-07-14 21:15:48 INFO Jobs still running: ['5157807', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:16:04 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157808: status FINISHED
2026-07-14 21:16:04 INFO None 5157810: status FINISHED
2026-07-14 21:16:04 INFO None 5157811: status FINISHED
2026-07-14 21:16:04 INFO None 5157812: status FINISHED
2026-07-14 21:16:04 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:16:04 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:16:04 INFO None 5157825: status FINISHED
2026-07-14 21:16:04 INFO Jobs still running: ['5157807', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:16:19 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:16:19 INFO None 5157808: status FINISHED
2026-07-14 21:16:20 INFO None 5157810: status FINISHED
2026-07-14 21:16:20 INFO None 5157811: status FINISHED
2026-07-14 21:16:20 INFO None 5157812: status FINISHED
2026-07-14 21:16:20 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157814: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157815: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157817: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157819: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:16:20 INFO None 5157821: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157823: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157824: status RUNNING/PENDING
2026-07-14 21:16:20 INFO None 5157825: status FINISHED
2026-07-14 21:16:20 INFO Jobs still running: ['5157807', '5157813', '5157814', '5157815', '5157817', '5157819', '5157821', '5157823', '5157824']. Waiting...
2026-07-14 21:16:35 INFO None 5157807: status RUNNING/PENDING
2026-07-14 21:18:33 INFO None 5157808: status FINISHED
2026-07-14 21:18:33 INFO None 5157810: status FINISHED
2026-07-14 21:18:33 INFO None 5157811: status FINISHED
2026-07-14 21:18:33 INFO None 5157812: status FINISHED
2026-07-14 21:18:33 INFO None 5157813: status RUNNING/PENDING
2026-07-14 21:18:33 INFO None 5157814: status FINISHED
2026-07-14 21:18:33 INFO None 5157815: status FINISHED
2026-07-14 21:18:33 INFO None 5157817: status FINISHED
2026-07-14 21:18:33 INFO None 5157819: status FINISHED
2026-07-14 21:18:33 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:18:33 INFO None 5157821: status FINISHED
2026-07-14 21:18:33 INFO None 5157823: status FINISHED
2026-07-14 21:18:33 INFO None 5157824: status FINISHED
2026-07-14 21:18:33 INFO None 5157825: status FINISHED
2026-07-14 21:18:33 INFO Jobs still running: ['5157807', '5157813']. Waiting...
2026-07-14 21:18:48 INFO None 5157807: status FINISHED
2026-07-14 21:18:48 INFO None 5157808: status FINISHED
2026-07-14 21:18:48 INFO None 5157810: status FINISHED
2026-07-14 21:18:48 INFO None 5157811: status FINISHED
2026-07-14 21:18:48 INFO None 5157812: status FINISHED
2026-07-14 21:18:48 INFO None 5157813: status FINISHED
2026-07-14 21:18:48 INFO None 5157814: status FINISHED
2026-07-14 21:18:48 INFO None 5157815: status FINISHED
2026-07-14 21:18:48 INFO None 5157817: status FINISHED
2026-07-14 21:18:48 INFO None 5157819: status FINISHED
2026-07-14 21:18:48 INFO None 5157820: status FINISHED (not in squeue)
2026-07-14 21:18:48 INFO None 5157821: status FINISHED
2026-07-14 21:18:48 INFO None 5157823: status FINISHED
2026-07-14 21:18:48 INFO None 5157824: status FINISHED
2026-07-14 21:18:48 INFO None 5157825: status FINISHED
2026-07-14 21:18:48 INFO Jobs ['5157807', '5157808', '5157810', '5157811', '5157812', '5157813', '5157814', '5157815', '5157817', '5157819', '5157820', '5157821', '5157823', '5157824', '5157825'] have finished
2026-07-14 21:18:48 INFO Checking restart files were created ...
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-14 21:18:48 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-14 21:18:48 INFO Check chimere log file at: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/ENS15_2020020614.out
2026-07-14 21:18:48 ERROR [PIPELINE] Error: The following chimere ENS run(s) failed (exit code 1): [11]
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 334, in run_model
    raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
pipeline_errors.ModelRunError: The following chimere ENS run(s) failed (exit code 1): [11]
+ exit 0
