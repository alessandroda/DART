+ SCRIPT_PID=1559787
+ /bin/bash -x /tmp/tmp.mCVP4X4NEJ
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
2026-07-15 16:35:22 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-15 16:35:22 INFO [PIPELINE] =======================================
2026-07-15 16:35:22 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-15 16:35:22 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-15 16:35:22 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-15 16:35:22 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260715_163522.log
2026-07-15 16:35:22 INFO [PIPELINE] =======================================
2026-07-15 16:35:22 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-15 16:35:22 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-15 16:35:22 INFO [STEP] ---- TIME LOOP START ----
2026-07-15 16:35:22 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:35:22 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-15 16:35:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:31 INFO Hourly dataset computed and listing created
2026-07-15 16:35:37 INFO Hourly dataset computed
2026-07-15 16:35:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:38 INFO Hourly dataset computed and listing created
2026-07-15 16:35:41 INFO Hourly dataset computed
2026-07-15 16:35:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:42 INFO Hourly dataset computed and listing created
2026-07-15 16:35:45 INFO Hourly dataset computed
2026-07-15 16:35:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:46 INFO Hourly dataset computed and listing created
2026-07-15 16:35:48 INFO Hourly dataset computed
2026-07-15 16:35:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:49 INFO Hourly dataset computed and listing created
2026-07-15 16:35:52 INFO Hourly dataset computed
2026-07-15 16:35:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:53 INFO Hourly dataset computed and listing created
2026-07-15 16:35:55 INFO Hourly dataset computed
2026-07-15 16:35:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:35:57 INFO Hourly dataset computed and listing created
2026-07-15 16:36:00 INFO Hourly dataset computed
2026-07-15 16:36:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:01 INFO Hourly dataset computed and listing created
2026-07-15 16:36:03 INFO Hourly dataset computed
2026-07-15 16:36:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:05 INFO Hourly dataset computed and listing created
2026-07-15 16:36:08 INFO Hourly dataset computed
2026-07-15 16:36:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:09 INFO Hourly dataset computed and listing created
2026-07-15 16:36:12 INFO Hourly dataset computed
2026-07-15 16:36:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:13 INFO Hourly dataset computed and listing created
2026-07-15 16:36:15 INFO Hourly dataset computed
2026-07-15 16:36:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:17 INFO Hourly dataset computed and listing created
2026-07-15 16:36:19 INFO Hourly dataset computed
2026-07-15 16:36:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:20 INFO Hourly dataset computed and listing created
2026-07-15 16:36:22 INFO Hourly dataset computed
2026-07-15 16:36:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:23 INFO Hourly dataset computed and listing created
2026-07-15 16:36:25 INFO Hourly dataset computed
2026-07-15 16:36:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:36:26 INFO Hourly dataset computed and listing created
2026-07-15 16:36:29 INFO Hourly dataset computed
2026-07-15 16:36:29 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-15 16:36:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:36:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 16:36:29 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-15 16:36:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 16:36:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:36:29 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 16:36:29 INFO Queuing job for member 1...
2026-07-15 16:36:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:36:29 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 16:36:30 INFO Found: ['5165568']
2026-07-15 16:36:35 INFO [TGCC-IRENE] Submitted job with ID:['5165568']
2026-07-15 16:36:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:36:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 16:36:35 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-15 16:36:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 16:36:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:36:35 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 16:36:35 INFO Queuing job for member 2...
2026-07-15 16:36:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:36:35 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 16:36:35 INFO Found: ['5165569']
2026-07-15 16:36:40 INFO [TGCC-IRENE] Submitted job with ID:['5165569']
2026-07-15 16:36:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:36:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 16:36:40 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-15 16:36:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 16:36:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:36:40 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 16:36:40 INFO Queuing job for member 3...
2026-07-15 16:36:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:36:40 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 16:36:41 INFO Found: ['5165570']
2026-07-15 16:36:46 INFO [TGCC-IRENE] Submitted job with ID:['5165570']
2026-07-15 16:36:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:36:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 16:36:46 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-15 16:36:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 16:36:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:36:46 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 16:36:46 INFO Queuing job for member 4...
2026-07-15 16:36:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:36:46 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 16:36:48 INFO Found: ['5165571']
2026-07-15 16:36:53 INFO [TGCC-IRENE] Submitted job with ID:['5165571']
2026-07-15 16:36:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:36:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 16:36:53 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-15 16:36:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 16:36:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:36:53 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 16:36:53 INFO Queuing job for member 5...
2026-07-15 16:36:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:36:53 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 16:36:55 INFO Found: ['5165572']
2026-07-15 16:37:00 INFO [TGCC-IRENE] Submitted job with ID:['5165572']
2026-07-15 16:37:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 16:37:00 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-15 16:37:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 16:37:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:00 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 16:37:00 INFO Queuing job for member 6...
2026-07-15 16:37:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:00 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 16:37:01 INFO Found: ['5165573']
2026-07-15 16:37:06 INFO [TGCC-IRENE] Submitted job with ID:['5165573']
2026-07-15 16:37:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 16:37:06 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-15 16:37:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 16:37:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:06 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 16:37:06 INFO Queuing job for member 7...
2026-07-15 16:37:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:06 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 16:37:07 INFO Found: ['5165575']
2026-07-15 16:37:12 INFO [TGCC-IRENE] Submitted job with ID:['5165575']
2026-07-15 16:37:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 16:37:12 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-15 16:37:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 16:37:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:12 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 16:37:12 INFO Queuing job for member 8...
2026-07-15 16:37:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:12 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 16:37:12 INFO Found: ['5165577']
2026-07-15 16:37:17 INFO [TGCC-IRENE] Submitted job with ID:['5165577']
2026-07-15 16:37:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 16:37:17 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-15 16:37:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 16:37:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:17 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 16:37:17 INFO Queuing job for member 9...
2026-07-15 16:37:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:17 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 16:37:18 INFO Found: ['5165578']
2026-07-15 16:37:23 INFO [TGCC-IRENE] Submitted job with ID:['5165578']
2026-07-15 16:37:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 16:37:23 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-15 16:37:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 16:37:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:23 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 16:37:23 INFO Queuing job for member 10...
2026-07-15 16:37:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:23 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 16:37:24 INFO Found: ['5165580']
2026-07-15 16:37:29 INFO [TGCC-IRENE] Submitted job with ID:['5165580']
2026-07-15 16:37:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 16:37:29 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-15 16:37:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 16:37:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:29 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 16:37:29 INFO Queuing job for member 11...
2026-07-15 16:37:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:29 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 16:37:30 INFO Found: ['5165581']
2026-07-15 16:37:35 INFO [TGCC-IRENE] Submitted job with ID:['5165581']
2026-07-15 16:37:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 16:37:35 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-15 16:37:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 16:37:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:35 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 16:37:35 INFO Queuing job for member 12...
2026-07-15 16:37:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:35 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 16:37:36 INFO Found: ['5165582']
2026-07-15 16:37:41 INFO [TGCC-IRENE] Submitted job with ID:['5165582']
2026-07-15 16:37:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 16:37:41 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-15 16:37:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 16:37:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:41 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 16:37:41 INFO Queuing job for member 13...
2026-07-15 16:37:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:41 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 16:37:43 INFO Found: ['5165583']
2026-07-15 16:37:48 INFO [TGCC-IRENE] Submitted job with ID:['5165583']
2026-07-15 16:37:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 16:37:48 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-15 16:37:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 16:37:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:48 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 16:37:48 INFO Queuing job for member 14...
2026-07-15 16:37:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:48 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 16:37:49 INFO Found: ['5165584']
2026-07-15 16:37:54 INFO [TGCC-IRENE] Submitted job with ID:['5165584']
2026-07-15 16:37:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:37:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 16:37:54 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-15 16:37:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 16:37:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:37:54 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 16:37:54 INFO Queuing job for member 15...
2026-07-15 16:37:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:37:54 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 16:37:55 INFO Found: ['5165586']
2026-07-15 16:38:00 INFO [TGCC-IRENE] Submitted job with ID:['5165586']
2026-07-15 16:38:00 INFO Checking job status ...
2026-07-15 16:38:00 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:38:00 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:38:00 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:38:15 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:38:15 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:38:15 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:38:30 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:38:30 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:38:31 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:38:31 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:38:46 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:38:46 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:38:46 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:39:01 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:39:01 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:39:01 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:39:16 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:39:16 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:39:16 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:39:16 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:39:17 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:39:17 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:39:32 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:39:32 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:39:32 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:39:32 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:39:32 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:39:34 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:39:34 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:39:49 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:39:49 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:39:49 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:40:04 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165570: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:40:04 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:40:05 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:40:05 INFO Jobs still running: ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:40:20 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165570: status FINISHED
2026-07-15 16:40:20 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:40:20 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:40:20 INFO Jobs still running: ['5165568', '5165569', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:40:37 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165569: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165570: status FINISHED
2026-07-15 16:40:37 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:40:37 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:40:37 INFO Jobs still running: ['5165568', '5165569', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:40:52 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165569: status FINISHED
2026-07-15 16:40:52 INFO None 5165570: status FINISHED
2026-07-15 16:40:52 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:40:52 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:40:52 INFO Jobs still running: ['5165568', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:41:07 INFO None 5165568: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165569: status FINISHED
2026-07-15 16:41:07 INFO None 5165570: status FINISHED
2026-07-15 16:41:07 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:41:07 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:41:08 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:41:08 INFO Jobs still running: ['5165568', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:41:23 INFO None 5165568: status FINISHED
2026-07-15 16:41:23 INFO None 5165569: status FINISHED
2026-07-15 16:41:23 INFO None 5165570: status FINISHED
2026-07-15 16:41:23 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165577: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165578: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:41:23 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:41:23 INFO Jobs still running: ['5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:41:38 INFO None 5165568: status FINISHED
2026-07-15 16:41:38 INFO None 5165569: status FINISHED
2026-07-15 16:41:38 INFO None 5165570: status FINISHED
2026-07-15 16:41:38 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165577: status FINISHED
2026-07-15 16:41:38 INFO None 5165578: status FINISHED
2026-07-15 16:41:38 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:41:38 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:41:38 INFO Jobs still running: ['5165571', '5165572', '5165573', '5165575', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:41:53 INFO None 5165568: status FINISHED
2026-07-15 16:41:53 INFO None 5165569: status FINISHED
2026-07-15 16:41:53 INFO None 5165570: status FINISHED
2026-07-15 16:41:53 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165573: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165577: status FINISHED
2026-07-15 16:41:53 INFO None 5165578: status FINISHED
2026-07-15 16:41:53 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:41:53 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:41:54 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:41:54 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:41:54 INFO Jobs still running: ['5165571', '5165572', '5165573', '5165575', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:42:09 INFO None 5165568: status FINISHED
2026-07-15 16:42:09 INFO None 5165569: status FINISHED
2026-07-15 16:42:09 INFO None 5165570: status FINISHED
2026-07-15 16:42:09 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165573: status FINISHED
2026-07-15 16:42:09 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165577: status FINISHED
2026-07-15 16:42:09 INFO None 5165578: status FINISHED
2026-07-15 16:42:09 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:42:09 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:42:09 INFO Jobs still running: ['5165571', '5165572', '5165575', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:42:25 INFO None 5165568: status FINISHED
2026-07-15 16:42:25 INFO None 5165569: status FINISHED
2026-07-15 16:42:25 INFO None 5165570: status FINISHED
2026-07-15 16:42:25 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165573: status FINISHED
2026-07-15 16:42:25 INFO None 5165575: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165577: status FINISHED
2026-07-15 16:42:25 INFO None 5165578: status FINISHED
2026-07-15 16:42:25 INFO None 5165580: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165581: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165582: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:42:25 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:42:25 INFO Jobs still running: ['5165571', '5165572', '5165575', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:42:40 INFO None 5165568: status FINISHED
2026-07-15 16:42:40 INFO None 5165569: status FINISHED
2026-07-15 16:42:40 INFO None 5165570: status FINISHED
2026-07-15 16:42:40 INFO None 5165571: status RUNNING/PENDING
2026-07-15 16:42:40 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:42:40 INFO None 5165573: status FINISHED
2026-07-15 16:42:40 INFO None 5165575: status FINISHED
2026-07-15 16:42:40 INFO None 5165577: status FINISHED
2026-07-15 16:42:40 INFO None 5165578: status FINISHED
2026-07-15 16:42:40 INFO None 5165580: status FINISHED
2026-07-15 16:42:40 INFO None 5165581: status FINISHED
2026-07-15 16:42:40 INFO None 5165582: status FINISHED
2026-07-15 16:42:40 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:42:41 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:42:41 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:42:41 INFO Jobs still running: ['5165571', '5165572', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:42:56 INFO None 5165568: status FINISHED
2026-07-15 16:42:56 INFO None 5165569: status FINISHED
2026-07-15 16:42:56 INFO None 5165570: status FINISHED
2026-07-15 16:42:56 INFO None 5165571: status FINISHED
2026-07-15 16:42:56 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:42:56 INFO None 5165573: status FINISHED
2026-07-15 16:42:56 INFO None 5165575: status FINISHED
2026-07-15 16:42:56 INFO None 5165577: status FINISHED
2026-07-15 16:42:56 INFO None 5165578: status FINISHED
2026-07-15 16:42:56 INFO None 5165580: status FINISHED
2026-07-15 16:42:56 INFO None 5165581: status FINISHED
2026-07-15 16:42:56 INFO None 5165582: status FINISHED
2026-07-15 16:42:56 INFO None 5165583: status RUNNING/PENDING
2026-07-15 16:42:56 INFO None 5165584: status RUNNING/PENDING
2026-07-15 16:42:56 INFO None 5165586: status RUNNING/PENDING
2026-07-15 16:42:56 INFO Jobs still running: ['5165572', '5165583', '5165584', '5165586']. Waiting...
2026-07-15 16:43:11 INFO None 5165568: status FINISHED
2026-07-15 16:43:11 INFO None 5165569: status FINISHED
2026-07-15 16:43:11 INFO None 5165570: status FINISHED
2026-07-15 16:43:11 INFO None 5165571: status FINISHED
2026-07-15 16:43:11 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:43:11 INFO None 5165573: status FINISHED
2026-07-15 16:43:11 INFO None 5165575: status FINISHED
2026-07-15 16:43:11 INFO None 5165577: status FINISHED
2026-07-15 16:43:11 INFO None 5165578: status FINISHED
2026-07-15 16:43:11 INFO None 5165580: status FINISHED
2026-07-15 16:43:11 INFO None 5165581: status FINISHED
2026-07-15 16:43:11 INFO None 5165582: status FINISHED
2026-07-15 16:43:11 INFO None 5165583: status FINISHED
2026-07-15 16:43:11 INFO None 5165584: status FINISHED
2026-07-15 16:43:11 INFO None 5165586: status FINISHED
2026-07-15 16:43:11 INFO Jobs still running: ['5165572']. Waiting...
2026-07-15 16:43:26 INFO None 5165568: status FINISHED
2026-07-15 16:43:26 INFO None 5165569: status FINISHED
2026-07-15 16:43:26 INFO None 5165570: status FINISHED
2026-07-15 16:43:26 INFO None 5165571: status FINISHED
2026-07-15 16:43:26 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:43:26 INFO None 5165573: status FINISHED
2026-07-15 16:43:26 INFO None 5165575: status FINISHED
2026-07-15 16:43:26 INFO None 5165577: status FINISHED
2026-07-15 16:43:26 INFO None 5165578: status FINISHED
2026-07-15 16:43:26 INFO None 5165580: status FINISHED
2026-07-15 16:43:26 INFO None 5165581: status FINISHED
2026-07-15 16:43:26 INFO None 5165582: status FINISHED
2026-07-15 16:43:26 INFO None 5165583: status FINISHED
2026-07-15 16:43:26 INFO None 5165584: status FINISHED
2026-07-15 16:43:27 INFO None 5165586: status FINISHED
2026-07-15 16:43:27 INFO Jobs still running: ['5165572']. Waiting...
2026-07-15 16:43:42 INFO None 5165568: status FINISHED
2026-07-15 16:43:42 INFO None 5165569: status FINISHED
2026-07-15 16:43:42 INFO None 5165570: status FINISHED
2026-07-15 16:43:42 INFO None 5165571: status FINISHED
2026-07-15 16:43:42 INFO None 5165572: status RUNNING/PENDING
2026-07-15 16:43:42 INFO None 5165573: status FINISHED
2026-07-15 16:43:42 INFO None 5165575: status FINISHED
2026-07-15 16:43:42 INFO None 5165577: status FINISHED
2026-07-15 16:43:42 INFO None 5165578: status FINISHED
2026-07-15 16:43:42 INFO None 5165580: status FINISHED
2026-07-15 16:43:42 INFO None 5165581: status FINISHED
2026-07-15 16:43:42 INFO None 5165582: status FINISHED
2026-07-15 16:43:42 INFO None 5165583: status FINISHED
2026-07-15 16:43:42 INFO None 5165584: status FINISHED
2026-07-15 16:43:42 INFO None 5165586: status FINISHED
2026-07-15 16:43:42 INFO Jobs still running: ['5165572']. Waiting...
2026-07-15 16:43:57 INFO None 5165568: status FINISHED
2026-07-15 16:43:57 INFO None 5165569: status FINISHED
2026-07-15 16:43:57 INFO None 5165570: status FINISHED
2026-07-15 16:43:57 INFO None 5165571: status FINISHED
2026-07-15 16:43:57 INFO None 5165572: status FINISHED
2026-07-15 16:43:57 INFO None 5165573: status FINISHED
2026-07-15 16:43:57 INFO None 5165575: status FINISHED
2026-07-15 16:43:57 INFO None 5165577: status FINISHED
2026-07-15 16:43:57 INFO None 5165578: status FINISHED
2026-07-15 16:43:57 INFO None 5165580: status FINISHED
2026-07-15 16:43:57 INFO None 5165581: status FINISHED
2026-07-15 16:43:57 INFO None 5165582: status FINISHED
2026-07-15 16:43:57 INFO None 5165583: status FINISHED
2026-07-15 16:43:57 INFO None 5165584: status FINISHED
2026-07-15 16:43:57 INFO None 5165586: status FINISHED
2026-07-15 16:43:57 INFO Jobs ['5165568', '5165569', '5165570', '5165571', '5165572', '5165573', '5165575', '5165577', '5165578', '5165580', '5165581', '5165582', '5165583', '5165584', '5165586'] have finished
2026-07-15 16:43:57 INFO Checking restart files were created ...
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-15 16:43:57 INFO  Run_model() completed successfully.
2026-07-15 16:43:57 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:43:57 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-15 16:43:57 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 16:43:57 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-15 16:43:57 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:43:57 INFO ---------->>> Running process_satellite_data()
2026-07-15 16:43:57 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-15 16:43:57 INFO ---------->>> Running run_obs_converter()
2026-07-15 16:43:57 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-15 16:43:57 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-15 16:43:57 INFO ---------->>> Running DART
2026-07-15 16:43:57 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 16:43:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-15 16:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-15 16:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-15 16:43:58 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-15 16:43:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-15 16:43:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-15 16:43:59 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-15 16:44:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-15 16:44:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-15 16:44:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-15 16:44:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-15 16:44:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-15 16:44:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-15 16:44:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-15 16:44:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-15 16:44:03 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 16:44:03 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 16:44:03 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 16:44:03 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 16:44:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 16:44:03 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 16:44:17 INFO Found: []
2026-07-15 16:44:17 INFO No job id returned by command ./run_filter.bsh
2026-07-15 16:44:17 INFO No monitoring will be performed
2026-07-15 16:44:17 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-15 16:44:17 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:17 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020613'
2026-07-15 16:44:18 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 16:44:21 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 16:44:21 INFO run_dart() is DONE.
2026-07-15 16:44:21 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 16:44:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-15 16:44:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-15 16:44:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-15 16:44:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-15 16:44:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-15 16:44:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-15 16:44:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:44:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-15 16:45:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-15 16:45:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-15 16:45:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-15 16:45:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-15 16:45:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-15 16:45:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-15 16:45:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-15 16:45:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-15 16:45:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:45:53 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 16:45:53 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:45:53 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:45:53 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-15 16:45:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:45:54 INFO Hourly dataset computed and listing created
2026-07-15 16:45:59 INFO Hourly dataset computed
2026-07-15 16:45:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:00 INFO Hourly dataset computed and listing created
2026-07-15 16:46:01 INFO Hourly dataset computed
2026-07-15 16:46:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:01 INFO Hourly dataset computed and listing created
2026-07-15 16:46:02 INFO Hourly dataset computed
2026-07-15 16:46:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:03 INFO Hourly dataset computed and listing created
2026-07-15 16:46:04 INFO Hourly dataset computed
2026-07-15 16:46:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:05 INFO Hourly dataset computed and listing created
2026-07-15 16:46:05 INFO Hourly dataset computed
2026-07-15 16:46:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:06 INFO Hourly dataset computed and listing created
2026-07-15 16:46:07 INFO Hourly dataset computed
2026-07-15 16:46:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:08 INFO Hourly dataset computed and listing created
2026-07-15 16:46:08 INFO Hourly dataset computed
2026-07-15 16:46:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:09 INFO Hourly dataset computed and listing created
2026-07-15 16:46:10 INFO Hourly dataset computed
2026-07-15 16:46:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:11 INFO Hourly dataset computed and listing created
2026-07-15 16:46:12 INFO Hourly dataset computed
2026-07-15 16:46:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:12 INFO Hourly dataset computed and listing created
2026-07-15 16:46:13 INFO Hourly dataset computed
2026-07-15 16:46:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:14 INFO Hourly dataset computed and listing created
2026-07-15 16:46:15 INFO Hourly dataset computed
2026-07-15 16:46:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:15 INFO Hourly dataset computed and listing created
2026-07-15 16:46:16 INFO Hourly dataset computed
2026-07-15 16:46:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:17 INFO Hourly dataset computed and listing created
2026-07-15 16:46:18 INFO Hourly dataset computed
2026-07-15 16:46:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:19 INFO Hourly dataset computed and listing created
2026-07-15 16:46:19 INFO Hourly dataset computed
2026-07-15 16:46:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:46:20 INFO Hourly dataset computed and listing created
2026-07-15 16:46:21 INFO Hourly dataset computed
2026-07-15 16:46:21 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-15 16:46:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 16:46:21 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-15 16:46:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 16:46:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:21 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 16:46:21 INFO Queuing job for member 1...
2026-07-15 16:46:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:21 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 16:46:22 INFO Found: ['5165643']
2026-07-15 16:46:27 INFO [TGCC-IRENE] Submitted job with ID:['5165643']
2026-07-15 16:46:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 16:46:27 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-15 16:46:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 16:46:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:27 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 16:46:27 INFO Queuing job for member 2...
2026-07-15 16:46:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:27 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 16:46:28 INFO Found: ['5165645']
2026-07-15 16:46:33 INFO [TGCC-IRENE] Submitted job with ID:['5165645']
2026-07-15 16:46:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 16:46:33 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-15 16:46:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 16:46:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:33 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 16:46:33 INFO Queuing job for member 3...
2026-07-15 16:46:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:33 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 16:46:33 INFO Found: ['5165646']
2026-07-15 16:46:38 INFO [TGCC-IRENE] Submitted job with ID:['5165646']
2026-07-15 16:46:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 16:46:38 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-15 16:46:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 16:46:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:38 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 16:46:38 INFO Queuing job for member 4...
2026-07-15 16:46:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:38 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 16:46:39 INFO Found: ['5165647']
2026-07-15 16:46:44 INFO [TGCC-IRENE] Submitted job with ID:['5165647']
2026-07-15 16:46:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 16:46:44 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-15 16:46:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 16:46:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:44 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 16:46:44 INFO Queuing job for member 5...
2026-07-15 16:46:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:44 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 16:46:45 INFO Found: ['5165648']
2026-07-15 16:46:50 INFO [TGCC-IRENE] Submitted job with ID:['5165648']
2026-07-15 16:46:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 16:46:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-15 16:46:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 16:46:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 16:46:50 INFO Queuing job for member 6...
2026-07-15 16:46:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 16:46:51 INFO Found: ['5165649']
2026-07-15 16:46:56 INFO [TGCC-IRENE] Submitted job with ID:['5165649']
2026-07-15 16:46:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:46:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 16:46:56 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-15 16:46:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 16:46:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:46:56 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 16:46:56 INFO Queuing job for member 7...
2026-07-15 16:46:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:46:56 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 16:46:56 INFO Found: ['5165650']
2026-07-15 16:47:01 INFO [TGCC-IRENE] Submitted job with ID:['5165650']
2026-07-15 16:47:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 16:47:01 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-15 16:47:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 16:47:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:01 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 16:47:01 INFO Queuing job for member 8...
2026-07-15 16:47:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:01 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 16:47:02 INFO Found: ['5165653']
2026-07-15 16:47:07 INFO [TGCC-IRENE] Submitted job with ID:['5165653']
2026-07-15 16:47:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 16:47:07 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-15 16:47:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 16:47:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:07 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 16:47:07 INFO Queuing job for member 9...
2026-07-15 16:47:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:07 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 16:47:08 INFO Found: ['5165655']
2026-07-15 16:47:13 INFO [TGCC-IRENE] Submitted job with ID:['5165655']
2026-07-15 16:47:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 16:47:13 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-15 16:47:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 16:47:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:13 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 16:47:13 INFO Queuing job for member 10...
2026-07-15 16:47:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:13 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 16:47:14 INFO Found: ['5165656']
2026-07-15 16:47:19 INFO [TGCC-IRENE] Submitted job with ID:['5165656']
2026-07-15 16:47:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 16:47:19 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-15 16:47:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 16:47:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:19 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 16:47:19 INFO Queuing job for member 11...
2026-07-15 16:47:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:19 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 16:47:20 INFO Found: ['5165658']
2026-07-15 16:47:25 INFO [TGCC-IRENE] Submitted job with ID:['5165658']
2026-07-15 16:47:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 16:47:25 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-15 16:47:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 16:47:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:25 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 16:47:25 INFO Queuing job for member 12...
2026-07-15 16:47:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:25 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 16:47:26 INFO Found: ['5165659']
2026-07-15 16:47:31 INFO [TGCC-IRENE] Submitted job with ID:['5165659']
2026-07-15 16:47:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 16:47:31 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-15 16:47:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 16:47:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:31 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 16:47:31 INFO Queuing job for member 13...
2026-07-15 16:47:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:31 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 16:47:32 INFO Found: ['5165660']
2026-07-15 16:47:37 INFO [TGCC-IRENE] Submitted job with ID:['5165660']
2026-07-15 16:47:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 16:47:37 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-15 16:47:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 16:47:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:37 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 16:47:37 INFO Queuing job for member 14...
2026-07-15 16:47:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:37 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 16:47:37 INFO Found: ['5165662']
2026-07-15 16:47:42 INFO [TGCC-IRENE] Submitted job with ID:['5165662']
2026-07-15 16:47:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 16:47:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 16:47:42 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-15 16:47:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 16:47:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 16:47:42 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 16:47:42 INFO Queuing job for member 15...
2026-07-15 16:47:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 16:47:42 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 16:47:43 INFO Found: ['5165663']
2026-07-15 16:47:48 INFO [TGCC-IRENE] Submitted job with ID:['5165663']
2026-07-15 16:47:48 INFO Checking job status ...
2026-07-15 16:47:48 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165650: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165653: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:47:48 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:47:48 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:48:04 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165650: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165653: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:48:04 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:48:04 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:48:19 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165650: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165653: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:48:19 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:48:19 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:48:34 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165650: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165653: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:48:34 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:48:35 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:48:35 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:48:50 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165650: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165653: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:48:50 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:48:50 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:49:06 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165650: status FINISHED
2026-07-15 16:49:06 INFO None 5165653: status FINISHED
2026-07-15 16:49:06 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:49:06 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:49:07 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:49:07 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:49:07 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:49:07 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:49:22 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165650: status FINISHED
2026-07-15 16:49:22 INFO None 5165653: status FINISHED
2026-07-15 16:49:22 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:49:22 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:49:22 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:49:37 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165650: status FINISHED
2026-07-15 16:49:37 INFO None 5165653: status FINISHED
2026-07-15 16:49:37 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165662: status RUNNING/PENDING
2026-07-15 16:49:37 INFO None 5165663: status RUNNING/PENDING
2026-07-15 16:49:37 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663']. Waiting...
2026-07-15 16:49:52 INFO None 5165643: status RUNNING/PENDING
2026-07-15 16:49:52 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:49:52 INFO None 5165646: status RUNNING/PENDING
2026-07-15 16:49:52 INFO None 5165647: status RUNNING/PENDING
2026-07-15 16:49:52 INFO None 5165648: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165650: status FINISHED
2026-07-15 16:49:54 INFO None 5165653: status FINISHED
2026-07-15 16:49:54 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:49:54 INFO None 5165662: status FINISHED
2026-07-15 16:49:54 INFO None 5165663: status FINISHED
2026-07-15 16:49:54 INFO Jobs still running: ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:50:09 INFO None 5165643: status FINISHED
2026-07-15 16:50:10 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165646: status FINISHED
2026-07-15 16:50:10 INFO None 5165647: status FINISHED
2026-07-15 16:50:10 INFO None 5165648: status FINISHED
2026-07-15 16:50:10 INFO None 5165649: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165650: status FINISHED
2026-07-15 16:50:10 INFO None 5165653: status FINISHED
2026-07-15 16:50:10 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:50:10 INFO None 5165662: status FINISHED
2026-07-15 16:50:10 INFO None 5165663: status FINISHED
2026-07-15 16:50:10 INFO Jobs still running: ['5165645', '5165649', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:50:25 INFO None 5165643: status FINISHED
2026-07-15 16:50:25 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165646: status FINISHED
2026-07-15 16:50:25 INFO None 5165647: status FINISHED
2026-07-15 16:50:25 INFO None 5165648: status FINISHED
2026-07-15 16:50:25 INFO None 5165649: status FINISHED
2026-07-15 16:50:25 INFO None 5165650: status FINISHED
2026-07-15 16:50:25 INFO None 5165653: status FINISHED
2026-07-15 16:50:25 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:50:25 INFO None 5165662: status FINISHED
2026-07-15 16:50:25 INFO None 5165663: status FINISHED
2026-07-15 16:50:25 INFO Jobs still running: ['5165645', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:50:40 INFO None 5165643: status FINISHED
2026-07-15 16:50:40 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165646: status FINISHED
2026-07-15 16:50:40 INFO None 5165647: status FINISHED
2026-07-15 16:50:40 INFO None 5165648: status FINISHED
2026-07-15 16:50:40 INFO None 5165649: status FINISHED
2026-07-15 16:50:40 INFO None 5165650: status FINISHED
2026-07-15 16:50:40 INFO None 5165653: status FINISHED
2026-07-15 16:50:40 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:50:40 INFO None 5165662: status FINISHED
2026-07-15 16:50:40 INFO None 5165663: status FINISHED
2026-07-15 16:50:40 INFO Jobs still running: ['5165645', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:50:57 INFO None 5165643: status FINISHED
2026-07-15 16:50:57 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165646: status FINISHED
2026-07-15 16:50:57 INFO None 5165647: status FINISHED
2026-07-15 16:50:57 INFO None 5165648: status FINISHED
2026-07-15 16:50:57 INFO None 5165649: status FINISHED
2026-07-15 16:50:57 INFO None 5165650: status FINISHED
2026-07-15 16:50:57 INFO None 5165653: status FINISHED
2026-07-15 16:50:57 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:50:57 INFO None 5165662: status FINISHED
2026-07-15 16:50:57 INFO None 5165663: status FINISHED
2026-07-15 16:50:57 INFO Jobs still running: ['5165645', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:51:12 INFO None 5165643: status FINISHED
2026-07-15 16:51:12 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:51:12 INFO None 5165646: status FINISHED
2026-07-15 16:51:13 INFO None 5165647: status FINISHED
2026-07-15 16:51:13 INFO None 5165648: status FINISHED
2026-07-15 16:51:13 INFO None 5165649: status FINISHED
2026-07-15 16:51:13 INFO None 5165650: status FINISHED
2026-07-15 16:51:13 INFO None 5165653: status FINISHED
2026-07-15 16:51:13 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:51:13 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:51:13 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:51:13 INFO None 5165659: status RUNNING/PENDING
2026-07-15 16:51:13 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:51:13 INFO None 5165662: status FINISHED
2026-07-15 16:51:13 INFO None 5165663: status FINISHED
2026-07-15 16:51:13 INFO Jobs still running: ['5165645', '5165655', '5165656', '5165658', '5165659', '5165660']. Waiting...
2026-07-15 16:51:28 INFO None 5165643: status FINISHED
2026-07-15 16:51:28 INFO None 5165645: status RUNNING/PENDING
2026-07-15 16:51:28 INFO None 5165646: status FINISHED
2026-07-15 16:51:28 INFO None 5165647: status FINISHED
2026-07-15 16:51:28 INFO None 5165648: status FINISHED
2026-07-15 16:51:28 INFO None 5165649: status FINISHED
2026-07-15 16:51:28 INFO None 5165650: status FINISHED
2026-07-15 16:51:28 INFO None 5165653: status FINISHED
2026-07-15 16:51:28 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:51:28 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:51:28 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:51:28 INFO None 5165659: status FINISHED
2026-07-15 16:51:28 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:51:28 INFO None 5165662: status FINISHED
2026-07-15 16:51:28 INFO None 5165663: status FINISHED
2026-07-15 16:51:28 INFO Jobs still running: ['5165645', '5165655', '5165656', '5165658', '5165660']. Waiting...
2026-07-15 16:51:43 INFO None 5165643: status FINISHED
2026-07-15 16:51:43 INFO None 5165645: status FINISHED
2026-07-15 16:51:43 INFO None 5165646: status FINISHED
2026-07-15 16:51:43 INFO None 5165647: status FINISHED
2026-07-15 16:51:43 INFO None 5165648: status FINISHED
2026-07-15 16:51:43 INFO None 5165649: status FINISHED
2026-07-15 16:51:43 INFO None 5165650: status FINISHED
2026-07-15 16:51:43 INFO None 5165653: status FINISHED
2026-07-15 16:51:43 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:51:43 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:51:43 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:51:43 INFO None 5165659: status FINISHED
2026-07-15 16:51:43 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:51:43 INFO None 5165662: status FINISHED
2026-07-15 16:51:43 INFO None 5165663: status FINISHED
2026-07-15 16:51:43 INFO Jobs still running: ['5165655', '5165656', '5165658', '5165660']. Waiting...
2026-07-15 16:51:58 INFO None 5165643: status FINISHED
2026-07-15 16:51:58 INFO None 5165645: status FINISHED
2026-07-15 16:51:58 INFO None 5165646: status FINISHED
2026-07-15 16:51:58 INFO None 5165647: status FINISHED
2026-07-15 16:51:58 INFO None 5165648: status FINISHED
2026-07-15 16:51:58 INFO None 5165649: status FINISHED
2026-07-15 16:51:58 INFO None 5165650: status FINISHED
2026-07-15 16:51:58 INFO None 5165653: status FINISHED
2026-07-15 16:51:58 INFO None 5165655: status RUNNING/PENDING
2026-07-15 16:51:58 INFO None 5165656: status RUNNING/PENDING
2026-07-15 16:51:59 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:51:59 INFO None 5165659: status FINISHED
2026-07-15 16:51:59 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:51:59 INFO None 5165662: status FINISHED
2026-07-15 16:51:59 INFO None 5165663: status FINISHED
2026-07-15 16:51:59 INFO Jobs still running: ['5165655', '5165656', '5165658', '5165660']. Waiting...
2026-07-15 16:52:14 INFO None 5165643: status FINISHED
2026-07-15 16:52:14 INFO None 5165645: status FINISHED
2026-07-15 16:52:14 INFO None 5165646: status FINISHED
2026-07-15 16:52:14 INFO None 5165647: status FINISHED
2026-07-15 16:52:14 INFO None 5165648: status FINISHED
2026-07-15 16:52:14 INFO None 5165649: status FINISHED
2026-07-15 16:52:14 INFO None 5165650: status FINISHED
2026-07-15 16:52:14 INFO None 5165653: status FINISHED
2026-07-15 16:52:14 INFO None 5165655: status FINISHED
2026-07-15 16:52:14 INFO None 5165656: status FINISHED
2026-07-15 16:52:14 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:52:14 INFO None 5165659: status FINISHED
2026-07-15 16:52:14 INFO None 5165660: status RUNNING/PENDING
2026-07-15 16:52:14 INFO None 5165662: status FINISHED
2026-07-15 16:52:14 INFO None 5165663: status FINISHED
2026-07-15 16:52:14 INFO Jobs still running: ['5165658', '5165660']. Waiting...
2026-07-15 16:52:29 INFO None 5165643: status FINISHED
2026-07-15 16:52:29 INFO None 5165645: status FINISHED
2026-07-15 16:52:29 INFO None 5165646: status FINISHED
2026-07-15 16:52:29 INFO None 5165647: status FINISHED
2026-07-15 16:52:29 INFO None 5165648: status FINISHED
2026-07-15 16:52:29 INFO None 5165649: status FINISHED
2026-07-15 16:52:29 INFO None 5165650: status FINISHED
2026-07-15 16:52:29 INFO None 5165653: status FINISHED
2026-07-15 16:52:29 INFO None 5165655: status FINISHED
2026-07-15 16:52:29 INFO None 5165656: status FINISHED
2026-07-15 16:52:29 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:52:29 INFO None 5165659: status FINISHED
2026-07-15 16:52:29 INFO None 5165660: status FINISHED
2026-07-15 16:52:29 INFO None 5165662: status FINISHED
2026-07-15 16:52:29 INFO None 5165663: status FINISHED
2026-07-15 16:52:29 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:52:46 INFO None 5165643: status FINISHED
2026-07-15 16:52:46 INFO None 5165645: status FINISHED
2026-07-15 16:52:46 INFO None 5165646: status FINISHED
2026-07-15 16:52:46 INFO None 5165647: status FINISHED
2026-07-15 16:52:46 INFO None 5165648: status FINISHED
2026-07-15 16:52:46 INFO None 5165649: status FINISHED
2026-07-15 16:52:46 INFO None 5165650: status FINISHED
2026-07-15 16:52:46 INFO None 5165653: status FINISHED
2026-07-15 16:52:46 INFO None 5165655: status FINISHED
2026-07-15 16:52:46 INFO None 5165656: status FINISHED
2026-07-15 16:52:46 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:52:46 INFO None 5165659: status FINISHED
2026-07-15 16:52:46 INFO None 5165660: status FINISHED
2026-07-15 16:52:46 INFO None 5165662: status FINISHED
2026-07-15 16:52:46 INFO None 5165663: status FINISHED
2026-07-15 16:52:46 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:53:01 INFO None 5165643: status FINISHED
2026-07-15 16:53:01 INFO None 5165645: status FINISHED
2026-07-15 16:53:01 INFO None 5165646: status FINISHED
2026-07-15 16:53:01 INFO None 5165647: status FINISHED
2026-07-15 16:53:01 INFO None 5165648: status FINISHED
2026-07-15 16:53:01 INFO None 5165649: status FINISHED
2026-07-15 16:53:01 INFO None 5165650: status FINISHED
2026-07-15 16:53:01 INFO None 5165653: status FINISHED
2026-07-15 16:53:01 INFO None 5165655: status FINISHED
2026-07-15 16:53:01 INFO None 5165656: status FINISHED
2026-07-15 16:53:01 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:53:01 INFO None 5165659: status FINISHED
2026-07-15 16:53:01 INFO None 5165660: status FINISHED
2026-07-15 16:53:01 INFO None 5165662: status FINISHED
2026-07-15 16:53:01 INFO None 5165663: status FINISHED
2026-07-15 16:53:01 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:53:16 INFO None 5165643: status FINISHED
2026-07-15 16:53:17 INFO None 5165645: status FINISHED
2026-07-15 16:53:17 INFO None 5165646: status FINISHED
2026-07-15 16:53:17 INFO None 5165647: status FINISHED
2026-07-15 16:53:17 INFO None 5165648: status FINISHED
2026-07-15 16:53:17 INFO None 5165649: status FINISHED
2026-07-15 16:53:17 INFO None 5165650: status FINISHED
2026-07-15 16:53:17 INFO None 5165653: status FINISHED
2026-07-15 16:53:17 INFO None 5165655: status FINISHED
2026-07-15 16:53:17 INFO None 5165656: status FINISHED
2026-07-15 16:53:17 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:53:17 INFO None 5165659: status FINISHED
2026-07-15 16:53:17 INFO None 5165660: status FINISHED
2026-07-15 16:53:17 INFO None 5165662: status FINISHED
2026-07-15 16:53:17 INFO None 5165663: status FINISHED
2026-07-15 16:53:17 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:53:32 INFO None 5165643: status FINISHED
2026-07-15 16:53:32 INFO None 5165645: status FINISHED
2026-07-15 16:53:32 INFO None 5165646: status FINISHED
2026-07-15 16:53:32 INFO None 5165647: status FINISHED
2026-07-15 16:53:32 INFO None 5165648: status FINISHED
2026-07-15 16:53:32 INFO None 5165649: status FINISHED
2026-07-15 16:53:32 INFO None 5165650: status FINISHED
2026-07-15 16:53:32 INFO None 5165653: status FINISHED
2026-07-15 16:53:32 INFO None 5165655: status FINISHED
2026-07-15 16:53:32 INFO None 5165656: status FINISHED
2026-07-15 16:53:32 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:53:32 INFO None 5165659: status FINISHED
2026-07-15 16:53:32 INFO None 5165660: status FINISHED
2026-07-15 16:53:32 INFO None 5165662: status FINISHED
2026-07-15 16:53:32 INFO None 5165663: status FINISHED
2026-07-15 16:53:32 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:53:47 INFO None 5165643: status FINISHED
2026-07-15 16:53:47 INFO None 5165645: status FINISHED
2026-07-15 16:53:47 INFO None 5165646: status FINISHED
2026-07-15 16:53:47 INFO None 5165647: status FINISHED
2026-07-15 16:53:47 INFO None 5165648: status FINISHED
2026-07-15 16:53:47 INFO None 5165649: status FINISHED
2026-07-15 16:53:47 INFO None 5165650: status FINISHED
2026-07-15 16:53:47 INFO None 5165653: status FINISHED
2026-07-15 16:53:47 INFO None 5165655: status FINISHED
2026-07-15 16:53:47 INFO None 5165656: status FINISHED
2026-07-15 16:53:47 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:53:47 INFO None 5165659: status FINISHED
2026-07-15 16:53:47 INFO None 5165660: status FINISHED
2026-07-15 16:53:47 INFO None 5165662: status FINISHED
2026-07-15 16:53:47 INFO None 5165663: status FINISHED
2026-07-15 16:53:47 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:54:02 INFO None 5165643: status FINISHED
2026-07-15 16:54:02 INFO None 5165645: status FINISHED
2026-07-15 16:54:02 INFO None 5165646: status FINISHED
2026-07-15 16:54:02 INFO None 5165647: status FINISHED
2026-07-15 16:54:02 INFO None 5165648: status FINISHED
2026-07-15 16:54:02 INFO None 5165649: status FINISHED
2026-07-15 16:54:03 INFO None 5165650: status FINISHED
2026-07-15 16:54:03 INFO None 5165653: status FINISHED
2026-07-15 16:54:03 INFO None 5165655: status FINISHED
2026-07-15 16:54:03 INFO None 5165656: status FINISHED
2026-07-15 16:54:03 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:54:03 INFO None 5165659: status FINISHED
2026-07-15 16:54:03 INFO None 5165660: status FINISHED
2026-07-15 16:54:03 INFO None 5165662: status FINISHED
2026-07-15 16:54:03 INFO None 5165663: status FINISHED
2026-07-15 16:54:03 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:54:18 INFO None 5165643: status FINISHED
2026-07-15 16:54:18 INFO None 5165645: status FINISHED
2026-07-15 16:54:18 INFO None 5165646: status FINISHED
2026-07-15 16:54:18 INFO None 5165647: status FINISHED
2026-07-15 16:54:18 INFO None 5165648: status FINISHED
2026-07-15 16:54:18 INFO None 5165649: status FINISHED
2026-07-15 16:54:18 INFO None 5165650: status FINISHED
2026-07-15 16:54:18 INFO None 5165653: status FINISHED
2026-07-15 16:54:18 INFO None 5165655: status FINISHED
2026-07-15 16:54:18 INFO None 5165656: status FINISHED
2026-07-15 16:54:18 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:54:18 INFO None 5165659: status FINISHED
2026-07-15 16:54:18 INFO None 5165660: status FINISHED
2026-07-15 16:54:18 INFO None 5165662: status FINISHED
2026-07-15 16:54:18 INFO None 5165663: status FINISHED
2026-07-15 16:54:18 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:54:33 INFO None 5165643: status FINISHED
2026-07-15 16:54:33 INFO None 5165645: status FINISHED
2026-07-15 16:54:33 INFO None 5165646: status FINISHED
2026-07-15 16:54:33 INFO None 5165647: status FINISHED
2026-07-15 16:54:33 INFO None 5165648: status FINISHED
2026-07-15 16:54:33 INFO None 5165649: status FINISHED
2026-07-15 16:54:33 INFO None 5165650: status FINISHED
2026-07-15 16:54:33 INFO None 5165653: status FINISHED
2026-07-15 16:54:33 INFO None 5165655: status FINISHED
2026-07-15 16:54:33 INFO None 5165656: status FINISHED
2026-07-15 16:54:33 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:54:33 INFO None 5165659: status FINISHED
2026-07-15 16:54:33 INFO None 5165660: status FINISHED
2026-07-15 16:54:33 INFO None 5165662: status FINISHED
2026-07-15 16:54:33 INFO None 5165663: status FINISHED
2026-07-15 16:54:33 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:54:48 INFO None 5165643: status FINISHED
2026-07-15 16:54:48 INFO None 5165645: status FINISHED
2026-07-15 16:54:48 INFO None 5165646: status FINISHED
2026-07-15 16:54:48 INFO None 5165647: status FINISHED
2026-07-15 16:54:48 INFO None 5165648: status FINISHED
2026-07-15 16:54:48 INFO None 5165649: status FINISHED
2026-07-15 16:54:48 INFO None 5165650: status FINISHED
2026-07-15 16:54:48 INFO None 5165653: status FINISHED
2026-07-15 16:54:48 INFO None 5165655: status FINISHED
2026-07-15 16:54:48 INFO None 5165656: status FINISHED
2026-07-15 16:54:48 INFO None 5165658: status RUNNING/PENDING
2026-07-15 16:54:48 INFO None 5165659: status FINISHED
2026-07-15 16:54:48 INFO None 5165660: status FINISHED
2026-07-15 16:54:48 INFO None 5165662: status FINISHED
2026-07-15 16:54:49 INFO None 5165663: status FINISHED
2026-07-15 16:54:49 INFO Jobs still running: ['5165658']. Waiting...
2026-07-15 16:55:04 INFO None 5165643: status FINISHED
2026-07-15 16:55:04 INFO None 5165645: status FINISHED
2026-07-15 16:55:04 INFO None 5165646: status FINISHED
2026-07-15 16:55:04 INFO None 5165647: status FINISHED
2026-07-15 16:55:04 INFO None 5165648: status FINISHED
2026-07-15 16:55:04 INFO None 5165649: status FINISHED
2026-07-15 16:55:04 INFO None 5165650: status FINISHED
2026-07-15 16:55:04 INFO None 5165653: status FINISHED
2026-07-15 16:55:04 INFO None 5165655: status FINISHED
2026-07-15 16:55:04 INFO None 5165656: status FINISHED
2026-07-15 16:55:04 INFO None 5165658: status FINISHED
2026-07-15 16:55:04 INFO None 5165659: status FINISHED
2026-07-15 16:55:04 INFO None 5165660: status FINISHED
2026-07-15 16:55:04 INFO None 5165662: status FINISHED
2026-07-15 16:55:04 INFO None 5165663: status FINISHED
2026-07-15 16:55:04 INFO Jobs ['5165643', '5165645', '5165646', '5165647', '5165648', '5165649', '5165650', '5165653', '5165655', '5165656', '5165658', '5165659', '5165660', '5165662', '5165663'] have finished
2026-07-15 16:55:04 INFO Checking restart files were created ...
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-15 16:55:04 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-15 16:55:04 INFO  Run_model() completed successfully.
2026-07-15 16:55:04 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:55:04 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-15 16:55:04 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 16:55:04 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-15 16:55:04 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:55:04 INFO ---------->>> Running process_satellite_data()
2026-07-15 16:55:04 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-15 16:55:04 INFO ---------->>> Running run_obs_converter()
2026-07-15 16:55:04 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-15 16:55:04 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-15 16:55:04 INFO ---------->>> Running DART
2026-07-15 16:55:04 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 16:55:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-15 16:55:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-15 16:55:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-15 16:55:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-15 16:55:05 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-15 16:55:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-15 16:55:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-15 16:55:06 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-15 16:55:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-15 16:55:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-15 16:55:07 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-15 16:55:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-15 16:55:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-15 16:55:08 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-15 16:55:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-15 16:55:09 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 16:55:09 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 16:55:09 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 16:55:09 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 16:55:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 16:55:09 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 16:55:17 INFO Found: []
2026-07-15 16:55:17 INFO No job id returned by command ./run_filter.bsh
2026-07-15 16:55:17 INFO No monitoring will be performed
2026-07-15 16:55:17 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-15 16:55:17 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:17 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:17 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:17 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:17 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020614'
2026-07-15 16:55:18 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 16:55:18 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 16:55:18 INFO run_dart() is DONE.
2026-07-15 16:55:18 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 16:55:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-15 16:55:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-15 16:55:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-15 16:55:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-15 16:55:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-15 16:55:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-15 16:55:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-15 16:55:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-15 16:55:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-15 16:55:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-15 16:55:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:55:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-15 16:56:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:56:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-15 16:56:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:56:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-15 16:56:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:56:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-15 16:56:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:56:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-15 16:56:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 16:56:19 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 16:56:19 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:56:19 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 16:56:19 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-15 16:56:19 INFO Copying EMIS of next day ...
2026-07-15 16:56:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-15 16:56:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:56:22 INFO Hourly dataset computed and listing created
2026-07-15 16:56:42 INFO Hourly dataset computed
2026-07-15 16:56:42 INFO Copying EMIS of next day ...
2026-07-15 16:56:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-15 16:56:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:56:45 INFO Hourly dataset computed and listing created
2026-07-15 16:57:05 INFO Hourly dataset computed
2026-07-15 16:57:05 INFO Copying EMIS of next day ...
2026-07-15 16:57:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-15 16:57:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:57:08 INFO Hourly dataset computed and listing created
2026-07-15 16:57:27 INFO Hourly dataset computed
2026-07-15 16:57:27 INFO Copying EMIS of next day ...
2026-07-15 16:57:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-15 16:57:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:57:29 INFO Hourly dataset computed and listing created
2026-07-15 16:57:45 INFO Hourly dataset computed
2026-07-15 16:57:45 INFO Copying EMIS of next day ...
2026-07-15 16:57:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-15 16:57:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:57:47 INFO Hourly dataset computed and listing created
2026-07-15 16:58:07 INFO Hourly dataset computed
2026-07-15 16:58:07 INFO Copying EMIS of next day ...
2026-07-15 16:58:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-15 16:58:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:58:09 INFO Hourly dataset computed and listing created
2026-07-15 16:58:25 INFO Hourly dataset computed
2026-07-15 16:58:25 INFO Copying EMIS of next day ...
2026-07-15 16:58:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-15 16:58:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:58:27 INFO Hourly dataset computed and listing created
2026-07-15 16:58:45 INFO Hourly dataset computed
2026-07-15 16:58:45 INFO Copying EMIS of next day ...
2026-07-15 16:58:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-15 16:58:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:58:47 INFO Hourly dataset computed and listing created
2026-07-15 16:59:03 INFO Hourly dataset computed
2026-07-15 16:59:03 INFO Copying EMIS of next day ...
2026-07-15 16:59:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-15 16:59:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:59:06 INFO Hourly dataset computed and listing created
2026-07-15 16:59:24 INFO Hourly dataset computed
2026-07-15 16:59:24 INFO Copying EMIS of next day ...
2026-07-15 16:59:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-15 16:59:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:59:26 INFO Hourly dataset computed and listing created
2026-07-15 16:59:44 INFO Hourly dataset computed
2026-07-15 16:59:44 INFO Copying EMIS of next day ...
2026-07-15 16:59:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-15 16:59:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 16:59:46 INFO Hourly dataset computed and listing created
2026-07-15 17:00:04 INFO Hourly dataset computed
2026-07-15 17:00:04 INFO Copying EMIS of next day ...
2026-07-15 17:00:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-15 17:00:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:00:06 INFO Hourly dataset computed and listing created
2026-07-15 17:00:23 INFO Hourly dataset computed
2026-07-15 17:00:23 INFO Copying EMIS of next day ...
2026-07-15 17:00:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-15 17:00:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:00:25 INFO Hourly dataset computed and listing created
2026-07-15 17:00:39 INFO Hourly dataset computed
2026-07-15 17:00:39 INFO Copying EMIS of next day ...
2026-07-15 17:00:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-15 17:00:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:00:41 INFO Hourly dataset computed and listing created
2026-07-15 17:00:55 INFO Hourly dataset computed
2026-07-15 17:00:55 INFO Copying EMIS of next day ...
2026-07-15 17:00:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-15 17:00:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:00:57 INFO Hourly dataset computed and listing created
2026-07-15 17:01:15 INFO Hourly dataset computed
2026-07-15 17:01:15 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-15 17:01:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 17:01:15 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-15 17:01:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 17:01:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:15 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 17:01:15 INFO Queuing job for member 1...
2026-07-15 17:01:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:15 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 17:01:16 INFO Found: ['5165751']
2026-07-15 17:01:21 INFO [TGCC-IRENE] Submitted job with ID:['5165751']
2026-07-15 17:01:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 17:01:21 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-15 17:01:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 17:01:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:21 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 17:01:21 INFO Queuing job for member 2...
2026-07-15 17:01:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:21 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 17:01:24 INFO Found: ['5165752']
2026-07-15 17:01:29 INFO [TGCC-IRENE] Submitted job with ID:['5165752']
2026-07-15 17:01:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 17:01:29 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-15 17:01:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 17:01:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:29 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 17:01:29 INFO Queuing job for member 3...
2026-07-15 17:01:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:29 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 17:01:29 INFO Found: ['5165753']
2026-07-15 17:01:34 INFO [TGCC-IRENE] Submitted job with ID:['5165753']
2026-07-15 17:01:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 17:01:34 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-15 17:01:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 17:01:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:34 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 17:01:35 INFO Queuing job for member 4...
2026-07-15 17:01:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:35 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 17:01:35 INFO Found: ['5165754']
2026-07-15 17:01:40 INFO [TGCC-IRENE] Submitted job with ID:['5165754']
2026-07-15 17:01:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 17:01:40 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-15 17:01:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 17:01:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:40 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 17:01:40 INFO Queuing job for member 5...
2026-07-15 17:01:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:40 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 17:01:41 INFO Found: ['5165756']
2026-07-15 17:01:46 INFO [TGCC-IRENE] Submitted job with ID:['5165756']
2026-07-15 17:01:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 17:01:46 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-15 17:01:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 17:01:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:46 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 17:01:46 INFO Queuing job for member 6...
2026-07-15 17:01:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:46 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 17:01:47 INFO Found: ['5165757']
2026-07-15 17:01:52 INFO [TGCC-IRENE] Submitted job with ID:['5165757']
2026-07-15 17:01:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 17:01:52 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-15 17:01:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 17:01:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:52 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 17:01:52 INFO Queuing job for member 7...
2026-07-15 17:01:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:52 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 17:01:52 INFO Found: ['5165758']
2026-07-15 17:01:57 INFO [TGCC-IRENE] Submitted job with ID:['5165758']
2026-07-15 17:01:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:01:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 17:01:57 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-15 17:01:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 17:01:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:01:58 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 17:01:58 INFO Queuing job for member 8...
2026-07-15 17:01:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:01:58 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 17:01:58 INFO Found: ['5165759']
2026-07-15 17:02:03 INFO [TGCC-IRENE] Submitted job with ID:['5165759']
2026-07-15 17:02:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 17:02:03 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-15 17:02:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 17:02:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:03 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 17:02:03 INFO Queuing job for member 9...
2026-07-15 17:02:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:03 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 17:02:04 INFO Found: ['5165761']
2026-07-15 17:02:09 INFO [TGCC-IRENE] Submitted job with ID:['5165761']
2026-07-15 17:02:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 17:02:09 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-15 17:02:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 17:02:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:09 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 17:02:09 INFO Queuing job for member 10...
2026-07-15 17:02:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:09 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 17:02:11 INFO Found: ['5165762']
2026-07-15 17:02:16 INFO [TGCC-IRENE] Submitted job with ID:['5165762']
2026-07-15 17:02:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 17:02:16 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-15 17:02:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 17:02:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:16 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 17:02:16 INFO Queuing job for member 11...
2026-07-15 17:02:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:16 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 17:02:18 INFO Found: ['5165763']
2026-07-15 17:02:23 INFO [TGCC-IRENE] Submitted job with ID:['5165763']
2026-07-15 17:02:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 17:02:23 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-15 17:02:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 17:02:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:23 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 17:02:23 INFO Queuing job for member 12...
2026-07-15 17:02:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:23 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 17:02:24 INFO Found: ['5165764']
2026-07-15 17:02:29 INFO [TGCC-IRENE] Submitted job with ID:['5165764']
2026-07-15 17:02:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 17:02:29 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-15 17:02:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 17:02:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:29 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 17:02:29 INFO Queuing job for member 13...
2026-07-15 17:02:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:29 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 17:02:30 INFO Found: ['5165766']
2026-07-15 17:02:35 INFO [TGCC-IRENE] Submitted job with ID:['5165766']
2026-07-15 17:02:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 17:02:35 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-15 17:02:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 17:02:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:35 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 17:02:35 INFO Queuing job for member 14...
2026-07-15 17:02:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:35 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 17:02:36 INFO Found: ['5165767']
2026-07-15 17:02:41 INFO [TGCC-IRENE] Submitted job with ID:['5165767']
2026-07-15 17:02:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:02:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 17:02:41 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-15 17:02:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 17:02:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:02:41 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 17:02:41 INFO Queuing job for member 15...
2026-07-15 17:02:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:02:41 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 17:02:41 INFO Found: ['5165769']
2026-07-15 17:02:46 INFO [TGCC-IRENE] Submitted job with ID:['5165769']
2026-07-15 17:02:46 INFO Checking job status ...
2026-07-15 17:02:46 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:02:46 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:02:46 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:02:46 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:02:47 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:02:47 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:03:02 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:03:02 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:03:02 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:03:18 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:03:18 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:03:18 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:03:33 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:03:33 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:03:33 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:03:48 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:03:48 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:03:48 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:04:05 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:04:05 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:04:05 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:04:20 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:04:20 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:04:20 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:04:35 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:04:35 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:04:36 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:04:36 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:04:51 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:04:51 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:04:51 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:05:06 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:05:06 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:05:06 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:05:21 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:05:21 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:05:22 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:05:22 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:05:22 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:05:37 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:05:37 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:05:37 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:05:52 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:05:52 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:05:52 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:06:07 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:06:07 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:06:07 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:06:23 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:06:23 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:06:23 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:06:38 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:06:38 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:06:38 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:06:55 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:06:55 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:06:55 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:07:10 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:07:10 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:07:10 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:07:25 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:07:25 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:07:26 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:07:26 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:07:41 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:07:41 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:07:41 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:07:58 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:07:58 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:07:58 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:08:13 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:08:13 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:08:13 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:08:28 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:08:28 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:08:29 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:08:29 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:08:29 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:08:29 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:08:44 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:08:44 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:08:44 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:08:59 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:08:59 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:08:59 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:09:14 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:09:14 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:09:14 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:09:29 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:09:30 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:09:30 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:09:45 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:09:45 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:09:45 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:10:00 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:10:00 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:10:01 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:10:01 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:10:16 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:10:16 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:10:16 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:10:31 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:10:31 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:10:31 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:10:47 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:10:47 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:10:47 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:11:02 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:11:02 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:11:02 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:11:17 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:11:17 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:11:17 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:11:17 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:11:17 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:11:17 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:11:18 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:11:18 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:11:33 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165759: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:11:33 INFO None 5165769: status RUNNING/PENDING
2026-07-15 17:11:33 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769']. Waiting...
2026-07-15 17:11:48 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165758: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165759: status FINISHED
2026-07-15 17:11:48 INFO None 5165761: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:11:48 INFO None 5165769: status FINISHED
2026-07-15 17:11:48 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767']. Waiting...
2026-07-15 17:12:03 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165752: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165754: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:12:03 INFO None 5165758: status FINISHED
2026-07-15 17:12:03 INFO None 5165759: status FINISHED
2026-07-15 17:12:03 INFO None 5165761: status FINISHED
2026-07-15 17:12:03 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:12:04 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:12:04 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:12:04 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:12:04 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:12:04 INFO None 5165769: status FINISHED
2026-07-15 17:12:04 INFO Jobs still running: ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165762', '5165763', '5165764', '5165766', '5165767']. Waiting...
2026-07-15 17:12:19 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165752: status FINISHED
2026-07-15 17:12:19 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165754: status FINISHED
2026-07-15 17:12:19 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165758: status FINISHED
2026-07-15 17:12:19 INFO None 5165759: status FINISHED
2026-07-15 17:12:19 INFO None 5165761: status FINISHED
2026-07-15 17:12:19 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:12:19 INFO None 5165769: status FINISHED
2026-07-15 17:12:19 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165757', '5165762', '5165763', '5165764', '5165766', '5165767']. Waiting...
2026-07-15 17:12:34 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165752: status FINISHED
2026-07-15 17:12:34 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165754: status FINISHED
2026-07-15 17:12:34 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165758: status FINISHED
2026-07-15 17:12:34 INFO None 5165759: status FINISHED
2026-07-15 17:12:34 INFO None 5165761: status FINISHED
2026-07-15 17:12:34 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165763: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165764: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:12:34 INFO None 5165769: status FINISHED
2026-07-15 17:12:34 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165757', '5165762', '5165763', '5165764', '5165766', '5165767']. Waiting...
2026-07-15 17:12:50 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165752: status FINISHED
2026-07-15 17:12:50 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165754: status FINISHED
2026-07-15 17:12:50 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165758: status FINISHED
2026-07-15 17:12:50 INFO None 5165759: status FINISHED
2026-07-15 17:12:50 INFO None 5165761: status FINISHED
2026-07-15 17:12:50 INFO None 5165762: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165763: status FINISHED
2026-07-15 17:12:50 INFO None 5165764: status FINISHED
2026-07-15 17:12:50 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:12:50 INFO None 5165769: status FINISHED
2026-07-15 17:12:50 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165757', '5165762', '5165766', '5165767']. Waiting...
2026-07-15 17:13:05 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165752: status FINISHED
2026-07-15 17:13:05 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165754: status FINISHED
2026-07-15 17:13:05 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165758: status FINISHED
2026-07-15 17:13:05 INFO None 5165759: status FINISHED
2026-07-15 17:13:05 INFO None 5165761: status FINISHED
2026-07-15 17:13:05 INFO None 5165762: status FINISHED
2026-07-15 17:13:05 INFO None 5165763: status FINISHED
2026-07-15 17:13:05 INFO None 5165764: status FINISHED
2026-07-15 17:13:05 INFO None 5165766: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:13:05 INFO None 5165769: status FINISHED
2026-07-15 17:13:05 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165757', '5165766', '5165767']. Waiting...
2026-07-15 17:13:20 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:13:20 INFO None 5165752: status FINISHED
2026-07-15 17:13:20 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:13:20 INFO None 5165754: status FINISHED
2026-07-15 17:13:21 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:13:21 INFO None 5165757: status RUNNING/PENDING
2026-07-15 17:13:21 INFO None 5165758: status FINISHED
2026-07-15 17:13:21 INFO None 5165759: status FINISHED
2026-07-15 17:13:21 INFO None 5165761: status FINISHED
2026-07-15 17:13:21 INFO None 5165762: status FINISHED
2026-07-15 17:13:21 INFO None 5165763: status FINISHED
2026-07-15 17:13:21 INFO None 5165764: status FINISHED
2026-07-15 17:13:21 INFO None 5165766: status FINISHED
2026-07-15 17:13:21 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:13:21 INFO None 5165769: status FINISHED
2026-07-15 17:13:21 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165757', '5165767']. Waiting...
2026-07-15 17:13:37 INFO None 5165751: status RUNNING/PENDING
2026-07-15 17:13:37 INFO None 5165752: status FINISHED
2026-07-15 17:13:37 INFO None 5165753: status RUNNING/PENDING
2026-07-15 17:13:37 INFO None 5165754: status FINISHED
2026-07-15 17:13:37 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:13:37 INFO None 5165757: status FINISHED
2026-07-15 17:13:37 INFO None 5165758: status FINISHED
2026-07-15 17:13:37 INFO None 5165759: status FINISHED
2026-07-15 17:13:37 INFO None 5165761: status FINISHED
2026-07-15 17:13:37 INFO None 5165762: status FINISHED
2026-07-15 17:13:37 INFO None 5165763: status FINISHED
2026-07-15 17:13:37 INFO None 5165764: status FINISHED
2026-07-15 17:13:37 INFO None 5165766: status FINISHED
2026-07-15 17:13:37 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:13:37 INFO None 5165769: status FINISHED
2026-07-15 17:13:37 INFO Jobs still running: ['5165751', '5165753', '5165756', '5165767']. Waiting...
2026-07-15 17:13:52 INFO None 5165751: status FINISHED
2026-07-15 17:13:52 INFO None 5165752: status FINISHED
2026-07-15 17:13:52 INFO None 5165753: status FINISHED
2026-07-15 17:13:52 INFO None 5165754: status FINISHED
2026-07-15 17:13:52 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:13:52 INFO None 5165757: status FINISHED
2026-07-15 17:13:52 INFO None 5165758: status FINISHED
2026-07-15 17:13:52 INFO None 5165759: status FINISHED
2026-07-15 17:13:52 INFO None 5165761: status FINISHED
2026-07-15 17:13:52 INFO None 5165762: status FINISHED
2026-07-15 17:13:52 INFO None 5165763: status FINISHED
2026-07-15 17:13:52 INFO None 5165764: status FINISHED
2026-07-15 17:13:52 INFO None 5165766: status FINISHED
2026-07-15 17:13:52 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:13:52 INFO None 5165769: status FINISHED
2026-07-15 17:13:52 INFO Jobs still running: ['5165756', '5165767']. Waiting...
2026-07-15 17:14:07 INFO None 5165751: status FINISHED
2026-07-15 17:14:07 INFO None 5165752: status FINISHED
2026-07-15 17:14:07 INFO None 5165753: status FINISHED
2026-07-15 17:14:07 INFO None 5165754: status FINISHED
2026-07-15 17:14:08 INFO None 5165756: status RUNNING/PENDING
2026-07-15 17:14:08 INFO None 5165757: status FINISHED
2026-07-15 17:14:08 INFO None 5165758: status FINISHED
2026-07-15 17:14:08 INFO None 5165759: status FINISHED
2026-07-15 17:14:08 INFO None 5165761: status FINISHED
2026-07-15 17:14:08 INFO None 5165762: status FINISHED
2026-07-15 17:14:08 INFO None 5165763: status FINISHED
2026-07-15 17:14:08 INFO None 5165764: status FINISHED
2026-07-15 17:14:08 INFO None 5165766: status FINISHED
2026-07-15 17:14:08 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:14:08 INFO None 5165769: status FINISHED
2026-07-15 17:14:08 INFO Jobs still running: ['5165756', '5165767']. Waiting...
2026-07-15 17:14:23 INFO None 5165751: status FINISHED
2026-07-15 17:14:23 INFO None 5165752: status FINISHED
2026-07-15 17:14:23 INFO None 5165753: status FINISHED
2026-07-15 17:14:23 INFO None 5165754: status FINISHED
2026-07-15 17:14:23 INFO None 5165756: status FINISHED
2026-07-15 17:14:23 INFO None 5165757: status FINISHED
2026-07-15 17:14:23 INFO None 5165758: status FINISHED
2026-07-15 17:14:23 INFO None 5165759: status FINISHED
2026-07-15 17:14:23 INFO None 5165761: status FINISHED
2026-07-15 17:14:23 INFO None 5165762: status FINISHED
2026-07-15 17:14:23 INFO None 5165763: status FINISHED
2026-07-15 17:14:23 INFO None 5165764: status FINISHED
2026-07-15 17:14:23 INFO None 5165766: status FINISHED
2026-07-15 17:14:23 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:14:23 INFO None 5165769: status FINISHED
2026-07-15 17:14:23 INFO Jobs still running: ['5165767']. Waiting...
2026-07-15 17:14:38 INFO None 5165751: status FINISHED
2026-07-15 17:14:38 INFO None 5165752: status FINISHED
2026-07-15 17:14:38 INFO None 5165753: status FINISHED
2026-07-15 17:14:38 INFO None 5165754: status FINISHED
2026-07-15 17:14:38 INFO None 5165756: status FINISHED
2026-07-15 17:14:38 INFO None 5165757: status FINISHED
2026-07-15 17:14:38 INFO None 5165758: status FINISHED
2026-07-15 17:14:38 INFO None 5165759: status FINISHED
2026-07-15 17:14:40 INFO None 5165761: status FINISHED
2026-07-15 17:14:40 INFO None 5165762: status FINISHED
2026-07-15 17:14:40 INFO None 5165763: status FINISHED
2026-07-15 17:14:40 INFO None 5165764: status FINISHED
2026-07-15 17:14:40 INFO None 5165766: status FINISHED
2026-07-15 17:14:40 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:14:40 INFO None 5165769: status FINISHED
2026-07-15 17:14:40 INFO Jobs still running: ['5165767']. Waiting...
2026-07-15 17:14:55 INFO None 5165751: status FINISHED
2026-07-15 17:14:55 INFO None 5165752: status FINISHED
2026-07-15 17:14:55 INFO None 5165753: status FINISHED
2026-07-15 17:14:55 INFO None 5165754: status FINISHED
2026-07-15 17:14:55 INFO None 5165756: status FINISHED
2026-07-15 17:14:55 INFO None 5165757: status FINISHED
2026-07-15 17:14:55 INFO None 5165758: status FINISHED
2026-07-15 17:14:55 INFO None 5165759: status FINISHED
2026-07-15 17:14:55 INFO None 5165761: status FINISHED
2026-07-15 17:14:55 INFO None 5165762: status FINISHED
2026-07-15 17:14:56 INFO None 5165763: status FINISHED
2026-07-15 17:14:56 INFO None 5165764: status FINISHED
2026-07-15 17:14:56 INFO None 5165766: status FINISHED
2026-07-15 17:14:56 INFO None 5165767: status RUNNING/PENDING
2026-07-15 17:14:56 INFO None 5165769: status FINISHED
2026-07-15 17:14:56 INFO Jobs still running: ['5165767']. Waiting...
2026-07-15 17:15:11 INFO None 5165751: status FINISHED
2026-07-15 17:15:11 INFO None 5165752: status FINISHED
2026-07-15 17:15:11 INFO None 5165753: status FINISHED
2026-07-15 17:15:11 INFO None 5165754: status FINISHED
2026-07-15 17:15:11 INFO None 5165756: status FINISHED
2026-07-15 17:15:11 INFO None 5165757: status FINISHED
2026-07-15 17:15:11 INFO None 5165758: status FINISHED
2026-07-15 17:15:11 INFO None 5165759: status FINISHED
2026-07-15 17:15:11 INFO None 5165761: status FINISHED
2026-07-15 17:15:11 INFO None 5165762: status FINISHED
2026-07-15 17:15:11 INFO None 5165763: status FINISHED
2026-07-15 17:15:11 INFO None 5165764: status FINISHED
2026-07-15 17:15:11 INFO None 5165766: status FINISHED
2026-07-15 17:15:11 INFO None 5165767: status FINISHED
2026-07-15 17:15:11 INFO None 5165769: status FINISHED
2026-07-15 17:15:11 INFO Jobs ['5165751', '5165752', '5165753', '5165754', '5165756', '5165757', '5165758', '5165759', '5165761', '5165762', '5165763', '5165764', '5165766', '5165767', '5165769'] have finished
2026-07-15 17:15:11 INFO Checking restart files were created ...
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-15 17:15:11 INFO  Run_model() completed successfully.
2026-07-15 17:15:11 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:15:11 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-15 17:15:11 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 17:15:11 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-15 17:15:11 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:15:11 INFO ---------->>> Running process_satellite_data()
2026-07-15 17:15:11 INFO [DART] No satellite data found, skipping assimilation
2026-07-15 17:15:11 INFO after_assimilation() skipped
2026-07-15 17:15:11 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 17:15:11 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:15:11 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:15:11 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-15 17:15:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:12 INFO Hourly dataset computed and listing created
2026-07-15 17:15:17 INFO Hourly dataset computed
2026-07-15 17:15:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:18 INFO Hourly dataset computed and listing created
2026-07-15 17:15:18 INFO Hourly dataset computed
2026-07-15 17:15:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:19 INFO Hourly dataset computed and listing created
2026-07-15 17:15:20 INFO Hourly dataset computed
2026-07-15 17:15:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:21 INFO Hourly dataset computed and listing created
2026-07-15 17:15:21 INFO Hourly dataset computed
2026-07-15 17:15:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:22 INFO Hourly dataset computed and listing created
2026-07-15 17:15:23 INFO Hourly dataset computed
2026-07-15 17:15:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:24 INFO Hourly dataset computed and listing created
2026-07-15 17:15:24 INFO Hourly dataset computed
2026-07-15 17:15:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:25 INFO Hourly dataset computed and listing created
2026-07-15 17:15:26 INFO Hourly dataset computed
2026-07-15 17:15:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:27 INFO Hourly dataset computed and listing created
2026-07-15 17:15:28 INFO Hourly dataset computed
2026-07-15 17:15:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:29 INFO Hourly dataset computed and listing created
2026-07-15 17:15:29 INFO Hourly dataset computed
2026-07-15 17:15:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:30 INFO Hourly dataset computed and listing created
2026-07-15 17:15:31 INFO Hourly dataset computed
2026-07-15 17:15:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:32 INFO Hourly dataset computed and listing created
2026-07-15 17:15:32 INFO Hourly dataset computed
2026-07-15 17:15:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:33 INFO Hourly dataset computed and listing created
2026-07-15 17:15:34 INFO Hourly dataset computed
2026-07-15 17:15:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:35 INFO Hourly dataset computed and listing created
2026-07-15 17:15:35 INFO Hourly dataset computed
2026-07-15 17:15:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:36 INFO Hourly dataset computed and listing created
2026-07-15 17:15:37 INFO Hourly dataset computed
2026-07-15 17:15:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:15:38 INFO Hourly dataset computed and listing created
2026-07-15 17:15:38 INFO Hourly dataset computed
2026-07-15 17:15:39 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-15 17:15:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:15:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 17:15:39 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-15 17:15:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 17:15:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:15:39 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 17:15:39 INFO Queuing job for member 1...
2026-07-15 17:15:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:15:39 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 17:15:41 INFO Found: ['5166020']
2026-07-15 17:15:46 INFO [TGCC-IRENE] Submitted job with ID:['5166020']
2026-07-15 17:15:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:15:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 17:15:46 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-15 17:15:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 17:15:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:15:46 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 17:15:46 INFO Queuing job for member 2...
2026-07-15 17:15:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:15:46 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 17:15:47 INFO Found: ['5166022']
2026-07-15 17:15:52 INFO [TGCC-IRENE] Submitted job with ID:['5166022']
2026-07-15 17:15:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:15:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 17:15:52 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-15 17:15:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 17:15:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:15:52 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 17:15:52 INFO Queuing job for member 3...
2026-07-15 17:15:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:15:52 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 17:15:53 INFO Found: ['5166023']
2026-07-15 17:15:58 INFO [TGCC-IRENE] Submitted job with ID:['5166023']
2026-07-15 17:15:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:15:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 17:15:58 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-15 17:15:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 17:15:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:15:58 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 17:15:58 INFO Queuing job for member 4...
2026-07-15 17:15:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:15:58 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 17:15:58 INFO Found: ['5166024']
2026-07-15 17:16:03 INFO [TGCC-IRENE] Submitted job with ID:['5166024']
2026-07-15 17:16:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 17:16:03 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-15 17:16:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 17:16:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:03 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 17:16:03 INFO Queuing job for member 5...
2026-07-15 17:16:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:03 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 17:16:04 INFO Found: ['5166026']
2026-07-15 17:16:09 INFO [TGCC-IRENE] Submitted job with ID:['5166026']
2026-07-15 17:16:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 17:16:09 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-15 17:16:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 17:16:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:09 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 17:16:09 INFO Queuing job for member 6...
2026-07-15 17:16:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:09 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 17:16:10 INFO Found: ['5166027']
2026-07-15 17:16:15 INFO [TGCC-IRENE] Submitted job with ID:['5166027']
2026-07-15 17:16:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 17:16:15 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-15 17:16:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 17:16:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:15 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 17:16:15 INFO Queuing job for member 7...
2026-07-15 17:16:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:15 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 17:16:15 INFO Found: ['5166028']
2026-07-15 17:16:20 INFO [TGCC-IRENE] Submitted job with ID:['5166028']
2026-07-15 17:16:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 17:16:20 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-15 17:16:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 17:16:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:20 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 17:16:20 INFO Queuing job for member 8...
2026-07-15 17:16:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:20 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 17:16:21 INFO Found: ['5166030']
2026-07-15 17:16:26 INFO [TGCC-IRENE] Submitted job with ID:['5166030']
2026-07-15 17:16:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 17:16:26 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-15 17:16:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 17:16:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:26 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 17:16:26 INFO Queuing job for member 9...
2026-07-15 17:16:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:26 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 17:16:27 INFO Found: ['5166031']
2026-07-15 17:16:32 INFO [TGCC-IRENE] Submitted job with ID:['5166031']
2026-07-15 17:16:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 17:16:32 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-15 17:16:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 17:16:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:32 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 17:16:32 INFO Queuing job for member 10...
2026-07-15 17:16:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:32 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 17:16:34 INFO Found: ['5166032']
2026-07-15 17:16:39 INFO [TGCC-IRENE] Submitted job with ID:['5166032']
2026-07-15 17:16:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 17:16:39 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-15 17:16:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 17:16:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:39 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 17:16:39 INFO Queuing job for member 11...
2026-07-15 17:16:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:39 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 17:16:42 INFO Found: ['5166033']
2026-07-15 17:16:47 INFO [TGCC-IRENE] Submitted job with ID:['5166033']
2026-07-15 17:16:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 17:16:47 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-15 17:16:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 17:16:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:47 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 17:16:47 INFO Queuing job for member 12...
2026-07-15 17:16:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:47 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 17:16:48 INFO Found: ['5166035']
2026-07-15 17:16:53 INFO [TGCC-IRENE] Submitted job with ID:['5166035']
2026-07-15 17:16:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 17:16:53 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-15 17:16:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 17:16:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:53 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 17:16:53 INFO Queuing job for member 13...
2026-07-15 17:16:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:53 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 17:16:53 INFO Found: ['5166036']
2026-07-15 17:16:58 INFO [TGCC-IRENE] Submitted job with ID:['5166036']
2026-07-15 17:16:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:16:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 17:16:58 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-15 17:16:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 17:16:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:16:58 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 17:16:58 INFO Queuing job for member 14...
2026-07-15 17:16:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:16:58 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 17:16:59 INFO Found: ['5166037']
2026-07-15 17:17:04 INFO [TGCC-IRENE] Submitted job with ID:['5166037']
2026-07-15 17:17:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:17:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 17:17:04 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-15 17:17:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 17:17:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:17:04 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 17:17:04 INFO Queuing job for member 15...
2026-07-15 17:17:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:17:04 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 17:17:05 INFO Found: ['5166040']
2026-07-15 17:17:10 INFO [TGCC-IRENE] Submitted job with ID:['5166040']
2026-07-15 17:17:10 INFO Checking job status ...
2026-07-15 17:17:10 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:17:10 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:17:10 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:17:25 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:17:25 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:17:26 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:17:26 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:17:41 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:17:41 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:17:41 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:17:56 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:17:56 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:17:56 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:18:11 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:18:11 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:18:12 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:18:12 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:18:28 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:18:28 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:18:28 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:18:43 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166023: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:18:43 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:18:44 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:18:44 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:18:44 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:18:47 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:18:47 INFO Jobs still running: ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:19:02 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:19:02 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:19:02 INFO None 5166023: status FINISHED
2026-07-15 17:19:02 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:19:02 INFO None 5166026: status RUNNING/PENDING
2026-07-15 17:19:02 INFO None 5166027: status RUNNING/PENDING
2026-07-15 17:19:05 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166035: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:19:06 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:19:06 INFO Jobs still running: ['5166020', '5166022', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:19:21 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166023: status FINISHED
2026-07-15 17:19:21 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166026: status FINISHED
2026-07-15 17:19:21 INFO None 5166027: status FINISHED
2026-07-15 17:19:21 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166035: status FINISHED
2026-07-15 17:19:21 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:19:21 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:19:21 INFO Jobs still running: ['5166020', '5166022', '5166024', '5166028', '5166030', '5166031', '5166032', '5166033', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:19:36 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166023: status FINISHED
2026-07-15 17:19:36 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166026: status FINISHED
2026-07-15 17:19:36 INFO None 5166027: status FINISHED
2026-07-15 17:19:36 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166035: status FINISHED
2026-07-15 17:19:36 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:19:36 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:19:36 INFO Jobs still running: ['5166020', '5166022', '5166024', '5166028', '5166030', '5166031', '5166032', '5166033', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:19:51 INFO None 5166020: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166022: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166023: status FINISHED
2026-07-15 17:19:51 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166026: status FINISHED
2026-07-15 17:19:51 INFO None 5166027: status FINISHED
2026-07-15 17:19:51 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166030: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166035: status FINISHED
2026-07-15 17:19:51 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:19:51 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:19:52 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:19:52 INFO Jobs still running: ['5166020', '5166022', '5166024', '5166028', '5166030', '5166031', '5166032', '5166033', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:20:07 INFO None 5166020: status FINISHED
2026-07-15 17:20:07 INFO None 5166022: status FINISHED
2026-07-15 17:20:07 INFO None 5166023: status FINISHED
2026-07-15 17:20:07 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166026: status FINISHED
2026-07-15 17:20:07 INFO None 5166027: status FINISHED
2026-07-15 17:20:07 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166030: status FINISHED
2026-07-15 17:20:07 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166035: status FINISHED
2026-07-15 17:20:07 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:20:07 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:20:07 INFO Jobs still running: ['5166024', '5166028', '5166031', '5166032', '5166033', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:20:22 INFO None 5166020: status FINISHED
2026-07-15 17:20:22 INFO None 5166022: status FINISHED
2026-07-15 17:20:22 INFO None 5166023: status FINISHED
2026-07-15 17:20:22 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166026: status FINISHED
2026-07-15 17:20:24 INFO None 5166027: status FINISHED
2026-07-15 17:20:24 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166030: status FINISHED
2026-07-15 17:20:24 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166035: status FINISHED
2026-07-15 17:20:24 INFO None 5166036: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166037: status RUNNING/PENDING
2026-07-15 17:20:24 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:20:24 INFO Jobs still running: ['5166024', '5166028', '5166031', '5166032', '5166033', '5166036', '5166037', '5166040']. Waiting...
2026-07-15 17:20:39 INFO None 5166020: status FINISHED
2026-07-15 17:20:39 INFO None 5166022: status FINISHED
2026-07-15 17:20:39 INFO None 5166023: status FINISHED
2026-07-15 17:20:39 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:20:39 INFO None 5166026: status FINISHED
2026-07-15 17:20:39 INFO None 5166027: status FINISHED
2026-07-15 17:20:39 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:20:39 INFO None 5166030: status FINISHED
2026-07-15 17:20:39 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:20:39 INFO None 5166032: status RUNNING/PENDING
2026-07-15 17:20:39 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:20:39 INFO None 5166035: status FINISHED
2026-07-15 17:20:39 INFO None 5166036: status FINISHED
2026-07-15 17:20:39 INFO None 5166037: status FINISHED
2026-07-15 17:20:39 INFO None 5166040: status RUNNING/PENDING
2026-07-15 17:20:39 INFO Jobs still running: ['5166024', '5166028', '5166031', '5166032', '5166033', '5166040']. Waiting...
2026-07-15 17:20:54 INFO None 5166020: status FINISHED
2026-07-15 17:20:54 INFO None 5166022: status FINISHED
2026-07-15 17:20:54 INFO None 5166023: status FINISHED
2026-07-15 17:20:54 INFO None 5166024: status RUNNING/PENDING
2026-07-15 17:20:55 INFO None 5166026: status FINISHED
2026-07-15 17:20:55 INFO None 5166027: status FINISHED
2026-07-15 17:20:55 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:20:55 INFO None 5166030: status FINISHED
2026-07-15 17:20:55 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:20:55 INFO None 5166032: status FINISHED
2026-07-15 17:20:55 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:20:55 INFO None 5166035: status FINISHED
2026-07-15 17:20:55 INFO None 5166036: status FINISHED
2026-07-15 17:20:55 INFO None 5166037: status FINISHED
2026-07-15 17:20:55 INFO None 5166040: status FINISHED
2026-07-15 17:20:55 INFO Jobs still running: ['5166024', '5166028', '5166031', '5166033']. Waiting...
2026-07-15 17:21:10 INFO None 5166020: status FINISHED
2026-07-15 17:21:10 INFO None 5166022: status FINISHED
2026-07-15 17:21:10 INFO None 5166023: status FINISHED
2026-07-15 17:21:10 INFO None 5166024: status FINISHED
2026-07-15 17:21:10 INFO None 5166026: status FINISHED
2026-07-15 17:21:10 INFO None 5166027: status FINISHED
2026-07-15 17:21:10 INFO None 5166028: status RUNNING/PENDING
2026-07-15 17:21:10 INFO None 5166030: status FINISHED
2026-07-15 17:21:10 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:21:10 INFO None 5166032: status FINISHED
2026-07-15 17:21:10 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:21:10 INFO None 5166035: status FINISHED
2026-07-15 17:21:10 INFO None 5166036: status FINISHED
2026-07-15 17:21:10 INFO None 5166037: status FINISHED
2026-07-15 17:21:10 INFO None 5166040: status FINISHED
2026-07-15 17:21:10 INFO Jobs still running: ['5166028', '5166031', '5166033']. Waiting...
2026-07-15 17:21:26 INFO None 5166020: status FINISHED
2026-07-15 17:21:26 INFO None 5166022: status FINISHED
2026-07-15 17:21:26 INFO None 5166023: status FINISHED
2026-07-15 17:21:26 INFO None 5166024: status FINISHED
2026-07-15 17:21:26 INFO None 5166026: status FINISHED
2026-07-15 17:21:26 INFO None 5166027: status FINISHED
2026-07-15 17:21:26 INFO None 5166028: status FINISHED
2026-07-15 17:21:26 INFO None 5166030: status FINISHED
2026-07-15 17:21:26 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:21:26 INFO None 5166032: status FINISHED
2026-07-15 17:21:26 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:21:26 INFO None 5166035: status FINISHED
2026-07-15 17:21:26 INFO None 5166036: status FINISHED
2026-07-15 17:21:26 INFO None 5166037: status FINISHED
2026-07-15 17:21:26 INFO None 5166040: status FINISHED
2026-07-15 17:21:26 INFO Jobs still running: ['5166031', '5166033']. Waiting...
2026-07-15 17:21:41 INFO None 5166020: status FINISHED
2026-07-15 17:21:41 INFO None 5166022: status FINISHED
2026-07-15 17:21:41 INFO None 5166023: status FINISHED
2026-07-15 17:21:41 INFO None 5166024: status FINISHED
2026-07-15 17:21:41 INFO None 5166026: status FINISHED
2026-07-15 17:21:41 INFO None 5166027: status FINISHED
2026-07-15 17:21:41 INFO None 5166028: status FINISHED
2026-07-15 17:21:41 INFO None 5166030: status FINISHED
2026-07-15 17:21:41 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:21:41 INFO None 5166032: status FINISHED
2026-07-15 17:21:41 INFO None 5166033: status RUNNING/PENDING
2026-07-15 17:21:41 INFO None 5166035: status FINISHED
2026-07-15 17:21:41 INFO None 5166036: status FINISHED
2026-07-15 17:21:41 INFO None 5166037: status FINISHED
2026-07-15 17:21:41 INFO None 5166040: status FINISHED
2026-07-15 17:21:41 INFO Jobs still running: ['5166031', '5166033']. Waiting...
2026-07-15 17:21:56 INFO None 5166020: status FINISHED
2026-07-15 17:21:56 INFO None 5166022: status FINISHED
2026-07-15 17:21:56 INFO None 5166023: status FINISHED
2026-07-15 17:21:56 INFO None 5166024: status FINISHED
2026-07-15 17:21:56 INFO None 5166026: status FINISHED
2026-07-15 17:21:56 INFO None 5166027: status FINISHED
2026-07-15 17:21:56 INFO None 5166028: status FINISHED
2026-07-15 17:21:56 INFO None 5166030: status FINISHED
2026-07-15 17:21:56 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:21:56 INFO None 5166032: status FINISHED
2026-07-15 17:21:56 INFO None 5166033: status FINISHED
2026-07-15 17:21:56 INFO None 5166035: status FINISHED
2026-07-15 17:21:57 INFO None 5166036: status FINISHED
2026-07-15 17:21:57 INFO None 5166037: status FINISHED
2026-07-15 17:21:57 INFO None 5166040: status FINISHED
2026-07-15 17:21:57 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:22:12 INFO None 5166020: status FINISHED
2026-07-15 17:22:12 INFO None 5166022: status FINISHED
2026-07-15 17:22:12 INFO None 5166023: status FINISHED
2026-07-15 17:22:12 INFO None 5166024: status FINISHED
2026-07-15 17:22:12 INFO None 5166026: status FINISHED
2026-07-15 17:22:12 INFO None 5166027: status FINISHED
2026-07-15 17:22:12 INFO None 5166028: status FINISHED
2026-07-15 17:22:12 INFO None 5166030: status FINISHED
2026-07-15 17:22:12 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:22:12 INFO None 5166032: status FINISHED
2026-07-15 17:22:12 INFO None 5166033: status FINISHED
2026-07-15 17:22:12 INFO None 5166035: status FINISHED
2026-07-15 17:22:12 INFO None 5166036: status FINISHED
2026-07-15 17:22:12 INFO None 5166037: status FINISHED
2026-07-15 17:22:12 INFO None 5166040: status FINISHED
2026-07-15 17:22:12 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:22:27 INFO None 5166020: status FINISHED
2026-07-15 17:22:27 INFO None 5166022: status FINISHED
2026-07-15 17:22:27 INFO None 5166023: status FINISHED
2026-07-15 17:22:27 INFO None 5166024: status FINISHED
2026-07-15 17:22:27 INFO None 5166026: status FINISHED
2026-07-15 17:22:27 INFO None 5166027: status FINISHED
2026-07-15 17:22:27 INFO None 5166028: status FINISHED
2026-07-15 17:22:28 INFO None 5166030: status FINISHED
2026-07-15 17:22:28 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:22:28 INFO None 5166032: status FINISHED
2026-07-15 17:22:28 INFO None 5166033: status FINISHED
2026-07-15 17:22:28 INFO None 5166035: status FINISHED
2026-07-15 17:22:28 INFO None 5166036: status FINISHED
2026-07-15 17:22:28 INFO None 5166037: status FINISHED
2026-07-15 17:22:28 INFO None 5166040: status FINISHED
2026-07-15 17:22:28 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:22:43 INFO None 5166020: status FINISHED
2026-07-15 17:22:43 INFO None 5166022: status FINISHED
2026-07-15 17:22:43 INFO None 5166023: status FINISHED
2026-07-15 17:22:43 INFO None 5166024: status FINISHED
2026-07-15 17:22:43 INFO None 5166026: status FINISHED
2026-07-15 17:22:43 INFO None 5166027: status FINISHED
2026-07-15 17:22:43 INFO None 5166028: status FINISHED
2026-07-15 17:22:43 INFO None 5166030: status FINISHED
2026-07-15 17:22:43 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:22:43 INFO None 5166032: status FINISHED
2026-07-15 17:22:43 INFO None 5166033: status FINISHED
2026-07-15 17:22:43 INFO None 5166035: status FINISHED
2026-07-15 17:22:43 INFO None 5166036: status FINISHED
2026-07-15 17:22:43 INFO None 5166037: status FINISHED
2026-07-15 17:22:43 INFO None 5166040: status FINISHED
2026-07-15 17:22:43 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:22:58 INFO None 5166020: status FINISHED
2026-07-15 17:22:58 INFO None 5166022: status FINISHED
2026-07-15 17:22:58 INFO None 5166023: status FINISHED
2026-07-15 17:22:58 INFO None 5166024: status FINISHED
2026-07-15 17:22:58 INFO None 5166026: status FINISHED
2026-07-15 17:22:58 INFO None 5166027: status FINISHED
2026-07-15 17:22:58 INFO None 5166028: status FINISHED
2026-07-15 17:22:58 INFO None 5166030: status FINISHED
2026-07-15 17:22:58 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:22:58 INFO None 5166032: status FINISHED
2026-07-15 17:22:58 INFO None 5166033: status FINISHED
2026-07-15 17:22:58 INFO None 5166035: status FINISHED
2026-07-15 17:22:58 INFO None 5166036: status FINISHED
2026-07-15 17:22:58 INFO None 5166037: status FINISHED
2026-07-15 17:22:58 INFO None 5166040: status FINISHED
2026-07-15 17:22:58 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:23:13 INFO None 5166020: status FINISHED
2026-07-15 17:23:13 INFO None 5166022: status FINISHED
2026-07-15 17:23:13 INFO None 5166023: status FINISHED
2026-07-15 17:23:13 INFO None 5166024: status FINISHED
2026-07-15 17:23:13 INFO None 5166026: status FINISHED
2026-07-15 17:23:13 INFO None 5166027: status FINISHED
2026-07-15 17:23:13 INFO None 5166028: status FINISHED
2026-07-15 17:23:13 INFO None 5166030: status FINISHED
2026-07-15 17:23:13 INFO None 5166031: status RUNNING/PENDING
2026-07-15 17:23:13 INFO None 5166032: status FINISHED
2026-07-15 17:23:13 INFO None 5166033: status FINISHED
2026-07-15 17:23:13 INFO None 5166035: status FINISHED
2026-07-15 17:23:13 INFO None 5166036: status FINISHED
2026-07-15 17:23:13 INFO None 5166037: status FINISHED
2026-07-15 17:23:13 INFO None 5166040: status FINISHED
2026-07-15 17:23:13 INFO Jobs still running: ['5166031']. Waiting...
2026-07-15 17:23:29 INFO None 5166020: status FINISHED
2026-07-15 17:23:29 INFO None 5166022: status FINISHED
2026-07-15 17:23:29 INFO None 5166023: status FINISHED
2026-07-15 17:23:29 INFO None 5166024: status FINISHED
2026-07-15 17:23:29 INFO None 5166026: status FINISHED
2026-07-15 17:23:29 INFO None 5166027: status FINISHED
2026-07-15 17:23:29 INFO None 5166028: status FINISHED
2026-07-15 17:23:29 INFO None 5166030: status FINISHED
2026-07-15 17:23:29 INFO None 5166031: status FINISHED
2026-07-15 17:23:29 INFO None 5166032: status FINISHED
2026-07-15 17:23:29 INFO None 5166033: status FINISHED
2026-07-15 17:23:29 INFO None 5166035: status FINISHED
2026-07-15 17:23:29 INFO None 5166036: status FINISHED
2026-07-15 17:23:29 INFO None 5166037: status FINISHED
2026-07-15 17:23:29 INFO None 5166040: status FINISHED
2026-07-15 17:23:29 INFO Jobs ['5166020', '5166022', '5166023', '5166024', '5166026', '5166027', '5166028', '5166030', '5166031', '5166032', '5166033', '5166035', '5166036', '5166037', '5166040'] have finished
2026-07-15 17:23:29 INFO Checking restart files were created ...
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-15 17:23:29 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-15 17:23:29 INFO  Run_model() completed successfully.
2026-07-15 17:23:29 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:23:29 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-15 17:23:29 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 17:23:29 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-15 17:23:29 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:23:29 INFO ---------->>> Running process_satellite_data()
2026-07-15 17:23:29 INFO [DART] No satellite data found, skipping assimilation
2026-07-15 17:23:29 INFO after_assimilation() skipped
2026-07-15 17:23:29 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 17:23:29 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:23:29 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:23:29 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-15 17:23:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:31 INFO Hourly dataset computed and listing created
2026-07-15 17:23:43 INFO Hourly dataset computed
2026-07-15 17:23:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:44 INFO Hourly dataset computed and listing created
2026-07-15 17:23:47 INFO Hourly dataset computed
2026-07-15 17:23:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:48 INFO Hourly dataset computed and listing created
2026-07-15 17:23:51 INFO Hourly dataset computed
2026-07-15 17:23:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:52 INFO Hourly dataset computed and listing created
2026-07-15 17:23:54 INFO Hourly dataset computed
2026-07-15 17:23:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:55 INFO Hourly dataset computed and listing created
2026-07-15 17:23:57 INFO Hourly dataset computed
2026-07-15 17:23:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:23:58 INFO Hourly dataset computed and listing created
2026-07-15 17:24:00 INFO Hourly dataset computed
2026-07-15 17:24:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:02 INFO Hourly dataset computed and listing created
2026-07-15 17:24:04 INFO Hourly dataset computed
2026-07-15 17:24:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:05 INFO Hourly dataset computed and listing created
2026-07-15 17:24:07 INFO Hourly dataset computed
2026-07-15 17:24:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:08 INFO Hourly dataset computed and listing created
2026-07-15 17:24:10 INFO Hourly dataset computed
2026-07-15 17:24:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:11 INFO Hourly dataset computed and listing created
2026-07-15 17:24:13 INFO Hourly dataset computed
2026-07-15 17:24:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:14 INFO Hourly dataset computed and listing created
2026-07-15 17:24:17 INFO Hourly dataset computed
2026-07-15 17:24:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:18 INFO Hourly dataset computed and listing created
2026-07-15 17:24:20 INFO Hourly dataset computed
2026-07-15 17:24:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:21 INFO Hourly dataset computed and listing created
2026-07-15 17:24:23 INFO Hourly dataset computed
2026-07-15 17:24:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:24 INFO Hourly dataset computed and listing created
2026-07-15 17:24:26 INFO Hourly dataset computed
2026-07-15 17:24:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:24:27 INFO Hourly dataset computed and listing created
2026-07-15 17:24:29 INFO Hourly dataset computed
2026-07-15 17:24:29 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-15 17:24:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 17:24:29 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-15 17:24:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 17:24:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:29 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 17:24:29 INFO Queuing job for member 1...
2026-07-15 17:24:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:29 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 17:24:30 INFO Found: ['5166133']
2026-07-15 17:24:35 INFO [TGCC-IRENE] Submitted job with ID:['5166133']
2026-07-15 17:24:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 17:24:35 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-15 17:24:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 17:24:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:35 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 17:24:35 INFO Queuing job for member 2...
2026-07-15 17:24:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:35 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 17:24:36 INFO Found: ['5166138']
2026-07-15 17:24:41 INFO [TGCC-IRENE] Submitted job with ID:['5166138']
2026-07-15 17:24:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 17:24:41 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-15 17:24:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 17:24:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:41 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 17:24:41 INFO Queuing job for member 3...
2026-07-15 17:24:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:41 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 17:24:42 INFO Found: ['5166139']
2026-07-15 17:24:47 INFO [TGCC-IRENE] Submitted job with ID:['5166139']
2026-07-15 17:24:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 17:24:47 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-15 17:24:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 17:24:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:47 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 17:24:47 INFO Queuing job for member 4...
2026-07-15 17:24:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:47 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 17:24:47 INFO Found: ['5166140']
2026-07-15 17:24:52 INFO [TGCC-IRENE] Submitted job with ID:['5166140']
2026-07-15 17:24:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 17:24:52 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-15 17:24:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 17:24:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:52 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 17:24:52 INFO Queuing job for member 5...
2026-07-15 17:24:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:52 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 17:24:53 INFO Found: ['5166141']
2026-07-15 17:24:58 INFO [TGCC-IRENE] Submitted job with ID:['5166141']
2026-07-15 17:24:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:24:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 17:24:58 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-15 17:24:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 17:24:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:24:58 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 17:24:58 INFO Queuing job for member 6...
2026-07-15 17:24:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:24:58 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 17:24:59 INFO Found: ['5166142']
2026-07-15 17:25:04 INFO [TGCC-IRENE] Submitted job with ID:['5166142']
2026-07-15 17:25:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 17:25:04 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-15 17:25:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 17:25:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:04 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 17:25:04 INFO Queuing job for member 7...
2026-07-15 17:25:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:04 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 17:25:05 INFO Found: ['5166144']
2026-07-15 17:25:10 INFO [TGCC-IRENE] Submitted job with ID:['5166144']
2026-07-15 17:25:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 17:25:10 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-15 17:25:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 17:25:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:10 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 17:25:10 INFO Queuing job for member 8...
2026-07-15 17:25:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:10 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 17:25:10 INFO Found: ['5166146']
2026-07-15 17:25:15 INFO [TGCC-IRENE] Submitted job with ID:['5166146']
2026-07-15 17:25:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 17:25:15 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-15 17:25:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 17:25:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:15 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 17:25:15 INFO Queuing job for member 9...
2026-07-15 17:25:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:15 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 17:25:16 INFO Found: ['5166147']
2026-07-15 17:25:21 INFO [TGCC-IRENE] Submitted job with ID:['5166147']
2026-07-15 17:25:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 17:25:21 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-15 17:25:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 17:25:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:21 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 17:25:21 INFO Queuing job for member 10...
2026-07-15 17:25:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:21 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 17:25:23 INFO Found: ['5166148']
2026-07-15 17:25:28 INFO [TGCC-IRENE] Submitted job with ID:['5166148']
2026-07-15 17:25:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 17:25:28 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-15 17:25:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 17:25:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:28 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 17:25:28 INFO Queuing job for member 11...
2026-07-15 17:25:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:28 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 17:25:29 INFO Found: ['5166149']
2026-07-15 17:25:34 INFO [TGCC-IRENE] Submitted job with ID:['5166149']
2026-07-15 17:25:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 17:25:34 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-15 17:25:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 17:25:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:34 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 17:25:34 INFO Queuing job for member 12...
2026-07-15 17:25:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:34 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 17:25:35 INFO Found: ['5166150']
2026-07-15 17:25:40 INFO [TGCC-IRENE] Submitted job with ID:['5166150']
2026-07-15 17:25:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 17:25:40 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-15 17:25:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 17:25:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:40 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 17:25:40 INFO Queuing job for member 13...
2026-07-15 17:25:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:40 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 17:25:40 INFO Found: ['5166151']
2026-07-15 17:25:45 INFO [TGCC-IRENE] Submitted job with ID:['5166151']
2026-07-15 17:25:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 17:25:45 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-15 17:25:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 17:25:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:45 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 17:25:45 INFO Queuing job for member 14...
2026-07-15 17:25:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:45 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 17:25:46 INFO Found: ['5166152']
2026-07-15 17:25:51 INFO [TGCC-IRENE] Submitted job with ID:['5166152']
2026-07-15 17:25:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:25:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 17:25:51 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-15 17:25:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 17:25:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:25:51 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 17:25:51 INFO Queuing job for member 15...
2026-07-15 17:25:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:25:51 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 17:25:52 INFO Found: ['5166153']
2026-07-15 17:25:57 INFO [TGCC-IRENE] Submitted job with ID:['5166153']
2026-07-15 17:25:57 INFO Checking job status ...
2026-07-15 17:25:57 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:25:57 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:25:57 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:26:12 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:26:12 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:26:12 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:26:12 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:26:12 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:26:12 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:26:13 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:26:13 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:26:28 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:26:28 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:26:28 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:26:43 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:26:43 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:26:43 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:26:58 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:26:58 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:26:59 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:26:59 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:27:15 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:27:15 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:27:16 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:27:16 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:27:16 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:27:31 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:27:31 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:27:31 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:27:46 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:27:46 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:27:46 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:28:01 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:28:01 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:28:01 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:28:17 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:28:17 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:28:19 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:28:19 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:28:19 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:28:19 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:28:19 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:28:34 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:28:34 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:28:35 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:28:35 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:28:50 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:28:50 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:28:50 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:29:05 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:29:05 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:29:05 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:29:21 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:29:21 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:29:21 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:29:36 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:29:36 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:29:36 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:29:52 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:29:52 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:29:52 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:30:08 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:30:08 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:30:08 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:30:23 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:30:23 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:30:23 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:30:38 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:30:39 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:30:39 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:30:54 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:30:54 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:30:54 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:31:09 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:31:09 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:31:10 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:31:10 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:31:25 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:31:25 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:31:25 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:31:40 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:31:40 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:31:40 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:31:55 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:31:55 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:31:56 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:31:56 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:32:11 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:32:11 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:32:11 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:32:11 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:32:11 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:32:11 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:32:12 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:32:12 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:32:27 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:32:27 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:32:27 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:32:42 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:32:42 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:32:42 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:32:57 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166138: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:32:57 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:32:58 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:32:58 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:32:58 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:32:58 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:32:58 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:32:58 INFO Jobs still running: ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:33:13 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166138: status FINISHED
2026-07-15 17:33:13 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:33:13 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:33:13 INFO Jobs still running: ['5166133', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:33:28 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166138: status FINISHED
2026-07-15 17:33:28 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:33:28 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:33:28 INFO Jobs still running: ['5166133', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:33:44 INFO None 5166133: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166138: status FINISHED
2026-07-15 17:33:44 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:33:44 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:33:44 INFO Jobs still running: ['5166133', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:33:59 INFO None 5166133: status FINISHED
2026-07-15 17:33:59 INFO None 5166138: status FINISHED
2026-07-15 17:33:59 INFO None 5166139: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166146: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:33:59 INFO None 5166153: status RUNNING/PENDING
2026-07-15 17:33:59 INFO Jobs still running: ['5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153']. Waiting...
2026-07-15 17:34:16 INFO None 5166133: status FINISHED
2026-07-15 17:34:16 INFO None 5166138: status FINISHED
2026-07-15 17:34:16 INFO None 5166139: status FINISHED
2026-07-15 17:34:16 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166146: status FINISHED
2026-07-15 17:34:16 INFO None 5166147: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:34:16 INFO None 5166153: status FINISHED
2026-07-15 17:34:16 INFO Jobs still running: ['5166140', '5166141', '5166142', '5166144', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152']. Waiting...
2026-07-15 17:34:31 INFO None 5166133: status FINISHED
2026-07-15 17:34:31 INFO None 5166138: status FINISHED
2026-07-15 17:34:31 INFO None 5166139: status FINISHED
2026-07-15 17:34:31 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166141: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166146: status FINISHED
2026-07-15 17:34:31 INFO None 5166147: status FINISHED
2026-07-15 17:34:31 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166150: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166151: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:34:31 INFO None 5166153: status FINISHED
2026-07-15 17:34:31 INFO Jobs still running: ['5166140', '5166141', '5166142', '5166144', '5166148', '5166149', '5166150', '5166151', '5166152']. Waiting...
2026-07-15 17:34:46 INFO None 5166133: status FINISHED
2026-07-15 17:34:46 INFO None 5166138: status FINISHED
2026-07-15 17:34:46 INFO None 5166139: status FINISHED
2026-07-15 17:34:46 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:34:46 INFO None 5166141: status FINISHED
2026-07-15 17:34:46 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:34:47 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:34:47 INFO None 5166146: status FINISHED
2026-07-15 17:34:47 INFO None 5166147: status FINISHED
2026-07-15 17:34:47 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:34:47 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:34:47 INFO None 5166150: status FINISHED
2026-07-15 17:34:47 INFO None 5166151: status FINISHED
2026-07-15 17:34:47 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:34:47 INFO None 5166153: status FINISHED
2026-07-15 17:34:47 INFO Jobs still running: ['5166140', '5166142', '5166144', '5166148', '5166149', '5166152']. Waiting...
2026-07-15 17:35:02 INFO None 5166133: status FINISHED
2026-07-15 17:35:02 INFO None 5166138: status FINISHED
2026-07-15 17:35:02 INFO None 5166139: status FINISHED
2026-07-15 17:35:02 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166141: status FINISHED
2026-07-15 17:35:02 INFO None 5166142: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166146: status FINISHED
2026-07-15 17:35:02 INFO None 5166147: status FINISHED
2026-07-15 17:35:02 INFO None 5166148: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166150: status FINISHED
2026-07-15 17:35:02 INFO None 5166151: status FINISHED
2026-07-15 17:35:02 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:35:02 INFO None 5166153: status FINISHED
2026-07-15 17:35:02 INFO Jobs still running: ['5166140', '5166142', '5166144', '5166148', '5166149', '5166152']. Waiting...
2026-07-15 17:35:18 INFO None 5166133: status FINISHED
2026-07-15 17:35:18 INFO None 5166138: status FINISHED
2026-07-15 17:35:19 INFO None 5166139: status FINISHED
2026-07-15 17:35:19 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:35:19 INFO None 5166141: status FINISHED
2026-07-15 17:35:19 INFO None 5166142: status FINISHED
2026-07-15 17:35:19 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:35:19 INFO None 5166146: status FINISHED
2026-07-15 17:35:19 INFO None 5166147: status FINISHED
2026-07-15 17:35:19 INFO None 5166148: status FINISHED
2026-07-15 17:35:19 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:35:19 INFO None 5166150: status FINISHED
2026-07-15 17:35:19 INFO None 5166151: status FINISHED
2026-07-15 17:35:19 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:35:19 INFO None 5166153: status FINISHED
2026-07-15 17:35:19 INFO Jobs still running: ['5166140', '5166144', '5166149', '5166152']. Waiting...
2026-07-15 17:35:34 INFO None 5166133: status FINISHED
2026-07-15 17:35:34 INFO None 5166138: status FINISHED
2026-07-15 17:35:34 INFO None 5166139: status FINISHED
2026-07-15 17:35:34 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:35:34 INFO None 5166141: status FINISHED
2026-07-15 17:35:34 INFO None 5166142: status FINISHED
2026-07-15 17:35:34 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:35:34 INFO None 5166146: status FINISHED
2026-07-15 17:35:34 INFO None 5166147: status FINISHED
2026-07-15 17:35:34 INFO None 5166148: status FINISHED
2026-07-15 17:35:34 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:35:34 INFO None 5166150: status FINISHED
2026-07-15 17:35:34 INFO None 5166151: status FINISHED
2026-07-15 17:35:34 INFO None 5166152: status RUNNING/PENDING
2026-07-15 17:35:34 INFO None 5166153: status FINISHED
2026-07-15 17:35:34 INFO Jobs still running: ['5166140', '5166144', '5166149', '5166152']. Waiting...
2026-07-15 17:35:49 INFO None 5166133: status FINISHED
2026-07-15 17:35:49 INFO None 5166138: status FINISHED
2026-07-15 17:35:49 INFO None 5166139: status FINISHED
2026-07-15 17:35:49 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:35:49 INFO None 5166141: status FINISHED
2026-07-15 17:35:49 INFO None 5166142: status FINISHED
2026-07-15 17:35:49 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:35:49 INFO None 5166146: status FINISHED
2026-07-15 17:35:49 INFO None 5166147: status FINISHED
2026-07-15 17:35:49 INFO None 5166148: status FINISHED
2026-07-15 17:35:49 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:35:49 INFO None 5166150: status FINISHED
2026-07-15 17:35:49 INFO None 5166151: status FINISHED
2026-07-15 17:35:49 INFO None 5166152: status FINISHED
2026-07-15 17:35:49 INFO None 5166153: status FINISHED
2026-07-15 17:35:49 INFO Jobs still running: ['5166140', '5166144', '5166149']. Waiting...
2026-07-15 17:36:04 INFO None 5166133: status FINISHED
2026-07-15 17:36:04 INFO None 5166138: status FINISHED
2026-07-15 17:36:04 INFO None 5166139: status FINISHED
2026-07-15 17:36:04 INFO None 5166140: status RUNNING/PENDING
2026-07-15 17:36:04 INFO None 5166141: status FINISHED
2026-07-15 17:36:04 INFO None 5166142: status FINISHED
2026-07-15 17:36:05 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:36:05 INFO None 5166146: status FINISHED
2026-07-15 17:36:05 INFO None 5166147: status FINISHED
2026-07-15 17:36:05 INFO None 5166148: status FINISHED
2026-07-15 17:36:05 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:36:05 INFO None 5166150: status FINISHED
2026-07-15 17:36:05 INFO None 5166151: status FINISHED
2026-07-15 17:36:05 INFO None 5166152: status FINISHED
2026-07-15 17:36:05 INFO None 5166153: status FINISHED
2026-07-15 17:36:05 INFO Jobs still running: ['5166140', '5166144', '5166149']. Waiting...
2026-07-15 17:36:20 INFO None 5166133: status FINISHED
2026-07-15 17:36:20 INFO None 5166138: status FINISHED
2026-07-15 17:36:20 INFO None 5166139: status FINISHED
2026-07-15 17:36:20 INFO None 5166140: status FINISHED
2026-07-15 17:36:20 INFO None 5166141: status FINISHED
2026-07-15 17:36:20 INFO None 5166142: status FINISHED
2026-07-15 17:36:20 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:36:20 INFO None 5166146: status FINISHED
2026-07-15 17:36:20 INFO None 5166147: status FINISHED
2026-07-15 17:36:20 INFO None 5166148: status FINISHED
2026-07-15 17:36:20 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:36:20 INFO None 5166150: status FINISHED
2026-07-15 17:36:20 INFO None 5166151: status FINISHED
2026-07-15 17:36:20 INFO None 5166152: status FINISHED
2026-07-15 17:36:20 INFO None 5166153: status FINISHED
2026-07-15 17:36:20 INFO Jobs still running: ['5166144', '5166149']. Waiting...
2026-07-15 17:36:35 INFO None 5166133: status FINISHED
2026-07-15 17:36:35 INFO None 5166138: status FINISHED
2026-07-15 17:36:35 INFO None 5166139: status FINISHED
2026-07-15 17:36:35 INFO None 5166140: status FINISHED
2026-07-15 17:36:35 INFO None 5166141: status FINISHED
2026-07-15 17:36:35 INFO None 5166142: status FINISHED
2026-07-15 17:36:35 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:36:35 INFO None 5166146: status FINISHED
2026-07-15 17:36:35 INFO None 5166147: status FINISHED
2026-07-15 17:36:35 INFO None 5166148: status FINISHED
2026-07-15 17:36:35 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:36:35 INFO None 5166150: status FINISHED
2026-07-15 17:36:36 INFO None 5166151: status FINISHED
2026-07-15 17:36:36 INFO None 5166152: status FINISHED
2026-07-15 17:36:36 INFO None 5166153: status FINISHED
2026-07-15 17:36:36 INFO Jobs still running: ['5166144', '5166149']. Waiting...
2026-07-15 17:36:51 INFO None 5166133: status FINISHED
2026-07-15 17:36:51 INFO None 5166138: status FINISHED
2026-07-15 17:36:51 INFO None 5166139: status FINISHED
2026-07-15 17:36:51 INFO None 5166140: status FINISHED
2026-07-15 17:36:51 INFO None 5166141: status FINISHED
2026-07-15 17:36:51 INFO None 5166142: status FINISHED
2026-07-15 17:36:51 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:36:51 INFO None 5166146: status FINISHED
2026-07-15 17:36:51 INFO None 5166147: status FINISHED
2026-07-15 17:36:51 INFO None 5166148: status FINISHED
2026-07-15 17:36:51 INFO None 5166149: status RUNNING/PENDING
2026-07-15 17:36:51 INFO None 5166150: status FINISHED
2026-07-15 17:36:51 INFO None 5166151: status FINISHED
2026-07-15 17:36:51 INFO None 5166152: status FINISHED
2026-07-15 17:36:51 INFO None 5166153: status FINISHED
2026-07-15 17:36:51 INFO Jobs still running: ['5166144', '5166149']. Waiting...
2026-07-15 17:37:06 INFO None 5166133: status FINISHED
2026-07-15 17:37:06 INFO None 5166138: status FINISHED
2026-07-15 17:37:06 INFO None 5166139: status FINISHED
2026-07-15 17:37:06 INFO None 5166140: status FINISHED
2026-07-15 17:37:06 INFO None 5166141: status FINISHED
2026-07-15 17:37:06 INFO None 5166142: status FINISHED
2026-07-15 17:37:06 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:37:06 INFO None 5166146: status FINISHED
2026-07-15 17:37:06 INFO None 5166147: status FINISHED
2026-07-15 17:37:06 INFO None 5166148: status FINISHED
2026-07-15 17:37:06 INFO None 5166149: status FINISHED
2026-07-15 17:37:06 INFO None 5166150: status FINISHED
2026-07-15 17:37:06 INFO None 5166151: status FINISHED
2026-07-15 17:37:06 INFO None 5166152: status FINISHED
2026-07-15 17:37:06 INFO None 5166153: status FINISHED
2026-07-15 17:37:06 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:37:23 INFO None 5166133: status FINISHED
2026-07-15 17:37:23 INFO None 5166138: status FINISHED
2026-07-15 17:37:23 INFO None 5166139: status FINISHED
2026-07-15 17:37:23 INFO None 5166140: status FINISHED
2026-07-15 17:37:23 INFO None 5166141: status FINISHED
2026-07-15 17:37:23 INFO None 5166142: status FINISHED
2026-07-15 17:37:23 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:37:23 INFO None 5166146: status FINISHED
2026-07-15 17:37:23 INFO None 5166147: status FINISHED
2026-07-15 17:37:23 INFO None 5166148: status FINISHED
2026-07-15 17:37:23 INFO None 5166149: status FINISHED
2026-07-15 17:37:23 INFO None 5166150: status FINISHED
2026-07-15 17:37:23 INFO None 5166151: status FINISHED
2026-07-15 17:37:23 INFO None 5166152: status FINISHED
2026-07-15 17:37:23 INFO None 5166153: status FINISHED
2026-07-15 17:37:23 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:37:38 INFO None 5166133: status FINISHED
2026-07-15 17:37:38 INFO None 5166138: status FINISHED
2026-07-15 17:37:38 INFO None 5166139: status FINISHED
2026-07-15 17:37:38 INFO None 5166140: status FINISHED
2026-07-15 17:37:38 INFO None 5166141: status FINISHED
2026-07-15 17:37:38 INFO None 5166142: status FINISHED
2026-07-15 17:37:38 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:37:38 INFO None 5166146: status FINISHED
2026-07-15 17:37:38 INFO None 5166147: status FINISHED
2026-07-15 17:37:38 INFO None 5166148: status FINISHED
2026-07-15 17:37:38 INFO None 5166149: status FINISHED
2026-07-15 17:37:38 INFO None 5166150: status FINISHED
2026-07-15 17:37:38 INFO None 5166151: status FINISHED
2026-07-15 17:37:38 INFO None 5166152: status FINISHED
2026-07-15 17:37:38 INFO None 5166153: status FINISHED
2026-07-15 17:37:38 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:37:53 INFO None 5166133: status FINISHED
2026-07-15 17:37:53 INFO None 5166138: status FINISHED
2026-07-15 17:37:53 INFO None 5166139: status FINISHED
2026-07-15 17:37:53 INFO None 5166140: status FINISHED
2026-07-15 17:37:53 INFO None 5166141: status FINISHED
2026-07-15 17:37:53 INFO None 5166142: status FINISHED
2026-07-15 17:37:53 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:37:53 INFO None 5166146: status FINISHED
2026-07-15 17:37:54 INFO None 5166147: status FINISHED
2026-07-15 17:37:54 INFO None 5166148: status FINISHED
2026-07-15 17:37:54 INFO None 5166149: status FINISHED
2026-07-15 17:37:54 INFO None 5166150: status FINISHED
2026-07-15 17:37:54 INFO None 5166151: status FINISHED
2026-07-15 17:37:54 INFO None 5166152: status FINISHED
2026-07-15 17:37:54 INFO None 5166153: status FINISHED
2026-07-15 17:37:54 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:38:09 INFO None 5166133: status FINISHED
2026-07-15 17:38:09 INFO None 5166138: status FINISHED
2026-07-15 17:38:09 INFO None 5166139: status FINISHED
2026-07-15 17:38:09 INFO None 5166140: status FINISHED
2026-07-15 17:38:09 INFO None 5166141: status FINISHED
2026-07-15 17:38:09 INFO None 5166142: status FINISHED
2026-07-15 17:38:09 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:38:09 INFO None 5166146: status FINISHED
2026-07-15 17:38:09 INFO None 5166147: status FINISHED
2026-07-15 17:38:09 INFO None 5166148: status FINISHED
2026-07-15 17:38:09 INFO None 5166149: status FINISHED
2026-07-15 17:38:09 INFO None 5166150: status FINISHED
2026-07-15 17:38:09 INFO None 5166151: status FINISHED
2026-07-15 17:38:09 INFO None 5166152: status FINISHED
2026-07-15 17:38:09 INFO None 5166153: status FINISHED
2026-07-15 17:38:09 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:38:24 INFO None 5166133: status FINISHED
2026-07-15 17:38:24 INFO None 5166138: status FINISHED
2026-07-15 17:38:24 INFO None 5166139: status FINISHED
2026-07-15 17:38:24 INFO None 5166140: status FINISHED
2026-07-15 17:38:24 INFO None 5166141: status FINISHED
2026-07-15 17:38:24 INFO None 5166142: status FINISHED
2026-07-15 17:38:24 INFO None 5166144: status RUNNING/PENDING
2026-07-15 17:38:24 INFO None 5166146: status FINISHED
2026-07-15 17:38:24 INFO None 5166147: status FINISHED
2026-07-15 17:38:24 INFO None 5166148: status FINISHED
2026-07-15 17:38:24 INFO None 5166149: status FINISHED
2026-07-15 17:38:24 INFO None 5166150: status FINISHED
2026-07-15 17:38:24 INFO None 5166151: status FINISHED
2026-07-15 17:38:24 INFO None 5166152: status FINISHED
2026-07-15 17:38:24 INFO None 5166153: status FINISHED
2026-07-15 17:38:24 INFO Jobs still running: ['5166144']. Waiting...
2026-07-15 17:38:39 INFO None 5166133: status FINISHED
2026-07-15 17:38:39 INFO None 5166138: status FINISHED
2026-07-15 17:38:39 INFO None 5166139: status FINISHED
2026-07-15 17:38:39 INFO None 5166140: status FINISHED
2026-07-15 17:38:39 INFO None 5166141: status FINISHED
2026-07-15 17:38:39 INFO None 5166142: status FINISHED
2026-07-15 17:38:39 INFO None 5166144: status FINISHED
2026-07-15 17:38:39 INFO None 5166146: status FINISHED
2026-07-15 17:38:39 INFO None 5166147: status FINISHED
2026-07-15 17:38:39 INFO None 5166148: status FINISHED
2026-07-15 17:38:39 INFO None 5166149: status FINISHED
2026-07-15 17:38:40 INFO None 5166150: status FINISHED
2026-07-15 17:38:40 INFO None 5166151: status FINISHED
2026-07-15 17:38:40 INFO None 5166152: status FINISHED
2026-07-15 17:38:40 INFO None 5166153: status FINISHED
2026-07-15 17:38:40 INFO Jobs ['5166133', '5166138', '5166139', '5166140', '5166141', '5166142', '5166144', '5166146', '5166147', '5166148', '5166149', '5166150', '5166151', '5166152', '5166153'] have finished
2026-07-15 17:38:40 INFO Checking restart files were created ...
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-15 17:38:40 INFO  Run_model() completed successfully.
2026-07-15 17:38:40 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:38:40 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-15 17:38:40 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 17:38:40 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-15 17:38:40 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:38:40 INFO ---------->>> Running process_satellite_data()
2026-07-15 17:38:40 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-15 17:38:40 INFO ---------->>> Running run_obs_converter()
2026-07-15 17:38:40 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-15 17:38:40 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-15 17:38:40 INFO ---------->>> Running DART
2026-07-15 17:38:40 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 17:38:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-15 17:38:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-15 17:38:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-15 17:38:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-15 17:38:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-15 17:38:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-15 17:38:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-15 17:38:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-15 17:38:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-15 17:38:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-15 17:38:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-15 17:38:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-15 17:38:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-15 17:38:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-15 17:38:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-15 17:38:45 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 17:38:46 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 17:38:46 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 17:38:46 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 17:38:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 17:38:46 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 17:38:52 INFO Found: []
2026-07-15 17:38:52 INFO No job id returned by command ./run_filter.bsh
2026-07-15 17:38:52 INFO No monitoring will be performed
2026-07-15 17:38:52 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-15 17:38:52 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020709'
2026-07-15 17:38:52 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 17:38:53 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 17:38:53 INFO run_dart() is DONE.
2026-07-15 17:38:53 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 17:38:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-15 17:39:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:39:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-15 17:39:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:39:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-15 17:39:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:39:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-15 17:39:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:39:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-15 17:40:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:40:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-15 17:40:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:40:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-15 17:40:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:40:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-15 17:40:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:40:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-15 17:40:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:40:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-15 17:41:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:41:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-15 17:41:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:41:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-15 17:41:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:41:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-15 17:41:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:41:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-15 17:42:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:42:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-15 17:42:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:42:23 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 17:42:23 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:42:23 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:42:23 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-15 17:42:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:24 INFO Hourly dataset computed and listing created
2026-07-15 17:42:30 INFO Hourly dataset computed
2026-07-15 17:42:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:31 INFO Hourly dataset computed and listing created
2026-07-15 17:42:32 INFO Hourly dataset computed
2026-07-15 17:42:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:33 INFO Hourly dataset computed and listing created
2026-07-15 17:42:33 INFO Hourly dataset computed
2026-07-15 17:42:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:34 INFO Hourly dataset computed and listing created
2026-07-15 17:42:35 INFO Hourly dataset computed
2026-07-15 17:42:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:36 INFO Hourly dataset computed and listing created
2026-07-15 17:42:37 INFO Hourly dataset computed
2026-07-15 17:42:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:38 INFO Hourly dataset computed and listing created
2026-07-15 17:42:39 INFO Hourly dataset computed
2026-07-15 17:42:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:40 INFO Hourly dataset computed and listing created
2026-07-15 17:42:41 INFO Hourly dataset computed
2026-07-15 17:42:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:42 INFO Hourly dataset computed and listing created
2026-07-15 17:42:43 INFO Hourly dataset computed
2026-07-15 17:42:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:44 INFO Hourly dataset computed and listing created
2026-07-15 17:42:45 INFO Hourly dataset computed
2026-07-15 17:42:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:46 INFO Hourly dataset computed and listing created
2026-07-15 17:42:46 INFO Hourly dataset computed
2026-07-15 17:42:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:47 INFO Hourly dataset computed and listing created
2026-07-15 17:42:48 INFO Hourly dataset computed
2026-07-15 17:42:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:49 INFO Hourly dataset computed and listing created
2026-07-15 17:42:50 INFO Hourly dataset computed
2026-07-15 17:42:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:51 INFO Hourly dataset computed and listing created
2026-07-15 17:42:52 INFO Hourly dataset computed
2026-07-15 17:42:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:53 INFO Hourly dataset computed and listing created
2026-07-15 17:42:54 INFO Hourly dataset computed
2026-07-15 17:42:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:42:55 INFO Hourly dataset computed and listing created
2026-07-15 17:42:56 INFO Hourly dataset computed
2026-07-15 17:42:56 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-15 17:42:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:42:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 17:42:56 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-15 17:42:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 17:42:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:42:56 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 17:42:56 INFO Queuing job for member 1...
2026-07-15 17:42:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:42:56 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 17:42:57 INFO Found: ['5166297']
2026-07-15 17:43:02 INFO [TGCC-IRENE] Submitted job with ID:['5166297']
2026-07-15 17:43:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 17:43:02 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-15 17:43:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 17:43:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:02 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 17:43:02 INFO Queuing job for member 2...
2026-07-15 17:43:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:02 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 17:43:02 INFO Found: ['5166300']
2026-07-15 17:43:07 INFO [TGCC-IRENE] Submitted job with ID:['5166300']
2026-07-15 17:43:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 17:43:07 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-15 17:43:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 17:43:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:07 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 17:43:07 INFO Queuing job for member 3...
2026-07-15 17:43:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:07 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 17:43:08 INFO Found: ['5166301']
2026-07-15 17:43:13 INFO [TGCC-IRENE] Submitted job with ID:['5166301']
2026-07-15 17:43:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 17:43:13 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-15 17:43:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 17:43:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:13 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 17:43:13 INFO Queuing job for member 4...
2026-07-15 17:43:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:13 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 17:43:14 INFO Found: ['5166302']
2026-07-15 17:43:19 INFO [TGCC-IRENE] Submitted job with ID:['5166302']
2026-07-15 17:43:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 17:43:19 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-15 17:43:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 17:43:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:19 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 17:43:19 INFO Queuing job for member 5...
2026-07-15 17:43:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:19 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 17:43:20 INFO Found: ['5166303']
2026-07-15 17:43:25 INFO [TGCC-IRENE] Submitted job with ID:['5166303']
2026-07-15 17:43:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 17:43:25 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-15 17:43:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 17:43:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:25 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 17:43:25 INFO Queuing job for member 6...
2026-07-15 17:43:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:25 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 17:43:27 INFO Found: ['5166304']
2026-07-15 17:43:32 INFO [TGCC-IRENE] Submitted job with ID:['5166304']
2026-07-15 17:43:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 17:43:32 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-15 17:43:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 17:43:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:32 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 17:43:32 INFO Queuing job for member 7...
2026-07-15 17:43:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:32 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 17:43:34 INFO Found: ['5166305']
2026-07-15 17:43:39 INFO [TGCC-IRENE] Submitted job with ID:['5166305']
2026-07-15 17:43:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 17:43:39 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-15 17:43:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 17:43:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:39 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 17:43:39 INFO Queuing job for member 8...
2026-07-15 17:43:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:39 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 17:43:40 INFO Found: ['5166306']
2026-07-15 17:43:45 INFO [TGCC-IRENE] Submitted job with ID:['5166306']
2026-07-15 17:43:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 17:43:45 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-15 17:43:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 17:43:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:45 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 17:43:45 INFO Queuing job for member 9...
2026-07-15 17:43:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:45 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 17:43:46 INFO Found: ['5166307']
2026-07-15 17:43:51 INFO [TGCC-IRENE] Submitted job with ID:['5166307']
2026-07-15 17:43:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 17:43:51 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-15 17:43:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 17:43:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:51 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 17:43:51 INFO Queuing job for member 10...
2026-07-15 17:43:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:51 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 17:43:52 INFO Found: ['5166308']
2026-07-15 17:43:57 INFO [TGCC-IRENE] Submitted job with ID:['5166308']
2026-07-15 17:43:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:43:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 17:43:57 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-15 17:43:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 17:43:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:43:57 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 17:43:57 INFO Queuing job for member 11...
2026-07-15 17:43:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:43:57 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 17:43:58 INFO Found: ['5166311']
2026-07-15 17:44:03 INFO [TGCC-IRENE] Submitted job with ID:['5166311']
2026-07-15 17:44:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:44:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 17:44:03 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-15 17:44:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 17:44:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:44:03 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 17:44:03 INFO Queuing job for member 12...
2026-07-15 17:44:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:44:03 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 17:44:03 INFO Found: ['5166314']
2026-07-15 17:44:08 INFO [TGCC-IRENE] Submitted job with ID:['5166314']
2026-07-15 17:44:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:44:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 17:44:08 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-15 17:44:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 17:44:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:44:08 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 17:44:08 INFO Queuing job for member 13...
2026-07-15 17:44:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:44:08 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 17:44:09 INFO Found: ['5166315']
2026-07-15 17:44:14 INFO [TGCC-IRENE] Submitted job with ID:['5166315']
2026-07-15 17:44:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:44:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 17:44:14 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-15 17:44:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 17:44:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:44:14 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 17:44:14 INFO Queuing job for member 14...
2026-07-15 17:44:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:44:14 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 17:44:15 INFO Found: ['5166318']
2026-07-15 17:44:20 INFO [TGCC-IRENE] Submitted job with ID:['5166318']
2026-07-15 17:44:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:44:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 17:44:20 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-15 17:44:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 17:44:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:44:20 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 17:44:20 INFO Queuing job for member 15...
2026-07-15 17:44:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:44:20 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 17:44:21 INFO Found: ['5166320']
2026-07-15 17:44:26 INFO [TGCC-IRENE] Submitted job with ID:['5166320']
2026-07-15 17:44:26 INFO Checking job status ...
2026-07-15 17:44:27 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:44:27 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:44:27 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:44:27 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:44:28 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:44:28 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:44:43 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:44:43 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:44:43 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:44:58 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:44:58 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:44:58 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:45:13 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:45:13 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:45:14 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:45:14 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:45:14 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:45:14 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:45:14 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:45:14 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:45:29 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:45:29 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:45:29 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:45:44 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:45:44 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:45:44 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:45:59 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:45:59 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:46:00 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:46:00 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:46:15 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166301: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:46:15 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:46:15 INFO Jobs still running: ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:46:30 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166301: status FINISHED
2026-07-15 17:46:30 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:46:30 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:46:30 INFO Jobs still running: ['5166297', '5166300', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:46:45 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:46:45 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:46:45 INFO None 5166301: status FINISHED
2026-07-15 17:46:46 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:46:46 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:46:46 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:46:46 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:46:48 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:46:48 INFO Jobs still running: ['5166297', '5166300', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:47:03 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166301: status FINISHED
2026-07-15 17:47:03 INFO None 5166302: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166303: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166306: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166307: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:47:03 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:47:03 INFO Jobs still running: ['5166297', '5166300', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:47:18 INFO None 5166297: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166300: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166301: status FINISHED
2026-07-15 17:47:18 INFO None 5166302: status FINISHED
2026-07-15 17:47:18 INFO None 5166303: status FINISHED
2026-07-15 17:47:18 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166306: status FINISHED
2026-07-15 17:47:18 INFO None 5166307: status FINISHED
2026-07-15 17:47:18 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:47:18 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:47:18 INFO Jobs still running: ['5166297', '5166300', '5166304', '5166305', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:47:33 INFO None 5166297: status FINISHED
2026-07-15 17:47:33 INFO None 5166300: status FINISHED
2026-07-15 17:47:35 INFO None 5166301: status FINISHED
2026-07-15 17:47:35 INFO None 5166302: status FINISHED
2026-07-15 17:47:35 INFO None 5166303: status FINISHED
2026-07-15 17:47:35 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:47:35 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:47:35 INFO None 5166306: status FINISHED
2026-07-15 17:47:36 INFO None 5166307: status FINISHED
2026-07-15 17:47:36 INFO None 5166308: status RUNNING/PENDING
2026-07-15 17:47:36 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:47:36 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:47:36 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:47:36 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:47:36 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:47:36 INFO Jobs still running: ['5166304', '5166305', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:47:51 INFO None 5166297: status FINISHED
2026-07-15 17:47:51 INFO None 5166300: status FINISHED
2026-07-15 17:47:51 INFO None 5166301: status FINISHED
2026-07-15 17:47:51 INFO None 5166302: status FINISHED
2026-07-15 17:47:51 INFO None 5166303: status FINISHED
2026-07-15 17:47:51 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166306: status FINISHED
2026-07-15 17:47:51 INFO None 5166307: status FINISHED
2026-07-15 17:47:51 INFO None 5166308: status FINISHED
2026-07-15 17:47:51 INFO None 5166311: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:47:51 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:47:51 INFO Jobs still running: ['5166304', '5166305', '5166311', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:48:06 INFO None 5166297: status FINISHED
2026-07-15 17:48:06 INFO None 5166300: status FINISHED
2026-07-15 17:48:06 INFO None 5166301: status FINISHED
2026-07-15 17:48:06 INFO None 5166302: status FINISHED
2026-07-15 17:48:06 INFO None 5166303: status FINISHED
2026-07-15 17:48:06 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:48:06 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:48:06 INFO None 5166306: status FINISHED
2026-07-15 17:48:06 INFO None 5166307: status FINISHED
2026-07-15 17:48:06 INFO None 5166308: status FINISHED
2026-07-15 17:48:06 INFO None 5166311: status FINISHED
2026-07-15 17:48:06 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:48:06 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:48:06 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:48:06 INFO None 5166320: status RUNNING/PENDING
2026-07-15 17:48:06 INFO Jobs still running: ['5166304', '5166305', '5166314', '5166315', '5166318', '5166320']. Waiting...
2026-07-15 17:48:21 INFO None 5166297: status FINISHED
2026-07-15 17:48:21 INFO None 5166300: status FINISHED
2026-07-15 17:48:21 INFO None 5166301: status FINISHED
2026-07-15 17:48:21 INFO None 5166302: status FINISHED
2026-07-15 17:48:21 INFO None 5166303: status FINISHED
2026-07-15 17:48:21 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:48:21 INFO None 5166305: status RUNNING/PENDING
2026-07-15 17:48:21 INFO None 5166306: status FINISHED
2026-07-15 17:48:21 INFO None 5166307: status FINISHED
2026-07-15 17:48:21 INFO None 5166308: status FINISHED
2026-07-15 17:48:21 INFO None 5166311: status FINISHED
2026-07-15 17:48:21 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:48:21 INFO None 5166315: status RUNNING/PENDING
2026-07-15 17:48:21 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:48:22 INFO None 5166320: status FINISHED
2026-07-15 17:48:22 INFO Jobs still running: ['5166304', '5166305', '5166314', '5166315', '5166318']. Waiting...
2026-07-15 17:48:37 INFO None 5166297: status FINISHED
2026-07-15 17:48:37 INFO None 5166300: status FINISHED
2026-07-15 17:48:37 INFO None 5166301: status FINISHED
2026-07-15 17:48:37 INFO None 5166302: status FINISHED
2026-07-15 17:48:37 INFO None 5166303: status FINISHED
2026-07-15 17:48:37 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:48:37 INFO None 5166305: status FINISHED
2026-07-15 17:48:37 INFO None 5166306: status FINISHED
2026-07-15 17:48:37 INFO None 5166307: status FINISHED
2026-07-15 17:48:37 INFO None 5166308: status FINISHED
2026-07-15 17:48:37 INFO None 5166311: status FINISHED
2026-07-15 17:48:37 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:48:37 INFO None 5166315: status FINISHED
2026-07-15 17:48:37 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:48:37 INFO None 5166320: status FINISHED
2026-07-15 17:48:37 INFO Jobs still running: ['5166304', '5166314', '5166318']. Waiting...
2026-07-15 17:48:53 INFO None 5166297: status FINISHED
2026-07-15 17:48:53 INFO None 5166300: status FINISHED
2026-07-15 17:48:53 INFO None 5166301: status FINISHED
2026-07-15 17:48:53 INFO None 5166302: status FINISHED
2026-07-15 17:48:54 INFO None 5166303: status FINISHED
2026-07-15 17:48:54 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:48:54 INFO None 5166305: status FINISHED
2026-07-15 17:48:54 INFO None 5166306: status FINISHED
2026-07-15 17:48:54 INFO None 5166307: status FINISHED
2026-07-15 17:48:54 INFO None 5166308: status FINISHED
2026-07-15 17:48:54 INFO None 5166311: status FINISHED
2026-07-15 17:48:54 INFO None 5166314: status RUNNING/PENDING
2026-07-15 17:48:54 INFO None 5166315: status FINISHED
2026-07-15 17:48:54 INFO None 5166318: status RUNNING/PENDING
2026-07-15 17:48:54 INFO None 5166320: status FINISHED
2026-07-15 17:48:54 INFO Jobs still running: ['5166304', '5166314', '5166318']. Waiting...
2026-07-15 17:49:09 INFO None 5166297: status FINISHED
2026-07-15 17:49:09 INFO None 5166300: status FINISHED
2026-07-15 17:49:09 INFO None 5166301: status FINISHED
2026-07-15 17:49:09 INFO None 5166302: status FINISHED
2026-07-15 17:49:09 INFO None 5166303: status FINISHED
2026-07-15 17:49:09 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:49:09 INFO None 5166305: status FINISHED
2026-07-15 17:49:09 INFO None 5166306: status FINISHED
2026-07-15 17:49:09 INFO None 5166307: status FINISHED
2026-07-15 17:49:09 INFO None 5166308: status FINISHED
2026-07-15 17:49:09 INFO None 5166311: status FINISHED
2026-07-15 17:49:09 INFO None 5166314: status FINISHED
2026-07-15 17:49:09 INFO None 5166315: status FINISHED
2026-07-15 17:49:09 INFO None 5166318: status FINISHED
2026-07-15 17:49:09 INFO None 5166320: status FINISHED
2026-07-15 17:49:09 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:49:24 INFO None 5166297: status FINISHED
2026-07-15 17:49:24 INFO None 5166300: status FINISHED
2026-07-15 17:49:24 INFO None 5166301: status FINISHED
2026-07-15 17:49:24 INFO None 5166302: status FINISHED
2026-07-15 17:49:24 INFO None 5166303: status FINISHED
2026-07-15 17:49:24 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:49:24 INFO None 5166305: status FINISHED
2026-07-15 17:49:24 INFO None 5166306: status FINISHED
2026-07-15 17:49:24 INFO None 5166307: status FINISHED
2026-07-15 17:49:24 INFO None 5166308: status FINISHED
2026-07-15 17:49:24 INFO None 5166311: status FINISHED
2026-07-15 17:49:24 INFO None 5166314: status FINISHED
2026-07-15 17:49:24 INFO None 5166315: status FINISHED
2026-07-15 17:49:24 INFO None 5166318: status FINISHED
2026-07-15 17:49:24 INFO None 5166320: status FINISHED
2026-07-15 17:49:24 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:49:39 INFO None 5166297: status FINISHED
2026-07-15 17:49:39 INFO None 5166300: status FINISHED
2026-07-15 17:49:39 INFO None 5166301: status FINISHED
2026-07-15 17:49:39 INFO None 5166302: status FINISHED
2026-07-15 17:49:39 INFO None 5166303: status FINISHED
2026-07-15 17:49:39 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:49:39 INFO None 5166305: status FINISHED
2026-07-15 17:49:39 INFO None 5166306: status FINISHED
2026-07-15 17:49:40 INFO None 5166307: status FINISHED
2026-07-15 17:49:40 INFO None 5166308: status FINISHED
2026-07-15 17:49:40 INFO None 5166311: status FINISHED
2026-07-15 17:49:40 INFO None 5166314: status FINISHED
2026-07-15 17:49:40 INFO None 5166315: status FINISHED
2026-07-15 17:49:40 INFO None 5166318: status FINISHED
2026-07-15 17:49:40 INFO None 5166320: status FINISHED
2026-07-15 17:49:40 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:49:55 INFO None 5166297: status FINISHED
2026-07-15 17:49:55 INFO None 5166300: status FINISHED
2026-07-15 17:49:55 INFO None 5166301: status FINISHED
2026-07-15 17:49:55 INFO None 5166302: status FINISHED
2026-07-15 17:49:55 INFO None 5166303: status FINISHED
2026-07-15 17:49:55 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:49:55 INFO None 5166305: status FINISHED
2026-07-15 17:49:55 INFO None 5166306: status FINISHED
2026-07-15 17:49:55 INFO None 5166307: status FINISHED
2026-07-15 17:49:55 INFO None 5166308: status FINISHED
2026-07-15 17:49:55 INFO None 5166311: status FINISHED
2026-07-15 17:49:55 INFO None 5166314: status FINISHED
2026-07-15 17:49:55 INFO None 5166315: status FINISHED
2026-07-15 17:49:55 INFO None 5166318: status FINISHED
2026-07-15 17:49:55 INFO None 5166320: status FINISHED
2026-07-15 17:49:55 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:50:10 INFO None 5166297: status FINISHED
2026-07-15 17:50:10 INFO None 5166300: status FINISHED
2026-07-15 17:50:10 INFO None 5166301: status FINISHED
2026-07-15 17:50:10 INFO None 5166302: status FINISHED
2026-07-15 17:50:10 INFO None 5166303: status FINISHED
2026-07-15 17:50:10 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:50:10 INFO None 5166305: status FINISHED
2026-07-15 17:50:10 INFO None 5166306: status FINISHED
2026-07-15 17:50:10 INFO None 5166307: status FINISHED
2026-07-15 17:50:10 INFO None 5166308: status FINISHED
2026-07-15 17:50:10 INFO None 5166311: status FINISHED
2026-07-15 17:50:10 INFO None 5166314: status FINISHED
2026-07-15 17:50:10 INFO None 5166315: status FINISHED
2026-07-15 17:50:10 INFO None 5166318: status FINISHED
2026-07-15 17:50:10 INFO None 5166320: status FINISHED
2026-07-15 17:50:10 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:50:26 INFO None 5166297: status FINISHED
2026-07-15 17:50:26 INFO None 5166300: status FINISHED
2026-07-15 17:50:26 INFO None 5166301: status FINISHED
2026-07-15 17:50:26 INFO None 5166302: status FINISHED
2026-07-15 17:50:26 INFO None 5166303: status FINISHED
2026-07-15 17:50:26 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:50:26 INFO None 5166305: status FINISHED
2026-07-15 17:50:26 INFO None 5166306: status FINISHED
2026-07-15 17:50:26 INFO None 5166307: status FINISHED
2026-07-15 17:50:26 INFO None 5166308: status FINISHED
2026-07-15 17:50:26 INFO None 5166311: status FINISHED
2026-07-15 17:50:26 INFO None 5166314: status FINISHED
2026-07-15 17:50:26 INFO None 5166315: status FINISHED
2026-07-15 17:50:26 INFO None 5166318: status FINISHED
2026-07-15 17:50:26 INFO None 5166320: status FINISHED
2026-07-15 17:50:26 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:50:41 INFO None 5166297: status FINISHED
2026-07-15 17:50:41 INFO None 5166300: status FINISHED
2026-07-15 17:50:41 INFO None 5166301: status FINISHED
2026-07-15 17:50:41 INFO None 5166302: status FINISHED
2026-07-15 17:50:41 INFO None 5166303: status FINISHED
2026-07-15 17:50:41 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:50:41 INFO None 5166305: status FINISHED
2026-07-15 17:50:41 INFO None 5166306: status FINISHED
2026-07-15 17:50:41 INFO None 5166307: status FINISHED
2026-07-15 17:50:41 INFO None 5166308: status FINISHED
2026-07-15 17:50:41 INFO None 5166311: status FINISHED
2026-07-15 17:50:41 INFO None 5166314: status FINISHED
2026-07-15 17:50:41 INFO None 5166315: status FINISHED
2026-07-15 17:50:41 INFO None 5166318: status FINISHED
2026-07-15 17:50:41 INFO None 5166320: status FINISHED
2026-07-15 17:50:41 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:50:58 INFO None 5166297: status FINISHED
2026-07-15 17:50:58 INFO None 5166300: status FINISHED
2026-07-15 17:50:58 INFO None 5166301: status FINISHED
2026-07-15 17:50:58 INFO None 5166302: status FINISHED
2026-07-15 17:50:58 INFO None 5166303: status FINISHED
2026-07-15 17:50:58 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:50:58 INFO None 5166305: status FINISHED
2026-07-15 17:50:58 INFO None 5166306: status FINISHED
2026-07-15 17:50:58 INFO None 5166307: status FINISHED
2026-07-15 17:50:58 INFO None 5166308: status FINISHED
2026-07-15 17:50:58 INFO None 5166311: status FINISHED
2026-07-15 17:50:58 INFO None 5166314: status FINISHED
2026-07-15 17:50:58 INFO None 5166315: status FINISHED
2026-07-15 17:50:58 INFO None 5166318: status FINISHED
2026-07-15 17:50:58 INFO None 5166320: status FINISHED
2026-07-15 17:50:58 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:51:13 INFO None 5166297: status FINISHED
2026-07-15 17:51:13 INFO None 5166300: status FINISHED
2026-07-15 17:51:13 INFO None 5166301: status FINISHED
2026-07-15 17:51:13 INFO None 5166302: status FINISHED
2026-07-15 17:51:13 INFO None 5166303: status FINISHED
2026-07-15 17:51:13 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:51:13 INFO None 5166305: status FINISHED
2026-07-15 17:51:13 INFO None 5166306: status FINISHED
2026-07-15 17:51:13 INFO None 5166307: status FINISHED
2026-07-15 17:51:13 INFO None 5166308: status FINISHED
2026-07-15 17:51:13 INFO None 5166311: status FINISHED
2026-07-15 17:51:13 INFO None 5166314: status FINISHED
2026-07-15 17:51:13 INFO None 5166315: status FINISHED
2026-07-15 17:51:13 INFO None 5166318: status FINISHED
2026-07-15 17:51:13 INFO None 5166320: status FINISHED
2026-07-15 17:51:13 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:51:28 INFO None 5166297: status FINISHED
2026-07-15 17:51:28 INFO None 5166300: status FINISHED
2026-07-15 17:51:28 INFO None 5166301: status FINISHED
2026-07-15 17:51:28 INFO None 5166302: status FINISHED
2026-07-15 17:51:29 INFO None 5166303: status FINISHED
2026-07-15 17:51:29 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:51:29 INFO None 5166305: status FINISHED
2026-07-15 17:51:29 INFO None 5166306: status FINISHED
2026-07-15 17:51:29 INFO None 5166307: status FINISHED
2026-07-15 17:51:29 INFO None 5166308: status FINISHED
2026-07-15 17:51:29 INFO None 5166311: status FINISHED
2026-07-15 17:51:29 INFO None 5166314: status FINISHED
2026-07-15 17:51:29 INFO None 5166315: status FINISHED
2026-07-15 17:51:29 INFO None 5166318: status FINISHED
2026-07-15 17:51:29 INFO None 5166320: status FINISHED
2026-07-15 17:51:29 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:51:44 INFO None 5166297: status FINISHED
2026-07-15 17:51:44 INFO None 5166300: status FINISHED
2026-07-15 17:51:44 INFO None 5166301: status FINISHED
2026-07-15 17:51:44 INFO None 5166302: status FINISHED
2026-07-15 17:51:44 INFO None 5166303: status FINISHED
2026-07-15 17:51:44 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:51:44 INFO None 5166305: status FINISHED
2026-07-15 17:51:44 INFO None 5166306: status FINISHED
2026-07-15 17:51:44 INFO None 5166307: status FINISHED
2026-07-15 17:51:44 INFO None 5166308: status FINISHED
2026-07-15 17:51:44 INFO None 5166311: status FINISHED
2026-07-15 17:51:44 INFO None 5166314: status FINISHED
2026-07-15 17:51:44 INFO None 5166315: status FINISHED
2026-07-15 17:51:44 INFO None 5166318: status FINISHED
2026-07-15 17:51:44 INFO None 5166320: status FINISHED
2026-07-15 17:51:44 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:52:00 INFO None 5166297: status FINISHED
2026-07-15 17:52:00 INFO None 5166300: status FINISHED
2026-07-15 17:52:00 INFO None 5166301: status FINISHED
2026-07-15 17:52:00 INFO None 5166302: status FINISHED
2026-07-15 17:52:00 INFO None 5166303: status FINISHED
2026-07-15 17:52:00 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:52:00 INFO None 5166305: status FINISHED
2026-07-15 17:52:00 INFO None 5166306: status FINISHED
2026-07-15 17:52:00 INFO None 5166307: status FINISHED
2026-07-15 17:52:00 INFO None 5166308: status FINISHED
2026-07-15 17:52:00 INFO None 5166311: status FINISHED
2026-07-15 17:52:00 INFO None 5166314: status FINISHED
2026-07-15 17:52:00 INFO None 5166315: status FINISHED
2026-07-15 17:52:00 INFO None 5166318: status FINISHED
2026-07-15 17:52:00 INFO None 5166320: status FINISHED
2026-07-15 17:52:00 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:52:15 INFO None 5166297: status FINISHED
2026-07-15 17:52:15 INFO None 5166300: status FINISHED
2026-07-15 17:52:15 INFO None 5166301: status FINISHED
2026-07-15 17:52:15 INFO None 5166302: status FINISHED
2026-07-15 17:52:15 INFO None 5166303: status FINISHED
2026-07-15 17:52:15 INFO None 5166304: status RUNNING/PENDING
2026-07-15 17:52:15 INFO None 5166305: status FINISHED
2026-07-15 17:52:15 INFO None 5166306: status FINISHED
2026-07-15 17:52:15 INFO None 5166307: status FINISHED
2026-07-15 17:52:15 INFO None 5166308: status FINISHED
2026-07-15 17:52:15 INFO None 5166311: status FINISHED
2026-07-15 17:52:15 INFO None 5166314: status FINISHED
2026-07-15 17:52:15 INFO None 5166315: status FINISHED
2026-07-15 17:52:15 INFO None 5166318: status FINISHED
2026-07-15 17:52:15 INFO None 5166320: status FINISHED
2026-07-15 17:52:15 INFO Jobs still running: ['5166304']. Waiting...
2026-07-15 17:52:30 INFO None 5166297: status FINISHED
2026-07-15 17:52:30 INFO None 5166300: status FINISHED
2026-07-15 17:52:30 INFO None 5166301: status FINISHED
2026-07-15 17:52:30 INFO None 5166302: status FINISHED
2026-07-15 17:52:30 INFO None 5166303: status FINISHED
2026-07-15 17:52:30 INFO None 5166304: status FINISHED
2026-07-15 17:52:30 INFO None 5166305: status FINISHED
2026-07-15 17:52:30 INFO None 5166306: status FINISHED
2026-07-15 17:52:31 INFO None 5166307: status FINISHED
2026-07-15 17:52:31 INFO None 5166308: status FINISHED
2026-07-15 17:52:31 INFO None 5166311: status FINISHED
2026-07-15 17:52:31 INFO None 5166314: status FINISHED
2026-07-15 17:52:31 INFO None 5166315: status FINISHED
2026-07-15 17:52:31 INFO None 5166318: status FINISHED
2026-07-15 17:52:31 INFO None 5166320: status FINISHED
2026-07-15 17:52:31 INFO Jobs ['5166297', '5166300', '5166301', '5166302', '5166303', '5166304', '5166305', '5166306', '5166307', '5166308', '5166311', '5166314', '5166315', '5166318', '5166320'] have finished
2026-07-15 17:52:31 INFO Checking restart files were created ...
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-15 17:52:31 INFO  Run_model() completed successfully.
2026-07-15 17:52:31 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:52:31 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-15 17:52:31 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 17:52:31 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-15 17:52:31 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:52:31 INFO ---------->>> Running process_satellite_data()
2026-07-15 17:52:31 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-15 17:52:31 INFO ---------->>> Running run_obs_converter()
2026-07-15 17:52:31 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-15 17:52:31 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-15 17:52:31 INFO ---------->>> Running DART
2026-07-15 17:52:31 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 17:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-15 17:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-15 17:52:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-15 17:52:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-15 17:52:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-15 17:52:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-15 17:52:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-15 17:52:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-15 17:52:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-15 17:52:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-15 17:52:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-15 17:52:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-15 17:52:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-15 17:52:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-15 17:52:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-15 17:52:35 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 17:52:35 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 17:52:35 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 17:52:35 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 17:52:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 17:52:35 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 17:52:46 INFO Found: []
2026-07-15 17:52:46 INFO No job id returned by command ./run_filter.bsh
2026-07-15 17:52:46 INFO No monitoring will be performed
2026-07-15 17:52:46 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-15 17:52:46 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020711'
2026-07-15 17:52:46 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 17:52:46 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 17:52:46 INFO run_dart() is DONE.
2026-07-15 17:52:46 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 17:52:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-15 17:52:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:52:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-15 17:52:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:52:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-15 17:53:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-15 17:53:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-15 17:53:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-15 17:53:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-15 17:53:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-15 17:53:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-15 17:53:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-15 17:53:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-15 17:53:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-15 17:53:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-15 17:53:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:53:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-15 17:54:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:54:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-15 17:54:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 17:54:10 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 17:54:10 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:54:10 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 17:54:10 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-15 17:54:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:11 INFO Hourly dataset computed and listing created
2026-07-15 17:54:14 INFO Hourly dataset computed
2026-07-15 17:54:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:15 INFO Hourly dataset computed and listing created
2026-07-15 17:54:16 INFO Hourly dataset computed
2026-07-15 17:54:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:17 INFO Hourly dataset computed and listing created
2026-07-15 17:54:17 INFO Hourly dataset computed
2026-07-15 17:54:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:18 INFO Hourly dataset computed and listing created
2026-07-15 17:54:19 INFO Hourly dataset computed
2026-07-15 17:54:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:20 INFO Hourly dataset computed and listing created
2026-07-15 17:54:20 INFO Hourly dataset computed
2026-07-15 17:54:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:21 INFO Hourly dataset computed and listing created
2026-07-15 17:54:22 INFO Hourly dataset computed
2026-07-15 17:54:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:23 INFO Hourly dataset computed and listing created
2026-07-15 17:54:23 INFO Hourly dataset computed
2026-07-15 17:54:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:24 INFO Hourly dataset computed and listing created
2026-07-15 17:54:25 INFO Hourly dataset computed
2026-07-15 17:54:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:26 INFO Hourly dataset computed and listing created
2026-07-15 17:54:27 INFO Hourly dataset computed
2026-07-15 17:54:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:28 INFO Hourly dataset computed and listing created
2026-07-15 17:54:28 INFO Hourly dataset computed
2026-07-15 17:54:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:29 INFO Hourly dataset computed and listing created
2026-07-15 17:54:30 INFO Hourly dataset computed
2026-07-15 17:54:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:31 INFO Hourly dataset computed and listing created
2026-07-15 17:54:31 INFO Hourly dataset computed
2026-07-15 17:54:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:32 INFO Hourly dataset computed and listing created
2026-07-15 17:54:33 INFO Hourly dataset computed
2026-07-15 17:54:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:34 INFO Hourly dataset computed and listing created
2026-07-15 17:54:34 INFO Hourly dataset computed
2026-07-15 17:54:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 17:54:35 INFO Hourly dataset computed and listing created
2026-07-15 17:54:40 INFO Hourly dataset computed
2026-07-15 17:54:40 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-15 17:54:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:54:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 17:54:40 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-15 17:54:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 17:54:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:54:40 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 17:54:40 INFO Queuing job for member 1...
2026-07-15 17:54:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:54:40 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 17:54:41 INFO Found: ['5166473']
2026-07-15 17:54:46 INFO [TGCC-IRENE] Submitted job with ID:['5166473']
2026-07-15 17:54:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:54:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 17:54:46 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-15 17:54:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 17:54:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:54:46 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 17:54:46 INFO Queuing job for member 2...
2026-07-15 17:54:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:54:46 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 17:54:46 INFO Found: ['5166474']
2026-07-15 17:54:51 INFO [TGCC-IRENE] Submitted job with ID:['5166474']
2026-07-15 17:54:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:54:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 17:54:51 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-15 17:54:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 17:54:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:54:51 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 17:54:51 INFO Queuing job for member 3...
2026-07-15 17:54:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:54:51 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 17:54:52 INFO Found: ['5166475']
2026-07-15 17:54:57 INFO [TGCC-IRENE] Submitted job with ID:['5166475']
2026-07-15 17:54:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:54:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 17:54:57 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-15 17:54:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 17:54:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:54:57 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 17:54:57 INFO Queuing job for member 4...
2026-07-15 17:54:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:54:57 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 17:54:58 INFO Found: ['5166478']
2026-07-15 17:55:03 INFO [TGCC-IRENE] Submitted job with ID:['5166478']
2026-07-15 17:55:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 17:55:03 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-15 17:55:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 17:55:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:03 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 17:55:03 INFO Queuing job for member 5...
2026-07-15 17:55:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:03 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 17:55:04 INFO Found: ['5166481']
2026-07-15 17:55:09 INFO [TGCC-IRENE] Submitted job with ID:['5166481']
2026-07-15 17:55:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 17:55:09 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-15 17:55:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 17:55:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:09 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 17:55:09 INFO Queuing job for member 6...
2026-07-15 17:55:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:09 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 17:55:09 INFO Found: ['5166482']
2026-07-15 17:55:14 INFO [TGCC-IRENE] Submitted job with ID:['5166482']
2026-07-15 17:55:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 17:55:14 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-15 17:55:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 17:55:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:14 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 17:55:14 INFO Queuing job for member 7...
2026-07-15 17:55:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:14 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 17:55:15 INFO Found: ['5166483']
2026-07-15 17:55:20 INFO [TGCC-IRENE] Submitted job with ID:['5166483']
2026-07-15 17:55:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 17:55:20 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-15 17:55:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 17:55:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:20 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 17:55:20 INFO Queuing job for member 8...
2026-07-15 17:55:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:20 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 17:55:21 INFO Found: ['5166485']
2026-07-15 17:55:26 INFO [TGCC-IRENE] Submitted job with ID:['5166485']
2026-07-15 17:55:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 17:55:26 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-15 17:55:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 17:55:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:26 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 17:55:26 INFO Queuing job for member 9...
2026-07-15 17:55:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:26 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 17:55:29 INFO Found: ['5166487']
2026-07-15 17:55:34 INFO [TGCC-IRENE] Submitted job with ID:['5166487']
2026-07-15 17:55:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 17:55:34 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-15 17:55:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 17:55:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:34 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 17:55:34 INFO Queuing job for member 10...
2026-07-15 17:55:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:34 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 17:55:36 INFO Found: ['5166489']
2026-07-15 17:55:41 INFO [TGCC-IRENE] Submitted job with ID:['5166489']
2026-07-15 17:55:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 17:55:41 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-15 17:55:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 17:55:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:41 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 17:55:41 INFO Queuing job for member 11...
2026-07-15 17:55:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:41 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 17:55:42 INFO Found: ['5166490']
2026-07-15 17:55:47 INFO [TGCC-IRENE] Submitted job with ID:['5166490']
2026-07-15 17:55:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 17:55:47 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-15 17:55:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 17:55:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:47 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 17:55:47 INFO Queuing job for member 12...
2026-07-15 17:55:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:47 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 17:55:48 INFO Found: ['5166491']
2026-07-15 17:55:53 INFO [TGCC-IRENE] Submitted job with ID:['5166491']
2026-07-15 17:55:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 17:55:53 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-15 17:55:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 17:55:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:53 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 17:55:53 INFO Queuing job for member 13...
2026-07-15 17:55:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:53 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 17:55:54 INFO Found: ['5166492']
2026-07-15 17:55:59 INFO [TGCC-IRENE] Submitted job with ID:['5166492']
2026-07-15 17:55:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:55:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 17:55:59 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-15 17:55:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 17:55:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:55:59 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 17:55:59 INFO Queuing job for member 14...
2026-07-15 17:55:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:55:59 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 17:55:59 INFO Found: ['5166493']
2026-07-15 17:56:04 INFO [TGCC-IRENE] Submitted job with ID:['5166493']
2026-07-15 17:56:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 17:56:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 17:56:04 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-15 17:56:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 17:56:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 17:56:04 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 17:56:04 INFO Queuing job for member 15...
2026-07-15 17:56:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 17:56:04 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 17:56:05 INFO Found: ['5166495']
2026-07-15 17:56:10 INFO [TGCC-IRENE] Submitted job with ID:['5166495']
2026-07-15 17:56:10 INFO Checking job status ...
2026-07-15 17:56:10 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:56:10 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:56:10 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:56:27 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:56:27 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:56:27 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:56:28 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:56:28 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:56:43 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:56:43 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:56:44 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:56:44 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:56:44 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:56:44 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:56:44 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:56:44 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:56:59 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:57:00 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:57:00 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:57:15 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:57:15 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:57:15 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:57:31 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:57:31 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:57:32 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:57:32 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:57:32 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:57:32 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:57:32 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:57:47 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:57:47 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:57:49 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:57:49 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:57:49 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:58:04 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166475: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:58:04 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:58:04 INFO Jobs still running: ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:58:19 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166475: status FINISHED
2026-07-15 17:58:19 INFO None 5166478: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:58:19 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:58:19 INFO Jobs still running: ['5166473', '5166474', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:58:34 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166474: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166475: status FINISHED
2026-07-15 17:58:35 INFO None 5166478: status FINISHED
2026-07-15 17:58:35 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166483: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:58:35 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:58:35 INFO Jobs still running: ['5166473', '5166474', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:58:50 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166474: status FINISHED
2026-07-15 17:58:50 INFO None 5166475: status FINISHED
2026-07-15 17:58:50 INFO None 5166478: status FINISHED
2026-07-15 17:58:50 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166482: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166483: status FINISHED
2026-07-15 17:58:50 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:58:50 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:58:50 INFO Jobs still running: ['5166473', '5166481', '5166482', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:59:05 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166474: status FINISHED
2026-07-15 17:59:05 INFO None 5166475: status FINISHED
2026-07-15 17:59:05 INFO None 5166478: status FINISHED
2026-07-15 17:59:05 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166482: status FINISHED
2026-07-15 17:59:05 INFO None 5166483: status FINISHED
2026-07-15 17:59:05 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:59:05 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:59:05 INFO Jobs still running: ['5166473', '5166481', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:59:20 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:59:20 INFO None 5166474: status FINISHED
2026-07-15 17:59:20 INFO None 5166475: status FINISHED
2026-07-15 17:59:20 INFO None 5166478: status FINISHED
2026-07-15 17:59:20 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166482: status FINISHED
2026-07-15 17:59:21 INFO None 5166483: status FINISHED
2026-07-15 17:59:21 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166491: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166492: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:59:21 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:59:21 INFO Jobs still running: ['5166473', '5166481', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495']. Waiting...
2026-07-15 17:59:36 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166474: status FINISHED
2026-07-15 17:59:36 INFO None 5166475: status FINISHED
2026-07-15 17:59:36 INFO None 5166478: status FINISHED
2026-07-15 17:59:36 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166482: status FINISHED
2026-07-15 17:59:36 INFO None 5166483: status FINISHED
2026-07-15 17:59:36 INFO None 5166485: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166487: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166489: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166490: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166491: status FINISHED
2026-07-15 17:59:36 INFO None 5166492: status FINISHED
2026-07-15 17:59:36 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:59:36 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:59:36 INFO Jobs still running: ['5166473', '5166481', '5166485', '5166487', '5166489', '5166490', '5166493', '5166495']. Waiting...
2026-07-15 17:59:51 INFO None 5166473: status RUNNING/PENDING
2026-07-15 17:59:51 INFO None 5166474: status FINISHED
2026-07-15 17:59:51 INFO None 5166475: status FINISHED
2026-07-15 17:59:51 INFO None 5166478: status FINISHED
2026-07-15 17:59:51 INFO None 5166481: status RUNNING/PENDING
2026-07-15 17:59:51 INFO None 5166482: status FINISHED
2026-07-15 17:59:51 INFO None 5166483: status FINISHED
2026-07-15 17:59:52 INFO None 5166485: status FINISHED
2026-07-15 17:59:52 INFO None 5166487: status FINISHED
2026-07-15 17:59:52 INFO None 5166489: status FINISHED
2026-07-15 17:59:52 INFO None 5166490: status FINISHED
2026-07-15 17:59:52 INFO None 5166491: status FINISHED
2026-07-15 17:59:52 INFO None 5166492: status FINISHED
2026-07-15 17:59:52 INFO None 5166493: status RUNNING/PENDING
2026-07-15 17:59:52 INFO None 5166495: status RUNNING/PENDING
2026-07-15 17:59:52 INFO Jobs still running: ['5166473', '5166481', '5166493', '5166495']. Waiting...
2026-07-15 18:00:07 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:00:07 INFO None 5166474: status FINISHED
2026-07-15 18:00:07 INFO None 5166475: status FINISHED
2026-07-15 18:00:07 INFO None 5166478: status FINISHED
2026-07-15 18:00:07 INFO None 5166481: status FINISHED
2026-07-15 18:00:07 INFO None 5166482: status FINISHED
2026-07-15 18:00:07 INFO None 5166483: status FINISHED
2026-07-15 18:00:07 INFO None 5166485: status FINISHED
2026-07-15 18:00:07 INFO None 5166487: status FINISHED
2026-07-15 18:00:07 INFO None 5166489: status FINISHED
2026-07-15 18:00:07 INFO None 5166490: status FINISHED
2026-07-15 18:00:07 INFO None 5166491: status FINISHED
2026-07-15 18:00:07 INFO None 5166492: status FINISHED
2026-07-15 18:00:07 INFO None 5166493: status FINISHED
2026-07-15 18:00:07 INFO None 5166495: status FINISHED
2026-07-15 18:00:07 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:00:22 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:00:22 INFO None 5166474: status FINISHED
2026-07-15 18:00:22 INFO None 5166475: status FINISHED
2026-07-15 18:00:22 INFO None 5166478: status FINISHED
2026-07-15 18:00:22 INFO None 5166481: status FINISHED
2026-07-15 18:00:22 INFO None 5166482: status FINISHED
2026-07-15 18:00:22 INFO None 5166483: status FINISHED
2026-07-15 18:00:22 INFO None 5166485: status FINISHED
2026-07-15 18:00:22 INFO None 5166487: status FINISHED
2026-07-15 18:00:22 INFO None 5166489: status FINISHED
2026-07-15 18:00:22 INFO None 5166490: status FINISHED
2026-07-15 18:00:22 INFO None 5166491: status FINISHED
2026-07-15 18:00:22 INFO None 5166492: status FINISHED
2026-07-15 18:00:22 INFO None 5166493: status FINISHED
2026-07-15 18:00:22 INFO None 5166495: status FINISHED
2026-07-15 18:00:22 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:00:37 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:00:37 INFO None 5166474: status FINISHED
2026-07-15 18:00:37 INFO None 5166475: status FINISHED
2026-07-15 18:00:37 INFO None 5166478: status FINISHED
2026-07-15 18:00:37 INFO None 5166481: status FINISHED
2026-07-15 18:00:37 INFO None 5166482: status FINISHED
2026-07-15 18:00:37 INFO None 5166483: status FINISHED
2026-07-15 18:00:37 INFO None 5166485: status FINISHED
2026-07-15 18:00:37 INFO None 5166487: status FINISHED
2026-07-15 18:00:37 INFO None 5166489: status FINISHED
2026-07-15 18:00:37 INFO None 5166490: status FINISHED
2026-07-15 18:00:37 INFO None 5166491: status FINISHED
2026-07-15 18:00:38 INFO None 5166492: status FINISHED
2026-07-15 18:00:38 INFO None 5166493: status FINISHED
2026-07-15 18:00:38 INFO None 5166495: status FINISHED
2026-07-15 18:00:38 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:00:53 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:00:53 INFO None 5166474: status FINISHED
2026-07-15 18:00:53 INFO None 5166475: status FINISHED
2026-07-15 18:00:53 INFO None 5166478: status FINISHED
2026-07-15 18:00:54 INFO None 5166481: status FINISHED
2026-07-15 18:00:54 INFO None 5166482: status FINISHED
2026-07-15 18:00:54 INFO None 5166483: status FINISHED
2026-07-15 18:00:54 INFO None 5166485: status FINISHED
2026-07-15 18:00:54 INFO None 5166487: status FINISHED
2026-07-15 18:00:54 INFO None 5166489: status FINISHED
2026-07-15 18:00:54 INFO None 5166490: status FINISHED
2026-07-15 18:00:54 INFO None 5166491: status FINISHED
2026-07-15 18:00:54 INFO None 5166492: status FINISHED
2026-07-15 18:00:54 INFO None 5166493: status FINISHED
2026-07-15 18:00:54 INFO None 5166495: status FINISHED
2026-07-15 18:00:54 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:01:10 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:01:10 INFO None 5166474: status FINISHED
2026-07-15 18:01:10 INFO None 5166475: status FINISHED
2026-07-15 18:01:10 INFO None 5166478: status FINISHED
2026-07-15 18:01:10 INFO None 5166481: status FINISHED
2026-07-15 18:01:10 INFO None 5166482: status FINISHED
2026-07-15 18:01:10 INFO None 5166483: status FINISHED
2026-07-15 18:01:10 INFO None 5166485: status FINISHED
2026-07-15 18:01:10 INFO None 5166487: status FINISHED
2026-07-15 18:01:10 INFO None 5166489: status FINISHED
2026-07-15 18:01:10 INFO None 5166490: status FINISHED
2026-07-15 18:01:10 INFO None 5166491: status FINISHED
2026-07-15 18:01:10 INFO None 5166492: status FINISHED
2026-07-15 18:01:10 INFO None 5166493: status FINISHED
2026-07-15 18:01:10 INFO None 5166495: status FINISHED
2026-07-15 18:01:10 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:01:25 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:01:25 INFO None 5166474: status FINISHED
2026-07-15 18:01:25 INFO None 5166475: status FINISHED
2026-07-15 18:01:25 INFO None 5166478: status FINISHED
2026-07-15 18:01:25 INFO None 5166481: status FINISHED
2026-07-15 18:01:25 INFO None 5166482: status FINISHED
2026-07-15 18:01:25 INFO None 5166483: status FINISHED
2026-07-15 18:01:25 INFO None 5166485: status FINISHED
2026-07-15 18:01:25 INFO None 5166487: status FINISHED
2026-07-15 18:01:25 INFO None 5166489: status FINISHED
2026-07-15 18:01:25 INFO None 5166490: status FINISHED
2026-07-15 18:01:25 INFO None 5166491: status FINISHED
2026-07-15 18:01:25 INFO None 5166492: status FINISHED
2026-07-15 18:01:25 INFO None 5166493: status FINISHED
2026-07-15 18:01:25 INFO None 5166495: status FINISHED
2026-07-15 18:01:25 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:01:40 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:01:40 INFO None 5166474: status FINISHED
2026-07-15 18:01:40 INFO None 5166475: status FINISHED
2026-07-15 18:01:40 INFO None 5166478: status FINISHED
2026-07-15 18:01:40 INFO None 5166481: status FINISHED
2026-07-15 18:01:40 INFO None 5166482: status FINISHED
2026-07-15 18:01:40 INFO None 5166483: status FINISHED
2026-07-15 18:01:40 INFO None 5166485: status FINISHED
2026-07-15 18:01:40 INFO None 5166487: status FINISHED
2026-07-15 18:01:40 INFO None 5166489: status FINISHED
2026-07-15 18:01:40 INFO None 5166490: status FINISHED
2026-07-15 18:01:40 INFO None 5166491: status FINISHED
2026-07-15 18:01:40 INFO None 5166492: status FINISHED
2026-07-15 18:01:40 INFO None 5166493: status FINISHED
2026-07-15 18:01:40 INFO None 5166495: status FINISHED
2026-07-15 18:01:40 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:01:57 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:01:57 INFO None 5166474: status FINISHED
2026-07-15 18:01:57 INFO None 5166475: status FINISHED
2026-07-15 18:01:57 INFO None 5166478: status FINISHED
2026-07-15 18:01:57 INFO None 5166481: status FINISHED
2026-07-15 18:01:57 INFO None 5166482: status FINISHED
2026-07-15 18:01:57 INFO None 5166483: status FINISHED
2026-07-15 18:01:57 INFO None 5166485: status FINISHED
2026-07-15 18:01:57 INFO None 5166487: status FINISHED
2026-07-15 18:01:57 INFO None 5166489: status FINISHED
2026-07-15 18:01:57 INFO None 5166490: status FINISHED
2026-07-15 18:01:57 INFO None 5166491: status FINISHED
2026-07-15 18:01:57 INFO None 5166492: status FINISHED
2026-07-15 18:01:57 INFO None 5166493: status FINISHED
2026-07-15 18:01:57 INFO None 5166495: status FINISHED
2026-07-15 18:01:57 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:02:12 INFO None 5166473: status RUNNING/PENDING
2026-07-15 18:02:12 INFO None 5166474: status FINISHED
2026-07-15 18:02:12 INFO None 5166475: status FINISHED
2026-07-15 18:02:12 INFO None 5166478: status FINISHED
2026-07-15 18:02:12 INFO None 5166481: status FINISHED
2026-07-15 18:02:12 INFO None 5166482: status FINISHED
2026-07-15 18:02:12 INFO None 5166483: status FINISHED
2026-07-15 18:02:12 INFO None 5166485: status FINISHED
2026-07-15 18:02:12 INFO None 5166487: status FINISHED
2026-07-15 18:02:12 INFO None 5166489: status FINISHED
2026-07-15 18:02:12 INFO None 5166490: status FINISHED
2026-07-15 18:02:12 INFO None 5166491: status FINISHED
2026-07-15 18:02:12 INFO None 5166492: status FINISHED
2026-07-15 18:02:12 INFO None 5166493: status FINISHED
2026-07-15 18:02:12 INFO None 5166495: status FINISHED
2026-07-15 18:02:12 INFO Jobs still running: ['5166473']. Waiting...
2026-07-15 18:02:27 INFO None 5166473: status FINISHED
2026-07-15 18:02:27 INFO None 5166474: status FINISHED
2026-07-15 18:02:28 INFO None 5166475: status FINISHED
2026-07-15 18:02:28 INFO None 5166478: status FINISHED
2026-07-15 18:02:28 INFO None 5166481: status FINISHED
2026-07-15 18:02:28 INFO None 5166482: status FINISHED
2026-07-15 18:02:28 INFO None 5166483: status FINISHED
2026-07-15 18:02:28 INFO None 5166485: status FINISHED
2026-07-15 18:02:28 INFO None 5166487: status FINISHED
2026-07-15 18:02:28 INFO None 5166489: status FINISHED
2026-07-15 18:02:28 INFO None 5166490: status FINISHED
2026-07-15 18:02:28 INFO None 5166491: status FINISHED
2026-07-15 18:02:28 INFO None 5166492: status FINISHED
2026-07-15 18:02:28 INFO None 5166493: status FINISHED
2026-07-15 18:02:28 INFO None 5166495: status FINISHED
2026-07-15 18:02:28 INFO Jobs ['5166473', '5166474', '5166475', '5166478', '5166481', '5166482', '5166483', '5166485', '5166487', '5166489', '5166490', '5166491', '5166492', '5166493', '5166495'] have finished
2026-07-15 18:02:28 INFO Checking restart files were created ...
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-15 18:02:28 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-15 18:02:28 INFO  Run_model() completed successfully.
2026-07-15 18:02:28 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:02:28 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-15 18:02:28 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 18:02:28 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-15 18:02:28 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:02:28 INFO ---------->>> Running process_satellite_data()
2026-07-15 18:02:28 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-15 18:02:28 INFO ---------->>> Running run_obs_converter()
2026-07-15 18:02:28 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-15 18:02:28 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-15 18:02:28 INFO ---------->>> Running DART
2026-07-15 18:02:28 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 18:02:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-15 18:02:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-15 18:02:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-15 18:02:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-15 18:02:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-15 18:02:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-15 18:02:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-15 18:02:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-15 18:02:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-15 18:02:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-15 18:02:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-15 18:02:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-15 18:02:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-15 18:02:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-15 18:02:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-15 18:02:33 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 18:02:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 18:02:33 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 18:02:33 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 18:02:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 18:02:33 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 18:02:47 INFO Found: []
2026-07-15 18:02:47 INFO No job id returned by command ./run_filter.bsh
2026-07-15 18:02:47 INFO No monitoring will be performed
2026-07-15 18:02:47 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-15 18:02:47 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:47 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:47 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:47 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020712'
2026-07-15 18:02:48 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 18:02:48 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 18:02:48 INFO run_dart() is DONE.
2026-07-15 18:02:48 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 18:02:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-15 18:02:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:02:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-15 18:02:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:02:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-15 18:02:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-15 18:03:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-15 18:03:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-15 18:03:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-15 18:03:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-15 18:03:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-15 18:03:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-15 18:03:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-15 18:03:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-15 18:03:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-15 18:03:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-15 18:03:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-15 18:03:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:03:50 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 18:03:50 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:03:50 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:03:50 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-15 18:03:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:03:52 INFO Hourly dataset computed and listing created
2026-07-15 18:03:57 INFO Hourly dataset computed
2026-07-15 18:03:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:03:58 INFO Hourly dataset computed and listing created
2026-07-15 18:03:59 INFO Hourly dataset computed
2026-07-15 18:03:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:00 INFO Hourly dataset computed and listing created
2026-07-15 18:04:00 INFO Hourly dataset computed
2026-07-15 18:04:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:01 INFO Hourly dataset computed and listing created
2026-07-15 18:04:02 INFO Hourly dataset computed
2026-07-15 18:04:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:03 INFO Hourly dataset computed and listing created
2026-07-15 18:04:04 INFO Hourly dataset computed
2026-07-15 18:04:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:05 INFO Hourly dataset computed and listing created
2026-07-15 18:04:06 INFO Hourly dataset computed
2026-07-15 18:04:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:07 INFO Hourly dataset computed and listing created
2026-07-15 18:04:08 INFO Hourly dataset computed
2026-07-15 18:04:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:09 INFO Hourly dataset computed and listing created
2026-07-15 18:04:09 INFO Hourly dataset computed
2026-07-15 18:04:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:10 INFO Hourly dataset computed and listing created
2026-07-15 18:04:11 INFO Hourly dataset computed
2026-07-15 18:04:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:12 INFO Hourly dataset computed and listing created
2026-07-15 18:04:13 INFO Hourly dataset computed
2026-07-15 18:04:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:14 INFO Hourly dataset computed and listing created
2026-07-15 18:04:15 INFO Hourly dataset computed
2026-07-15 18:04:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:16 INFO Hourly dataset computed and listing created
2026-07-15 18:04:17 INFO Hourly dataset computed
2026-07-15 18:04:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:18 INFO Hourly dataset computed and listing created
2026-07-15 18:04:19 INFO Hourly dataset computed
2026-07-15 18:04:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:20 INFO Hourly dataset computed and listing created
2026-07-15 18:04:20 INFO Hourly dataset computed
2026-07-15 18:04:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:04:21 INFO Hourly dataset computed and listing created
2026-07-15 18:04:22 INFO Hourly dataset computed
2026-07-15 18:04:22 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-15 18:04:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 18:04:22 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-15 18:04:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 18:04:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:22 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 18:04:22 INFO Queuing job for member 1...
2026-07-15 18:04:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:22 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 18:04:23 INFO Found: ['5166563']
2026-07-15 18:04:28 INFO [TGCC-IRENE] Submitted job with ID:['5166563']
2026-07-15 18:04:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 18:04:28 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-15 18:04:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 18:04:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:28 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 18:04:28 INFO Queuing job for member 2...
2026-07-15 18:04:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:28 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 18:04:29 INFO Found: ['5166565']
2026-07-15 18:04:34 INFO [TGCC-IRENE] Submitted job with ID:['5166565']
2026-07-15 18:04:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 18:04:34 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-15 18:04:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 18:04:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:34 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 18:04:34 INFO Queuing job for member 3...
2026-07-15 18:04:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:34 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 18:04:34 INFO Found: ['5166580']
2026-07-15 18:04:39 INFO [TGCC-IRENE] Submitted job with ID:['5166580']
2026-07-15 18:04:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 18:04:39 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-15 18:04:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 18:04:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:39 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 18:04:40 INFO Queuing job for member 4...
2026-07-15 18:04:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:40 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 18:04:40 INFO Found: ['5166581']
2026-07-15 18:04:45 INFO [TGCC-IRENE] Submitted job with ID:['5166581']
2026-07-15 18:04:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 18:04:45 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-15 18:04:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 18:04:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:45 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 18:04:45 INFO Queuing job for member 5...
2026-07-15 18:04:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:45 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 18:04:46 INFO Found: ['5166582']
2026-07-15 18:04:51 INFO [TGCC-IRENE] Submitted job with ID:['5166582']
2026-07-15 18:04:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 18:04:51 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-15 18:04:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 18:04:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:51 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 18:04:51 INFO Queuing job for member 6...
2026-07-15 18:04:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:51 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 18:04:52 INFO Found: ['5166583']
2026-07-15 18:04:57 INFO [TGCC-IRENE] Submitted job with ID:['5166583']
2026-07-15 18:04:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:04:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 18:04:57 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-15 18:04:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 18:04:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:04:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 18:04:57 INFO Queuing job for member 7...
2026-07-15 18:04:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:04:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 18:04:57 INFO Found: ['5166584']
2026-07-15 18:05:02 INFO [TGCC-IRENE] Submitted job with ID:['5166584']
2026-07-15 18:05:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 18:05:02 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-15 18:05:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 18:05:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:03 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 18:05:03 INFO Queuing job for member 8...
2026-07-15 18:05:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:03 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 18:05:03 INFO Found: ['5166585']
2026-07-15 18:05:08 INFO [TGCC-IRENE] Submitted job with ID:['5166585']
2026-07-15 18:05:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 18:05:08 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-15 18:05:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 18:05:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:08 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 18:05:08 INFO Queuing job for member 9...
2026-07-15 18:05:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:08 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 18:05:11 INFO Found: ['5166587']
2026-07-15 18:05:16 INFO [TGCC-IRENE] Submitted job with ID:['5166587']
2026-07-15 18:05:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 18:05:16 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-15 18:05:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 18:05:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:16 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 18:05:16 INFO Queuing job for member 10...
2026-07-15 18:05:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:16 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 18:05:17 INFO Found: ['5166588']
2026-07-15 18:05:22 INFO [TGCC-IRENE] Submitted job with ID:['5166588']
2026-07-15 18:05:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 18:05:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-15 18:05:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 18:05:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 18:05:23 INFO Queuing job for member 11...
2026-07-15 18:05:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:23 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 18:05:23 INFO Found: ['5166589']
2026-07-15 18:05:28 INFO [TGCC-IRENE] Submitted job with ID:['5166589']
2026-07-15 18:05:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 18:05:28 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-15 18:05:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 18:05:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:28 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 18:05:28 INFO Queuing job for member 12...
2026-07-15 18:05:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:28 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 18:05:29 INFO Found: ['5166591']
2026-07-15 18:05:34 INFO [TGCC-IRENE] Submitted job with ID:['5166591']
2026-07-15 18:05:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 18:05:34 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-15 18:05:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 18:05:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:34 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 18:05:34 INFO Queuing job for member 13...
2026-07-15 18:05:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:34 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 18:05:35 INFO Found: ['5166592']
2026-07-15 18:05:40 INFO [TGCC-IRENE] Submitted job with ID:['5166592']
2026-07-15 18:05:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 18:05:40 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-15 18:05:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 18:05:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:40 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 18:05:40 INFO Queuing job for member 14...
2026-07-15 18:05:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:40 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 18:05:41 INFO Found: ['5166593']
2026-07-15 18:05:46 INFO [TGCC-IRENE] Submitted job with ID:['5166593']
2026-07-15 18:05:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:05:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 18:05:46 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-15 18:05:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 18:05:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:05:46 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 18:05:46 INFO Queuing job for member 15...
2026-07-15 18:05:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:05:46 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 18:05:46 INFO Found: ['5166594']
2026-07-15 18:05:51 INFO [TGCC-IRENE] Submitted job with ID:['5166594']
2026-07-15 18:05:51 INFO Checking job status ...
2026-07-15 18:05:51 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:05:51 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:05:52 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:05:52 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:06:07 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:06:07 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:06:07 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:06:22 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:06:22 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:06:22 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:06:37 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:06:38 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:06:38 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:06:53 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:06:53 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:06:53 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:07:08 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:07:08 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:07:09 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:07:09 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:07:24 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:07:24 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:07:24 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:07:39 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:07:39 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:07:39 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:07:54 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:07:54 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:07:54 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:08:10 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:08:10 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:08:10 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:08:25 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:08:25 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:08:25 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:08:40 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166585: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:08:41 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:08:41 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:08:56 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166565: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166580: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166581: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166582: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166583: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166584: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166585: status FINISHED
2026-07-15 18:08:56 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:08:56 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:08:56 INFO Jobs still running: ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:09:13 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166565: status FINISHED
2026-07-15 18:09:13 INFO None 5166580: status FINISHED
2026-07-15 18:09:13 INFO None 5166581: status FINISHED
2026-07-15 18:09:13 INFO None 5166582: status FINISHED
2026-07-15 18:09:13 INFO None 5166583: status FINISHED
2026-07-15 18:09:13 INFO None 5166584: status FINISHED
2026-07-15 18:09:13 INFO None 5166585: status FINISHED
2026-07-15 18:09:13 INFO None 5166587: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:09:13 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:09:13 INFO Jobs still running: ['5166563', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:09:28 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166565: status FINISHED
2026-07-15 18:09:28 INFO None 5166580: status FINISHED
2026-07-15 18:09:28 INFO None 5166581: status FINISHED
2026-07-15 18:09:28 INFO None 5166582: status FINISHED
2026-07-15 18:09:28 INFO None 5166583: status FINISHED
2026-07-15 18:09:28 INFO None 5166584: status FINISHED
2026-07-15 18:09:28 INFO None 5166585: status FINISHED
2026-07-15 18:09:28 INFO None 5166587: status FINISHED
2026-07-15 18:09:28 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:09:28 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:09:28 INFO Jobs still running: ['5166563', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:09:43 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:09:43 INFO None 5166565: status FINISHED
2026-07-15 18:09:44 INFO None 5166580: status FINISHED
2026-07-15 18:09:44 INFO None 5166581: status FINISHED
2026-07-15 18:09:44 INFO None 5166582: status FINISHED
2026-07-15 18:09:44 INFO None 5166583: status FINISHED
2026-07-15 18:09:44 INFO None 5166584: status FINISHED
2026-07-15 18:09:44 INFO None 5166585: status FINISHED
2026-07-15 18:09:44 INFO None 5166587: status FINISHED
2026-07-15 18:09:44 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:09:44 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:09:44 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:09:44 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:09:44 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:09:44 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:09:44 INFO Jobs still running: ['5166563', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:09:59 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166565: status FINISHED
2026-07-15 18:09:59 INFO None 5166580: status FINISHED
2026-07-15 18:09:59 INFO None 5166581: status FINISHED
2026-07-15 18:09:59 INFO None 5166582: status FINISHED
2026-07-15 18:09:59 INFO None 5166583: status FINISHED
2026-07-15 18:09:59 INFO None 5166584: status FINISHED
2026-07-15 18:09:59 INFO None 5166585: status FINISHED
2026-07-15 18:09:59 INFO None 5166587: status FINISHED
2026-07-15 18:09:59 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166593: status RUNNING/PENDING
2026-07-15 18:09:59 INFO None 5166594: status RUNNING/PENDING
2026-07-15 18:09:59 INFO Jobs still running: ['5166563', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594']. Waiting...
2026-07-15 18:10:15 INFO None 5166563: status RUNNING/PENDING
2026-07-15 18:10:15 INFO None 5166565: status FINISHED
2026-07-15 18:10:15 INFO None 5166580: status FINISHED
2026-07-15 18:10:15 INFO None 5166581: status FINISHED
2026-07-15 18:10:15 INFO None 5166582: status FINISHED
2026-07-15 18:10:15 INFO None 5166583: status FINISHED
2026-07-15 18:10:15 INFO None 5166584: status FINISHED
2026-07-15 18:10:15 INFO None 5166585: status FINISHED
2026-07-15 18:10:15 INFO None 5166587: status FINISHED
2026-07-15 18:10:15 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:10:15 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:10:15 INFO None 5166591: status RUNNING/PENDING
2026-07-15 18:10:15 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:10:15 INFO None 5166593: status FINISHED
2026-07-15 18:10:15 INFO None 5166594: status FINISHED
2026-07-15 18:10:15 INFO Jobs still running: ['5166563', '5166588', '5166589', '5166591', '5166592']. Waiting...
2026-07-15 18:10:30 INFO None 5166563: status FINISHED
2026-07-15 18:10:30 INFO None 5166565: status FINISHED
2026-07-15 18:10:30 INFO None 5166580: status FINISHED
2026-07-15 18:10:30 INFO None 5166581: status FINISHED
2026-07-15 18:10:30 INFO None 5166582: status FINISHED
2026-07-15 18:10:30 INFO None 5166583: status FINISHED
2026-07-15 18:10:30 INFO None 5166584: status FINISHED
2026-07-15 18:10:30 INFO None 5166585: status FINISHED
2026-07-15 18:10:30 INFO None 5166587: status FINISHED
2026-07-15 18:10:30 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:10:30 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:10:30 INFO None 5166591: status FINISHED
2026-07-15 18:10:30 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:10:30 INFO None 5166593: status FINISHED
2026-07-15 18:10:30 INFO None 5166594: status FINISHED
2026-07-15 18:10:30 INFO Jobs still running: ['5166588', '5166589', '5166592']. Waiting...
2026-07-15 18:10:45 INFO None 5166563: status FINISHED
2026-07-15 18:10:45 INFO None 5166565: status FINISHED
2026-07-15 18:10:45 INFO None 5166580: status FINISHED
2026-07-15 18:10:45 INFO None 5166581: status FINISHED
2026-07-15 18:10:45 INFO None 5166582: status FINISHED
2026-07-15 18:10:45 INFO None 5166583: status FINISHED
2026-07-15 18:10:45 INFO None 5166584: status FINISHED
2026-07-15 18:10:46 INFO None 5166585: status FINISHED
2026-07-15 18:10:46 INFO None 5166587: status FINISHED
2026-07-15 18:10:46 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:10:46 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:10:46 INFO None 5166591: status FINISHED
2026-07-15 18:10:46 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:10:46 INFO None 5166593: status FINISHED
2026-07-15 18:10:46 INFO None 5166594: status FINISHED
2026-07-15 18:10:46 INFO Jobs still running: ['5166588', '5166589', '5166592']. Waiting...
2026-07-15 18:11:01 INFO None 5166563: status FINISHED
2026-07-15 18:11:01 INFO None 5166565: status FINISHED
2026-07-15 18:11:01 INFO None 5166580: status FINISHED
2026-07-15 18:11:01 INFO None 5166581: status FINISHED
2026-07-15 18:11:01 INFO None 5166582: status FINISHED
2026-07-15 18:11:01 INFO None 5166583: status FINISHED
2026-07-15 18:11:01 INFO None 5166584: status FINISHED
2026-07-15 18:11:01 INFO None 5166585: status FINISHED
2026-07-15 18:11:01 INFO None 5166587: status FINISHED
2026-07-15 18:11:01 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:11:01 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:11:01 INFO None 5166591: status FINISHED
2026-07-15 18:11:01 INFO None 5166592: status RUNNING/PENDING
2026-07-15 18:11:01 INFO None 5166593: status FINISHED
2026-07-15 18:11:01 INFO None 5166594: status FINISHED
2026-07-15 18:11:01 INFO Jobs still running: ['5166588', '5166589', '5166592']. Waiting...
2026-07-15 18:11:16 INFO None 5166563: status FINISHED
2026-07-15 18:11:16 INFO None 5166565: status FINISHED
2026-07-15 18:11:16 INFO None 5166580: status FINISHED
2026-07-15 18:11:16 INFO None 5166581: status FINISHED
2026-07-15 18:11:16 INFO None 5166582: status FINISHED
2026-07-15 18:11:16 INFO None 5166583: status FINISHED
2026-07-15 18:11:18 INFO None 5166584: status FINISHED
2026-07-15 18:11:18 INFO None 5166585: status FINISHED
2026-07-15 18:11:18 INFO None 5166587: status FINISHED
2026-07-15 18:11:18 INFO None 5166588: status RUNNING/PENDING
2026-07-15 18:11:18 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:11:18 INFO None 5166591: status FINISHED
2026-07-15 18:11:18 INFO None 5166592: status FINISHED
2026-07-15 18:11:18 INFO None 5166593: status FINISHED
2026-07-15 18:11:18 INFO None 5166594: status FINISHED
2026-07-15 18:11:18 INFO Jobs still running: ['5166588', '5166589']. Waiting...
2026-07-15 18:11:33 INFO None 5166563: status FINISHED
2026-07-15 18:11:33 INFO None 5166565: status FINISHED
2026-07-15 18:11:33 INFO None 5166580: status FINISHED
2026-07-15 18:11:33 INFO None 5166581: status FINISHED
2026-07-15 18:11:33 INFO None 5166582: status FINISHED
2026-07-15 18:11:33 INFO None 5166583: status FINISHED
2026-07-15 18:11:33 INFO None 5166584: status FINISHED
2026-07-15 18:11:33 INFO None 5166585: status FINISHED
2026-07-15 18:11:33 INFO None 5166587: status FINISHED
2026-07-15 18:11:33 INFO None 5166588: status FINISHED
2026-07-15 18:11:33 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:11:33 INFO None 5166591: status FINISHED
2026-07-15 18:11:33 INFO None 5166592: status FINISHED
2026-07-15 18:11:33 INFO None 5166593: status FINISHED
2026-07-15 18:11:33 INFO None 5166594: status FINISHED
2026-07-15 18:11:33 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:11:49 INFO None 5166563: status FINISHED
2026-07-15 18:11:49 INFO None 5166565: status FINISHED
2026-07-15 18:11:49 INFO None 5166580: status FINISHED
2026-07-15 18:11:49 INFO None 5166581: status FINISHED
2026-07-15 18:11:49 INFO None 5166582: status FINISHED
2026-07-15 18:11:49 INFO None 5166583: status FINISHED
2026-07-15 18:11:49 INFO None 5166584: status FINISHED
2026-07-15 18:11:49 INFO None 5166585: status FINISHED
2026-07-15 18:11:49 INFO None 5166587: status FINISHED
2026-07-15 18:11:49 INFO None 5166588: status FINISHED
2026-07-15 18:11:49 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:11:49 INFO None 5166591: status FINISHED
2026-07-15 18:11:49 INFO None 5166592: status FINISHED
2026-07-15 18:11:49 INFO None 5166593: status FINISHED
2026-07-15 18:11:49 INFO None 5166594: status FINISHED
2026-07-15 18:11:49 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:12:04 INFO None 5166563: status FINISHED
2026-07-15 18:12:04 INFO None 5166565: status FINISHED
2026-07-15 18:12:04 INFO None 5166580: status FINISHED
2026-07-15 18:12:04 INFO None 5166581: status FINISHED
2026-07-15 18:12:04 INFO None 5166582: status FINISHED
2026-07-15 18:12:04 INFO None 5166583: status FINISHED
2026-07-15 18:12:04 INFO None 5166584: status FINISHED
2026-07-15 18:12:04 INFO None 5166585: status FINISHED
2026-07-15 18:12:04 INFO None 5166587: status FINISHED
2026-07-15 18:12:04 INFO None 5166588: status FINISHED
2026-07-15 18:12:04 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:12:04 INFO None 5166591: status FINISHED
2026-07-15 18:12:04 INFO None 5166592: status FINISHED
2026-07-15 18:12:04 INFO None 5166593: status FINISHED
2026-07-15 18:12:04 INFO None 5166594: status FINISHED
2026-07-15 18:12:04 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:12:19 INFO None 5166563: status FINISHED
2026-07-15 18:12:19 INFO None 5166565: status FINISHED
2026-07-15 18:12:19 INFO None 5166580: status FINISHED
2026-07-15 18:12:19 INFO None 5166581: status FINISHED
2026-07-15 18:12:19 INFO None 5166582: status FINISHED
2026-07-15 18:12:19 INFO None 5166583: status FINISHED
2026-07-15 18:12:19 INFO None 5166584: status FINISHED
2026-07-15 18:12:19 INFO None 5166585: status FINISHED
2026-07-15 18:12:19 INFO None 5166587: status FINISHED
2026-07-15 18:12:21 INFO None 5166588: status FINISHED
2026-07-15 18:12:21 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:12:21 INFO None 5166591: status FINISHED
2026-07-15 18:12:21 INFO None 5166592: status FINISHED
2026-07-15 18:12:21 INFO None 5166593: status FINISHED
2026-07-15 18:12:21 INFO None 5166594: status FINISHED
2026-07-15 18:12:21 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:12:36 INFO None 5166563: status FINISHED
2026-07-15 18:12:36 INFO None 5166565: status FINISHED
2026-07-15 18:12:36 INFO None 5166580: status FINISHED
2026-07-15 18:12:36 INFO None 5166581: status FINISHED
2026-07-15 18:12:36 INFO None 5166582: status FINISHED
2026-07-15 18:12:36 INFO None 5166583: status FINISHED
2026-07-15 18:12:36 INFO None 5166584: status FINISHED
2026-07-15 18:12:37 INFO None 5166585: status FINISHED
2026-07-15 18:12:37 INFO None 5166587: status FINISHED
2026-07-15 18:12:37 INFO None 5166588: status FINISHED
2026-07-15 18:12:37 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:12:37 INFO None 5166591: status FINISHED
2026-07-15 18:12:37 INFO None 5166592: status FINISHED
2026-07-15 18:12:37 INFO None 5166593: status FINISHED
2026-07-15 18:12:37 INFO None 5166594: status FINISHED
2026-07-15 18:12:37 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:12:52 INFO None 5166563: status FINISHED
2026-07-15 18:12:52 INFO None 5166565: status FINISHED
2026-07-15 18:12:52 INFO None 5166580: status FINISHED
2026-07-15 18:12:52 INFO None 5166581: status FINISHED
2026-07-15 18:12:52 INFO None 5166582: status FINISHED
2026-07-15 18:12:52 INFO None 5166583: status FINISHED
2026-07-15 18:12:52 INFO None 5166584: status FINISHED
2026-07-15 18:12:52 INFO None 5166585: status FINISHED
2026-07-15 18:12:52 INFO None 5166587: status FINISHED
2026-07-15 18:12:52 INFO None 5166588: status FINISHED
2026-07-15 18:12:52 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:12:52 INFO None 5166591: status FINISHED
2026-07-15 18:12:52 INFO None 5166592: status FINISHED
2026-07-15 18:12:52 INFO None 5166593: status FINISHED
2026-07-15 18:12:52 INFO None 5166594: status FINISHED
2026-07-15 18:12:52 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:13:07 INFO None 5166563: status FINISHED
2026-07-15 18:13:07 INFO None 5166565: status FINISHED
2026-07-15 18:13:07 INFO None 5166580: status FINISHED
2026-07-15 18:13:07 INFO None 5166581: status FINISHED
2026-07-15 18:13:07 INFO None 5166582: status FINISHED
2026-07-15 18:13:07 INFO None 5166583: status FINISHED
2026-07-15 18:13:07 INFO None 5166584: status FINISHED
2026-07-15 18:13:07 INFO None 5166585: status FINISHED
2026-07-15 18:13:07 INFO None 5166587: status FINISHED
2026-07-15 18:13:07 INFO None 5166588: status FINISHED
2026-07-15 18:13:07 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:13:07 INFO None 5166591: status FINISHED
2026-07-15 18:13:07 INFO None 5166592: status FINISHED
2026-07-15 18:13:07 INFO None 5166593: status FINISHED
2026-07-15 18:13:07 INFO None 5166594: status FINISHED
2026-07-15 18:13:07 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:13:23 INFO None 5166563: status FINISHED
2026-07-15 18:13:23 INFO None 5166565: status FINISHED
2026-07-15 18:13:23 INFO None 5166580: status FINISHED
2026-07-15 18:13:23 INFO None 5166581: status FINISHED
2026-07-15 18:13:23 INFO None 5166582: status FINISHED
2026-07-15 18:13:23 INFO None 5166583: status FINISHED
2026-07-15 18:13:23 INFO None 5166584: status FINISHED
2026-07-15 18:13:23 INFO None 5166585: status FINISHED
2026-07-15 18:13:23 INFO None 5166587: status FINISHED
2026-07-15 18:13:23 INFO None 5166588: status FINISHED
2026-07-15 18:13:23 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:13:23 INFO None 5166591: status FINISHED
2026-07-15 18:13:23 INFO None 5166592: status FINISHED
2026-07-15 18:13:23 INFO None 5166593: status FINISHED
2026-07-15 18:13:23 INFO None 5166594: status FINISHED
2026-07-15 18:13:23 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:13:38 INFO None 5166563: status FINISHED
2026-07-15 18:13:38 INFO None 5166565: status FINISHED
2026-07-15 18:13:38 INFO None 5166580: status FINISHED
2026-07-15 18:13:38 INFO None 5166581: status FINISHED
2026-07-15 18:13:38 INFO None 5166582: status FINISHED
2026-07-15 18:13:38 INFO None 5166583: status FINISHED
2026-07-15 18:13:38 INFO None 5166584: status FINISHED
2026-07-15 18:13:38 INFO None 5166585: status FINISHED
2026-07-15 18:13:38 INFO None 5166587: status FINISHED
2026-07-15 18:13:38 INFO None 5166588: status FINISHED
2026-07-15 18:13:38 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:13:38 INFO None 5166591: status FINISHED
2026-07-15 18:13:38 INFO None 5166592: status FINISHED
2026-07-15 18:13:38 INFO None 5166593: status FINISHED
2026-07-15 18:13:38 INFO None 5166594: status FINISHED
2026-07-15 18:13:38 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:13:53 INFO None 5166563: status FINISHED
2026-07-15 18:13:53 INFO None 5166565: status FINISHED
2026-07-15 18:13:53 INFO None 5166580: status FINISHED
2026-07-15 18:13:53 INFO None 5166581: status FINISHED
2026-07-15 18:13:53 INFO None 5166582: status FINISHED
2026-07-15 18:13:53 INFO None 5166583: status FINISHED
2026-07-15 18:13:53 INFO None 5166584: status FINISHED
2026-07-15 18:13:53 INFO None 5166585: status FINISHED
2026-07-15 18:13:53 INFO None 5166587: status FINISHED
2026-07-15 18:13:53 INFO None 5166588: status FINISHED
2026-07-15 18:13:53 INFO None 5166589: status RUNNING/PENDING
2026-07-15 18:13:54 INFO None 5166591: status FINISHED
2026-07-15 18:13:54 INFO None 5166592: status FINISHED
2026-07-15 18:13:54 INFO None 5166593: status FINISHED
2026-07-15 18:13:54 INFO None 5166594: status FINISHED
2026-07-15 18:13:54 INFO Jobs still running: ['5166589']. Waiting...
2026-07-15 18:14:09 INFO None 5166563: status FINISHED
2026-07-15 18:14:09 INFO None 5166565: status FINISHED
2026-07-15 18:14:09 INFO None 5166580: status FINISHED
2026-07-15 18:14:09 INFO None 5166581: status FINISHED
2026-07-15 18:14:09 INFO None 5166582: status FINISHED
2026-07-15 18:14:09 INFO None 5166583: status FINISHED
2026-07-15 18:14:09 INFO None 5166584: status FINISHED
2026-07-15 18:14:09 INFO None 5166585: status FINISHED
2026-07-15 18:14:09 INFO None 5166587: status FINISHED
2026-07-15 18:14:09 INFO None 5166588: status FINISHED
2026-07-15 18:14:09 INFO None 5166589: status FINISHED
2026-07-15 18:14:09 INFO None 5166591: status FINISHED
2026-07-15 18:14:09 INFO None 5166592: status FINISHED
2026-07-15 18:14:09 INFO None 5166593: status FINISHED
2026-07-15 18:14:09 INFO None 5166594: status FINISHED
2026-07-15 18:14:09 INFO Jobs ['5166563', '5166565', '5166580', '5166581', '5166582', '5166583', '5166584', '5166585', '5166587', '5166588', '5166589', '5166591', '5166592', '5166593', '5166594'] have finished
2026-07-15 18:14:09 INFO Checking restart files were created ...
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc(1002685915 bytes)
2026-07-15 18:14:09 INFO  Run_model() completed successfully.
2026-07-15 18:14:09 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:14:09 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 14:00:00 days=153073 seconds=50400
2026-07-15 18:14:09 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 18:14:09 INFO [TIME] increment current_time 2020-02-07 12:00:00 -> 2020-02-07 14:00:00
2026-07-15 18:14:09 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:14:09 INFO ---------->>> Running process_satellite_data()
2026-07-15 18:14:09 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12019.nc
2026-07-15 18:14:09 INFO ---------->>> Running run_obs_converter()
2026-07-15 18:14:09 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-15 18:14:09 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-15 18:14:09 INFO ---------->>> Running DART
2026-07-15 18:14:09 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-15 18:14:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020714_1_out_toDART.nc
2026-07-15 18:14:09 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020714_1_out_toDART.nc
2026-07-15 18:14:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020714_1_out_toDART.nc
2026-07-15 18:14:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020714_1_out_toDART.nc
2026-07-15 18:14:10 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020714_1_out_toDART.nc
2026-07-15 18:14:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020714_1_out_toDART.nc
2026-07-15 18:14:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020714_1_out_toDART.nc
2026-07-15 18:14:11 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020714_1_out_toDART.nc
2026-07-15 18:14:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020714_1_out_toDART.nc
2026-07-15 18:14:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020714_1_out_toDART.nc
2026-07-15 18:14:12 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020714_1_out_toDART.nc
2026-07-15 18:14:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020714_1_out_toDART.nc
2026-07-15 18:14:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020714_1_out_toDART.nc
2026-07-15 18:14:13 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020714_1_out_toDART.nc
2026-07-15 18:14:14 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020714_1_out_toDART.nc
2026-07-15 18:14:14 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-15 18:14:14 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-15 18:14:14 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-15 18:14:14 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-15 18:14:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-15 18:14:14 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-15 18:14:23 INFO Found: []
2026-07-15 18:14:23 INFO No job id returned by command ./run_filter.bsh
2026-07-15 18:14:23 INFO No monitoring will be performed
2026-07-15 18:14:23 INFO Moving DART output files to analysis and preassim directories for date 2020020714 if present ...
2026-07-15 18:14:23 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020714'
2026-07-15 18:14:23 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-15 18:14:23 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-15 18:14:23 INFO run_dart() is DONE.
2026-07-15 18:14:23 INFO ---------->>> Running update_pollutant_in_end()
2026-07-15 18:14:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-15 18:14:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-15 18:14:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-15 18:14:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-15 18:14:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-15 18:14:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-15 18:14:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:14:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-15 18:15:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-15 18:15:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-15 18:15:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-15 18:15:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-15 18:15:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-15 18:15:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-15 18:15:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-15 18:15:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-15 18:15:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-15 18:15:49 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 18:15:49 INFO [TIME] step_end current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:15:49 INFO [TIME] step_start current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:15:49 INFO [TIME] window start=2020-02-07 14:00:00 end=2020-02-08 00:00:00 run_hours=10 has_assimilation=False
2026-07-15 18:15:49 INFO Copying EMIS of next day ...
2026-07-15 18:15:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-15 18:15:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:15:51 INFO Hourly dataset computed and listing created
2026-07-15 18:16:10 INFO Hourly dataset computed
2026-07-15 18:16:10 INFO Copying EMIS of next day ...
2026-07-15 18:16:10 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-15 18:16:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:12 INFO Hourly dataset computed and listing created
2026-07-15 18:16:27 INFO Hourly dataset computed
2026-07-15 18:16:27 INFO Copying EMIS of next day ...
2026-07-15 18:16:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-15 18:16:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:29 INFO Hourly dataset computed and listing created
2026-07-15 18:16:41 INFO Hourly dataset computed
2026-07-15 18:16:41 INFO Copying EMIS of next day ...
2026-07-15 18:16:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-15 18:16:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:43 INFO Hourly dataset computed and listing created
2026-07-15 18:16:45 INFO Hourly dataset computed
2026-07-15 18:16:45 INFO Copying EMIS of next day ...
2026-07-15 18:16:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-15 18:16:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:47 INFO Hourly dataset computed and listing created
2026-07-15 18:16:49 INFO Hourly dataset computed
2026-07-15 18:16:49 INFO Copying EMIS of next day ...
2026-07-15 18:16:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-15 18:16:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:52 INFO Hourly dataset computed and listing created
2026-07-15 18:16:54 INFO Hourly dataset computed
2026-07-15 18:16:54 INFO Copying EMIS of next day ...
2026-07-15 18:16:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-15 18:16:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:16:56 INFO Hourly dataset computed and listing created
2026-07-15 18:16:59 INFO Hourly dataset computed
2026-07-15 18:16:59 INFO Copying EMIS of next day ...
2026-07-15 18:16:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-15 18:16:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:00 INFO Hourly dataset computed and listing created
2026-07-15 18:17:03 INFO Hourly dataset computed
2026-07-15 18:17:03 INFO Copying EMIS of next day ...
2026-07-15 18:17:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-15 18:17:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:04 INFO Hourly dataset computed and listing created
2026-07-15 18:17:07 INFO Hourly dataset computed
2026-07-15 18:17:07 INFO Copying EMIS of next day ...
2026-07-15 18:17:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-15 18:17:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:09 INFO Hourly dataset computed and listing created
2026-07-15 18:17:11 INFO Hourly dataset computed
2026-07-15 18:17:12 INFO Copying EMIS of next day ...
2026-07-15 18:17:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-15 18:17:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:13 INFO Hourly dataset computed and listing created
2026-07-15 18:17:16 INFO Hourly dataset computed
2026-07-15 18:17:16 INFO Copying EMIS of next day ...
2026-07-15 18:17:16 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-15 18:17:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:17 INFO Hourly dataset computed and listing created
2026-07-15 18:17:20 INFO Hourly dataset computed
2026-07-15 18:17:20 INFO Copying EMIS of next day ...
2026-07-15 18:17:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-15 18:17:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:21 INFO Hourly dataset computed and listing created
2026-07-15 18:17:24 INFO Hourly dataset computed
2026-07-15 18:17:24 INFO Copying EMIS of next day ...
2026-07-15 18:17:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-15 18:17:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:26 INFO Hourly dataset computed and listing created
2026-07-15 18:17:28 INFO Hourly dataset computed
2026-07-15 18:17:28 INFO Copying EMIS of next day ...
2026-07-15 18:17:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-15 18:17:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-15 18:17:30 INFO Hourly dataset computed and listing created
2026-07-15 18:17:32 INFO Hourly dataset computed
2026-07-15 18:17:32 INFO ---------->>> Running CHIMERE model from 2020-02-07 14:00:00 to 2020-02-08 00:00:00
2026-07-15 18:17:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:17:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-15 18:17:32 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-15 18:17:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-15 18:17:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:17:32 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-15 18:17:32 INFO Queuing job for member 1...
2026-07-15 18:17:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:17:32 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-15 18:17:33 INFO Found: ['5166680']
2026-07-15 18:17:38 INFO [TGCC-IRENE] Submitted job with ID:['5166680']
2026-07-15 18:17:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:17:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-15 18:17:38 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-15 18:17:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-15 18:17:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:17:38 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-15 18:17:38 INFO Queuing job for member 2...
2026-07-15 18:17:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:17:38 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-15 18:17:41 INFO Found: ['5166681']
2026-07-15 18:17:46 INFO [TGCC-IRENE] Submitted job with ID:['5166681']
2026-07-15 18:17:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:17:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-15 18:17:46 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-15 18:17:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-15 18:17:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:17:46 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-15 18:17:46 INFO Queuing job for member 3...
2026-07-15 18:17:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:17:46 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-15 18:17:46 INFO Found: ['5166683']
2026-07-15 18:17:51 INFO [TGCC-IRENE] Submitted job with ID:['5166683']
2026-07-15 18:17:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:17:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-15 18:17:51 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-15 18:17:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-15 18:17:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:17:51 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-15 18:17:51 INFO Queuing job for member 4...
2026-07-15 18:17:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:17:51 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-15 18:17:52 INFO Found: ['5166684']
2026-07-15 18:17:57 INFO [TGCC-IRENE] Submitted job with ID:['5166684']
2026-07-15 18:17:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:17:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-15 18:17:57 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-15 18:17:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-15 18:17:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:17:57 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-15 18:17:57 INFO Queuing job for member 5...
2026-07-15 18:17:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:17:57 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-15 18:17:58 INFO Found: ['5166685']
2026-07-15 18:18:03 INFO [TGCC-IRENE] Submitted job with ID:['5166685']
2026-07-15 18:18:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-15 18:18:03 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-15 18:18:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-15 18:18:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:03 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-15 18:18:03 INFO Queuing job for member 6...
2026-07-15 18:18:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:03 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-15 18:18:04 INFO Found: ['5166687']
2026-07-15 18:18:09 INFO [TGCC-IRENE] Submitted job with ID:['5166687']
2026-07-15 18:18:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-15 18:18:09 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-15 18:18:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-15 18:18:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:09 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-15 18:18:09 INFO Queuing job for member 7...
2026-07-15 18:18:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:09 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-15 18:18:09 INFO Found: ['5166688']
2026-07-15 18:18:14 INFO [TGCC-IRENE] Submitted job with ID:['5166688']
2026-07-15 18:18:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-15 18:18:14 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-15 18:18:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-15 18:18:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:14 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-15 18:18:14 INFO Queuing job for member 8...
2026-07-15 18:18:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:14 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-15 18:18:15 INFO Found: ['5166690']
2026-07-15 18:18:20 INFO [TGCC-IRENE] Submitted job with ID:['5166690']
2026-07-15 18:18:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-15 18:18:20 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-15 18:18:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-15 18:18:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:20 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-15 18:18:20 INFO Queuing job for member 9...
2026-07-15 18:18:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:20 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-15 18:18:21 INFO Found: ['5166691']
2026-07-15 18:18:26 INFO [TGCC-IRENE] Submitted job with ID:['5166691']
2026-07-15 18:18:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-15 18:18:26 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-15 18:18:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-15 18:18:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:26 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-15 18:18:26 INFO Queuing job for member 10...
2026-07-15 18:18:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:26 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-15 18:18:27 INFO Found: ['5166692']
2026-07-15 18:18:32 INFO [TGCC-IRENE] Submitted job with ID:['5166692']
2026-07-15 18:18:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-15 18:18:32 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-15 18:18:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-15 18:18:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:32 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-15 18:18:32 INFO Queuing job for member 11...
2026-07-15 18:18:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:32 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-15 18:18:32 INFO Found: ['5166694']
2026-07-15 18:18:37 INFO [TGCC-IRENE] Submitted job with ID:['5166694']
2026-07-15 18:18:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-15 18:18:37 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-15 18:18:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-15 18:18:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:37 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-15 18:18:37 INFO Queuing job for member 12...
2026-07-15 18:18:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:37 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-15 18:18:40 INFO Found: ['5166696']
2026-07-15 18:18:45 INFO [TGCC-IRENE] Submitted job with ID:['5166696']
2026-07-15 18:18:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-15 18:18:45 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-15 18:18:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-15 18:18:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:45 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-15 18:18:45 INFO Queuing job for member 13...
2026-07-15 18:18:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:45 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-15 18:18:46 INFO Found: ['5166698']
2026-07-15 18:18:51 INFO [TGCC-IRENE] Submitted job with ID:['5166698']
2026-07-15 18:18:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-15 18:18:51 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-15 18:18:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-15 18:18:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:51 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-15 18:18:51 INFO Queuing job for member 14...
2026-07-15 18:18:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:51 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-15 18:18:52 INFO Found: ['5166699']
2026-07-15 18:18:57 INFO [TGCC-IRENE] Submitted job with ID:['5166699']
2026-07-15 18:18:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-15 18:18:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-15 18:18:57 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-15 18:18:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-15 18:18:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-15 18:18:57 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-15 18:18:57 INFO Queuing job for member 15...
2026-07-15 18:18:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-15 18:18:57 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-15 18:18:57 INFO Found: ['5166700']
2026-07-15 18:19:02 INFO [TGCC-IRENE] Submitted job with ID:['5166700']
2026-07-15 18:19:02 INFO Checking job status ...
2026-07-15 18:19:03 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:19:03 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:19:03 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:19:18 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:19:18 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:19:18 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:19:33 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:19:33 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:19:35 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:19:35 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:19:35 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:19:35 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:19:50 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:19:50 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:19:50 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:19:50 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:19:51 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:19:51 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:20:06 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:20:06 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:20:06 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:20:21 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:20:21 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:20:21 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:20:36 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:20:36 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:20:37 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:20:37 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:20:37 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:20:37 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:20:37 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:20:39 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:20:39 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:20:54 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:20:54 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:20:54 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:21:09 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:21:09 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:21:09 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:21:24 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:21:24 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:21:25 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:21:25 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:21:25 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:21:40 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:21:40 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:21:40 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:21:55 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:21:55 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:21:56 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:21:56 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:22:11 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:22:11 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:22:12 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:22:12 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:22:12 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:22:27 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:22:27 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:22:27 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:22:43 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:22:43 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:22:43 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:22:58 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:22:58 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:22:58 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:22:58 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:22:58 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:22:58 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:22:59 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:22:59 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:23:14 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:23:14 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:23:14 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:23:29 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:23:29 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:23:30 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:23:30 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:23:46 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:23:46 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:23:46 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:24:01 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:24:01 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:24:03 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:24:03 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:24:03 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:24:03 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:24:18 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:24:18 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:24:18 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:24:18 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:24:19 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:24:19 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:24:34 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:24:34 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:24:34 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:24:50 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:24:50 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:24:51 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:24:51 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:25:06 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:25:06 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:25:08 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:25:08 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:25:08 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:25:23 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:25:23 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:25:23 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:25:38 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:25:38 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:25:39 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:25:39 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:25:39 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:25:39 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:25:39 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:25:39 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:25:55 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:25:55 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:25:55 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:26:10 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:26:10 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:26:12 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:26:12 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:26:12 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:26:12 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:26:27 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:26:27 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:26:27 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:26:42 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:26:42 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:26:42 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:26:42 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:26:43 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:26:43 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:26:59 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:26:59 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:26:59 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:27:14 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:27:14 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:27:14 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:27:14 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:27:14 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:27:14 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:27:15 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:27:16 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:27:16 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:27:16 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:27:31 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166684: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:27:31 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:27:31 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:27:46 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:27:46 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:27:46 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:27:46 INFO None 5166684: status FINISHED
2026-07-15 18:27:46 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:27:46 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:27:46 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166691: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166692: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:27:47 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:27:47 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:28:03 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166684: status FINISHED
2026-07-15 18:28:03 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166691: status FINISHED
2026-07-15 18:28:03 INFO None 5166692: status FINISHED
2026-07-15 18:28:03 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:28:03 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:28:03 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166685', '5166687', '5166688', '5166690', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:28:18 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:28:18 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:28:18 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:28:18 INFO None 5166684: status FINISHED
2026-07-15 18:28:18 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:28:19 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:28:19 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:28:19 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:28:19 INFO None 5166691: status FINISHED
2026-07-15 18:28:20 INFO None 5166692: status FINISHED
2026-07-15 18:28:20 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:28:20 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:28:20 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:28:20 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:28:20 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:28:20 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166685', '5166687', '5166688', '5166690', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:28:35 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166684: status FINISHED
2026-07-15 18:28:35 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166688: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166691: status FINISHED
2026-07-15 18:28:35 INFO None 5166692: status FINISHED
2026-07-15 18:28:35 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166698: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:28:35 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:28:35 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166685', '5166687', '5166688', '5166690', '5166694', '5166696', '5166698', '5166699', '5166700']. Waiting...
2026-07-15 18:28:50 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166684: status FINISHED
2026-07-15 18:28:50 INFO None 5166685: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166688: status FINISHED
2026-07-15 18:28:50 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166691: status FINISHED
2026-07-15 18:28:50 INFO None 5166692: status FINISHED
2026-07-15 18:28:50 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166698: status FINISHED
2026-07-15 18:28:50 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:28:50 INFO None 5166700: status RUNNING/PENDING
2026-07-15 18:28:50 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166685', '5166687', '5166690', '5166694', '5166696', '5166699', '5166700']. Waiting...
2026-07-15 18:29:07 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166681: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166684: status FINISHED
2026-07-15 18:29:07 INFO None 5166685: status FINISHED
2026-07-15 18:29:07 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166688: status FINISHED
2026-07-15 18:29:07 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166691: status FINISHED
2026-07-15 18:29:07 INFO None 5166692: status FINISHED
2026-07-15 18:29:07 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166698: status FINISHED
2026-07-15 18:29:07 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:29:07 INFO None 5166700: status FINISHED
2026-07-15 18:29:07 INFO Jobs still running: ['5166680', '5166681', '5166683', '5166687', '5166690', '5166694', '5166696', '5166699']. Waiting...
2026-07-15 18:29:22 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:29:22 INFO None 5166681: status FINISHED
2026-07-15 18:29:22 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:29:22 INFO None 5166684: status FINISHED
2026-07-15 18:29:22 INFO None 5166685: status FINISHED
2026-07-15 18:29:22 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:29:22 INFO None 5166688: status FINISHED
2026-07-15 18:29:22 INFO None 5166690: status RUNNING/PENDING
2026-07-15 18:29:22 INFO None 5166691: status FINISHED
2026-07-15 18:29:22 INFO None 5166692: status FINISHED
2026-07-15 18:29:22 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:29:22 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:29:23 INFO None 5166698: status FINISHED
2026-07-15 18:29:23 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:29:23 INFO None 5166700: status FINISHED
2026-07-15 18:29:23 INFO Jobs still running: ['5166680', '5166683', '5166687', '5166690', '5166694', '5166696', '5166699']. Waiting...
2026-07-15 18:29:38 INFO None 5166680: status RUNNING/PENDING
2026-07-15 18:29:38 INFO None 5166681: status FINISHED
2026-07-15 18:29:38 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:29:38 INFO None 5166684: status FINISHED
2026-07-15 18:29:38 INFO None 5166685: status FINISHED
2026-07-15 18:29:38 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:29:38 INFO None 5166688: status FINISHED
2026-07-15 18:29:38 INFO None 5166690: status FINISHED
2026-07-15 18:29:38 INFO None 5166691: status FINISHED
2026-07-15 18:29:38 INFO None 5166692: status FINISHED
2026-07-15 18:29:38 INFO None 5166694: status RUNNING/PENDING
2026-07-15 18:29:38 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:29:39 INFO None 5166698: status FINISHED
2026-07-15 18:29:39 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:29:39 INFO None 5166700: status FINISHED
2026-07-15 18:29:39 INFO Jobs still running: ['5166680', '5166683', '5166687', '5166694', '5166696', '5166699']. Waiting...
2026-07-15 18:29:54 INFO None 5166680: status FINISHED
2026-07-15 18:29:54 INFO None 5166681: status FINISHED
2026-07-15 18:29:54 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:29:54 INFO None 5166684: status FINISHED
2026-07-15 18:29:54 INFO None 5166685: status FINISHED
2026-07-15 18:29:54 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:29:54 INFO None 5166688: status FINISHED
2026-07-15 18:29:54 INFO None 5166690: status FINISHED
2026-07-15 18:29:54 INFO None 5166691: status FINISHED
2026-07-15 18:29:54 INFO None 5166692: status FINISHED
2026-07-15 18:29:54 INFO None 5166694: status FINISHED
2026-07-15 18:29:54 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:29:54 INFO None 5166698: status FINISHED
2026-07-15 18:29:54 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:29:54 INFO None 5166700: status FINISHED
2026-07-15 18:29:54 INFO Jobs still running: ['5166683', '5166687', '5166696', '5166699']. Waiting...
2026-07-15 18:30:10 INFO None 5166680: status FINISHED
2026-07-15 18:30:10 INFO None 5166681: status FINISHED
2026-07-15 18:30:10 INFO None 5166683: status RUNNING/PENDING
2026-07-15 18:30:10 INFO None 5166684: status FINISHED
2026-07-15 18:30:10 INFO None 5166685: status FINISHED
2026-07-15 18:30:10 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:30:10 INFO None 5166688: status FINISHED
2026-07-15 18:30:10 INFO None 5166690: status FINISHED
2026-07-15 18:30:10 INFO None 5166691: status FINISHED
2026-07-15 18:30:10 INFO None 5166692: status FINISHED
2026-07-15 18:30:10 INFO None 5166694: status FINISHED
2026-07-15 18:30:10 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:30:10 INFO None 5166698: status FINISHED
2026-07-15 18:30:10 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:30:10 INFO None 5166700: status FINISHED
2026-07-15 18:30:10 INFO Jobs still running: ['5166683', '5166687', '5166696', '5166699']. Waiting...
2026-07-15 18:30:25 INFO None 5166680: status FINISHED
2026-07-15 18:30:25 INFO None 5166681: status FINISHED
2026-07-15 18:30:25 INFO None 5166683: status FINISHED
2026-07-15 18:30:25 INFO None 5166684: status FINISHED
2026-07-15 18:30:25 INFO None 5166685: status FINISHED
2026-07-15 18:30:25 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:30:25 INFO None 5166688: status FINISHED
2026-07-15 18:30:25 INFO None 5166690: status FINISHED
2026-07-15 18:30:26 INFO None 5166691: status FINISHED
2026-07-15 18:30:26 INFO None 5166692: status FINISHED
2026-07-15 18:30:26 INFO None 5166694: status FINISHED
2026-07-15 18:30:26 INFO None 5166696: status RUNNING/PENDING
2026-07-15 18:30:26 INFO None 5166698: status FINISHED
2026-07-15 18:30:27 INFO None 5166699: status RUNNING/PENDING
2026-07-15 18:30:27 INFO None 5166700: status FINISHED
2026-07-15 18:30:27 INFO Jobs still running: ['5166687', '5166696', '5166699']. Waiting...
2026-07-15 18:30:42 INFO None 5166680: status FINISHED
2026-07-15 18:30:42 INFO None 5166681: status FINISHED
2026-07-15 18:30:42 INFO None 5166683: status FINISHED
2026-07-15 18:30:42 INFO None 5166684: status FINISHED
2026-07-15 18:30:42 INFO None 5166685: status FINISHED
2026-07-15 18:30:42 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:30:42 INFO None 5166688: status FINISHED
2026-07-15 18:30:42 INFO None 5166690: status FINISHED
2026-07-15 18:30:42 INFO None 5166691: status FINISHED
2026-07-15 18:30:42 INFO None 5166692: status FINISHED
2026-07-15 18:30:42 INFO None 5166694: status FINISHED
2026-07-15 18:30:42 INFO None 5166696: status FINISHED
2026-07-15 18:30:42 INFO None 5166698: status FINISHED
2026-07-15 18:30:42 INFO None 5166699: status FINISHED
2026-07-15 18:30:42 INFO None 5166700: status FINISHED
2026-07-15 18:30:42 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:30:57 INFO None 5166680: status FINISHED
2026-07-15 18:30:57 INFO None 5166681: status FINISHED
2026-07-15 18:30:57 INFO None 5166683: status FINISHED
2026-07-15 18:30:57 INFO None 5166684: status FINISHED
2026-07-15 18:30:57 INFO None 5166685: status FINISHED
2026-07-15 18:30:57 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:30:57 INFO None 5166688: status FINISHED
2026-07-15 18:30:57 INFO None 5166690: status FINISHED
2026-07-15 18:30:57 INFO None 5166691: status FINISHED
2026-07-15 18:30:57 INFO None 5166692: status FINISHED
2026-07-15 18:30:57 INFO None 5166694: status FINISHED
2026-07-15 18:30:57 INFO None 5166696: status FINISHED
2026-07-15 18:30:57 INFO None 5166698: status FINISHED
2026-07-15 18:30:57 INFO None 5166699: status FINISHED
2026-07-15 18:30:57 INFO None 5166700: status FINISHED
2026-07-15 18:30:57 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:31:14 INFO None 5166680: status FINISHED
2026-07-15 18:31:14 INFO None 5166681: status FINISHED
2026-07-15 18:31:14 INFO None 5166683: status FINISHED
2026-07-15 18:31:14 INFO None 5166684: status FINISHED
2026-07-15 18:31:14 INFO None 5166685: status FINISHED
2026-07-15 18:31:14 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:31:14 INFO None 5166688: status FINISHED
2026-07-15 18:31:14 INFO None 5166690: status FINISHED
2026-07-15 18:31:14 INFO None 5166691: status FINISHED
2026-07-15 18:31:14 INFO None 5166692: status FINISHED
2026-07-15 18:31:14 INFO None 5166694: status FINISHED
2026-07-15 18:31:14 INFO None 5166696: status FINISHED
2026-07-15 18:31:14 INFO None 5166698: status FINISHED
2026-07-15 18:31:14 INFO None 5166699: status FINISHED
2026-07-15 18:31:14 INFO None 5166700: status FINISHED
2026-07-15 18:31:14 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:31:29 INFO None 5166680: status FINISHED
2026-07-15 18:31:29 INFO None 5166681: status FINISHED
2026-07-15 18:31:29 INFO None 5166683: status FINISHED
2026-07-15 18:31:29 INFO None 5166684: status FINISHED
2026-07-15 18:31:29 INFO None 5166685: status FINISHED
2026-07-15 18:31:29 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:31:29 INFO None 5166688: status FINISHED
2026-07-15 18:31:29 INFO None 5166690: status FINISHED
2026-07-15 18:31:29 INFO None 5166691: status FINISHED
2026-07-15 18:31:29 INFO None 5166692: status FINISHED
2026-07-15 18:31:29 INFO None 5166694: status FINISHED
2026-07-15 18:31:29 INFO None 5166696: status FINISHED
2026-07-15 18:31:30 INFO None 5166698: status FINISHED
2026-07-15 18:31:30 INFO None 5166699: status FINISHED
2026-07-15 18:31:30 INFO None 5166700: status FINISHED
2026-07-15 18:31:30 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:31:45 INFO None 5166680: status FINISHED
2026-07-15 18:31:45 INFO None 5166681: status FINISHED
2026-07-15 18:31:45 INFO None 5166683: status FINISHED
2026-07-15 18:31:45 INFO None 5166684: status FINISHED
2026-07-15 18:31:45 INFO None 5166685: status FINISHED
2026-07-15 18:31:45 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:31:45 INFO None 5166688: status FINISHED
2026-07-15 18:31:45 INFO None 5166690: status FINISHED
2026-07-15 18:31:45 INFO None 5166691: status FINISHED
2026-07-15 18:31:45 INFO None 5166692: status FINISHED
2026-07-15 18:31:45 INFO None 5166694: status FINISHED
2026-07-15 18:31:45 INFO None 5166696: status FINISHED
2026-07-15 18:31:45 INFO None 5166698: status FINISHED
2026-07-15 18:31:45 INFO None 5166699: status FINISHED
2026-07-15 18:31:46 INFO None 5166700: status FINISHED
2026-07-15 18:31:46 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:32:01 INFO None 5166680: status FINISHED
2026-07-15 18:32:01 INFO None 5166681: status FINISHED
2026-07-15 18:32:01 INFO None 5166683: status FINISHED
2026-07-15 18:32:01 INFO None 5166684: status FINISHED
2026-07-15 18:32:01 INFO None 5166685: status FINISHED
2026-07-15 18:32:01 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:32:01 INFO None 5166688: status FINISHED
2026-07-15 18:32:01 INFO None 5166690: status FINISHED
2026-07-15 18:32:01 INFO None 5166691: status FINISHED
2026-07-15 18:32:01 INFO None 5166692: status FINISHED
2026-07-15 18:32:01 INFO None 5166694: status FINISHED
2026-07-15 18:32:01 INFO None 5166696: status FINISHED
2026-07-15 18:32:01 INFO None 5166698: status FINISHED
2026-07-15 18:32:01 INFO None 5166699: status FINISHED
2026-07-15 18:32:01 INFO None 5166700: status FINISHED
2026-07-15 18:32:01 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:32:17 INFO None 5166680: status FINISHED
2026-07-15 18:32:17 INFO None 5166681: status FINISHED
2026-07-15 18:32:17 INFO None 5166683: status FINISHED
2026-07-15 18:32:17 INFO None 5166684: status FINISHED
2026-07-15 18:32:17 INFO None 5166685: status FINISHED
2026-07-15 18:32:17 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:32:17 INFO None 5166688: status FINISHED
2026-07-15 18:32:17 INFO None 5166690: status FINISHED
2026-07-15 18:32:17 INFO None 5166691: status FINISHED
2026-07-15 18:32:17 INFO None 5166692: status FINISHED
2026-07-15 18:32:17 INFO None 5166694: status FINISHED
2026-07-15 18:32:17 INFO None 5166696: status FINISHED
2026-07-15 18:32:17 INFO None 5166698: status FINISHED
2026-07-15 18:32:17 INFO None 5166699: status FINISHED
2026-07-15 18:32:17 INFO None 5166700: status FINISHED
2026-07-15 18:32:17 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:32:32 INFO None 5166680: status FINISHED
2026-07-15 18:32:32 INFO None 5166681: status FINISHED
2026-07-15 18:32:33 INFO None 5166683: status FINISHED
2026-07-15 18:32:33 INFO None 5166684: status FINISHED
2026-07-15 18:32:33 INFO None 5166685: status FINISHED
2026-07-15 18:32:33 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:32:33 INFO None 5166688: status FINISHED
2026-07-15 18:32:33 INFO None 5166690: status FINISHED
2026-07-15 18:32:33 INFO None 5166691: status FINISHED
2026-07-15 18:32:33 INFO None 5166692: status FINISHED
2026-07-15 18:32:33 INFO None 5166694: status FINISHED
2026-07-15 18:32:33 INFO None 5166696: status FINISHED
2026-07-15 18:32:33 INFO None 5166698: status FINISHED
2026-07-15 18:32:33 INFO None 5166699: status FINISHED
2026-07-15 18:32:33 INFO None 5166700: status FINISHED
2026-07-15 18:32:33 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:32:48 INFO None 5166680: status FINISHED
2026-07-15 18:32:48 INFO None 5166681: status FINISHED
2026-07-15 18:32:48 INFO None 5166683: status FINISHED
2026-07-15 18:32:48 INFO None 5166684: status FINISHED
2026-07-15 18:32:48 INFO None 5166685: status FINISHED
2026-07-15 18:32:48 INFO None 5166687: status RUNNING/PENDING
2026-07-15 18:32:48 INFO None 5166688: status FINISHED
2026-07-15 18:32:49 INFO None 5166690: status FINISHED
2026-07-15 18:32:49 INFO None 5166691: status FINISHED
2026-07-15 18:32:49 INFO None 5166692: status FINISHED
2026-07-15 18:32:49 INFO None 5166694: status FINISHED
2026-07-15 18:32:49 INFO None 5166696: status FINISHED
2026-07-15 18:32:49 INFO None 5166698: status FINISHED
2026-07-15 18:32:49 INFO None 5166699: status FINISHED
2026-07-15 18:32:49 INFO None 5166700: status FINISHED
2026-07-15 18:32:49 INFO Jobs still running: ['5166687']. Waiting...
2026-07-15 18:33:04 INFO None 5166680: status FINISHED
2026-07-15 18:33:04 INFO None 5166681: status FINISHED
2026-07-15 18:33:04 INFO None 5166683: status FINISHED
2026-07-15 18:33:04 INFO None 5166684: status FINISHED
2026-07-15 18:33:04 INFO None 5166685: status FINISHED
2026-07-15 18:33:04 INFO None 5166687: status FINISHED
2026-07-15 18:33:04 INFO None 5166688: status FINISHED
2026-07-15 18:33:04 INFO None 5166690: status FINISHED
2026-07-15 18:33:04 INFO None 5166691: status FINISHED
2026-07-15 18:33:04 INFO None 5166692: status FINISHED
2026-07-15 18:33:04 INFO None 5166694: status FINISHED
2026-07-15 18:33:04 INFO None 5166696: status FINISHED
2026-07-15 18:33:04 INFO None 5166698: status FINISHED
2026-07-15 18:33:04 INFO None 5166699: status FINISHED
2026-07-15 18:33:04 INFO None 5166700: status FINISHED
2026-07-15 18:33:04 INFO Jobs ['5166680', '5166681', '5166683', '5166684', '5166685', '5166687', '5166688', '5166690', '5166691', '5166692', '5166694', '5166696', '5166698', '5166699', '5166700'] have finished
2026-07-15 18:33:04 INFO Checking restart files were created ...
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020714_10_ENS1.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020714_10_ENS2.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020714_10_ENS3.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020714_10_ENS4.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020714_10_ENS5.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020714_10_ENS6.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020714_10_ENS7.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020714_10_ENS8.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020714_10_ENS9.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020714_10_ENS10.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020714_10_ENS11.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020714_10_ENS12.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020714_10_ENS13.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020714_10_ENS14.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020714_10_ENS15.nc(3673513755 bytes)
2026-07-15 18:33:04 INFO  Run_model() completed successfully.
2026-07-15 18:33:04 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 14:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:33:04 INFO [TIME] gregorian_conversion simulated_time=2020-02-08 00:00:00 days=153074 seconds=0
2026-07-15 18:33:04 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-15 18:33:04 INFO [TIME] increment current_time 2020-02-07 14:00:00 -> 2020-02-08 00:00:00
2026-07-15 18:33:04 INFO [TIME] after_increment_before_assimilation current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:33:04 INFO ---------->>> Running process_satellite_data()
2026-07-15 18:33:04 INFO [DART] No satellite data found, skipping assimilation
2026-07-15 18:33:04 INFO after_assimilation() skipped
2026-07-15 18:33:04 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-15 18:33:04 INFO [TIME] step_end current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-15 18:33:04 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
