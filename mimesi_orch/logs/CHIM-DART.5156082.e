+ SCRIPT_PID=1046662
+ /bin/bash -x /tmp/tmp.XHIQRzIku2
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
2026-07-14 16:03:25 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-14 16:03:25 INFO [PIPELINE] =======================================
2026-07-14 16:03:25 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-14 16:03:25 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-14 16:03:25 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-14 16:03:25 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260714_160325.log
2026-07-14 16:03:25 INFO [PIPELINE] =======================================
2026-07-14 16:03:25 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-14 16:03:25 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-14 16:03:25 INFO [STEP] ---- TIME LOOP START ----
2026-07-14 16:03:25 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:03:25 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-14 16:03:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:34 INFO Hourly dataset computed and listing created
2026-07-14 16:03:40 INFO Hourly dataset computed
2026-07-14 16:03:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:41 INFO Hourly dataset computed and listing created
2026-07-14 16:03:42 INFO Hourly dataset computed
2026-07-14 16:03:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:43 INFO Hourly dataset computed and listing created
2026-07-14 16:03:44 INFO Hourly dataset computed
2026-07-14 16:03:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:45 INFO Hourly dataset computed and listing created
2026-07-14 16:03:46 INFO Hourly dataset computed
2026-07-14 16:03:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:47 INFO Hourly dataset computed and listing created
2026-07-14 16:03:48 INFO Hourly dataset computed
2026-07-14 16:03:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:49 INFO Hourly dataset computed and listing created
2026-07-14 16:03:50 INFO Hourly dataset computed
2026-07-14 16:03:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:51 INFO Hourly dataset computed and listing created
2026-07-14 16:03:52 INFO Hourly dataset computed
2026-07-14 16:03:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:53 INFO Hourly dataset computed and listing created
2026-07-14 16:03:53 INFO Hourly dataset computed
2026-07-14 16:03:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:55 INFO Hourly dataset computed and listing created
2026-07-14 16:03:55 INFO Hourly dataset computed
2026-07-14 16:03:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:56 INFO Hourly dataset computed and listing created
2026-07-14 16:03:57 INFO Hourly dataset computed
2026-07-14 16:03:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:03:58 INFO Hourly dataset computed and listing created
2026-07-14 16:03:59 INFO Hourly dataset computed
2026-07-14 16:03:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:04:00 INFO Hourly dataset computed and listing created
2026-07-14 16:04:01 INFO Hourly dataset computed
2026-07-14 16:04:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:04:02 INFO Hourly dataset computed and listing created
2026-07-14 16:04:03 INFO Hourly dataset computed
2026-07-14 16:04:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:04:04 INFO Hourly dataset computed and listing created
2026-07-14 16:04:05 INFO Hourly dataset computed
2026-07-14 16:04:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:04:06 INFO Hourly dataset computed and listing created
2026-07-14 16:04:07 INFO Hourly dataset computed
2026-07-14 16:04:07 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-14 16:04:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 16:04:07 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-14 16:04:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 16:04:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:07 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 16:04:07 INFO Queuing job for member 1...
2026-07-14 16:04:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:07 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 16:04:08 INFO Found: ['5156947']
2026-07-14 16:04:13 INFO [TGCC-IRENE] Submitted job with ID:['5156947']
2026-07-14 16:04:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 16:04:13 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-14 16:04:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 16:04:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:13 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 16:04:13 INFO Queuing job for member 2...
2026-07-14 16:04:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:13 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 16:04:14 INFO Found: ['5156948']
2026-07-14 16:04:19 INFO [TGCC-IRENE] Submitted job with ID:['5156948']
2026-07-14 16:04:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 16:04:19 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-14 16:04:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 16:04:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:19 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 16:04:19 INFO Queuing job for member 3...
2026-07-14 16:04:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:19 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 16:04:20 INFO Found: ['5156949']
2026-07-14 16:04:25 INFO [TGCC-IRENE] Submitted job with ID:['5156949']
2026-07-14 16:04:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 16:04:25 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-14 16:04:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 16:04:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:25 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 16:04:25 INFO Queuing job for member 4...
2026-07-14 16:04:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:25 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 16:04:25 INFO Found: ['5156950']
2026-07-14 16:04:30 INFO [TGCC-IRENE] Submitted job with ID:['5156950']
2026-07-14 16:04:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 16:04:30 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-14 16:04:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 16:04:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:30 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 16:04:30 INFO Queuing job for member 5...
2026-07-14 16:04:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:30 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 16:04:31 INFO Found: ['5156951']
2026-07-14 16:04:36 INFO [TGCC-IRENE] Submitted job with ID:['5156951']
2026-07-14 16:04:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 16:04:36 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-14 16:04:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 16:04:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:36 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 16:04:36 INFO Queuing job for member 6...
2026-07-14 16:04:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:36 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 16:04:37 INFO Found: ['5156952']
2026-07-14 16:04:42 INFO [TGCC-IRENE] Submitted job with ID:['5156952']
2026-07-14 16:04:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 16:04:42 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-14 16:04:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 16:04:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:42 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 16:04:42 INFO Queuing job for member 7...
2026-07-14 16:04:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:42 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 16:04:43 INFO Found: ['5156953']
2026-07-14 16:04:48 INFO [TGCC-IRENE] Submitted job with ID:['5156953']
2026-07-14 16:04:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 16:04:48 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-14 16:04:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 16:04:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:48 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 16:04:48 INFO Queuing job for member 8...
2026-07-14 16:04:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:48 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 16:04:48 INFO Found: ['5156954']
2026-07-14 16:04:53 INFO [TGCC-IRENE] Submitted job with ID:['5156954']
2026-07-14 16:04:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 16:04:53 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-14 16:04:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 16:04:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:53 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 16:04:53 INFO Queuing job for member 9...
2026-07-14 16:04:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:53 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 16:04:54 INFO Found: ['5156955']
2026-07-14 16:04:59 INFO [TGCC-IRENE] Submitted job with ID:['5156955']
2026-07-14 16:04:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:04:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 16:04:59 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-14 16:04:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 16:04:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:04:59 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 16:04:59 INFO Queuing job for member 10...
2026-07-14 16:04:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:04:59 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 16:05:01 INFO Found: ['5156956']
2026-07-14 16:05:06 INFO [TGCC-IRENE] Submitted job with ID:['5156956']
2026-07-14 16:05:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:05:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 16:05:06 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-14 16:05:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 16:05:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:05:06 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 16:05:06 INFO Queuing job for member 11...
2026-07-14 16:05:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:05:06 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 16:05:08 INFO Found: ['5156958']
2026-07-14 16:05:13 INFO [TGCC-IRENE] Submitted job with ID:['5156958']
2026-07-14 16:05:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:05:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 16:05:13 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-14 16:05:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 16:05:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:05:13 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 16:05:13 INFO Queuing job for member 12...
2026-07-14 16:05:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:05:13 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 16:05:16 INFO Found: ['5156959']
2026-07-14 16:05:21 INFO [TGCC-IRENE] Submitted job with ID:['5156959']
2026-07-14 16:05:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:05:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 16:05:21 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-14 16:05:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 16:05:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:05:21 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 16:05:21 INFO Queuing job for member 13...
2026-07-14 16:05:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:05:21 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 16:05:22 INFO Found: ['5156960']
2026-07-14 16:05:27 INFO [TGCC-IRENE] Submitted job with ID:['5156960']
2026-07-14 16:05:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:05:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 16:05:27 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-14 16:05:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 16:05:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:05:27 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 16:05:27 INFO Queuing job for member 14...
2026-07-14 16:05:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:05:27 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 16:05:28 INFO Found: ['5156961']
2026-07-14 16:05:33 INFO [TGCC-IRENE] Submitted job with ID:['5156961']
2026-07-14 16:05:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:05:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 16:05:33 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-14 16:05:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 16:05:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:05:33 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 16:05:33 INFO Queuing job for member 15...
2026-07-14 16:05:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:05:33 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 16:05:33 INFO Found: ['5156962']
2026-07-14 16:05:38 INFO [TGCC-IRENE] Submitted job with ID:['5156962']
2026-07-14 16:05:38 INFO Checking job status ...
2026-07-14 16:05:38 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:05:38 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:05:39 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:05:39 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:05:39 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:05:39 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:05:39 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:05:54 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:05:54 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:05:54 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:06:09 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:06:09 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:06:09 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:06:24 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:06:24 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:06:24 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:06:39 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:06:40 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:06:40 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:06:55 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:06:55 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:06:55 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:07:10 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:07:10 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:07:10 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:07:25 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:07:25 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:07:26 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:07:26 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:07:41 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:07:41 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:07:41 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:07:56 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:07:56 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:07:56 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:08:13 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:08:13 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:08:13 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:08:28 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:08:28 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:08:29 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:08:29 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:08:44 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:08:44 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:08:44 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:08:59 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:08:59 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:08:59 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:09:16 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:09:16 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:09:17 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:09:17 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:09:17 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:09:32 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:09:32 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:09:32 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:09:47 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156951: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:09:47 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:09:47 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:10:02 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156950: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156951: status FINISHED
2026-07-14 16:10:02 INFO None 5156952: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156953: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156954: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156955: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:10:02 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:10:02 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156950', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:10:18 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156950: status FINISHED
2026-07-14 16:10:18 INFO None 5156951: status FINISHED
2026-07-14 16:10:18 INFO None 5156952: status FINISHED
2026-07-14 16:10:18 INFO None 5156953: status FINISHED
2026-07-14 16:10:18 INFO None 5156954: status FINISHED
2026-07-14 16:10:18 INFO None 5156955: status FINISHED
2026-07-14 16:10:18 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:10:18 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:10:18 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:10:33 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156950: status FINISHED
2026-07-14 16:10:33 INFO None 5156951: status FINISHED
2026-07-14 16:10:33 INFO None 5156952: status FINISHED
2026-07-14 16:10:33 INFO None 5156953: status FINISHED
2026-07-14 16:10:33 INFO None 5156954: status FINISHED
2026-07-14 16:10:33 INFO None 5156955: status FINISHED
2026-07-14 16:10:33 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156960: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:10:33 INFO None 5156962: status RUNNING/PENDING
2026-07-14 16:10:33 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962']. Waiting...
2026-07-14 16:10:48 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:10:48 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:10:48 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:10:48 INFO None 5156950: status FINISHED
2026-07-14 16:10:48 INFO None 5156951: status FINISHED
2026-07-14 16:10:48 INFO None 5156952: status FINISHED
2026-07-14 16:10:48 INFO None 5156953: status FINISHED
2026-07-14 16:10:48 INFO None 5156954: status FINISHED
2026-07-14 16:10:48 INFO None 5156955: status FINISHED
2026-07-14 16:10:48 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:10:48 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:10:49 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:10:49 INFO None 5156960: status FINISHED
2026-07-14 16:10:49 INFO None 5156961: status RUNNING/PENDING
2026-07-14 16:10:49 INFO None 5156962: status FINISHED
2026-07-14 16:10:49 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156956', '5156958', '5156959', '5156961']. Waiting...
2026-07-14 16:11:04 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156949: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156950: status FINISHED
2026-07-14 16:11:04 INFO None 5156951: status FINISHED
2026-07-14 16:11:04 INFO None 5156952: status FINISHED
2026-07-14 16:11:04 INFO None 5156953: status FINISHED
2026-07-14 16:11:04 INFO None 5156954: status FINISHED
2026-07-14 16:11:04 INFO None 5156955: status FINISHED
2026-07-14 16:11:04 INFO None 5156956: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:11:04 INFO None 5156960: status FINISHED
2026-07-14 16:11:04 INFO None 5156961: status FINISHED
2026-07-14 16:11:04 INFO None 5156962: status FINISHED
2026-07-14 16:11:04 INFO Jobs still running: ['5156947', '5156948', '5156949', '5156956', '5156958', '5156959']. Waiting...
2026-07-14 16:11:20 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:11:20 INFO None 5156948: status RUNNING/PENDING
2026-07-14 16:11:20 INFO None 5156949: status FINISHED
2026-07-14 16:11:20 INFO None 5156950: status FINISHED
2026-07-14 16:11:20 INFO None 5156951: status FINISHED
2026-07-14 16:11:20 INFO None 5156952: status FINISHED
2026-07-14 16:11:20 INFO None 5156953: status FINISHED
2026-07-14 16:11:20 INFO None 5156954: status FINISHED
2026-07-14 16:11:20 INFO None 5156955: status FINISHED
2026-07-14 16:11:20 INFO None 5156956: status FINISHED
2026-07-14 16:11:20 INFO None 5156958: status RUNNING/PENDING
2026-07-14 16:11:20 INFO None 5156959: status RUNNING/PENDING
2026-07-14 16:11:21 INFO None 5156960: status FINISHED
2026-07-14 16:11:21 INFO None 5156961: status FINISHED
2026-07-14 16:11:21 INFO None 5156962: status FINISHED
2026-07-14 16:11:21 INFO Jobs still running: ['5156947', '5156948', '5156958', '5156959']. Waiting...
2026-07-14 16:11:36 INFO None 5156947: status RUNNING/PENDING
2026-07-14 16:13:24 INFO None 5156948: status FINISHED
2026-07-14 16:13:24 INFO None 5156949: status FINISHED
2026-07-14 16:13:24 INFO None 5156950: status FINISHED
2026-07-14 16:13:24 INFO None 5156951: status FINISHED
2026-07-14 16:13:24 INFO None 5156952: status FINISHED
2026-07-14 16:13:24 INFO None 5156953: status FINISHED
2026-07-14 16:13:24 INFO None 5156954: status FINISHED
2026-07-14 16:13:24 INFO None 5156955: status FINISHED
2026-07-14 16:13:24 INFO None 5156956: status FINISHED
2026-07-14 16:13:24 INFO None 5156958: status FINISHED
2026-07-14 16:13:24 INFO None 5156959: status FINISHED
2026-07-14 16:13:24 INFO None 5156960: status FINISHED
2026-07-14 16:13:24 INFO None 5156961: status FINISHED
2026-07-14 16:13:24 INFO None 5156962: status FINISHED
2026-07-14 16:13:24 INFO Jobs still running: ['5156947']. Waiting...
2026-07-14 16:13:39 INFO None 5156947: status FINISHED
2026-07-14 16:13:39 INFO None 5156948: status FINISHED
2026-07-14 16:13:39 INFO None 5156949: status FINISHED
2026-07-14 16:13:39 INFO None 5156950: status FINISHED
2026-07-14 16:13:39 INFO None 5156951: status FINISHED
2026-07-14 16:13:39 INFO None 5156952: status FINISHED
2026-07-14 16:13:39 INFO None 5156953: status FINISHED
2026-07-14 16:13:39 INFO None 5156954: status FINISHED
2026-07-14 16:13:39 INFO None 5156955: status FINISHED
2026-07-14 16:13:39 INFO None 5156956: status FINISHED
2026-07-14 16:13:39 INFO None 5156958: status FINISHED
2026-07-14 16:13:39 INFO None 5156959: status FINISHED
2026-07-14 16:13:39 INFO None 5156960: status FINISHED
2026-07-14 16:13:39 INFO None 5156961: status FINISHED
2026-07-14 16:13:39 INFO None 5156962: status FINISHED
2026-07-14 16:13:39 INFO Jobs ['5156947', '5156948', '5156949', '5156950', '5156951', '5156952', '5156953', '5156954', '5156955', '5156956', '5156958', '5156959', '5156960', '5156961', '5156962'] have finished
2026-07-14 16:13:39 INFO Checking restart files were created ...
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-14 16:13:39 INFO  Run_model() completed successfully.
2026-07-14 16:13:39 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:13:39 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-14 16:13:39 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 16:13:39 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-14 16:13:39 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:13:39 INFO ---------->>> Running process_satellite_data()
2026-07-14 16:13:39 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-14 16:13:39 INFO ---------->>> Running run_obs_converter()
2026-07-14 16:13:39 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 16:13:39 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-14 16:13:39 INFO ---------->>> Running DART
2026-07-14 16:13:39 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 16:13:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-14 16:13:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-14 16:13:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-14 16:13:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-14 16:13:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-14 16:13:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-14 16:13:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-14 16:13:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-14 16:13:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-14 16:13:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-14 16:13:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-14 16:13:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-14 16:13:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-14 16:13:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-14 16:13:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-14 16:13:45 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 16:13:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 16:13:45 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 16:13:45 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 16:13:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 16:13:45 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 16:13:59 INFO Found: []
2026-07-14 16:13:59 INFO No job id returned by command ./run_filter.bsh
2026-07-14 16:13:59 INFO No monitoring will be performed
2026-07-14 16:13:59 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-14 16:13:59 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/analysis/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyIC_0607_15m_low_v2/preassim/2020020611'
2026-07-14 16:13:59 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 16:14:02 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 16:14:02 INFO run_dart() is DONE.
2026-07-14 16:14:02 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 16:14:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-14 16:14:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-14 16:14:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:14 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-14 16:14:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-14 16:14:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-14 16:14:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-14 16:14:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-14 16:14:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-14 16:14:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-14 16:14:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:14:58 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-14 16:15:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:15:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-14 16:15:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:15:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-14 16:15:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:15:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-14 16:15:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:15:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-14 16:15:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
/ccc/work/cont003/gen7232/demoling/mimesi_orch/orchestrator_utils.py:1302: SerializationWarning: Unable to decode time axis into full numpy.datetime64[ns] objects, continuing using cftime.datetime objects instead, reason: dates out of range. To silence this warning use a coarser resolution 'time_unit' or specify 'use_cftime=True'.
  with xr.open_dataset(dart_file, decode_timedelta=True) as dart_ds, xr.open_dataset(end_file, decode_timedelta=True) as end_ds, xr.open_dataset(out_file, decode_timedelta=True) as out_ds, xr.open_dataset(emis_file, decode_timedelta=True) as emis_ds,  xr.open_dataset(dart_in_file, decode_timedelta=True) as dart_in_ds:
2026-07-14 16:15:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-14 16:15:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 16:15:34 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 16:15:34 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:15:34 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:15:34 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-14 16:15:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:36 INFO Hourly dataset computed and listing created
2026-07-14 16:15:40 INFO Hourly dataset computed
2026-07-14 16:15:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:41 INFO Hourly dataset computed and listing created
2026-07-14 16:15:42 INFO Hourly dataset computed
2026-07-14 16:15:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:43 INFO Hourly dataset computed and listing created
2026-07-14 16:15:44 INFO Hourly dataset computed
2026-07-14 16:15:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:45 INFO Hourly dataset computed and listing created
2026-07-14 16:15:46 INFO Hourly dataset computed
2026-07-14 16:15:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:47 INFO Hourly dataset computed and listing created
2026-07-14 16:15:48 INFO Hourly dataset computed
2026-07-14 16:15:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:49 INFO Hourly dataset computed and listing created
2026-07-14 16:15:50 INFO Hourly dataset computed
2026-07-14 16:15:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:51 INFO Hourly dataset computed and listing created
2026-07-14 16:15:52 INFO Hourly dataset computed
2026-07-14 16:15:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:53 INFO Hourly dataset computed and listing created
2026-07-14 16:15:54 INFO Hourly dataset computed
2026-07-14 16:15:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:55 INFO Hourly dataset computed and listing created
2026-07-14 16:15:56 INFO Hourly dataset computed
2026-07-14 16:15:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:57 INFO Hourly dataset computed and listing created
2026-07-14 16:15:58 INFO Hourly dataset computed
2026-07-14 16:15:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:15:59 INFO Hourly dataset computed and listing created
2026-07-14 16:16:00 INFO Hourly dataset computed
2026-07-14 16:16:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:16:01 INFO Hourly dataset computed and listing created
2026-07-14 16:18:09 INFO Hourly dataset computed
2026-07-14 16:19:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:19:02 INFO Hourly dataset computed and listing created
2026-07-14 16:19:03 INFO Hourly dataset computed
2026-07-14 16:19:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:19:04 INFO Hourly dataset computed and listing created
2026-07-14 16:19:05 INFO Hourly dataset computed
2026-07-14 16:19:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 16:19:06 INFO Hourly dataset computed and listing created
2026-07-14 16:19:07 INFO Hourly dataset computed
2026-07-14 16:19:07 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-14 16:19:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1
2026-07-14 16:19:07 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-14 16:19:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 16:19:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:07 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 16:19:07 INFO Queuing job for member 1...
2026-07-14 16:19:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:07 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 16:19:08 INFO Found: ['5156994']
2026-07-14 16:19:13 INFO [TGCC-IRENE] Submitted job with ID:['5156994']
2026-07-14 16:19:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2
2026-07-14 16:19:13 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-14 16:19:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 16:19:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:13 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 16:19:13 INFO Queuing job for member 2...
2026-07-14 16:19:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:13 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 16:19:14 INFO Found: ['5156995']
2026-07-14 16:19:19 INFO [TGCC-IRENE] Submitted job with ID:['5156995']
2026-07-14 16:19:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3
2026-07-14 16:19:19 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-14 16:19:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 16:19:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:19 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 16:19:19 INFO Queuing job for member 3...
2026-07-14 16:19:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:19 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 16:19:19 INFO Found: ['5156996']
2026-07-14 16:19:24 INFO [TGCC-IRENE] Submitted job with ID:['5156996']
2026-07-14 16:19:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4
2026-07-14 16:19:24 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-14 16:19:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 16:19:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:24 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 16:19:24 INFO Queuing job for member 4...
2026-07-14 16:19:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:24 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 16:19:25 INFO Found: ['5156997']
2026-07-14 16:19:30 INFO [TGCC-IRENE] Submitted job with ID:['5156997']
2026-07-14 16:19:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5
2026-07-14 16:19:30 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-14 16:19:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 16:19:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:30 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 16:19:30 INFO Queuing job for member 5...
2026-07-14 16:19:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:30 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 16:19:31 INFO Found: ['5156998']
2026-07-14 16:19:36 INFO [TGCC-IRENE] Submitted job with ID:['5156998']
2026-07-14 16:19:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6
2026-07-14 16:19:36 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-14 16:19:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 16:19:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:36 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 16:19:36 INFO Queuing job for member 6...
2026-07-14 16:19:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:36 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 16:19:38 INFO Found: ['5156999']
2026-07-14 16:19:43 INFO [TGCC-IRENE] Submitted job with ID:['5156999']
2026-07-14 16:19:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7
2026-07-14 16:19:43 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-14 16:19:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 16:19:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:43 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 16:19:43 INFO Queuing job for member 7...
2026-07-14 16:19:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:43 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 16:19:46 INFO Found: ['5157000']
2026-07-14 16:19:51 INFO [TGCC-IRENE] Submitted job with ID:['5157000']
2026-07-14 16:19:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8
2026-07-14 16:19:51 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-14 16:19:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 16:19:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:51 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 16:19:51 INFO Queuing job for member 8...
2026-07-14 16:19:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:51 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 16:19:51 INFO Found: ['5157001']
2026-07-14 16:19:56 INFO [TGCC-IRENE] Submitted job with ID:['5157001']
2026-07-14 16:19:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:19:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9
2026-07-14 16:19:56 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-14 16:19:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 16:19:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:19:56 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 16:19:56 INFO Queuing job for member 9...
2026-07-14 16:19:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:19:56 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 16:19:57 INFO Found: ['5157002']
2026-07-14 16:20:02 INFO [TGCC-IRENE] Submitted job with ID:['5157002']
2026-07-14 16:20:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10
2026-07-14 16:20:02 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-14 16:20:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 16:20:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:02 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 16:20:02 INFO Queuing job for member 10...
2026-07-14 16:20:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:02 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 16:20:03 INFO Found: ['5157004']
2026-07-14 16:20:08 INFO [TGCC-IRENE] Submitted job with ID:['5157004']
2026-07-14 16:20:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11
2026-07-14 16:20:08 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-14 16:20:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 16:20:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:08 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 16:20:08 INFO Queuing job for member 11...
2026-07-14 16:20:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:08 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 16:20:09 INFO Found: ['5157005']
2026-07-14 16:20:14 INFO [TGCC-IRENE] Submitted job with ID:['5157005']
2026-07-14 16:20:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12
2026-07-14 16:20:14 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-14 16:20:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 16:20:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:14 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 16:20:14 INFO Queuing job for member 12...
2026-07-14 16:20:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:14 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 16:20:14 INFO Found: ['5157006']
2026-07-14 16:20:19 INFO [TGCC-IRENE] Submitted job with ID:['5157006']
2026-07-14 16:20:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13
2026-07-14 16:20:19 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-14 16:20:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 16:20:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:19 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 16:20:20 INFO Queuing job for member 13...
2026-07-14 16:20:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:20 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 16:20:20 INFO Found: ['5157007']
2026-07-14 16:20:25 INFO [TGCC-IRENE] Submitted job with ID:['5157007']
2026-07-14 16:20:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14
2026-07-14 16:20:25 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-14 16:20:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 16:20:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:25 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 16:20:25 INFO Queuing job for member 14...
2026-07-14 16:20:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:25 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 16:20:26 INFO Found: ['5157008']
2026-07-14 16:20:31 INFO [TGCC-IRENE] Submitted job with ID:['5157008']
2026-07-14 16:20:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 16:20:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15
2026-07-14 16:20:31 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-14 16:20:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 16:20:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 16:20:31 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 16:20:31 INFO Queuing job for member 15...
2026-07-14 16:20:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 16:20:31 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 16:20:32 INFO Found: ['5157009']
2026-07-14 16:20:37 INFO [TGCC-IRENE] Submitted job with ID:['5157009']
2026-07-14 16:20:37 INFO Checking job status ...
2026-07-14 16:20:37 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157001: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157002: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:20:37 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:20:37 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157001', '5157002', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:20:52 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:20:52 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:20:52 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:20:52 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157001: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157002: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:20:53 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:20:53 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157001', '5157002', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:21:08 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157001: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157002: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:21:08 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:21:08 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157001', '5157002', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:21:23 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157001: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157002: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:23:21 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:23:21 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157001', '5157002', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:23:38 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157001: status FINISHED
2026-07-14 16:23:38 INFO None 5157002: status FINISHED
2026-07-14 16:23:38 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:23:38 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:23:38 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:23:53 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157001: status FINISHED
2026-07-14 16:23:53 INFO None 5157002: status FINISHED
2026-07-14 16:23:53 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:23:53 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:23:53 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:24:08 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157001: status FINISHED
2026-07-14 16:24:08 INFO None 5157002: status FINISHED
2026-07-14 16:24:08 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:24:08 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:24:08 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:24:24 INFO None 5156994: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157001: status FINISHED
2026-07-14 16:24:24 INFO None 5157002: status FINISHED
2026-07-14 16:24:24 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:24:24 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:24:24 INFO Jobs still running: ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:24:39 INFO None 5156994: status FINISHED
2026-07-14 16:24:39 INFO None 5156995: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5156996: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:24:39 INFO None 5157001: status FINISHED
2026-07-14 16:24:39 INFO None 5157002: status FINISHED
2026-07-14 16:24:41 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:24:41 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:24:41 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:24:41 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:24:41 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:24:41 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:24:41 INFO Jobs still running: ['5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:24:56 INFO None 5156994: status FINISHED
2026-07-14 16:24:56 INFO None 5156995: status FINISHED
2026-07-14 16:24:56 INFO None 5156996: status FINISHED
2026-07-14 16:24:56 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157001: status FINISHED
2026-07-14 16:24:56 INFO None 5157002: status FINISHED
2026-07-14 16:24:56 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:24:56 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:24:56 INFO Jobs still running: ['5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:25:11 INFO None 5156994: status FINISHED
2026-07-14 16:25:11 INFO None 5156995: status FINISHED
2026-07-14 16:25:11 INFO None 5156996: status FINISHED
2026-07-14 16:25:11 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157001: status FINISHED
2026-07-14 16:25:12 INFO None 5157002: status FINISHED
2026-07-14 16:25:12 INFO None 5157004: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157005: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157006: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157007: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:25:12 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:25:12 INFO Jobs still running: ['5156997', '5156998', '5156999', '5157000', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009']. Waiting...
2026-07-14 16:25:27 INFO None 5156994: status FINISHED
2026-07-14 16:25:27 INFO None 5156995: status FINISHED
2026-07-14 16:25:27 INFO None 5156996: status FINISHED
2026-07-14 16:25:27 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:25:27 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:25:27 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:25:27 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:25:27 INFO None 5157001: status FINISHED
2026-07-14 16:25:27 INFO None 5157002: status FINISHED
2026-07-14 16:25:27 INFO None 5157004: status FINISHED
2026-07-14 16:25:27 INFO None 5157005: status FINISHED
2026-07-14 16:25:27 INFO None 5157006: status FINISHED
2026-07-14 16:25:27 INFO None 5157007: status FINISHED
2026-07-14 16:25:27 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:25:27 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:25:27 INFO Jobs still running: ['5156997', '5156998', '5156999', '5157000', '5157008', '5157009']. Waiting...
2026-07-14 16:25:42 INFO None 5156994: status FINISHED
2026-07-14 16:25:42 INFO None 5156995: status FINISHED
2026-07-14 16:25:42 INFO None 5156996: status FINISHED
2026-07-14 16:25:42 INFO None 5156997: status RUNNING/PENDING
2026-07-14 16:25:42 INFO None 5156998: status RUNNING/PENDING
2026-07-14 16:25:42 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:25:42 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:25:42 INFO None 5157001: status FINISHED
2026-07-14 16:25:42 INFO None 5157002: status FINISHED
2026-07-14 16:25:42 INFO None 5157004: status FINISHED
2026-07-14 16:25:42 INFO None 5157005: status FINISHED
2026-07-14 16:25:42 INFO None 5157006: status FINISHED
2026-07-14 16:25:42 INFO None 5157007: status FINISHED
2026-07-14 16:25:42 INFO None 5157008: status RUNNING/PENDING
2026-07-14 16:25:42 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:25:42 INFO Jobs still running: ['5156997', '5156998', '5156999', '5157000', '5157008', '5157009']. Waiting...
2026-07-14 16:25:57 INFO None 5156994: status FINISHED
2026-07-14 16:25:57 INFO None 5156995: status FINISHED
2026-07-14 16:25:57 INFO None 5156996: status FINISHED
2026-07-14 16:25:57 INFO None 5156997: status FINISHED
2026-07-14 16:25:57 INFO None 5156998: status FINISHED
2026-07-14 16:25:57 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:25:57 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:25:57 INFO None 5157001: status FINISHED
2026-07-14 16:25:57 INFO None 5157002: status FINISHED
2026-07-14 16:25:57 INFO None 5157004: status FINISHED
2026-07-14 16:25:57 INFO None 5157005: status FINISHED
2026-07-14 16:25:58 INFO None 5157006: status FINISHED
2026-07-14 16:25:58 INFO None 5157007: status FINISHED
2026-07-14 16:25:58 INFO None 5157008: status FINISHED
2026-07-14 16:25:58 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:25:58 INFO Jobs still running: ['5156999', '5157000', '5157009']. Waiting...
2026-07-14 16:26:13 INFO None 5156994: status FINISHED
2026-07-14 16:26:13 INFO None 5156995: status FINISHED
2026-07-14 16:26:13 INFO None 5156996: status FINISHED
2026-07-14 16:26:13 INFO None 5156997: status FINISHED
2026-07-14 16:26:13 INFO None 5156998: status FINISHED
2026-07-14 16:26:13 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:26:13 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:26:13 INFO None 5157001: status FINISHED
2026-07-14 16:26:13 INFO None 5157002: status FINISHED
2026-07-14 16:26:13 INFO None 5157004: status FINISHED
2026-07-14 16:26:13 INFO None 5157005: status FINISHED
2026-07-14 16:26:13 INFO None 5157006: status FINISHED
2026-07-14 16:26:13 INFO None 5157007: status FINISHED
2026-07-14 16:26:13 INFO None 5157008: status FINISHED
2026-07-14 16:26:13 INFO None 5157009: status RUNNING/PENDING
2026-07-14 16:26:13 INFO Jobs still running: ['5156999', '5157000', '5157009']. Waiting...
2026-07-14 16:26:28 INFO None 5156994: status FINISHED
2026-07-14 16:26:28 INFO None 5156995: status FINISHED
2026-07-14 16:26:28 INFO None 5156996: status FINISHED
2026-07-14 16:26:28 INFO None 5156997: status FINISHED
2026-07-14 16:26:28 INFO None 5156998: status FINISHED
2026-07-14 16:26:28 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:26:28 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:26:28 INFO None 5157001: status FINISHED
2026-07-14 16:26:28 INFO None 5157002: status FINISHED
2026-07-14 16:26:28 INFO None 5157004: status FINISHED
2026-07-14 16:26:28 INFO None 5157005: status FINISHED
2026-07-14 16:26:28 INFO None 5157006: status FINISHED
2026-07-14 16:26:28 INFO None 5157007: status FINISHED
2026-07-14 16:26:28 INFO None 5157008: status FINISHED
2026-07-14 16:26:28 INFO None 5157009: status FINISHED
2026-07-14 16:26:28 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:26:43 INFO None 5156994: status FINISHED
2026-07-14 16:26:43 INFO None 5156995: status FINISHED
2026-07-14 16:26:43 INFO None 5156996: status FINISHED
2026-07-14 16:26:43 INFO None 5156997: status FINISHED
2026-07-14 16:26:43 INFO None 5156998: status FINISHED
2026-07-14 16:26:45 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:26:45 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:26:45 INFO None 5157001: status FINISHED
2026-07-14 16:26:45 INFO None 5157002: status FINISHED
2026-07-14 16:26:45 INFO None 5157004: status FINISHED
2026-07-14 16:26:45 INFO None 5157005: status FINISHED
2026-07-14 16:26:45 INFO None 5157006: status FINISHED
2026-07-14 16:26:45 INFO None 5157007: status FINISHED
2026-07-14 16:26:45 INFO None 5157008: status FINISHED
2026-07-14 16:26:45 INFO None 5157009: status FINISHED
2026-07-14 16:26:45 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:27:00 INFO None 5156994: status FINISHED
2026-07-14 16:27:01 INFO None 5156995: status FINISHED
2026-07-14 16:27:01 INFO None 5156996: status FINISHED
2026-07-14 16:27:01 INFO None 5156997: status FINISHED
2026-07-14 16:27:01 INFO None 5156998: status FINISHED
2026-07-14 16:27:01 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:27:01 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:27:01 INFO None 5157001: status FINISHED
2026-07-14 16:27:01 INFO None 5157002: status FINISHED
2026-07-14 16:27:01 INFO None 5157004: status FINISHED
2026-07-14 16:27:01 INFO None 5157005: status FINISHED
2026-07-14 16:27:01 INFO None 5157006: status FINISHED
2026-07-14 16:27:01 INFO None 5157007: status FINISHED
2026-07-14 16:27:01 INFO None 5157008: status FINISHED
2026-07-14 16:27:01 INFO None 5157009: status FINISHED
2026-07-14 16:27:01 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:27:16 INFO None 5156994: status FINISHED
2026-07-14 16:28:47 INFO None 5156995: status FINISHED
2026-07-14 16:28:47 INFO None 5156996: status FINISHED
2026-07-14 16:28:47 INFO None 5156997: status FINISHED
2026-07-14 16:28:47 INFO None 5156998: status FINISHED
2026-07-14 16:28:47 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:28:47 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:28:47 INFO None 5157001: status FINISHED
2026-07-14 16:28:47 INFO None 5157002: status FINISHED
2026-07-14 16:28:47 INFO None 5157004: status FINISHED
2026-07-14 16:28:47 INFO None 5157005: status FINISHED
2026-07-14 16:28:47 INFO None 5157006: status FINISHED
2026-07-14 16:28:47 INFO None 5157007: status FINISHED
2026-07-14 16:28:47 INFO None 5157008: status FINISHED
2026-07-14 16:28:47 INFO None 5157009: status FINISHED
2026-07-14 16:28:47 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:29:02 INFO None 5156994: status FINISHED
2026-07-14 16:29:02 INFO None 5156995: status FINISHED
2026-07-14 16:29:03 INFO None 5156996: status FINISHED
2026-07-14 16:29:03 INFO None 5156997: status FINISHED
2026-07-14 16:29:03 INFO None 5156998: status FINISHED
2026-07-14 16:29:03 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:29:03 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:29:03 INFO None 5157001: status FINISHED
2026-07-14 16:29:03 INFO None 5157002: status FINISHED
2026-07-14 16:29:03 INFO None 5157004: status FINISHED
2026-07-14 16:29:03 INFO None 5157005: status FINISHED
2026-07-14 16:29:03 INFO None 5157006: status FINISHED
2026-07-14 16:29:03 INFO None 5157007: status FINISHED
2026-07-14 16:29:03 INFO None 5157008: status FINISHED
2026-07-14 16:29:03 INFO None 5157009: status FINISHED
2026-07-14 16:29:03 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:29:18 INFO None 5156994: status FINISHED
2026-07-14 16:29:18 INFO None 5156995: status FINISHED
2026-07-14 16:29:18 INFO None 5156996: status FINISHED
2026-07-14 16:29:18 INFO None 5156997: status FINISHED
2026-07-14 16:29:18 INFO None 5156998: status FINISHED
2026-07-14 16:29:18 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:29:18 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:29:18 INFO None 5157001: status FINISHED
2026-07-14 16:29:18 INFO None 5157002: status FINISHED
2026-07-14 16:29:18 INFO None 5157004: status FINISHED
2026-07-14 16:29:18 INFO None 5157005: status FINISHED
2026-07-14 16:29:18 INFO None 5157006: status FINISHED
2026-07-14 16:29:18 INFO None 5157007: status FINISHED
2026-07-14 16:29:18 INFO None 5157008: status FINISHED
2026-07-14 16:29:18 INFO None 5157009: status FINISHED
2026-07-14 16:29:18 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:29:33 INFO None 5156994: status FINISHED
2026-07-14 16:29:33 INFO None 5156995: status FINISHED
2026-07-14 16:29:33 INFO None 5156996: status FINISHED
2026-07-14 16:29:33 INFO None 5156997: status FINISHED
2026-07-14 16:29:33 INFO None 5156998: status FINISHED
2026-07-14 16:29:33 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:29:33 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:29:33 INFO None 5157001: status FINISHED
2026-07-14 16:29:33 INFO None 5157002: status FINISHED
2026-07-14 16:29:33 INFO None 5157004: status FINISHED
2026-07-14 16:29:33 INFO None 5157005: status FINISHED
2026-07-14 16:29:33 INFO None 5157006: status FINISHED
2026-07-14 16:29:33 INFO None 5157007: status FINISHED
2026-07-14 16:29:33 INFO None 5157008: status FINISHED
2026-07-14 16:29:33 INFO None 5157009: status FINISHED
2026-07-14 16:29:33 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:29:48 INFO None 5156994: status FINISHED
2026-07-14 16:29:48 INFO None 5156995: status FINISHED
2026-07-14 16:29:48 INFO None 5156996: status FINISHED
2026-07-14 16:29:48 INFO None 5156997: status FINISHED
2026-07-14 16:29:48 INFO None 5156998: status FINISHED
2026-07-14 16:29:48 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:29:49 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:29:49 INFO None 5157001: status FINISHED
2026-07-14 16:29:51 INFO None 5157002: status FINISHED
2026-07-14 16:29:51 INFO None 5157004: status FINISHED
2026-07-14 16:29:51 INFO None 5157005: status FINISHED
2026-07-14 16:29:51 INFO None 5157006: status FINISHED
2026-07-14 16:29:51 INFO None 5157007: status FINISHED
2026-07-14 16:29:51 INFO None 5157008: status FINISHED
2026-07-14 16:29:51 INFO None 5157009: status FINISHED
2026-07-14 16:29:51 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:30:06 INFO None 5156994: status FINISHED
2026-07-14 16:30:06 INFO None 5156995: status FINISHED
2026-07-14 16:30:06 INFO None 5156996: status FINISHED
2026-07-14 16:30:06 INFO None 5156997: status FINISHED
2026-07-14 16:30:06 INFO None 5156998: status FINISHED
2026-07-14 16:30:06 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:30:06 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:30:06 INFO None 5157001: status FINISHED
2026-07-14 16:30:06 INFO None 5157002: status FINISHED
2026-07-14 16:30:06 INFO None 5157004: status FINISHED
2026-07-14 16:30:06 INFO None 5157005: status FINISHED
2026-07-14 16:30:06 INFO None 5157006: status FINISHED
2026-07-14 16:30:06 INFO None 5157007: status FINISHED
2026-07-14 16:30:06 INFO None 5157008: status FINISHED
2026-07-14 16:30:06 INFO None 5157009: status FINISHED
2026-07-14 16:30:06 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:30:21 INFO None 5156994: status FINISHED
2026-07-14 16:30:21 INFO None 5156995: status FINISHED
2026-07-14 16:30:21 INFO None 5156996: status FINISHED
2026-07-14 16:30:21 INFO None 5156997: status FINISHED
2026-07-14 16:30:21 INFO None 5156998: status FINISHED
2026-07-14 16:30:21 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:30:21 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:30:21 INFO None 5157001: status FINISHED
2026-07-14 16:30:21 INFO None 5157002: status FINISHED
2026-07-14 16:30:21 INFO None 5157004: status FINISHED
2026-07-14 16:30:21 INFO None 5157005: status FINISHED
2026-07-14 16:30:21 INFO None 5157006: status FINISHED
2026-07-14 16:30:21 INFO None 5157007: status FINISHED
2026-07-14 16:30:21 INFO None 5157008: status FINISHED
2026-07-14 16:30:21 INFO None 5157009: status FINISHED
2026-07-14 16:30:21 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:30:36 INFO None 5156994: status FINISHED
2026-07-14 16:30:36 INFO None 5156995: status FINISHED
2026-07-14 16:30:36 INFO None 5156996: status FINISHED
2026-07-14 16:30:36 INFO None 5156997: status FINISHED
2026-07-14 16:30:36 INFO None 5156998: status FINISHED
2026-07-14 16:30:36 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:30:36 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:30:36 INFO None 5157001: status FINISHED
2026-07-14 16:30:36 INFO None 5157002: status FINISHED
2026-07-14 16:30:37 INFO None 5157004: status FINISHED
2026-07-14 16:30:37 INFO None 5157005: status FINISHED
2026-07-14 16:30:37 INFO None 5157006: status FINISHED
2026-07-14 16:30:37 INFO None 5157007: status FINISHED
2026-07-14 16:30:37 INFO None 5157008: status FINISHED
2026-07-14 16:30:37 INFO None 5157009: status FINISHED
2026-07-14 16:30:37 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:30:52 INFO None 5156994: status FINISHED
2026-07-14 16:30:52 INFO None 5156995: status FINISHED
2026-07-14 16:30:52 INFO None 5156996: status FINISHED
2026-07-14 16:30:52 INFO None 5156997: status FINISHED
2026-07-14 16:30:52 INFO None 5156998: status FINISHED
2026-07-14 16:30:52 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:30:52 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:30:52 INFO None 5157001: status FINISHED
2026-07-14 16:30:52 INFO None 5157002: status FINISHED
2026-07-14 16:30:52 INFO None 5157004: status FINISHED
2026-07-14 16:30:52 INFO None 5157005: status FINISHED
2026-07-14 16:30:52 INFO None 5157006: status FINISHED
2026-07-14 16:30:52 INFO None 5157007: status FINISHED
2026-07-14 16:30:52 INFO None 5157008: status FINISHED
2026-07-14 16:30:52 INFO None 5157009: status FINISHED
2026-07-14 16:30:52 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:31:07 INFO None 5156994: status FINISHED
2026-07-14 16:31:07 INFO None 5156995: status FINISHED
2026-07-14 16:31:07 INFO None 5156996: status FINISHED
2026-07-14 16:31:07 INFO None 5156997: status FINISHED
2026-07-14 16:31:07 INFO None 5156998: status FINISHED
2026-07-14 16:31:07 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:31:07 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:31:07 INFO None 5157001: status FINISHED
2026-07-14 16:31:07 INFO None 5157002: status FINISHED
2026-07-14 16:31:07 INFO None 5157004: status FINISHED
2026-07-14 16:31:07 INFO None 5157005: status FINISHED
2026-07-14 16:31:07 INFO None 5157006: status FINISHED
2026-07-14 16:31:07 INFO None 5157007: status FINISHED
2026-07-14 16:31:07 INFO None 5157008: status FINISHED
2026-07-14 16:31:07 INFO None 5157009: status FINISHED
2026-07-14 16:31:07 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:31:22 INFO None 5156994: status FINISHED
2026-07-14 16:31:22 INFO None 5156995: status FINISHED
2026-07-14 16:31:22 INFO None 5156996: status FINISHED
2026-07-14 16:31:22 INFO None 5156997: status FINISHED
2026-07-14 16:31:22 INFO None 5156998: status FINISHED
2026-07-14 16:31:22 INFO None 5156999: status RUNNING/PENDING
2026-07-14 16:31:22 INFO None 5157000: status RUNNING/PENDING
2026-07-14 16:31:22 INFO None 5157001: status FINISHED
2026-07-14 16:31:22 INFO None 5157002: status FINISHED
2026-07-14 16:31:22 INFO None 5157004: status FINISHED
2026-07-14 16:31:22 INFO None 5157005: status FINISHED
2026-07-14 16:31:22 INFO None 5157006: status FINISHED
2026-07-14 16:31:22 INFO None 5157007: status FINISHED
2026-07-14 16:31:22 INFO None 5157008: status FINISHED
2026-07-14 16:31:23 INFO None 5157009: status FINISHED
2026-07-14 16:31:23 INFO Jobs still running: ['5156999', '5157000']. Waiting...
2026-07-14 16:31:38 INFO None 5156994: status FINISHED
2026-07-14 16:33:36 INFO None 5156995: status FINISHED
2026-07-14 16:33:36 INFO None 5156996: status FINISHED
2026-07-14 16:33:36 INFO None 5156997: status FINISHED
2026-07-14 16:33:36 INFO None 5156998: status FINISHED
2026-07-14 16:33:36 INFO None 5156999: status FINISHED
2026-07-14 16:33:36 INFO None 5157000: status FINISHED
2026-07-14 16:33:36 INFO None 5157001: status FINISHED
2026-07-14 16:33:36 INFO None 5157002: status FINISHED
2026-07-14 16:33:36 INFO None 5157004: status FINISHED
2026-07-14 16:33:36 INFO None 5157005: status FINISHED
2026-07-14 16:33:36 INFO None 5157006: status FINISHED
2026-07-14 16:33:36 INFO None 5157007: status FINISHED
2026-07-14 16:33:36 INFO None 5157008: status FINISHED
2026-07-14 16:33:36 INFO None 5157009: status FINISHED
2026-07-14 16:33:36 INFO Jobs ['5156994', '5156995', '5156996', '5156997', '5156998', '5156999', '5157000', '5157001', '5157002', '5157004', '5157005', '5157006', '5157007', '5157008', '5157009'] have finished
2026-07-14 16:33:36 INFO Checking restart files were created ...
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-14 16:33:36 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-14 16:33:36 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-14 16:33:36 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-14 16:33:36 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-14 16:33:36 INFO Check chimere log file at: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyIC_0607_15m_low_v2/ENS15/ENS15_2020020611.out
2026-07-14 16:33:36 ERROR [PIPELINE] Error: The following chimere ENS run(s) failed (exit code 1): [6, 7, 8, 9]
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 334, in run_model
    raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
pipeline_errors.ModelRunError: The following chimere ENS run(s) failed (exit code 1): [6, 7, 8, 9]
+ exit 0
