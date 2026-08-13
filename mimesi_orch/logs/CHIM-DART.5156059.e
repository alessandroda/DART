+ SCRIPT_PID=630485
+ /bin/bash -x /tmp/tmp.qxvxXHKweZ
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
2026-07-14 13:43:34 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-14 13:43:34 INFO [PIPELINE] =======================================
2026-07-14 13:43:34 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-14 13:43:34 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-14 13:43:34 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-14 13:43:34 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260714_134334.log
2026-07-14 13:43:34 INFO [PIPELINE] =======================================
2026-07-14 13:43:34 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-14 13:43:34 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-14 13:43:34 INFO [STEP] ---- TIME LOOP START ----
2026-07-14 13:43:34 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 13:43:34 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-14 13:43:34 INFO Copying EMIS of next day ...
2026-07-14 13:43:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-07-14 13:43:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:43:44 INFO Hourly dataset computed and listing created
2026-07-14 13:44:03 INFO Hourly dataset computed
2026-07-14 13:44:03 INFO Copying EMIS of next day ...
2026-07-14 13:44:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-14 13:44:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:44:06 INFO Hourly dataset computed and listing created
2026-07-14 13:44:18 INFO Hourly dataset computed
2026-07-14 13:44:19 INFO Copying EMIS of next day ...
2026-07-14 13:44:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-14 13:44:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:44:20 INFO Hourly dataset computed and listing created
2026-07-14 13:44:36 INFO Hourly dataset computed
2026-07-14 13:44:36 INFO Copying EMIS of next day ...
2026-07-14 13:44:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-14 13:44:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:44:38 INFO Hourly dataset computed and listing created
2026-07-14 13:44:55 INFO Hourly dataset computed
2026-07-14 13:44:55 INFO Copying EMIS of next day ...
2026-07-14 13:44:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-14 13:44:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:44:57 INFO Hourly dataset computed and listing created
2026-07-14 13:45:13 INFO Hourly dataset computed
2026-07-14 13:45:13 INFO Copying EMIS of next day ...
2026-07-14 13:45:15 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-14 13:45:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:45:18 INFO Hourly dataset computed and listing created
2026-07-14 13:45:38 INFO Hourly dataset computed
2026-07-14 13:45:38 INFO Copying EMIS of next day ...
2026-07-14 13:45:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-14 13:45:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:45:40 INFO Hourly dataset computed and listing created
2026-07-14 13:45:54 INFO Hourly dataset computed
2026-07-14 13:45:54 INFO Copying EMIS of next day ...
2026-07-14 13:45:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-14 13:45:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:45:56 INFO Hourly dataset computed and listing created
2026-07-14 13:48:32 INFO Hourly dataset computed
2026-07-14 13:48:32 INFO Copying EMIS of next day ...
2026-07-14 13:48:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-14 13:48:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:48:34 INFO Hourly dataset computed and listing created
2026-07-14 13:48:51 INFO Hourly dataset computed
2026-07-14 13:48:51 INFO Copying EMIS of next day ...
2026-07-14 13:48:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-14 13:48:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:48:53 INFO Hourly dataset computed and listing created
2026-07-14 13:49:09 INFO Hourly dataset computed
2026-07-14 13:49:09 INFO Copying EMIS of next day ...
2026-07-14 13:49:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-14 13:49:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:49:11 INFO Hourly dataset computed and listing created
2026-07-14 13:49:26 INFO Hourly dataset computed
2026-07-14 13:49:26 INFO Copying EMIS of next day ...
2026-07-14 13:49:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-14 13:49:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:49:28 INFO Hourly dataset computed and listing created
2026-07-14 13:49:44 INFO Hourly dataset computed
2026-07-14 13:49:44 INFO Copying EMIS of next day ...
2026-07-14 13:49:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-14 13:49:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:49:46 INFO Hourly dataset computed and listing created
2026-07-14 13:50:03 INFO Hourly dataset computed
2026-07-14 13:50:03 INFO Copying EMIS of next day ...
2026-07-14 13:50:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-14 13:50:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:50:05 INFO Hourly dataset computed and listing created
2026-07-14 13:50:20 INFO Hourly dataset computed
2026-07-14 13:50:20 INFO Copying EMIS of next day ...
2026-07-14 13:50:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-14 13:50:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 13:50:22 INFO Hourly dataset computed and listing created
2026-07-14 13:50:38 INFO Hourly dataset computed
2026-07-14 13:50:38 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-14 13:50:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:50:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 13:50:38 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-14 13:50:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 13:50:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:50:38 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 13:50:38 INFO Queuing job for member 1...
2026-07-14 13:50:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:50:38 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 13:50:38 INFO Found: ['5156111']
2026-07-14 13:50:44 INFO [TGCC-IRENE] Submitted job with ID:['5156111']
2026-07-14 13:50:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:50:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 13:50:44 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-14 13:50:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 13:50:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:50:44 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 13:50:44 INFO Queuing job for member 2...
2026-07-14 13:50:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:50:44 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 13:50:44 INFO Found: ['5156112']
2026-07-14 13:50:49 INFO [TGCC-IRENE] Submitted job with ID:['5156112']
2026-07-14 13:50:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:50:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 13:50:49 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-14 13:50:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 13:50:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:50:49 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 13:50:49 INFO Queuing job for member 3...
2026-07-14 13:50:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:50:49 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 13:50:50 INFO Found: ['5156113']
2026-07-14 13:50:55 INFO [TGCC-IRENE] Submitted job with ID:['5156113']
2026-07-14 13:50:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:50:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 13:50:55 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-14 13:50:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 13:50:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:50:55 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 13:50:55 INFO Queuing job for member 4...
2026-07-14 13:50:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:50:55 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 13:50:56 INFO Found: ['5156115']
2026-07-14 13:51:01 INFO [TGCC-IRENE] Submitted job with ID:['5156115']
2026-07-14 13:51:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:51:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 13:51:01 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-14 13:55:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 13:55:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:07 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 13:55:07 INFO Queuing job for member 5...
2026-07-14 13:55:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:07 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 13:55:09 INFO Found: ['5156134']
2026-07-14 13:55:14 INFO [TGCC-IRENE] Submitted job with ID:['5156134']
2026-07-14 13:55:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 13:55:14 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-14 13:55:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 13:55:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:14 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 13:55:14 INFO Queuing job for member 6...
2026-07-14 13:55:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:14 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 13:55:15 INFO Found: ['5156135']
2026-07-14 13:55:20 INFO [TGCC-IRENE] Submitted job with ID:['5156135']
2026-07-14 13:55:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 13:55:20 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-14 13:55:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 13:55:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:20 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 13:55:20 INFO Queuing job for member 7...
2026-07-14 13:55:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:20 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 13:55:21 INFO Found: ['5156136']
2026-07-14 13:55:26 INFO [TGCC-IRENE] Submitted job with ID:['5156136']
2026-07-14 13:55:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 13:55:26 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-14 13:55:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 13:55:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:26 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 13:55:26 INFO Queuing job for member 8...
2026-07-14 13:55:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:26 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 13:55:26 INFO Found: ['5156137']
2026-07-14 13:55:31 INFO [TGCC-IRENE] Submitted job with ID:['5156137']
2026-07-14 13:55:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 13:55:31 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-14 13:55:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 13:55:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:32 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 13:55:32 INFO Queuing job for member 9...
2026-07-14 13:55:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:32 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 13:55:32 INFO Found: ['5156139']
2026-07-14 13:55:37 INFO [TGCC-IRENE] Submitted job with ID:['5156139']
2026-07-14 13:55:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 13:55:37 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-14 13:55:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 13:55:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:37 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 13:55:37 INFO Queuing job for member 10...
2026-07-14 13:55:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:37 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 13:55:38 INFO Found: ['5156140']
2026-07-14 13:55:43 INFO [TGCC-IRENE] Submitted job with ID:['5156140']
2026-07-14 13:55:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 13:55:43 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-14 13:55:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 13:55:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:43 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 13:55:43 INFO Queuing job for member 11...
2026-07-14 13:55:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:43 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 13:55:44 INFO Found: ['5156141']
2026-07-14 13:55:49 INFO [TGCC-IRENE] Submitted job with ID:['5156141']
2026-07-14 13:55:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 13:55:49 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-14 13:55:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 13:55:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:49 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 13:55:49 INFO Queuing job for member 12...
2026-07-14 13:55:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:49 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 13:55:50 INFO Found: ['5156142']
2026-07-14 13:55:55 INFO [TGCC-IRENE] Submitted job with ID:['5156142']
2026-07-14 13:55:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:55:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 13:55:55 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-14 13:55:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 13:55:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:55:55 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 13:55:55 INFO Queuing job for member 13...
2026-07-14 13:55:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:55:55 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 13:55:55 INFO Found: ['5156144']
2026-07-14 13:56:00 INFO [TGCC-IRENE] Submitted job with ID:['5156144']
2026-07-14 13:56:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:56:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 13:56:00 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-14 13:58:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 13:58:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:58:08 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 13:58:08 INFO Queuing job for member 14...
2026-07-14 13:58:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:58:08 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 13:58:09 INFO Found: ['5156157']
2026-07-14 13:58:14 INFO [TGCC-IRENE] Submitted job with ID:['5156157']
2026-07-14 13:58:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 13:58:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 13:58:14 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-14 13:58:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 13:58:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 13:58:14 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 13:58:14 INFO Queuing job for member 15...
2026-07-14 13:58:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 13:58:14 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 13:58:15 INFO Found: ['5156158']
2026-07-14 13:58:20 INFO [TGCC-IRENE] Submitted job with ID:['5156158']
2026-07-14 13:58:20 INFO Checking job status ...
2026-07-14 13:58:20 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:58:20 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:58:20 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:58:35 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:58:35 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:58:35 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:58:50 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:58:50 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:58:51 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:58:51 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:59:06 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:59:06 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:59:06 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:59:21 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:59:21 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:59:21 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:59:36 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:59:36 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:59:36 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 13:59:51 INFO None 5156111: status RUNNING/PENDING
2026-07-14 13:59:52 INFO None 5156112: status RUNNING/PENDING
2026-07-14 13:59:52 INFO None 5156113: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156115: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156134: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156135: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156136: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156137: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156139: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156140: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156141: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156142: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156144: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156157: status RUNNING/PENDING
2026-07-14 13:59:54 INFO None 5156158: status RUNNING/PENDING
2026-07-14 13:59:54 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:00:09 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156113: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156115: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:00:09 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:00:09 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:00:24 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156113: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156115: status FINISHED
2026-07-14 14:00:24 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:00:24 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:00:24 INFO Jobs still running: ['5156111', '5156112', '5156113', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:00:39 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:00:39 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:00:39 INFO None 5156113: status FINISHED
2026-07-14 14:00:39 INFO None 5156115: status FINISHED
2026-07-14 14:00:39 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:00:40 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:00:40 INFO Jobs still running: ['5156111', '5156112', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:00:55 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156113: status FINISHED
2026-07-14 14:00:55 INFO None 5156115: status FINISHED
2026-07-14 14:00:55 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:00:55 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:00:55 INFO Jobs still running: ['5156111', '5156112', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:01:10 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156113: status FINISHED
2026-07-14 14:01:10 INFO None 5156115: status FINISHED
2026-07-14 14:01:10 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:01:10 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:01:10 INFO Jobs still running: ['5156111', '5156112', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:01:25 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:01:25 INFO None 5156112: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156113: status FINISHED
2026-07-14 14:03:10 INFO None 5156115: status FINISHED
2026-07-14 14:03:10 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:03:10 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:03:10 INFO Jobs still running: ['5156111', '5156112', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:03:25 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:03:25 INFO None 5156112: status FINISHED
2026-07-14 14:03:25 INFO None 5156113: status FINISHED
2026-07-14 14:03:25 INFO None 5156115: status FINISHED
2026-07-14 14:03:26 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:03:26 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:03:26 INFO Jobs still running: ['5156111', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:03:41 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156112: status FINISHED
2026-07-14 14:03:41 INFO None 5156113: status FINISHED
2026-07-14 14:03:41 INFO None 5156115: status FINISHED
2026-07-14 14:03:41 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:03:41 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:03:41 INFO Jobs still running: ['5156111', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:03:56 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156112: status FINISHED
2026-07-14 14:03:56 INFO None 5156113: status FINISHED
2026-07-14 14:03:56 INFO None 5156115: status FINISHED
2026-07-14 14:03:56 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:03:56 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:03:57 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:03:57 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:03:57 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:03:57 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:03:57 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:03:57 INFO Jobs still running: ['5156111', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:04:12 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:04:12 INFO None 5156112: status FINISHED
2026-07-14 14:04:12 INFO None 5156113: status FINISHED
2026-07-14 14:04:12 INFO None 5156115: status FINISHED
2026-07-14 14:04:12 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:04:13 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:04:13 INFO Jobs still running: ['5156111', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:04:28 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156112: status FINISHED
2026-07-14 14:04:28 INFO None 5156113: status FINISHED
2026-07-14 14:04:28 INFO None 5156115: status FINISHED
2026-07-14 14:04:28 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156135: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:04:28 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:04:28 INFO Jobs still running: ['5156111', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:04:44 INFO None 5156111: status RUNNING/PENDING
2026-07-14 14:04:44 INFO None 5156112: status FINISHED
2026-07-14 14:04:44 INFO None 5156113: status FINISHED
2026-07-14 14:04:44 INFO None 5156115: status FINISHED
2026-07-14 14:04:44 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:04:44 INFO None 5156135: status FINISHED
2026-07-14 14:04:44 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:04:45 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:04:45 INFO Jobs still running: ['5156111', '5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:05:00 INFO None 5156111: status FINISHED
2026-07-14 14:05:00 INFO None 5156112: status FINISHED
2026-07-14 14:05:00 INFO None 5156113: status FINISHED
2026-07-14 14:05:00 INFO None 5156115: status FINISHED
2026-07-14 14:05:00 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156135: status FINISHED
2026-07-14 14:05:00 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:05:00 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:05:00 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:05:15 INFO None 5156111: status FINISHED
2026-07-14 14:05:15 INFO None 5156112: status FINISHED
2026-07-14 14:05:15 INFO None 5156113: status FINISHED
2026-07-14 14:05:15 INFO None 5156115: status FINISHED
2026-07-14 14:05:15 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156135: status FINISHED
2026-07-14 14:05:15 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:05:15 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:05:15 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:05:30 INFO None 5156111: status FINISHED
2026-07-14 14:05:30 INFO None 5156112: status FINISHED
2026-07-14 14:05:30 INFO None 5156113: status FINISHED
2026-07-14 14:05:30 INFO None 5156115: status FINISHED
2026-07-14 14:05:30 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:05:30 INFO None 5156135: status FINISHED
2026-07-14 14:05:30 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:05:30 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:05:30 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:05:31 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:05:31 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:05:47 INFO None 5156111: status FINISHED
2026-07-14 14:05:47 INFO None 5156112: status FINISHED
2026-07-14 14:05:47 INFO None 5156113: status FINISHED
2026-07-14 14:05:47 INFO None 5156115: status FINISHED
2026-07-14 14:05:47 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156135: status FINISHED
2026-07-14 14:05:48 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:05:48 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:05:48 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:06:03 INFO None 5156111: status FINISHED
2026-07-14 14:06:03 INFO None 5156112: status FINISHED
2026-07-14 14:06:03 INFO None 5156113: status FINISHED
2026-07-14 14:06:03 INFO None 5156115: status FINISHED
2026-07-14 14:06:03 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156135: status FINISHED
2026-07-14 14:06:03 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:06:03 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:06:03 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:06:18 INFO None 5156111: status FINISHED
2026-07-14 14:06:18 INFO None 5156112: status FINISHED
2026-07-14 14:06:18 INFO None 5156113: status FINISHED
2026-07-14 14:06:18 INFO None 5156115: status FINISHED
2026-07-14 14:06:18 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156135: status FINISHED
2026-07-14 14:06:18 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:06:18 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:06:18 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:06:33 INFO None 5156111: status FINISHED
2026-07-14 14:06:33 INFO None 5156112: status FINISHED
2026-07-14 14:06:33 INFO None 5156113: status FINISHED
2026-07-14 14:06:33 INFO None 5156115: status FINISHED
2026-07-14 14:06:33 INFO None 5156134: status RUNNING/PENDING
2026-07-14 14:06:33 INFO None 5156135: status FINISHED
2026-07-14 14:06:33 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:06:33 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:06:33 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:06:33 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:06:33 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:06:34 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:06:34 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:06:34 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:06:34 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:06:34 INFO Jobs still running: ['5156134', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:06:50 INFO None 5156111: status FINISHED
2026-07-14 14:08:19 INFO None 5156112: status FINISHED
2026-07-14 14:08:19 INFO None 5156113: status FINISHED
2026-07-14 14:08:19 INFO None 5156115: status FINISHED
2026-07-14 14:08:19 INFO None 5156134: status FINISHED
2026-07-14 14:08:19 INFO None 5156135: status FINISHED
2026-07-14 14:08:19 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:08:19 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:08:19 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:08:34 INFO None 5156111: status FINISHED
2026-07-14 14:08:34 INFO None 5156112: status FINISHED
2026-07-14 14:08:34 INFO None 5156113: status FINISHED
2026-07-14 14:08:34 INFO None 5156115: status FINISHED
2026-07-14 14:08:34 INFO None 5156134: status FINISHED
2026-07-14 14:08:34 INFO None 5156135: status FINISHED
2026-07-14 14:08:34 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:08:34 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:08:34 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:08:34 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:08:35 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:08:35 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:08:35 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:08:35 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:08:37 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:08:37 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:08:52 INFO None 5156111: status FINISHED
2026-07-14 14:08:52 INFO None 5156112: status FINISHED
2026-07-14 14:08:52 INFO None 5156113: status FINISHED
2026-07-14 14:08:52 INFO None 5156115: status FINISHED
2026-07-14 14:08:52 INFO None 5156134: status FINISHED
2026-07-14 14:08:52 INFO None 5156135: status FINISHED
2026-07-14 14:08:52 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:08:52 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:08:52 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:09:07 INFO None 5156111: status FINISHED
2026-07-14 14:09:07 INFO None 5156112: status FINISHED
2026-07-14 14:09:07 INFO None 5156113: status FINISHED
2026-07-14 14:09:07 INFO None 5156115: status FINISHED
2026-07-14 14:09:07 INFO None 5156134: status FINISHED
2026-07-14 14:09:07 INFO None 5156135: status FINISHED
2026-07-14 14:09:07 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:09:07 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:09:07 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:09:22 INFO None 5156111: status FINISHED
2026-07-14 14:09:22 INFO None 5156112: status FINISHED
2026-07-14 14:09:22 INFO None 5156113: status FINISHED
2026-07-14 14:09:22 INFO None 5156115: status FINISHED
2026-07-14 14:09:22 INFO None 5156134: status FINISHED
2026-07-14 14:09:22 INFO None 5156135: status FINISHED
2026-07-14 14:09:22 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:09:22 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:09:23 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:09:23 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:09:23 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:09:39 INFO None 5156111: status FINISHED
2026-07-14 14:09:39 INFO None 5156112: status FINISHED
2026-07-14 14:09:39 INFO None 5156113: status FINISHED
2026-07-14 14:09:39 INFO None 5156115: status FINISHED
2026-07-14 14:09:40 INFO None 5156134: status FINISHED
2026-07-14 14:09:40 INFO None 5156135: status FINISHED
2026-07-14 14:09:40 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:09:40 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:09:40 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:09:55 INFO None 5156111: status FINISHED
2026-07-14 14:09:55 INFO None 5156112: status FINISHED
2026-07-14 14:09:55 INFO None 5156113: status FINISHED
2026-07-14 14:09:55 INFO None 5156115: status FINISHED
2026-07-14 14:09:55 INFO None 5156134: status FINISHED
2026-07-14 14:09:55 INFO None 5156135: status FINISHED
2026-07-14 14:09:55 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:09:55 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:09:55 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:10:10 INFO None 5156111: status FINISHED
2026-07-14 14:10:10 INFO None 5156112: status FINISHED
2026-07-14 14:10:10 INFO None 5156113: status FINISHED
2026-07-14 14:10:10 INFO None 5156115: status FINISHED
2026-07-14 14:10:10 INFO None 5156134: status FINISHED
2026-07-14 14:10:10 INFO None 5156135: status FINISHED
2026-07-14 14:10:10 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:10:10 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:10:10 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:10:25 INFO None 5156111: status FINISHED
2026-07-14 14:10:25 INFO None 5156112: status FINISHED
2026-07-14 14:10:25 INFO None 5156113: status FINISHED
2026-07-14 14:10:25 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:10:25 INFO None 5156134: status FINISHED
2026-07-14 14:10:26 INFO None 5156135: status FINISHED
2026-07-14 14:10:26 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:10:26 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:10:26 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:10:41 INFO None 5156111: status FINISHED
2026-07-14 14:10:41 INFO None 5156112: status FINISHED
2026-07-14 14:10:41 INFO None 5156113: status FINISHED
2026-07-14 14:10:41 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:10:41 INFO None 5156134: status FINISHED
2026-07-14 14:10:41 INFO None 5156135: status FINISHED
2026-07-14 14:10:41 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:10:41 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:10:43 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:10:43 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:10:58 INFO None 5156111: status FINISHED
2026-07-14 14:10:58 INFO None 5156112: status FINISHED
2026-07-14 14:10:58 INFO None 5156113: status FINISHED
2026-07-14 14:10:58 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:10:58 INFO None 5156134: status FINISHED
2026-07-14 14:10:58 INFO None 5156135: status FINISHED
2026-07-14 14:10:58 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156140: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156141: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156142: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156144: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156157: status RUNNING/PENDING
2026-07-14 14:10:58 INFO None 5156158: status RUNNING/PENDING
2026-07-14 14:10:58 INFO Jobs still running: ['5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158']. Waiting...
2026-07-14 14:13:08 INFO None 5156111: status FINISHED
2026-07-14 14:13:08 INFO None 5156112: status FINISHED
2026-07-14 14:13:08 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:13:08 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:13:08 INFO None 5156134: status FINISHED
2026-07-14 14:13:08 INFO None 5156135: status FINISHED
2026-07-14 14:13:08 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:13:08 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:13:08 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:13:08 INFO None 5156140: status FINISHED
2026-07-14 14:13:08 INFO None 5156141: status FINISHED
2026-07-14 14:13:08 INFO None 5156142: status FINISHED
2026-07-14 14:13:08 INFO None 5156144: status FINISHED
2026-07-14 14:13:08 INFO None 5156157: status FINISHED
2026-07-14 14:13:08 INFO None 5156158: status FINISHED
2026-07-14 14:13:08 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:13:23 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:13:24 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:13:24 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:13:24 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:13:24 INFO None 5156134: status FINISHED
2026-07-14 14:13:24 INFO None 5156135: status FINISHED
2026-07-14 14:13:24 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:13:24 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:13:24 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:13:24 INFO None 5156140: status FINISHED
2026-07-14 14:13:24 INFO None 5156141: status FINISHED
2026-07-14 14:13:24 INFO None 5156142: status FINISHED
2026-07-14 14:13:24 INFO None 5156144: status FINISHED
2026-07-14 14:13:24 INFO None 5156157: status FINISHED
2026-07-14 14:13:24 INFO None 5156158: status FINISHED
2026-07-14 14:13:24 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:13:39 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:13:39 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:13:39 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:13:39 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:13:39 INFO None 5156134: status FINISHED
2026-07-14 14:13:39 INFO None 5156135: status FINISHED
2026-07-14 14:13:39 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:13:39 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:13:39 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:13:40 INFO None 5156140: status FINISHED
2026-07-14 14:13:40 INFO None 5156141: status FINISHED
2026-07-14 14:13:40 INFO None 5156142: status FINISHED
2026-07-14 14:13:40 INFO None 5156144: status FINISHED
2026-07-14 14:13:40 INFO None 5156157: status FINISHED
2026-07-14 14:13:40 INFO None 5156158: status FINISHED
2026-07-14 14:13:40 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:13:55 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:13:55 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:13:55 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:13:55 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:13:55 INFO None 5156134: status FINISHED
2026-07-14 14:13:55 INFO None 5156135: status FINISHED
2026-07-14 14:13:55 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:13:55 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:13:55 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:13:55 INFO None 5156140: status FINISHED
2026-07-14 14:13:55 INFO None 5156141: status FINISHED
2026-07-14 14:13:55 INFO None 5156142: status FINISHED
2026-07-14 14:13:55 INFO None 5156144: status FINISHED
2026-07-14 14:13:55 INFO None 5156157: status FINISHED
2026-07-14 14:13:55 INFO None 5156158: status FINISHED
2026-07-14 14:13:55 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:14:10 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:14:10 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:14:10 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:14:10 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:14:10 INFO None 5156134: status FINISHED
2026-07-14 14:14:10 INFO None 5156135: status FINISHED
2026-07-14 14:14:11 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:14:11 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:14:11 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:14:11 INFO None 5156140: status FINISHED
2026-07-14 14:14:11 INFO None 5156141: status FINISHED
2026-07-14 14:14:11 INFO None 5156142: status FINISHED
2026-07-14 14:14:11 INFO None 5156144: status FINISHED
2026-07-14 14:14:11 INFO None 5156157: status FINISHED
2026-07-14 14:14:11 INFO None 5156158: status FINISHED
2026-07-14 14:14:11 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:14:27 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:14:27 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:14:27 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:14:27 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:14:27 INFO None 5156134: status FINISHED
2026-07-14 14:14:27 INFO None 5156135: status FINISHED
2026-07-14 14:14:27 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:14:27 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:14:27 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:14:27 INFO None 5156140: status FINISHED
2026-07-14 14:14:27 INFO None 5156141: status FINISHED
2026-07-14 14:14:27 INFO None 5156142: status FINISHED
2026-07-14 14:14:27 INFO None 5156144: status FINISHED
2026-07-14 14:14:27 INFO None 5156157: status FINISHED
2026-07-14 14:14:27 INFO None 5156158: status FINISHED
2026-07-14 14:14:27 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:14:42 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:14:42 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:14:42 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:14:42 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:14:42 INFO None 5156134: status FINISHED
2026-07-14 14:14:42 INFO None 5156135: status FINISHED
2026-07-14 14:14:42 INFO None 5156136: status RUNNING/PENDING
2026-07-14 14:14:42 INFO None 5156137: status RUNNING/PENDING
2026-07-14 14:14:42 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:14:42 INFO None 5156140: status FINISHED
2026-07-14 14:14:43 INFO None 5156141: status FINISHED
2026-07-14 14:14:43 INFO None 5156142: status FINISHED
2026-07-14 14:14:43 INFO None 5156144: status FINISHED
2026-07-14 14:14:43 INFO None 5156157: status FINISHED
2026-07-14 14:14:43 INFO None 5156158: status FINISHED
2026-07-14 14:14:43 INFO Jobs still running: ['5156136', '5156137', '5156139']. Waiting...
2026-07-14 14:14:58 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:14:58 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:14:58 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:14:58 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:14:58 INFO None 5156134: status FINISHED
2026-07-14 14:14:58 INFO None 5156135: status FINISHED
2026-07-14 14:14:58 INFO None 5156136: status FINISHED
2026-07-14 14:14:58 INFO None 5156137: status FINISHED
2026-07-14 14:14:58 INFO None 5156139: status RUNNING/PENDING
2026-07-14 14:14:58 INFO None 5156140: status FINISHED
2026-07-14 14:14:58 INFO None 5156141: status FINISHED
2026-07-14 14:14:58 INFO None 5156142: status FINISHED
2026-07-14 14:14:58 INFO None 5156144: status FINISHED
2026-07-14 14:14:58 INFO None 5156157: status FINISHED
2026-07-14 14:14:58 INFO None 5156158: status FINISHED
2026-07-14 14:14:58 INFO Jobs still running: ['5156139']. Waiting...
2026-07-14 14:15:13 INFO None 5156111: status FINISHED (not in squeue)
2026-07-14 14:15:13 INFO None 5156112: status FINISHED (not in squeue)
2026-07-14 14:15:13 INFO None 5156113: status FINISHED (not in squeue)
2026-07-14 14:15:13 INFO None 5156115: status FINISHED (not in squeue)
2026-07-14 14:15:13 INFO None 5156134: status FINISHED
2026-07-14 14:15:13 INFO None 5156135: status FINISHED
2026-07-14 14:15:13 INFO None 5156136: status FINISHED
2026-07-14 14:15:13 INFO None 5156137: status FINISHED
2026-07-14 14:15:13 INFO None 5156139: status FINISHED
2026-07-14 14:15:13 INFO None 5156140: status FINISHED
2026-07-14 14:15:13 INFO None 5156141: status FINISHED
2026-07-14 14:15:13 INFO None 5156142: status FINISHED
2026-07-14 14:15:13 INFO None 5156144: status FINISHED
2026-07-14 14:15:13 INFO None 5156157: status FINISHED
2026-07-14 14:15:13 INFO None 5156158: status FINISHED
2026-07-14 14:15:13 INFO Jobs ['5156111', '5156112', '5156113', '5156115', '5156134', '5156135', '5156136', '5156137', '5156139', '5156140', '5156141', '5156142', '5156144', '5156157', '5156158'] have finished
2026-07-14 14:15:13 INFO Checking restart files were created ...
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-14 14:15:13 INFO  Run_model() completed successfully.
2026-07-14 14:15:13 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:15:13 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-14 14:15:13 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 14:15:13 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-14 14:15:13 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:15:13 INFO ---------->>> Running process_satellite_data()
2026-07-14 14:15:13 INFO [DART] No satellite data found, skipping assimilation
2026-07-14 14:15:13 INFO after_assimilation() skipped
2026-07-14 14:15:13 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 14:15:13 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:15:13 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:15:13 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-14 14:15:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:15 INFO Hourly dataset computed and listing created
2026-07-14 14:15:22 INFO Hourly dataset computed
2026-07-14 14:15:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:25 INFO Hourly dataset computed and listing created
2026-07-14 14:15:30 INFO Hourly dataset computed
2026-07-14 14:15:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:31 INFO Hourly dataset computed and listing created
2026-07-14 14:15:32 INFO Hourly dataset computed
2026-07-14 14:15:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:33 INFO Hourly dataset computed and listing created
2026-07-14 14:15:35 INFO Hourly dataset computed
2026-07-14 14:15:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:36 INFO Hourly dataset computed and listing created
2026-07-14 14:15:37 INFO Hourly dataset computed
2026-07-14 14:15:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:39 INFO Hourly dataset computed and listing created
2026-07-14 14:15:40 INFO Hourly dataset computed
2026-07-14 14:15:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:41 INFO Hourly dataset computed and listing created
2026-07-14 14:15:42 INFO Hourly dataset computed
2026-07-14 14:15:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:43 INFO Hourly dataset computed and listing created
2026-07-14 14:15:45 INFO Hourly dataset computed
2026-07-14 14:15:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:46 INFO Hourly dataset computed and listing created
2026-07-14 14:15:47 INFO Hourly dataset computed
2026-07-14 14:15:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:48 INFO Hourly dataset computed and listing created
2026-07-14 14:15:50 INFO Hourly dataset computed
2026-07-14 14:15:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:51 INFO Hourly dataset computed and listing created
2026-07-14 14:15:52 INFO Hourly dataset computed
2026-07-14 14:15:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:53 INFO Hourly dataset computed and listing created
2026-07-14 14:15:54 INFO Hourly dataset computed
2026-07-14 14:15:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:56 INFO Hourly dataset computed and listing created
2026-07-14 14:15:57 INFO Hourly dataset computed
2026-07-14 14:15:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:15:58 INFO Hourly dataset computed and listing created
2026-07-14 14:15:59 INFO Hourly dataset computed
2026-07-14 14:15:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:18:08 INFO Hourly dataset computed and listing created
2026-07-14 14:19:03 INFO Hourly dataset computed
2026-07-14 14:19:03 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-14 14:19:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 14:19:03 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-14 14:19:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 14:19:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:03 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 14:19:03 INFO Queuing job for member 1...
2026-07-14 14:19:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:03 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 14:19:04 INFO Found: ['5156217']
2026-07-14 14:19:09 INFO [TGCC-IRENE] Submitted job with ID:['5156217']
2026-07-14 14:19:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 14:19:09 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-14 14:19:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 14:19:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:09 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 14:19:09 INFO Queuing job for member 2...
2026-07-14 14:19:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:09 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 14:19:10 INFO Found: ['5156219']
2026-07-14 14:19:15 INFO [TGCC-IRENE] Submitted job with ID:['5156219']
2026-07-14 14:19:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 14:19:15 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-14 14:19:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 14:19:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:15 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 14:19:15 INFO Queuing job for member 3...
2026-07-14 14:19:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:15 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 14:19:17 INFO Found: ['5156221']
2026-07-14 14:19:22 INFO [TGCC-IRENE] Submitted job with ID:['5156221']
2026-07-14 14:19:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 14:19:22 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-14 14:19:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 14:19:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:22 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 14:19:22 INFO Queuing job for member 4...
2026-07-14 14:19:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:22 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 14:19:25 INFO Found: ['5156222']
2026-07-14 14:19:30 INFO [TGCC-IRENE] Submitted job with ID:['5156222']
2026-07-14 14:19:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 14:19:30 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-14 14:19:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 14:19:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:30 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 14:19:30 INFO Queuing job for member 5...
2026-07-14 14:19:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:30 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 14:19:31 INFO Found: ['5156223']
2026-07-14 14:19:36 INFO [TGCC-IRENE] Submitted job with ID:['5156223']
2026-07-14 14:19:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 14:19:36 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-14 14:19:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 14:19:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:36 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 14:19:36 INFO Queuing job for member 6...
2026-07-14 14:19:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:36 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 14:19:36 INFO Found: ['5156224']
2026-07-14 14:19:41 INFO [TGCC-IRENE] Submitted job with ID:['5156224']
2026-07-14 14:19:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 14:19:41 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-14 14:19:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 14:19:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:41 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 14:19:41 INFO Queuing job for member 7...
2026-07-14 14:19:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:42 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 14:19:42 INFO Found: ['5156225']
2026-07-14 14:19:47 INFO [TGCC-IRENE] Submitted job with ID:['5156225']
2026-07-14 14:19:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 14:19:47 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-14 14:19:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 14:19:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:47 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 14:19:47 INFO Queuing job for member 8...
2026-07-14 14:19:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:47 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 14:19:48 INFO Found: ['5156226']
2026-07-14 14:19:53 INFO [TGCC-IRENE] Submitted job with ID:['5156226']
2026-07-14 14:19:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 14:19:53 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-14 14:19:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 14:19:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:53 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 14:19:53 INFO Queuing job for member 9...
2026-07-14 14:19:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:53 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 14:19:54 INFO Found: ['5156227']
2026-07-14 14:19:59 INFO [TGCC-IRENE] Submitted job with ID:['5156227']
2026-07-14 14:19:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:19:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 14:19:59 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-14 14:19:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 14:19:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:19:59 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 14:19:59 INFO Queuing job for member 10...
2026-07-14 14:19:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:19:59 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 14:20:00 INFO Found: ['5156228']
2026-07-14 14:20:05 INFO [TGCC-IRENE] Submitted job with ID:['5156228']
2026-07-14 14:20:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:20:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 14:20:05 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-14 14:20:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 14:20:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:20:05 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 14:20:05 INFO Queuing job for member 11...
2026-07-14 14:20:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:20:05 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 14:20:05 INFO Found: ['5156230']
2026-07-14 14:20:10 INFO [TGCC-IRENE] Submitted job with ID:['5156230']
2026-07-14 14:20:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:20:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 14:20:10 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-14 14:20:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 14:20:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:20:10 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 14:20:10 INFO Queuing job for member 12...
2026-07-14 14:20:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:20:10 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 14:20:11 INFO Found: ['5156231']
2026-07-14 14:20:16 INFO [TGCC-IRENE] Submitted job with ID:['5156231']
2026-07-14 14:20:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:20:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 14:20:16 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-14 14:20:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 14:20:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:20:16 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 14:20:16 INFO Queuing job for member 13...
2026-07-14 14:20:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:20:16 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 14:20:18 INFO Found: ['5156232']
2026-07-14 14:20:23 INFO [TGCC-IRENE] Submitted job with ID:['5156232']
2026-07-14 14:20:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:20:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 14:20:23 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-14 14:20:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 14:20:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:20:23 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 14:20:23 INFO Queuing job for member 14...
2026-07-14 14:20:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:20:23 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 14:20:25 INFO Found: ['5156233']
2026-07-14 14:20:30 INFO [TGCC-IRENE] Submitted job with ID:['5156233']
2026-07-14 14:20:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:20:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 14:20:30 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-14 14:20:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 14:20:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:20:30 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 14:20:30 INFO Queuing job for member 15...
2026-07-14 14:20:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:20:30 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 14:20:31 INFO Found: ['5156235']
2026-07-14 14:20:36 INFO [TGCC-IRENE] Submitted job with ID:['5156235']
2026-07-14 14:20:36 INFO Checking job status ...
2026-07-14 14:20:36 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156219: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156225: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156226: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156227: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156228: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156231: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:20:36 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:20:36 INFO Jobs still running: ['5156217', '5156219', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:20:51 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156219: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156225: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156226: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156227: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156228: status RUNNING/PENDING
2026-07-14 14:20:51 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:20:52 INFO None 5156231: status RUNNING/PENDING
2026-07-14 14:20:52 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:20:52 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:20:52 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:20:52 INFO Jobs still running: ['5156217', '5156219', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:23:10 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156219: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156225: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156226: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156227: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156228: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156231: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:23:10 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:23:10 INFO Jobs still running: ['5156217', '5156219', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:23:25 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156219: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156225: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156226: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156227: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156228: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156231: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:23:25 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:23:25 INFO Jobs still running: ['5156217', '5156219', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:23:40 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:23:40 INFO None 5156219: status FINISHED
2026-07-14 14:23:40 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:23:40 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:23:40 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:23:40 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156225: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156226: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156227: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156228: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156231: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:23:41 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:23:41 INFO Jobs still running: ['5156217', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:23:56 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156219: status FINISHED
2026-07-14 14:23:56 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156225: status FINISHED
2026-07-14 14:23:56 INFO None 5156226: status FINISHED
2026-07-14 14:23:56 INFO None 5156227: status FINISHED
2026-07-14 14:23:56 INFO None 5156228: status FINISHED
2026-07-14 14:23:56 INFO None 5156230: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156231: status FINISHED
2026-07-14 14:23:56 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:23:56 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:23:56 INFO Jobs still running: ['5156217', '5156221', '5156222', '5156223', '5156224', '5156230', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:24:11 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156219: status FINISHED
2026-07-14 14:24:13 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156225: status FINISHED
2026-07-14 14:24:13 INFO None 5156226: status FINISHED
2026-07-14 14:24:13 INFO None 5156227: status FINISHED
2026-07-14 14:24:13 INFO None 5156228: status FINISHED
2026-07-14 14:24:13 INFO None 5156230: status FINISHED
2026-07-14 14:24:13 INFO None 5156231: status FINISHED
2026-07-14 14:24:13 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:24:13 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:24:13 INFO Jobs still running: ['5156217', '5156221', '5156222', '5156223', '5156224', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:24:28 INFO None 5156217: status RUNNING/PENDING
2026-07-14 14:24:28 INFO None 5156219: status FINISHED
2026-07-14 14:24:28 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:24:28 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:24:28 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:24:28 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:24:28 INFO None 5156225: status FINISHED
2026-07-14 14:24:28 INFO None 5156226: status FINISHED
2026-07-14 14:24:28 INFO None 5156227: status FINISHED
2026-07-14 14:24:29 INFO None 5156228: status FINISHED
2026-07-14 14:24:29 INFO None 5156230: status FINISHED
2026-07-14 14:24:29 INFO None 5156231: status FINISHED
2026-07-14 14:24:29 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:24:29 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:24:29 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:24:29 INFO Jobs still running: ['5156217', '5156221', '5156222', '5156223', '5156224', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:24:44 INFO None 5156217: status FINISHED
2026-07-14 14:24:44 INFO None 5156219: status FINISHED
2026-07-14 14:24:44 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156222: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156223: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156224: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156225: status FINISHED
2026-07-14 14:24:44 INFO None 5156226: status FINISHED
2026-07-14 14:24:44 INFO None 5156227: status FINISHED
2026-07-14 14:24:44 INFO None 5156228: status FINISHED
2026-07-14 14:24:44 INFO None 5156230: status FINISHED
2026-07-14 14:24:44 INFO None 5156231: status FINISHED
2026-07-14 14:24:44 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:24:44 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:24:44 INFO Jobs still running: ['5156221', '5156222', '5156223', '5156224', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:24:59 INFO None 5156217: status FINISHED
2026-07-14 14:24:59 INFO None 5156219: status FINISHED
2026-07-14 14:24:59 INFO None 5156221: status RUNNING/PENDING
2026-07-14 14:24:59 INFO None 5156222: status FINISHED
2026-07-14 14:24:59 INFO None 5156223: status FINISHED
2026-07-14 14:24:59 INFO None 5156224: status FINISHED
2026-07-14 14:24:59 INFO None 5156225: status FINISHED
2026-07-14 14:24:59 INFO None 5156226: status FINISHED
2026-07-14 14:24:59 INFO None 5156227: status FINISHED
2026-07-14 14:24:59 INFO None 5156228: status FINISHED
2026-07-14 14:24:59 INFO None 5156230: status FINISHED
2026-07-14 14:24:59 INFO None 5156231: status FINISHED
2026-07-14 14:24:59 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:24:59 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:24:59 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:24:59 INFO Jobs still running: ['5156221', '5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:25:16 INFO None 5156217: status FINISHED
2026-07-14 14:25:16 INFO None 5156219: status FINISHED
2026-07-14 14:25:16 INFO None 5156221: status FINISHED
2026-07-14 14:25:16 INFO None 5156222: status FINISHED
2026-07-14 14:25:16 INFO None 5156223: status FINISHED
2026-07-14 14:25:16 INFO None 5156224: status FINISHED
2026-07-14 14:25:16 INFO None 5156225: status FINISHED
2026-07-14 14:25:16 INFO None 5156226: status FINISHED
2026-07-14 14:25:16 INFO None 5156227: status FINISHED
2026-07-14 14:25:16 INFO None 5156228: status FINISHED
2026-07-14 14:25:16 INFO None 5156230: status FINISHED
2026-07-14 14:25:16 INFO None 5156231: status FINISHED
2026-07-14 14:25:16 INFO None 5156232: status RUNNING/PENDING
2026-07-14 14:25:16 INFO None 5156233: status RUNNING/PENDING
2026-07-14 14:25:16 INFO None 5156235: status RUNNING/PENDING
2026-07-14 14:25:16 INFO Jobs still running: ['5156232', '5156233', '5156235']. Waiting...
2026-07-14 14:25:31 INFO None 5156217: status FINISHED
2026-07-14 14:25:31 INFO None 5156219: status FINISHED
2026-07-14 14:25:31 INFO None 5156221: status FINISHED
2026-07-14 14:25:31 INFO None 5156222: status FINISHED
2026-07-14 14:25:31 INFO None 5156223: status FINISHED
2026-07-14 14:25:31 INFO None 5156224: status FINISHED
2026-07-14 14:25:31 INFO None 5156225: status FINISHED
2026-07-14 14:25:32 INFO None 5156226: status FINISHED
2026-07-14 14:25:32 INFO None 5156227: status FINISHED
2026-07-14 14:25:32 INFO None 5156228: status FINISHED
2026-07-14 14:25:32 INFO None 5156230: status FINISHED
2026-07-14 14:25:32 INFO None 5156231: status FINISHED
2026-07-14 14:25:32 INFO None 5156232: status FINISHED
2026-07-14 14:25:32 INFO None 5156233: status FINISHED
2026-07-14 14:25:32 INFO None 5156235: status FINISHED
2026-07-14 14:25:32 INFO Jobs ['5156217', '5156219', '5156221', '5156222', '5156223', '5156224', '5156225', '5156226', '5156227', '5156228', '5156230', '5156231', '5156232', '5156233', '5156235'] have finished
2026-07-14 14:25:32 INFO Checking restart files were created ...
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-14 14:25:32 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-14 14:25:32 INFO  Run_model() completed successfully.
2026-07-14 14:25:32 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:25:32 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-14 14:25:32 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 14:25:32 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-14 14:25:32 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:25:32 INFO ---------->>> Running process_satellite_data()
2026-07-14 14:25:32 INFO [DART] No satellite data found, skipping assimilation
2026-07-14 14:25:32 INFO after_assimilation() skipped
2026-07-14 14:25:32 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 14:25:32 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:25:32 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:25:32 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-14 14:25:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:25:34 INFO Hourly dataset computed and listing created
2026-07-14 14:25:50 INFO Hourly dataset computed
2026-07-14 14:25:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:25:52 INFO Hourly dataset computed and listing created
2026-07-14 14:26:01 INFO Hourly dataset computed
2026-07-14 14:26:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:28:10 INFO Hourly dataset computed and listing created
2026-07-14 14:28:25 INFO Hourly dataset computed
2026-07-14 14:28:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:28:26 INFO Hourly dataset computed and listing created
2026-07-14 14:28:37 INFO Hourly dataset computed
2026-07-14 14:28:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:28:38 INFO Hourly dataset computed and listing created
2026-07-14 14:28:47 INFO Hourly dataset computed
2026-07-14 14:28:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:28:48 INFO Hourly dataset computed and listing created
2026-07-14 14:28:57 INFO Hourly dataset computed
2026-07-14 14:28:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:28:58 INFO Hourly dataset computed and listing created
2026-07-14 14:29:08 INFO Hourly dataset computed
2026-07-14 14:29:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:29:09 INFO Hourly dataset computed and listing created
2026-07-14 14:29:20 INFO Hourly dataset computed
2026-07-14 14:29:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:29:21 INFO Hourly dataset computed and listing created
2026-07-14 14:29:31 INFO Hourly dataset computed
2026-07-14 14:29:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:29:32 INFO Hourly dataset computed and listing created
2026-07-14 14:29:42 INFO Hourly dataset computed
2026-07-14 14:29:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:29:43 INFO Hourly dataset computed and listing created
2026-07-14 14:29:52 INFO Hourly dataset computed
2026-07-14 14:29:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:29:53 INFO Hourly dataset computed and listing created
2026-07-14 14:30:02 INFO Hourly dataset computed
2026-07-14 14:30:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:30:03 INFO Hourly dataset computed and listing created
2026-07-14 14:30:14 INFO Hourly dataset computed
2026-07-14 14:30:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:30:15 INFO Hourly dataset computed and listing created
2026-07-14 14:30:29 INFO Hourly dataset computed
2026-07-14 14:30:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:30:30 INFO Hourly dataset computed and listing created
2026-07-14 14:30:40 INFO Hourly dataset computed
2026-07-14 14:30:40 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-14 14:30:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:30:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 14:30:40 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-14 14:30:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 14:30:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:30:40 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 14:30:42 INFO Queuing job for member 1...
2026-07-14 14:30:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:30:42 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 14:30:43 INFO Found: ['5156253']
2026-07-14 14:30:48 INFO [TGCC-IRENE] Submitted job with ID:['5156253']
2026-07-14 14:30:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:30:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 14:30:48 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-14 14:30:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 14:30:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:30:48 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 14:30:48 INFO Queuing job for member 2...
2026-07-14 14:30:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:30:48 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 14:30:48 INFO Found: ['5156254']
2026-07-14 14:30:53 INFO [TGCC-IRENE] Submitted job with ID:['5156254']
2026-07-14 14:30:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:30:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 14:30:53 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-14 14:30:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 14:30:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:30:53 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 14:30:53 INFO Queuing job for member 3...
2026-07-14 14:30:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:30:53 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 14:30:55 INFO Found: ['5156255']
2026-07-14 14:31:00 INFO [TGCC-IRENE] Submitted job with ID:['5156255']
2026-07-14 14:31:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:31:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 14:31:00 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-14 14:31:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 14:31:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:31:00 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 14:31:00 INFO Queuing job for member 4...
2026-07-14 14:31:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:31:00 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 14:31:03 INFO Found: ['5156257']
2026-07-14 14:31:08 INFO [TGCC-IRENE] Submitted job with ID:['5156257']
2026-07-14 14:31:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:31:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 14:31:08 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-14 14:33:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 14:33:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:16 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 14:33:16 INFO Queuing job for member 5...
2026-07-14 14:33:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:16 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 14:33:17 INFO Found: ['5156264']
2026-07-14 14:33:22 INFO [TGCC-IRENE] Submitted job with ID:['5156264']
2026-07-14 14:33:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 14:33:22 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-14 14:33:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 14:33:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:22 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 14:33:22 INFO Queuing job for member 6...
2026-07-14 14:33:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:22 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 14:33:23 INFO Found: ['5156265']
2026-07-14 14:33:28 INFO [TGCC-IRENE] Submitted job with ID:['5156265']
2026-07-14 14:33:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 14:33:28 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-14 14:33:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 14:33:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:28 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 14:33:28 INFO Queuing job for member 7...
2026-07-14 14:33:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:28 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 14:33:29 INFO Found: ['5156266']
2026-07-14 14:33:34 INFO [TGCC-IRENE] Submitted job with ID:['5156266']
2026-07-14 14:33:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 14:33:34 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-14 14:33:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 14:33:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:34 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 14:33:34 INFO Queuing job for member 8...
2026-07-14 14:33:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:34 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 14:33:35 INFO Found: ['5156267']
2026-07-14 14:33:40 INFO [TGCC-IRENE] Submitted job with ID:['5156267']
2026-07-14 14:33:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 14:33:40 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-14 14:33:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 14:33:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:40 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 14:33:40 INFO Queuing job for member 9...
2026-07-14 14:33:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:40 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 14:33:40 INFO Found: ['5156268']
2026-07-14 14:33:45 INFO [TGCC-IRENE] Submitted job with ID:['5156268']
2026-07-14 14:33:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 14:33:45 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-14 14:33:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 14:33:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:45 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 14:33:45 INFO Queuing job for member 10...
2026-07-14 14:33:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:45 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 14:33:46 INFO Found: ['5156269']
2026-07-14 14:33:51 INFO [TGCC-IRENE] Submitted job with ID:['5156269']
2026-07-14 14:33:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 14:33:51 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-14 14:33:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 14:33:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:51 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 14:33:51 INFO Queuing job for member 11...
2026-07-14 14:33:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:51 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 14:33:54 INFO Found: ['5156270']
2026-07-14 14:33:59 INFO [TGCC-IRENE] Submitted job with ID:['5156270']
2026-07-14 14:33:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:33:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 14:33:59 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-14 14:33:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 14:33:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:33:59 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 14:33:59 INFO Queuing job for member 12...
2026-07-14 14:33:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:33:59 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 14:33:59 INFO Found: ['5156271']
2026-07-14 14:34:04 INFO [TGCC-IRENE] Submitted job with ID:['5156271']
2026-07-14 14:34:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:34:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 14:34:04 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-14 14:34:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 14:34:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:34:04 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 14:34:04 INFO Queuing job for member 13...
2026-07-14 14:34:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:34:04 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 14:34:05 INFO Found: ['5156273']
2026-07-14 14:34:10 INFO [TGCC-IRENE] Submitted job with ID:['5156273']
2026-07-14 14:34:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:34:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 14:34:10 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-14 14:34:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 14:34:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:34:10 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 14:34:10 INFO Queuing job for member 14...
2026-07-14 14:34:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:34:10 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 14:34:11 INFO Found: ['5156274']
2026-07-14 14:34:16 INFO [TGCC-IRENE] Submitted job with ID:['5156274']
2026-07-14 14:34:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:34:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 14:34:16 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-14 14:34:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 14:34:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:34:16 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 14:34:16 INFO Queuing job for member 15...
2026-07-14 14:34:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:34:16 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 14:34:17 INFO Found: ['5156275']
2026-07-14 14:34:22 INFO [TGCC-IRENE] Submitted job with ID:['5156275']
2026-07-14 14:34:22 INFO Checking job status ...
2026-07-14 14:34:22 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:34:22 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:34:22 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:34:37 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:34:37 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:34:37 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:34:55 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:34:55 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:34:55 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:35:10 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:35:10 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:35:10 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:35:25 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:35:25 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:35:25 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:35:40 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:35:40 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:35:40 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:35:40 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:35:41 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:35:41 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:35:56 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:35:56 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:35:56 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:36:11 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:36:11 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:36:11 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:36:26 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:36:26 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:36:27 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:36:27 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:36:27 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:36:27 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:36:27 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:36:27 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:36:42 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:36:42 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:36:42 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:36:57 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:36:57 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:36:57 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:37:12 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:37:12 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:37:12 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:37:27 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:37:27 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:37:27 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:37:28 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:37:28 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:37:43 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:37:45 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:37:45 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:38:00 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:38:00 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:38:00 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:38:15 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:38:15 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:38:15 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:38:16 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:38:16 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:38:31 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:38:31 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:38:31 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:38:46 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:38:46 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:38:46 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:39:01 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:39:01 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:39:01 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:39:01 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:39:01 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:39:01 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:39:02 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:39:02 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:39:17 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:39:17 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:39:17 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:39:32 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:39:32 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:39:32 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:39:47 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:39:47 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:39:47 INFO None 5156255: status RUNNING/PENDING
2026-07-14 14:39:47 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:39:47 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:39:48 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:39:48 INFO Jobs still running: ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:40:03 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156255: status FINISHED
2026-07-14 14:40:03 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:40:03 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:40:04 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:40:04 INFO Jobs still running: ['5156253', '5156254', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:40:19 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156255: status FINISHED
2026-07-14 14:40:19 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:40:19 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:40:19 INFO Jobs still running: ['5156253', '5156254', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:40:35 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156255: status FINISHED
2026-07-14 14:40:35 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:40:35 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:40:35 INFO Jobs still running: ['5156253', '5156254', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:40:51 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156255: status FINISHED
2026-07-14 14:40:51 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:40:51 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:40:51 INFO Jobs still running: ['5156253', '5156254', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:41:06 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156255: status FINISHED
2026-07-14 14:41:06 INFO None 5156257: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:41:06 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:41:06 INFO Jobs still running: ['5156253', '5156254', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:41:21 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:43:21 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156255: status FINISHED
2026-07-14 14:43:57 INFO None 5156257: status FINISHED
2026-07-14 14:43:57 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:43:57 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:43:57 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:44:12 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156255: status FINISHED
2026-07-14 14:44:12 INFO None 5156257: status FINISHED
2026-07-14 14:44:12 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:44:12 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:44:12 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:44:28 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156255: status FINISHED
2026-07-14 14:44:28 INFO None 5156257: status FINISHED
2026-07-14 14:44:28 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:44:28 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:44:28 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:44:43 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156255: status FINISHED
2026-07-14 14:44:43 INFO None 5156257: status FINISHED
2026-07-14 14:44:43 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:44:43 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:44:43 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:44:59 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156255: status FINISHED
2026-07-14 14:44:59 INFO None 5156257: status FINISHED
2026-07-14 14:44:59 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:44:59 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:44:59 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:45:14 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156255: status FINISHED
2026-07-14 14:45:14 INFO None 5156257: status FINISHED
2026-07-14 14:45:14 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:45:14 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:45:14 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:45:30 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156255: status FINISHED
2026-07-14 14:45:30 INFO None 5156257: status FINISHED
2026-07-14 14:45:30 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:45:30 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:45:30 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:45:45 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156255: status FINISHED
2026-07-14 14:45:45 INFO None 5156257: status FINISHED
2026-07-14 14:45:45 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:45:45 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:45:45 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:46:00 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156254: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156255: status FINISHED
2026-07-14 14:46:00 INFO None 5156257: status FINISHED
2026-07-14 14:46:00 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156266: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156267: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156268: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156269: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156270: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156271: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:46:00 INFO None 5156275: status RUNNING/PENDING
2026-07-14 14:46:00 INFO Jobs still running: ['5156253', '5156254', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275']. Waiting...
2026-07-14 14:46:15 INFO None 5156253: status RUNNING/PENDING
2026-07-14 14:48:14 INFO None 5156254: status FINISHED
2026-07-14 14:48:14 INFO None 5156255: status FINISHED
2026-07-14 14:48:14 INFO None 5156257: status FINISHED
2026-07-14 14:48:14 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:48:14 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:48:14 INFO None 5156266: status FINISHED
2026-07-14 14:48:14 INFO None 5156267: status FINISHED
2026-07-14 14:48:14 INFO None 5156268: status FINISHED
2026-07-14 14:48:14 INFO None 5156269: status FINISHED
2026-07-14 14:48:14 INFO None 5156270: status FINISHED
2026-07-14 14:48:14 INFO None 5156271: status FINISHED
2026-07-14 14:48:14 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:48:14 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:48:14 INFO None 5156275: status FINISHED
2026-07-14 14:48:14 INFO Jobs still running: ['5156253', '5156264', '5156265', '5156273', '5156274']. Waiting...
2026-07-14 14:48:29 INFO None 5156253: status FINISHED
2026-07-14 14:48:29 INFO None 5156254: status FINISHED
2026-07-14 14:48:30 INFO None 5156255: status FINISHED
2026-07-14 14:48:30 INFO None 5156257: status FINISHED
2026-07-14 14:48:30 INFO None 5156264: status RUNNING/PENDING
2026-07-14 14:48:30 INFO None 5156265: status RUNNING/PENDING
2026-07-14 14:48:30 INFO None 5156266: status FINISHED
2026-07-14 14:48:30 INFO None 5156267: status FINISHED
2026-07-14 14:48:30 INFO None 5156268: status FINISHED
2026-07-14 14:48:30 INFO None 5156269: status FINISHED
2026-07-14 14:48:30 INFO None 5156270: status FINISHED
2026-07-14 14:48:30 INFO None 5156271: status FINISHED
2026-07-14 14:48:30 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:48:30 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:48:30 INFO None 5156275: status FINISHED
2026-07-14 14:48:30 INFO Jobs still running: ['5156264', '5156265', '5156273', '5156274']. Waiting...
2026-07-14 14:48:45 INFO None 5156253: status FINISHED
2026-07-14 14:48:45 INFO None 5156254: status FINISHED
2026-07-14 14:48:45 INFO None 5156255: status FINISHED
2026-07-14 14:48:45 INFO None 5156257: status FINISHED
2026-07-14 14:48:45 INFO None 5156264: status FINISHED
2026-07-14 14:48:45 INFO None 5156265: status FINISHED
2026-07-14 14:48:45 INFO None 5156266: status FINISHED
2026-07-14 14:48:45 INFO None 5156267: status FINISHED
2026-07-14 14:48:45 INFO None 5156268: status FINISHED
2026-07-14 14:48:45 INFO None 5156269: status FINISHED
2026-07-14 14:48:45 INFO None 5156270: status FINISHED
2026-07-14 14:48:45 INFO None 5156271: status FINISHED
2026-07-14 14:48:45 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:48:45 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:48:45 INFO None 5156275: status FINISHED
2026-07-14 14:48:45 INFO Jobs still running: ['5156273', '5156274']. Waiting...
2026-07-14 14:49:00 INFO None 5156253: status FINISHED
2026-07-14 14:49:00 INFO None 5156254: status FINISHED
2026-07-14 14:49:00 INFO None 5156255: status FINISHED
2026-07-14 14:49:00 INFO None 5156257: status FINISHED
2026-07-14 14:49:00 INFO None 5156264: status FINISHED
2026-07-14 14:49:00 INFO None 5156265: status FINISHED
2026-07-14 14:49:00 INFO None 5156266: status FINISHED
2026-07-14 14:49:00 INFO None 5156267: status FINISHED
2026-07-14 14:49:00 INFO None 5156268: status FINISHED
2026-07-14 14:49:00 INFO None 5156269: status FINISHED
2026-07-14 14:49:00 INFO None 5156270: status FINISHED
2026-07-14 14:49:00 INFO None 5156271: status FINISHED
2026-07-14 14:49:00 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:49:00 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:49:00 INFO None 5156275: status FINISHED
2026-07-14 14:49:00 INFO Jobs still running: ['5156273', '5156274']. Waiting...
2026-07-14 14:49:16 INFO None 5156253: status FINISHED
2026-07-14 14:49:16 INFO None 5156254: status FINISHED
2026-07-14 14:49:16 INFO None 5156255: status FINISHED
2026-07-14 14:49:16 INFO None 5156257: status FINISHED
2026-07-14 14:49:16 INFO None 5156264: status FINISHED
2026-07-14 14:49:16 INFO None 5156265: status FINISHED
2026-07-14 14:49:16 INFO None 5156266: status FINISHED
2026-07-14 14:49:16 INFO None 5156267: status FINISHED
2026-07-14 14:49:16 INFO None 5156268: status FINISHED
2026-07-14 14:49:16 INFO None 5156269: status FINISHED
2026-07-14 14:49:16 INFO None 5156270: status FINISHED
2026-07-14 14:49:16 INFO None 5156271: status FINISHED
2026-07-14 14:49:16 INFO None 5156273: status RUNNING/PENDING
2026-07-14 14:49:16 INFO None 5156274: status RUNNING/PENDING
2026-07-14 14:49:16 INFO None 5156275: status FINISHED
2026-07-14 14:49:16 INFO Jobs still running: ['5156273', '5156274']. Waiting...
2026-07-14 14:49:31 INFO None 5156253: status FINISHED
2026-07-14 14:49:31 INFO None 5156254: status FINISHED
2026-07-14 14:49:31 INFO None 5156255: status FINISHED
2026-07-14 14:49:31 INFO None 5156257: status FINISHED
2026-07-14 14:49:31 INFO None 5156264: status FINISHED
2026-07-14 14:49:31 INFO None 5156265: status FINISHED
2026-07-14 14:49:31 INFO None 5156266: status FINISHED
2026-07-14 14:49:31 INFO None 5156267: status FINISHED
2026-07-14 14:49:31 INFO None 5156268: status FINISHED
2026-07-14 14:49:31 INFO None 5156269: status FINISHED
2026-07-14 14:49:31 INFO None 5156270: status FINISHED
2026-07-14 14:49:31 INFO None 5156271: status FINISHED
2026-07-14 14:49:31 INFO None 5156273: status FINISHED
2026-07-14 14:49:31 INFO None 5156274: status FINISHED
2026-07-14 14:49:31 INFO None 5156275: status FINISHED
2026-07-14 14:49:31 INFO Jobs ['5156253', '5156254', '5156255', '5156257', '5156264', '5156265', '5156266', '5156267', '5156268', '5156269', '5156270', '5156271', '5156273', '5156274', '5156275'] have finished
2026-07-14 14:49:31 INFO Checking restart files were created ...
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-14 14:49:31 INFO  Run_model() completed successfully.
2026-07-14 14:49:31 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:49:31 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-14 14:49:31 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 14:49:31 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-14 14:49:31 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:49:31 INFO ---------->>> Running process_satellite_data()
2026-07-14 14:49:31 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-14 14:49:31 INFO ---------->>> Running run_obs_converter()
2026-07-14 14:49:31 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-14 14:49:31 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-14 14:49:31 INFO ---------->>> Running DART
2026-07-14 14:49:31 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 14:49:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-14 14:49:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-14 14:49:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-14 14:49:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-14 14:49:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-14 14:49:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-14 14:49:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-14 14:49:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-14 14:49:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-14 14:49:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-14 14:49:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-14 14:49:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-14 14:49:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-14 14:49:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-14 14:49:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-14 14:49:38 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 14:49:38 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 14:49:38 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 14:49:38 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 14:49:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 14:49:38 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 14:49:45 INFO Found: []
2026-07-14 14:49:45 INFO No job id returned by command ./run_filter.bsh
2026-07-14 14:49:45 INFO No monitoring will be performed
2026-07-14 14:49:45 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-14 14:49:45 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020709'
2026-07-14 14:49:45 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 14:49:49 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 14:49:49 INFO run_dart() is DONE.
2026-07-14 14:49:49 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 14:49:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:49 INFO No previous orbit memory found.
2026-07-14 14:49:49 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020709.nc
2026-07-14 14:49:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:50 INFO No previous orbit memory found.
2026-07-14 14:49:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020709.nc
2026-07-14 14:49:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:49:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:51 INFO No previous orbit memory found.
2026-07-14 14:49:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020709.nc
2026-07-14 14:49:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:52 INFO No previous orbit memory found.
2026-07-14 14:49:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020709.nc
2026-07-14 14:49:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:49:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:53 INFO No previous orbit memory found.
2026-07-14 14:49:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020709.nc
2026-07-14 14:49:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:54 INFO No previous orbit memory found.
2026-07-14 14:49:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020709.nc
2026-07-14 14:49:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:49:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:55 INFO No previous orbit memory found.
2026-07-14 14:49:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020709.nc
2026-07-14 14:49:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:56 INFO No previous orbit memory found.
2026-07-14 14:49:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020709.nc
2026-07-14 14:49:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:49:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:57 INFO No previous orbit memory found.
2026-07-14 14:49:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020709.nc
2026-07-14 14:49:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:58 INFO No previous orbit memory found.
2026-07-14 14:49:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020709.nc
2026-07-14 14:49:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:49:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:49:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:49:59 INFO No previous orbit memory found.
2026-07-14 14:49:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:49:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020709.nc
2026-07-14 14:50:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:00 INFO No previous orbit memory found.
2026-07-14 14:50:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020709.nc
2026-07-14 14:50:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:01 INFO No previous orbit memory found.
2026-07-14 14:50:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020709.nc
2026-07-14 14:50:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:02 INFO No previous orbit memory found.
2026-07-14 14:50:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020709.nc
2026-07-14 14:50:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:03 INFO No previous orbit memory found.
2026-07-14 14:50:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020709.nc
2026-07-14 14:50:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:04 INFO No previous orbit memory found.
2026-07-14 14:50:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020709.nc
2026-07-14 14:50:05 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:05 INFO No previous orbit memory found.
2026-07-14 14:50:05 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020709.nc
2026-07-14 14:50:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:06 INFO No previous orbit memory found.
2026-07-14 14:50:06 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020709.nc
2026-07-14 14:50:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:07 INFO No previous orbit memory found.
2026-07-14 14:50:07 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020709.nc
2026-07-14 14:50:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:08 INFO No previous orbit memory found.
2026-07-14 14:50:08 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:08 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020709.nc
2026-07-14 14:50:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:09 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:09 INFO No previous orbit memory found.
2026-07-14 14:50:09 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020709.nc
2026-07-14 14:50:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:10 INFO No previous orbit memory found.
2026-07-14 14:50:10 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020709.nc
2026-07-14 14:50:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:11 INFO No previous orbit memory found.
2026-07-14 14:50:11 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020709.nc
2026-07-14 14:50:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:12 INFO No previous orbit memory found.
2026-07-14 14:50:12 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020709.nc
2026-07-14 14:50:13 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:13 INFO No previous orbit memory found.
2026-07-14 14:50:13 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020709.nc
2026-07-14 14:50:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:14 INFO No previous orbit memory found.
2026-07-14 14:50:14 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020709.nc
2026-07-14 14:50:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:15 INFO No previous orbit memory found.
2026-07-14 14:50:15 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020709.nc
2026-07-14 14:50:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:16 INFO No previous orbit memory found.
2026-07-14 14:50:16 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020709.nc
2026-07-14 14:50:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:17 INFO No previous orbit memory found.
2026-07-14 14:50:17 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020709.nc
2026-07-14 14:50:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 14:50:18 INFO No previous orbit memory found.
2026-07-14 14:50:18 INFO Emission correction applied with pixel-based damping.
2026-07-14 14:50:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020709.nc
2026-07-14 14:50:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 14:50:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 14:50:19 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 14:50:19 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:50:19 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 14:50:19 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-14 14:50:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:20 INFO Hourly dataset computed and listing created
2026-07-14 14:50:25 INFO Hourly dataset computed
2026-07-14 14:50:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:26 INFO Hourly dataset computed and listing created
2026-07-14 14:50:28 INFO Hourly dataset computed
2026-07-14 14:50:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:29 INFO Hourly dataset computed and listing created
2026-07-14 14:50:30 INFO Hourly dataset computed
2026-07-14 14:50:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:31 INFO Hourly dataset computed and listing created
2026-07-14 14:50:32 INFO Hourly dataset computed
2026-07-14 14:50:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:33 INFO Hourly dataset computed and listing created
2026-07-14 14:50:34 INFO Hourly dataset computed
2026-07-14 14:50:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:35 INFO Hourly dataset computed and listing created
2026-07-14 14:50:36 INFO Hourly dataset computed
2026-07-14 14:50:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:37 INFO Hourly dataset computed and listing created
2026-07-14 14:50:38 INFO Hourly dataset computed
2026-07-14 14:50:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:39 INFO Hourly dataset computed and listing created
2026-07-14 14:50:40 INFO Hourly dataset computed
2026-07-14 14:50:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:41 INFO Hourly dataset computed and listing created
2026-07-14 14:50:42 INFO Hourly dataset computed
2026-07-14 14:50:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:43 INFO Hourly dataset computed and listing created
2026-07-14 14:50:43 INFO Hourly dataset computed
2026-07-14 14:50:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:44 INFO Hourly dataset computed and listing created
2026-07-14 14:50:45 INFO Hourly dataset computed
2026-07-14 14:50:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:46 INFO Hourly dataset computed and listing created
2026-07-14 14:50:47 INFO Hourly dataset computed
2026-07-14 14:50:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:48 INFO Hourly dataset computed and listing created
2026-07-14 14:50:49 INFO Hourly dataset computed
2026-07-14 14:50:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:50 INFO Hourly dataset computed and listing created
2026-07-14 14:50:50 INFO Hourly dataset computed
2026-07-14 14:50:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 14:50:51 INFO Hourly dataset computed and listing created
2026-07-14 14:50:52 INFO Hourly dataset computed
2026-07-14 14:50:52 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-14 14:50:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:50:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 14:50:52 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-14 14:50:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 14:50:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:50:52 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 14:50:52 INFO Queuing job for member 1...
2026-07-14 14:50:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:50:52 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 14:50:53 INFO Found: ['5156601']
2026-07-14 14:50:58 INFO [TGCC-IRENE] Submitted job with ID:['5156601']
2026-07-14 14:50:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:50:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 14:50:58 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-14 14:50:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 14:50:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:50:58 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 14:50:58 INFO Queuing job for member 2...
2026-07-14 14:50:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:50:58 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 14:50:59 INFO Found: ['5156603']
2026-07-14 14:51:04 INFO [TGCC-IRENE] Submitted job with ID:['5156603']
2026-07-14 14:51:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:51:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 14:51:04 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-14 14:54:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 14:54:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:54:42 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 14:54:42 INFO Queuing job for member 3...
2026-07-14 14:54:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:54:42 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 14:54:42 INFO Found: ['5156617']
2026-07-14 14:54:47 INFO [TGCC-IRENE] Submitted job with ID:['5156617']
2026-07-14 14:54:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:54:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 14:54:47 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-14 14:54:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 14:54:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:54:47 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 14:54:47 INFO Queuing job for member 4...
2026-07-14 14:54:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:54:47 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 14:54:48 INFO Found: ['5156619']
2026-07-14 14:54:53 INFO [TGCC-IRENE] Submitted job with ID:['5156619']
2026-07-14 14:54:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:54:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 14:54:53 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-14 14:54:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 14:54:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:54:53 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 14:54:53 INFO Queuing job for member 5...
2026-07-14 14:54:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:54:53 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 14:54:55 INFO Found: ['5156621']
2026-07-14 14:55:00 INFO [TGCC-IRENE] Submitted job with ID:['5156621']
2026-07-14 14:55:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 14:55:00 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-14 14:55:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 14:55:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:00 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 14:55:00 INFO Queuing job for member 6...
2026-07-14 14:55:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:00 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 14:55:02 INFO Found: ['5156624']
2026-07-14 14:55:07 INFO [TGCC-IRENE] Submitted job with ID:['5156624']
2026-07-14 14:55:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 14:55:07 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-14 14:55:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 14:55:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:07 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 14:55:07 INFO Queuing job for member 7...
2026-07-14 14:55:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:07 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 14:55:08 INFO Found: ['5156627']
2026-07-14 14:55:13 INFO [TGCC-IRENE] Submitted job with ID:['5156627']
2026-07-14 14:55:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 14:55:13 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-14 14:55:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 14:55:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:13 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 14:55:13 INFO Queuing job for member 8...
2026-07-14 14:55:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:13 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 14:55:14 INFO Found: ['5156629']
2026-07-14 14:55:19 INFO [TGCC-IRENE] Submitted job with ID:['5156629']
2026-07-14 14:55:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 14:55:19 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-14 14:55:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 14:55:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:19 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 14:55:19 INFO Queuing job for member 9...
2026-07-14 14:55:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:19 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 14:55:20 INFO Found: ['5156630']
2026-07-14 14:55:25 INFO [TGCC-IRENE] Submitted job with ID:['5156630']
2026-07-14 14:55:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 14:55:25 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-14 14:55:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 14:55:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:25 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 14:55:25 INFO Queuing job for member 10...
2026-07-14 14:55:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:25 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 14:55:26 INFO Found: ['5156631']
2026-07-14 14:55:31 INFO [TGCC-IRENE] Submitted job with ID:['5156631']
2026-07-14 14:55:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 14:55:31 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-14 14:55:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 14:55:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:31 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 14:55:31 INFO Queuing job for member 11...
2026-07-14 14:55:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:31 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 14:55:31 INFO Found: ['5156632']
2026-07-14 14:55:36 INFO [TGCC-IRENE] Submitted job with ID:['5156632']
2026-07-14 14:55:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 14:55:36 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-14 14:55:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 14:55:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:36 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 14:55:36 INFO Queuing job for member 12...
2026-07-14 14:55:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:36 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 14:55:37 INFO Found: ['5156633']
2026-07-14 14:55:42 INFO [TGCC-IRENE] Submitted job with ID:['5156633']
2026-07-14 14:55:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 14:55:42 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-14 14:55:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 14:55:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:42 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 14:55:42 INFO Queuing job for member 13...
2026-07-14 14:55:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:42 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 14:55:43 INFO Found: ['5156634']
2026-07-14 14:55:48 INFO [TGCC-IRENE] Submitted job with ID:['5156634']
2026-07-14 14:55:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 14:55:48 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-14 14:55:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 14:55:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:48 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 14:55:48 INFO Queuing job for member 14...
2026-07-14 14:55:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:48 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 14:55:50 INFO Found: ['5156635']
2026-07-14 14:55:55 INFO [TGCC-IRENE] Submitted job with ID:['5156635']
2026-07-14 14:55:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 14:55:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 14:55:55 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-14 14:55:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 14:55:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 14:55:55 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 14:55:55 INFO Queuing job for member 15...
2026-07-14 14:55:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 14:55:55 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 14:55:58 INFO Found: ['5156636']
2026-07-14 14:56:03 INFO [TGCC-IRENE] Submitted job with ID:['5156636']
2026-07-14 14:56:03 INFO Checking job status ...
2026-07-14 14:56:03 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:56:03 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:56:03 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:56:18 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:56:18 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:56:18 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:56:33 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:56:33 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:56:34 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:56:34 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:56:51 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:56:51 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:56:51 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:57:06 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:57:06 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:57:06 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:57:21 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:57:21 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:57:21 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:57:36 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:57:36 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:57:36 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:57:37 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:57:37 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:57:54 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156603: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:57:54 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:57:54 INFO Jobs still running: ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:58:09 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156603: status FINISHED
2026-07-14 14:58:09 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:58:09 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:58:09 INFO Jobs still running: ['5156601', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:58:24 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156603: status FINISHED
2026-07-14 14:58:24 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:58:24 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:58:25 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:58:25 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:58:25 INFO Jobs still running: ['5156601', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:58:40 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156603: status FINISHED
2026-07-14 14:58:40 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:58:40 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:58:40 INFO Jobs still running: ['5156601', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:58:55 INFO None 5156601: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156603: status FINISHED
2026-07-14 14:58:55 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:58:55 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:58:56 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:58:56 INFO Jobs still running: ['5156601', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:59:11 INFO None 5156601: status FINISHED
2026-07-14 14:59:11 INFO None 5156603: status FINISHED
2026-07-14 14:59:11 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:59:11 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:59:11 INFO Jobs still running: ['5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:59:26 INFO None 5156601: status FINISHED
2026-07-14 14:59:26 INFO None 5156603: status FINISHED
2026-07-14 14:59:26 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:59:26 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:59:27 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:59:27 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:59:27 INFO Jobs still running: ['5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:59:43 INFO None 5156601: status FINISHED
2026-07-14 14:59:43 INFO None 5156603: status FINISHED
2026-07-14 14:59:43 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:59:43 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:59:43 INFO Jobs still running: ['5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 14:59:58 INFO None 5156601: status FINISHED
2026-07-14 14:59:58 INFO None 5156603: status FINISHED
2026-07-14 14:59:58 INFO None 5156617: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156619: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156621: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156624: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156627: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156629: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156630: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156631: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156632: status RUNNING/PENDING
2026-07-14 14:59:58 INFO None 5156633: status RUNNING/PENDING
2026-07-14 14:59:59 INFO None 5156634: status RUNNING/PENDING
2026-07-14 14:59:59 INFO None 5156635: status RUNNING/PENDING
2026-07-14 14:59:59 INFO None 5156636: status RUNNING/PENDING
2026-07-14 14:59:59 INFO Jobs still running: ['5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 15:00:14 INFO None 5156601: status FINISHED
2026-07-14 15:00:14 INFO None 5156603: status FINISHED
2026-07-14 15:00:14 INFO None 5156617: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156619: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156621: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156624: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156627: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156629: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156630: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156631: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156632: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156633: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156634: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:00:14 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:00:14 INFO Jobs still running: ['5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 15:00:29 INFO None 5156601: status FINISHED
2026-07-14 15:00:29 INFO None 5156603: status FINISHED
2026-07-14 15:00:29 INFO None 5156617: status FINISHED
2026-07-14 15:00:29 INFO None 5156619: status FINISHED
2026-07-14 15:00:29 INFO None 5156621: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156624: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156627: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156629: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156630: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156631: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156632: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156633: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156634: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:00:29 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:00:29 INFO Jobs still running: ['5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 15:00:46 INFO None 5156601: status FINISHED
2026-07-14 15:00:46 INFO None 5156603: status FINISHED
2026-07-14 15:00:46 INFO None 5156617: status FINISHED
2026-07-14 15:00:46 INFO None 5156619: status FINISHED
2026-07-14 15:00:46 INFO None 5156621: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156624: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156627: status FINISHED
2026-07-14 15:00:46 INFO None 5156629: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156630: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156631: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156632: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156633: status FINISHED
2026-07-14 15:00:46 INFO None 5156634: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:00:46 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:00:46 INFO Jobs still running: ['5156621', '5156624', '5156629', '5156630', '5156631', '5156632', '5156634', '5156635', '5156636']. Waiting...
2026-07-14 15:01:01 INFO None 5156601: status FINISHED
2026-07-14 15:03:10 INFO None 5156603: status FINISHED
2026-07-14 15:03:10 INFO None 5156617: status FINISHED
2026-07-14 15:03:10 INFO None 5156619: status FINISHED
2026-07-14 15:03:10 INFO None 5156621: status FINISHED
2026-07-14 15:03:10 INFO None 5156624: status FINISHED
2026-07-14 15:03:10 INFO None 5156627: status FINISHED
2026-07-14 15:03:10 INFO None 5156629: status FINISHED
2026-07-14 15:03:10 INFO None 5156630: status FINISHED
2026-07-14 15:03:10 INFO None 5156631: status FINISHED
2026-07-14 15:03:10 INFO None 5156632: status FINISHED
2026-07-14 15:03:10 INFO None 5156633: status FINISHED
2026-07-14 15:03:10 INFO None 5156634: status FINISHED
2026-07-14 15:03:10 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:03:10 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:03:10 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:03:25 INFO None 5156601: status FINISHED
2026-07-14 15:03:25 INFO None 5156603: status FINISHED
2026-07-14 15:03:25 INFO None 5156617: status FINISHED
2026-07-14 15:03:25 INFO None 5156619: status FINISHED
2026-07-14 15:03:25 INFO None 5156621: status FINISHED
2026-07-14 15:03:25 INFO None 5156624: status FINISHED
2026-07-14 15:03:25 INFO None 5156627: status FINISHED
2026-07-14 15:03:25 INFO None 5156629: status FINISHED
2026-07-14 15:03:25 INFO None 5156630: status FINISHED
2026-07-14 15:03:25 INFO None 5156631: status FINISHED
2026-07-14 15:03:25 INFO None 5156632: status FINISHED
2026-07-14 15:03:25 INFO None 5156633: status FINISHED
2026-07-14 15:03:25 INFO None 5156634: status FINISHED
2026-07-14 15:03:25 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:03:25 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:03:25 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:03:42 INFO None 5156601: status FINISHED
2026-07-14 15:03:42 INFO None 5156603: status FINISHED
2026-07-14 15:03:42 INFO None 5156617: status FINISHED
2026-07-14 15:03:42 INFO None 5156619: status FINISHED
2026-07-14 15:03:42 INFO None 5156621: status FINISHED
2026-07-14 15:03:42 INFO None 5156624: status FINISHED
2026-07-14 15:03:42 INFO None 5156627: status FINISHED
2026-07-14 15:03:42 INFO None 5156629: status FINISHED
2026-07-14 15:03:42 INFO None 5156630: status FINISHED
2026-07-14 15:03:42 INFO None 5156631: status FINISHED
2026-07-14 15:03:42 INFO None 5156632: status FINISHED
2026-07-14 15:03:42 INFO None 5156633: status FINISHED
2026-07-14 15:03:42 INFO None 5156634: status FINISHED
2026-07-14 15:03:42 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:03:44 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:03:44 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:03:59 INFO None 5156601: status FINISHED
2026-07-14 15:03:59 INFO None 5156603: status FINISHED
2026-07-14 15:03:59 INFO None 5156617: status FINISHED
2026-07-14 15:03:59 INFO None 5156619: status FINISHED
2026-07-14 15:03:59 INFO None 5156621: status FINISHED
2026-07-14 15:03:59 INFO None 5156624: status FINISHED
2026-07-14 15:03:59 INFO None 5156627: status FINISHED
2026-07-14 15:03:59 INFO None 5156629: status FINISHED
2026-07-14 15:03:59 INFO None 5156630: status FINISHED
2026-07-14 15:03:59 INFO None 5156631: status FINISHED
2026-07-14 15:03:59 INFO None 5156632: status FINISHED
2026-07-14 15:03:59 INFO None 5156633: status FINISHED
2026-07-14 15:03:59 INFO None 5156634: status FINISHED
2026-07-14 15:03:59 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:03:59 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:03:59 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:04:14 INFO None 5156601: status FINISHED
2026-07-14 15:04:14 INFO None 5156603: status FINISHED
2026-07-14 15:04:14 INFO None 5156617: status FINISHED
2026-07-14 15:04:14 INFO None 5156619: status FINISHED
2026-07-14 15:04:14 INFO None 5156621: status FINISHED
2026-07-14 15:04:14 INFO None 5156624: status FINISHED
2026-07-14 15:04:14 INFO None 5156627: status FINISHED
2026-07-14 15:04:14 INFO None 5156629: status FINISHED
2026-07-14 15:04:14 INFO None 5156630: status FINISHED
2026-07-14 15:04:14 INFO None 5156631: status FINISHED
2026-07-14 15:04:14 INFO None 5156632: status FINISHED
2026-07-14 15:04:14 INFO None 5156633: status FINISHED
2026-07-14 15:04:15 INFO None 5156634: status FINISHED
2026-07-14 15:04:15 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:04:15 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:04:15 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:04:31 INFO None 5156601: status FINISHED
2026-07-14 15:04:31 INFO None 5156603: status FINISHED
2026-07-14 15:04:31 INFO None 5156617: status FINISHED
2026-07-14 15:04:31 INFO None 5156619: status FINISHED
2026-07-14 15:04:31 INFO None 5156621: status FINISHED
2026-07-14 15:04:31 INFO None 5156624: status FINISHED
2026-07-14 15:04:31 INFO None 5156627: status FINISHED
2026-07-14 15:04:31 INFO None 5156629: status FINISHED
2026-07-14 15:04:31 INFO None 5156630: status FINISHED
2026-07-14 15:04:31 INFO None 5156631: status FINISHED
2026-07-14 15:04:31 INFO None 5156632: status FINISHED
2026-07-14 15:04:31 INFO None 5156633: status FINISHED
2026-07-14 15:04:31 INFO None 5156634: status FINISHED
2026-07-14 15:04:31 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:04:31 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:04:31 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:04:46 INFO None 5156601: status FINISHED
2026-07-14 15:04:46 INFO None 5156603: status FINISHED
2026-07-14 15:04:46 INFO None 5156617: status FINISHED
2026-07-14 15:04:46 INFO None 5156619: status FINISHED
2026-07-14 15:04:46 INFO None 5156621: status FINISHED
2026-07-14 15:04:46 INFO None 5156624: status FINISHED
2026-07-14 15:04:46 INFO None 5156627: status FINISHED
2026-07-14 15:04:46 INFO None 5156629: status FINISHED
2026-07-14 15:04:46 INFO None 5156630: status FINISHED
2026-07-14 15:04:46 INFO None 5156631: status FINISHED
2026-07-14 15:04:46 INFO None 5156632: status FINISHED
2026-07-14 15:04:46 INFO None 5156633: status FINISHED
2026-07-14 15:04:46 INFO None 5156634: status FINISHED
2026-07-14 15:04:46 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:04:46 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:04:46 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:05:01 INFO None 5156601: status FINISHED
2026-07-14 15:05:01 INFO None 5156603: status FINISHED
2026-07-14 15:05:01 INFO None 5156617: status FINISHED
2026-07-14 15:05:01 INFO None 5156619: status FINISHED
2026-07-14 15:05:01 INFO None 5156621: status FINISHED
2026-07-14 15:05:01 INFO None 5156624: status FINISHED
2026-07-14 15:05:02 INFO None 5156627: status FINISHED
2026-07-14 15:05:02 INFO None 5156629: status FINISHED
2026-07-14 15:05:02 INFO None 5156630: status FINISHED
2026-07-14 15:05:02 INFO None 5156631: status FINISHED
2026-07-14 15:05:02 INFO None 5156632: status FINISHED
2026-07-14 15:05:02 INFO None 5156633: status FINISHED
2026-07-14 15:05:02 INFO None 5156634: status FINISHED
2026-07-14 15:05:02 INFO None 5156635: status RUNNING/PENDING
2026-07-14 15:05:02 INFO None 5156636: status RUNNING/PENDING
2026-07-14 15:05:02 INFO Jobs still running: ['5156635', '5156636']. Waiting...
2026-07-14 15:05:17 INFO None 5156601: status FINISHED
2026-07-14 15:05:17 INFO None 5156603: status FINISHED
2026-07-14 15:05:17 INFO None 5156617: status FINISHED
2026-07-14 15:05:17 INFO None 5156619: status FINISHED
2026-07-14 15:05:17 INFO None 5156621: status FINISHED
2026-07-14 15:05:17 INFO None 5156624: status FINISHED
2026-07-14 15:05:17 INFO None 5156627: status FINISHED
2026-07-14 15:05:17 INFO None 5156629: status FINISHED
2026-07-14 15:05:17 INFO None 5156630: status FINISHED
2026-07-14 15:05:17 INFO None 5156631: status FINISHED
2026-07-14 15:05:17 INFO None 5156632: status FINISHED
2026-07-14 15:05:17 INFO None 5156633: status FINISHED
2026-07-14 15:05:17 INFO None 5156634: status FINISHED
2026-07-14 15:05:17 INFO None 5156635: status FINISHED
2026-07-14 15:05:17 INFO None 5156636: status FINISHED
2026-07-14 15:05:17 INFO Jobs ['5156601', '5156603', '5156617', '5156619', '5156621', '5156624', '5156627', '5156629', '5156630', '5156631', '5156632', '5156633', '5156634', '5156635', '5156636'] have finished
2026-07-14 15:05:17 INFO Checking restart files were created ...
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-14 15:05:17 INFO  Run_model() completed successfully.
2026-07-14 15:05:17 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:05:17 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-14 15:05:17 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 15:05:17 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-14 15:05:17 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:05:17 INFO ---------->>> Running process_satellite_data()
2026-07-14 15:05:17 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-14 15:05:17 INFO ---------->>> Running run_obs_converter()
2026-07-14 15:05:17 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-14 15:05:17 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-14 15:05:17 INFO ---------->>> Running DART
2026-07-14 15:05:17 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 15:05:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-14 15:05:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-14 15:05:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-14 15:05:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-14 15:05:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-14 15:05:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-14 15:05:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-14 15:05:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-14 15:05:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-14 15:05:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-14 15:05:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-14 15:05:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-14 15:05:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-14 15:05:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-14 15:05:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-14 15:05:23 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 15:05:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 15:05:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 15:05:23 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 15:05:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 15:05:23 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 15:05:35 INFO Found: []
2026-07-14 15:05:35 INFO No job id returned by command ./run_filter.bsh
2026-07-14 15:05:35 INFO No monitoring will be performed
2026-07-14 15:05:35 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-14 15:05:35 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:35 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:35 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:35 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020711'
2026-07-14 15:05:36 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 15:05:36 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 15:05:36 INFO run_dart() is DONE.
2026-07-14 15:05:36 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 15:05:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:36 INFO No previous orbit memory found.
2026-07-14 15:05:36 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020711.nc
2026-07-14 15:05:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:37 INFO No previous orbit memory found.
2026-07-14 15:05:37 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020711.nc
2026-07-14 15:05:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:38 INFO No previous orbit memory found.
2026-07-14 15:05:38 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020711.nc
2026-07-14 15:05:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:39 INFO No previous orbit memory found.
2026-07-14 15:05:39 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020711.nc
2026-07-14 15:05:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:40 INFO No previous orbit memory found.
2026-07-14 15:05:40 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020711.nc
2026-07-14 15:05:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:41 INFO No previous orbit memory found.
2026-07-14 15:05:41 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020711.nc
2026-07-14 15:05:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:42 INFO No previous orbit memory found.
2026-07-14 15:05:42 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020711.nc
2026-07-14 15:05:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:43 INFO No previous orbit memory found.
2026-07-14 15:05:43 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020711.nc
2026-07-14 15:05:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:44 INFO No previous orbit memory found.
2026-07-14 15:05:44 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020711.nc
2026-07-14 15:05:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:45 INFO No previous orbit memory found.
2026-07-14 15:05:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020711.nc
2026-07-14 15:05:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:46 INFO No previous orbit memory found.
2026-07-14 15:05:46 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020711.nc
2026-07-14 15:05:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:47 INFO No previous orbit memory found.
2026-07-14 15:05:47 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020711.nc
2026-07-14 15:05:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:48 INFO No previous orbit memory found.
2026-07-14 15:05:48 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020711.nc
2026-07-14 15:05:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:50 INFO No previous orbit memory found.
2026-07-14 15:05:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020711.nc
2026-07-14 15:05:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:51 INFO No previous orbit memory found.
2026-07-14 15:05:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020711.nc
2026-07-14 15:05:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:52 INFO No previous orbit memory found.
2026-07-14 15:05:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020711.nc
2026-07-14 15:05:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:53 INFO No previous orbit memory found.
2026-07-14 15:05:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020711.nc
2026-07-14 15:05:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:54 INFO No previous orbit memory found.
2026-07-14 15:05:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020711.nc
2026-07-14 15:05:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:55 INFO No previous orbit memory found.
2026-07-14 15:05:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020711.nc
2026-07-14 15:05:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:56 INFO No previous orbit memory found.
2026-07-14 15:05:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020711.nc
2026-07-14 15:05:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:57 INFO No previous orbit memory found.
2026-07-14 15:05:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020711.nc
2026-07-14 15:05:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:58 INFO No previous orbit memory found.
2026-07-14 15:05:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020711.nc
2026-07-14 15:05:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:05:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:05:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:05:59 INFO No previous orbit memory found.
2026-07-14 15:05:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:05:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020711.nc
2026-07-14 15:08:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:14 INFO No previous orbit memory found.
2026-07-14 15:08:14 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020711.nc
2026-07-14 15:08:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:08:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:15 INFO No previous orbit memory found.
2026-07-14 15:08:15 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020711.nc
2026-07-14 15:08:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:16 INFO No previous orbit memory found.
2026-07-14 15:08:16 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020711.nc
2026-07-14 15:08:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:08:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:17 INFO No previous orbit memory found.
2026-07-14 15:08:17 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020711.nc
2026-07-14 15:08:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:18 INFO No previous orbit memory found.
2026-07-14 15:08:18 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020711.nc
2026-07-14 15:08:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:08:19 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:19 INFO No previous orbit memory found.
2026-07-14 15:08:19 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:19 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020711.nc
2026-07-14 15:08:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:08:20 INFO No previous orbit memory found.
2026-07-14 15:08:20 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:08:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020711.nc
2026-07-14 15:08:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:08:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:08:21 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 15:08:21 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:08:21 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:08:21 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-14 15:08:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:22 INFO Hourly dataset computed and listing created
2026-07-14 15:08:26 INFO Hourly dataset computed
2026-07-14 15:08:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:27 INFO Hourly dataset computed and listing created
2026-07-14 15:08:29 INFO Hourly dataset computed
2026-07-14 15:08:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:30 INFO Hourly dataset computed and listing created
2026-07-14 15:08:31 INFO Hourly dataset computed
2026-07-14 15:08:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:32 INFO Hourly dataset computed and listing created
2026-07-14 15:08:33 INFO Hourly dataset computed
2026-07-14 15:08:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:34 INFO Hourly dataset computed and listing created
2026-07-14 15:08:35 INFO Hourly dataset computed
2026-07-14 15:08:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:36 INFO Hourly dataset computed and listing created
2026-07-14 15:08:37 INFO Hourly dataset computed
2026-07-14 15:08:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:38 INFO Hourly dataset computed and listing created
2026-07-14 15:08:38 INFO Hourly dataset computed
2026-07-14 15:08:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:39 INFO Hourly dataset computed and listing created
2026-07-14 15:08:40 INFO Hourly dataset computed
2026-07-14 15:08:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:41 INFO Hourly dataset computed and listing created
2026-07-14 15:08:41 INFO Hourly dataset computed
2026-07-14 15:08:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:42 INFO Hourly dataset computed and listing created
2026-07-14 15:08:43 INFO Hourly dataset computed
2026-07-14 15:08:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:44 INFO Hourly dataset computed and listing created
2026-07-14 15:08:44 INFO Hourly dataset computed
2026-07-14 15:08:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:45 INFO Hourly dataset computed and listing created
2026-07-14 15:08:46 INFO Hourly dataset computed
2026-07-14 15:08:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:47 INFO Hourly dataset computed and listing created
2026-07-14 15:08:47 INFO Hourly dataset computed
2026-07-14 15:08:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:48 INFO Hourly dataset computed and listing created
2026-07-14 15:08:48 INFO Hourly dataset computed
2026-07-14 15:08:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:08:49 INFO Hourly dataset computed and listing created
2026-07-14 15:08:50 INFO Hourly dataset computed
2026-07-14 15:08:50 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-14 15:08:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:08:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 15:08:50 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-14 15:08:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 15:08:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:08:50 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 15:08:50 INFO Queuing job for member 1...
2026-07-14 15:08:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:08:50 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 15:08:51 INFO Found: ['5156665']
2026-07-14 15:08:56 INFO [TGCC-IRENE] Submitted job with ID:['5156665']
2026-07-14 15:08:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:08:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 15:08:56 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-14 15:08:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 15:08:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:08:56 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 15:08:56 INFO Queuing job for member 2...
2026-07-14 15:08:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:08:56 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 15:08:56 INFO Found: ['5156666']
2026-07-14 15:09:01 INFO [TGCC-IRENE] Submitted job with ID:['5156666']
2026-07-14 15:09:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 15:09:01 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-14 15:09:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 15:09:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:01 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 15:09:01 INFO Queuing job for member 3...
2026-07-14 15:09:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:01 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 15:09:02 INFO Found: ['5156668']
2026-07-14 15:09:07 INFO [TGCC-IRENE] Submitted job with ID:['5156668']
2026-07-14 15:09:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 15:09:07 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-14 15:09:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 15:09:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:07 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 15:09:07 INFO Queuing job for member 4...
2026-07-14 15:09:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:07 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 15:09:08 INFO Found: ['5156669']
2026-07-14 15:09:13 INFO [TGCC-IRENE] Submitted job with ID:['5156669']
2026-07-14 15:09:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 15:09:13 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-14 15:09:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 15:09:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:13 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 15:09:13 INFO Queuing job for member 5...
2026-07-14 15:09:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:13 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 15:09:14 INFO Found: ['5156670']
2026-07-14 15:09:19 INFO [TGCC-IRENE] Submitted job with ID:['5156670']
2026-07-14 15:09:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 15:09:19 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-14 15:09:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 15:09:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:19 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 15:09:19 INFO Queuing job for member 6...
2026-07-14 15:09:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:19 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 15:09:20 INFO Found: ['5156671']
2026-07-14 15:09:25 INFO [TGCC-IRENE] Submitted job with ID:['5156671']
2026-07-14 15:09:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 15:09:25 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-14 15:09:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 15:09:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:25 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 15:09:25 INFO Queuing job for member 7...
2026-07-14 15:09:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:25 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 15:09:28 INFO Found: ['5156673']
2026-07-14 15:09:33 INFO [TGCC-IRENE] Submitted job with ID:['5156673']
2026-07-14 15:09:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 15:09:33 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-14 15:09:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 15:09:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:33 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 15:09:33 INFO Queuing job for member 8...
2026-07-14 15:09:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:33 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 15:09:33 INFO Found: ['5156674']
2026-07-14 15:09:38 INFO [TGCC-IRENE] Submitted job with ID:['5156674']
2026-07-14 15:09:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 15:09:38 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-14 15:09:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 15:09:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:38 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 15:09:38 INFO Queuing job for member 9...
2026-07-14 15:09:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:38 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 15:09:39 INFO Found: ['5156675']
2026-07-14 15:09:44 INFO [TGCC-IRENE] Submitted job with ID:['5156675']
2026-07-14 15:09:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 15:09:44 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-14 15:09:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 15:09:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:44 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 15:09:44 INFO Queuing job for member 10...
2026-07-14 15:09:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:44 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 15:09:45 INFO Found: ['5156676']
2026-07-14 15:09:50 INFO [TGCC-IRENE] Submitted job with ID:['5156676']
2026-07-14 15:09:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 15:09:50 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-14 15:09:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 15:09:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:50 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 15:09:50 INFO Queuing job for member 11...
2026-07-14 15:09:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:50 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 15:09:51 INFO Found: ['5156677']
2026-07-14 15:09:56 INFO [TGCC-IRENE] Submitted job with ID:['5156677']
2026-07-14 15:09:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:09:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 15:09:56 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-14 15:09:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 15:09:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:09:56 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 15:09:56 INFO Queuing job for member 12...
2026-07-14 15:09:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:09:56 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 15:09:56 INFO Found: ['5156678']
2026-07-14 15:10:01 INFO [TGCC-IRENE] Submitted job with ID:['5156678']
2026-07-14 15:10:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:10:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 15:10:01 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-14 15:10:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 15:10:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:10:01 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 15:10:01 INFO Queuing job for member 13...
2026-07-14 15:10:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:10:01 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 15:10:02 INFO Found: ['5156680']
2026-07-14 15:10:07 INFO [TGCC-IRENE] Submitted job with ID:['5156680']
2026-07-14 15:10:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:10:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 15:10:07 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-14 15:10:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 15:10:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:10:07 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 15:10:07 INFO Queuing job for member 14...
2026-07-14 15:10:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:10:07 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 15:10:08 INFO Found: ['5156681']
2026-07-14 15:10:13 INFO [TGCC-IRENE] Submitted job with ID:['5156681']
2026-07-14 15:10:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:10:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 15:10:13 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-14 15:10:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 15:10:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:10:13 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 15:10:13 INFO Queuing job for member 15...
2026-07-14 15:10:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:10:13 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 15:10:14 INFO Found: ['5156682']
2026-07-14 15:10:19 INFO [TGCC-IRENE] Submitted job with ID:['5156682']
2026-07-14 15:10:19 INFO Checking job status ...
2026-07-14 15:10:19 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:10:19 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:10:19 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:10:19 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:10:20 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:10:20 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:10:35 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:10:35 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:10:35 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:10:50 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:10:50 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:10:50 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:13:10 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:13:10 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:13:10 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:13:25 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:13:25 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:13:26 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:13:26 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:13:26 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:13:26 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:13:26 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:13:26 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:13:41 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156675: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:13:41 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:13:41 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:13:56 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156675: status FINISHED
2026-07-14 15:13:56 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:13:56 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:13:56 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:14:13 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156671: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156673: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156675: status FINISHED
2026-07-14 15:14:13 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:14:13 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:14:13 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:14:28 INFO None 5156665: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156666: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156668: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156669: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156670: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156671: status FINISHED
2026-07-14 15:14:28 INFO None 5156673: status FINISHED
2026-07-14 15:14:28 INFO None 5156674: status RUNNING/PENDING
2026-07-14 15:14:28 INFO None 5156675: status FINISHED
2026-07-14 15:14:29 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:14:29 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:14:29 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:14:29 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:14:29 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:14:29 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:14:29 INFO Jobs still running: ['5156665', '5156666', '5156668', '5156669', '5156670', '5156674', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:14:44 INFO None 5156665: status FINISHED
2026-07-14 15:14:44 INFO None 5156666: status FINISHED
2026-07-14 15:14:44 INFO None 5156668: status FINISHED
2026-07-14 15:14:44 INFO None 5156669: status FINISHED
2026-07-14 15:14:44 INFO None 5156670: status FINISHED
2026-07-14 15:14:44 INFO None 5156671: status FINISHED
2026-07-14 15:14:44 INFO None 5156673: status FINISHED
2026-07-14 15:14:44 INFO None 5156674: status FINISHED
2026-07-14 15:14:44 INFO None 5156675: status FINISHED
2026-07-14 15:14:44 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:14:44 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:14:44 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:14:44 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:14:44 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:14:44 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:14:44 INFO Jobs still running: ['5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:14:59 INFO None 5156665: status FINISHED
2026-07-14 15:14:59 INFO None 5156666: status FINISHED
2026-07-14 15:14:59 INFO None 5156668: status FINISHED
2026-07-14 15:14:59 INFO None 5156669: status FINISHED
2026-07-14 15:14:59 INFO None 5156670: status FINISHED
2026-07-14 15:14:59 INFO None 5156671: status FINISHED
2026-07-14 15:14:59 INFO None 5156673: status FINISHED
2026-07-14 15:14:59 INFO None 5156674: status FINISHED
2026-07-14 15:14:59 INFO None 5156675: status FINISHED
2026-07-14 15:14:59 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:14:59 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:14:59 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:14:59 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:14:59 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:14:59 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:14:59 INFO Jobs still running: ['5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:15:16 INFO None 5156665: status FINISHED
2026-07-14 15:15:16 INFO None 5156666: status FINISHED
2026-07-14 15:15:16 INFO None 5156668: status FINISHED
2026-07-14 15:15:16 INFO None 5156669: status FINISHED
2026-07-14 15:15:16 INFO None 5156670: status FINISHED
2026-07-14 15:15:16 INFO None 5156671: status FINISHED
2026-07-14 15:15:16 INFO None 5156673: status FINISHED
2026-07-14 15:15:16 INFO None 5156674: status FINISHED
2026-07-14 15:15:16 INFO None 5156675: status FINISHED
2026-07-14 15:15:16 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:15:16 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:15:16 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:15:16 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:15:16 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:15:16 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:15:16 INFO Jobs still running: ['5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:15:31 INFO None 5156665: status FINISHED
2026-07-14 15:15:31 INFO None 5156666: status FINISHED
2026-07-14 15:15:31 INFO None 5156668: status FINISHED
2026-07-14 15:15:31 INFO None 5156669: status FINISHED
2026-07-14 15:15:31 INFO None 5156670: status FINISHED
2026-07-14 15:15:31 INFO None 5156671: status FINISHED
2026-07-14 15:15:31 INFO None 5156673: status FINISHED
2026-07-14 15:15:31 INFO None 5156674: status FINISHED
2026-07-14 15:15:31 INFO None 5156675: status FINISHED
2026-07-14 15:15:31 INFO None 5156676: status RUNNING/PENDING
2026-07-14 15:15:31 INFO None 5156677: status RUNNING/PENDING
2026-07-14 15:15:31 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:15:31 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:15:31 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:15:31 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:15:31 INFO Jobs still running: ['5156676', '5156677', '5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:15:46 INFO None 5156665: status FINISHED
2026-07-14 15:15:46 INFO None 5156666: status FINISHED
2026-07-14 15:15:46 INFO None 5156668: status FINISHED
2026-07-14 15:15:46 INFO None 5156669: status FINISHED
2026-07-14 15:15:46 INFO None 5156670: status FINISHED
2026-07-14 15:15:46 INFO None 5156671: status FINISHED
2026-07-14 15:15:47 INFO None 5156673: status FINISHED
2026-07-14 15:15:47 INFO None 5156674: status FINISHED
2026-07-14 15:15:47 INFO None 5156675: status FINISHED
2026-07-14 15:15:47 INFO None 5156676: status FINISHED
2026-07-14 15:15:47 INFO None 5156677: status FINISHED
2026-07-14 15:15:47 INFO None 5156678: status RUNNING/PENDING
2026-07-14 15:15:47 INFO None 5156680: status RUNNING/PENDING
2026-07-14 15:15:47 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:15:47 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:15:47 INFO Jobs still running: ['5156678', '5156680', '5156681', '5156682']. Waiting...
2026-07-14 15:16:02 INFO None 5156665: status FINISHED
2026-07-14 15:16:02 INFO None 5156666: status FINISHED
2026-07-14 15:16:02 INFO None 5156668: status FINISHED
2026-07-14 15:16:02 INFO None 5156669: status FINISHED
2026-07-14 15:16:02 INFO None 5156670: status FINISHED
2026-07-14 15:16:02 INFO None 5156671: status FINISHED
2026-07-14 15:16:02 INFO None 5156673: status FINISHED
2026-07-14 15:16:02 INFO None 5156674: status FINISHED
2026-07-14 15:16:02 INFO None 5156675: status FINISHED
2026-07-14 15:16:02 INFO None 5156676: status FINISHED
2026-07-14 15:16:02 INFO None 5156677: status FINISHED
2026-07-14 15:16:02 INFO None 5156678: status FINISHED
2026-07-14 15:16:02 INFO None 5156680: status FINISHED
2026-07-14 15:16:02 INFO None 5156681: status RUNNING/PENDING
2026-07-14 15:16:02 INFO None 5156682: status RUNNING/PENDING
2026-07-14 15:16:02 INFO Jobs still running: ['5156681', '5156682']. Waiting...
2026-07-14 15:16:17 INFO None 5156665: status FINISHED
2026-07-14 15:18:16 INFO None 5156666: status FINISHED
2026-07-14 15:18:16 INFO None 5156668: status FINISHED
2026-07-14 15:18:16 INFO None 5156669: status FINISHED
2026-07-14 15:18:16 INFO None 5156670: status FINISHED
2026-07-14 15:18:16 INFO None 5156671: status FINISHED
2026-07-14 15:18:16 INFO None 5156673: status FINISHED
2026-07-14 15:18:16 INFO None 5156674: status FINISHED
2026-07-14 15:18:16 INFO None 5156675: status FINISHED
2026-07-14 15:18:16 INFO None 5156676: status FINISHED
2026-07-14 15:18:16 INFO None 5156677: status FINISHED
2026-07-14 15:18:16 INFO None 5156678: status FINISHED
2026-07-14 15:18:16 INFO None 5156680: status FINISHED
2026-07-14 15:18:16 INFO None 5156681: status FINISHED
2026-07-14 15:18:16 INFO None 5156682: status FINISHED
2026-07-14 15:18:16 INFO Jobs ['5156665', '5156666', '5156668', '5156669', '5156670', '5156671', '5156673', '5156674', '5156675', '5156676', '5156677', '5156678', '5156680', '5156681', '5156682'] have finished
2026-07-14 15:18:16 INFO Checking restart files were created ...
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-14 15:18:16 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-14 15:18:17 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-14 15:18:17 INFO  Run_model() completed successfully.
2026-07-14 15:18:17 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:18:17 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-14 15:18:17 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 15:18:17 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-14 15:18:17 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:18:17 INFO ---------->>> Running process_satellite_data()
2026-07-14 15:18:17 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-14 15:18:17 INFO ---------->>> Running run_obs_converter()
2026-07-14 15:18:17 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-14 15:18:17 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-14 15:18:17 INFO ---------->>> Running DART
2026-07-14 15:18:17 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 15:18:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-14 15:18:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-14 15:18:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-14 15:18:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-14 15:18:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-14 15:18:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-14 15:18:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-14 15:18:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-14 15:18:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-14 15:18:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-14 15:18:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-14 15:18:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-14 15:18:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-14 15:18:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-14 15:18:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-14 15:18:21 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 15:18:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 15:18:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 15:18:21 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 15:18:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 15:18:21 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 15:18:35 INFO Found: []
2026-07-14 15:18:35 INFO No job id returned by command ./run_filter.bsh
2026-07-14 15:18:35 INFO No monitoring will be performed
2026-07-14 15:18:35 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-14 15:18:35 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:35 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:35 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:35 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:35 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020712'
2026-07-14 15:18:36 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 15:18:36 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 15:18:36 INFO run_dart() is DONE.
2026-07-14 15:18:36 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 15:18:36 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:36 INFO No previous orbit memory found.
2026-07-14 15:18:36 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:36 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020712.nc
2026-07-14 15:18:37 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:37 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:37 INFO No previous orbit memory found.
2026-07-14 15:18:37 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:37 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020712.nc
2026-07-14 15:18:38 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:38 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:38 INFO No previous orbit memory found.
2026-07-14 15:18:38 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:38 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020712.nc
2026-07-14 15:18:39 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:39 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:39 INFO No previous orbit memory found.
2026-07-14 15:18:39 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:39 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020712.nc
2026-07-14 15:18:40 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:40 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:40 INFO No previous orbit memory found.
2026-07-14 15:18:40 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:40 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020712.nc
2026-07-14 15:18:41 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:41 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:41 INFO No previous orbit memory found.
2026-07-14 15:18:41 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:41 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020712.nc
2026-07-14 15:18:42 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:42 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:42 INFO No previous orbit memory found.
2026-07-14 15:18:42 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:42 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020712.nc
2026-07-14 15:18:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:43 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:43 INFO No previous orbit memory found.
2026-07-14 15:18:43 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:43 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020712.nc
2026-07-14 15:18:43 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:44 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:44 INFO No previous orbit memory found.
2026-07-14 15:18:44 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:44 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020712.nc
2026-07-14 15:18:44 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:45 INFO No previous orbit memory found.
2026-07-14 15:18:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020712.nc
2026-07-14 15:18:45 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:45 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:45 INFO No previous orbit memory found.
2026-07-14 15:18:45 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:45 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020712.nc
2026-07-14 15:18:46 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:46 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:46 INFO No previous orbit memory found.
2026-07-14 15:18:46 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:46 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020712.nc
2026-07-14 15:18:47 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:47 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:47 INFO No previous orbit memory found.
2026-07-14 15:18:47 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:47 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020712.nc
2026-07-14 15:18:48 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:48 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:48 INFO No previous orbit memory found.
2026-07-14 15:18:48 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:48 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020712.nc
2026-07-14 15:18:49 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:49 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:49 INFO No previous orbit memory found.
2026-07-14 15:18:49 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:49 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020712.nc
2026-07-14 15:18:50 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:50 INFO No previous orbit memory found.
2026-07-14 15:18:50 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020712.nc
2026-07-14 15:18:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:51 INFO No previous orbit memory found.
2026-07-14 15:18:51 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020712.nc
2026-07-14 15:18:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:52 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:52 INFO No previous orbit memory found.
2026-07-14 15:18:52 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:52 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020712.nc
2026-07-14 15:18:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:53 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:53 INFO No previous orbit memory found.
2026-07-14 15:18:53 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:53 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020712.nc
2026-07-14 15:18:53 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:54 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:54 INFO No previous orbit memory found.
2026-07-14 15:18:54 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:54 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020712.nc
2026-07-14 15:18:54 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:55 INFO No previous orbit memory found.
2026-07-14 15:18:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:55 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020712.nc
2026-07-14 15:18:55 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:55 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:55 INFO No previous orbit memory found.
2026-07-14 15:18:55 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020712.nc
2026-07-14 15:18:56 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:56 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:56 INFO No previous orbit memory found.
2026-07-14 15:18:56 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:56 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020712.nc
2026-07-14 15:18:57 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:57 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:57 INFO No previous orbit memory found.
2026-07-14 15:18:57 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:57 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020712.nc
2026-07-14 15:18:58 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:18:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:58 INFO No previous orbit memory found.
2026-07-14 15:18:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020712.nc
2026-07-14 15:18:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:18:59 INFO No previous orbit memory found.
2026-07-14 15:18:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:18:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020712.nc
2026-07-14 15:18:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:18:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:19:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:19:00 INFO No previous orbit memory found.
2026-07-14 15:19:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:19:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020712.nc
2026-07-14 15:19:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:19:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:19:01 INFO No previous orbit memory found.
2026-07-14 15:19:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:19:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020712.nc
2026-07-14 15:19:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:19:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:19:02 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:19:02 INFO No previous orbit memory found.
2026-07-14 15:19:02 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:19:02 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020712.nc
2026-07-14 15:19:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:19:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:19:03 INFO No previous orbit memory found.
2026-07-14 15:19:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:19:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020712.nc
2026-07-14 15:19:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:19:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:19:03 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 15:19:03 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:19:03 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:19:03 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-14 15:19:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:04 INFO Hourly dataset computed and listing created
2026-07-14 15:19:08 INFO Hourly dataset computed
2026-07-14 15:19:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:09 INFO Hourly dataset computed and listing created
2026-07-14 15:19:10 INFO Hourly dataset computed
2026-07-14 15:19:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:11 INFO Hourly dataset computed and listing created
2026-07-14 15:19:12 INFO Hourly dataset computed
2026-07-14 15:19:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:13 INFO Hourly dataset computed and listing created
2026-07-14 15:19:13 INFO Hourly dataset computed
2026-07-14 15:19:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:14 INFO Hourly dataset computed and listing created
2026-07-14 15:19:15 INFO Hourly dataset computed
2026-07-14 15:19:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:16 INFO Hourly dataset computed and listing created
2026-07-14 15:19:17 INFO Hourly dataset computed
2026-07-14 15:19:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:18 INFO Hourly dataset computed and listing created
2026-07-14 15:19:19 INFO Hourly dataset computed
2026-07-14 15:19:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:19 INFO Hourly dataset computed and listing created
2026-07-14 15:19:20 INFO Hourly dataset computed
2026-07-14 15:19:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:21 INFO Hourly dataset computed and listing created
2026-07-14 15:19:22 INFO Hourly dataset computed
2026-07-14 15:19:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:23 INFO Hourly dataset computed and listing created
2026-07-14 15:19:24 INFO Hourly dataset computed
2026-07-14 15:19:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:25 INFO Hourly dataset computed and listing created
2026-07-14 15:19:25 INFO Hourly dataset computed
2026-07-14 15:19:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:26 INFO Hourly dataset computed and listing created
2026-07-14 15:19:27 INFO Hourly dataset computed
2026-07-14 15:19:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:28 INFO Hourly dataset computed and listing created
2026-07-14 15:19:29 INFO Hourly dataset computed
2026-07-14 15:19:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:30 INFO Hourly dataset computed and listing created
2026-07-14 15:19:30 INFO Hourly dataset computed
2026-07-14 15:19:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:19:31 INFO Hourly dataset computed and listing created
2026-07-14 15:19:32 INFO Hourly dataset computed
2026-07-14 15:19:32 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-14 15:19:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:19:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 15:19:32 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-14 15:19:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 15:19:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:19:32 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 15:19:32 INFO Queuing job for member 1...
2026-07-14 15:19:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:19:32 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 15:19:33 INFO Found: ['5156702']
2026-07-14 15:19:38 INFO [TGCC-IRENE] Submitted job with ID:['5156702']
2026-07-14 15:19:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:19:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 15:19:38 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-14 15:19:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 15:19:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:19:38 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 15:19:38 INFO Queuing job for member 2...
2026-07-14 15:19:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:19:38 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 15:19:39 INFO Found: ['5156703']
2026-07-14 15:19:44 INFO [TGCC-IRENE] Submitted job with ID:['5156703']
2026-07-14 15:19:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:19:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 15:19:44 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-14 15:19:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 15:19:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:19:44 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 15:19:44 INFO Queuing job for member 3...
2026-07-14 15:19:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:19:44 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 15:19:44 INFO Found: ['5156704']
2026-07-14 15:19:49 INFO [TGCC-IRENE] Submitted job with ID:['5156704']
2026-07-14 15:19:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:19:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 15:19:49 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-14 15:19:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 15:19:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:19:49 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 15:19:49 INFO Queuing job for member 4...
2026-07-14 15:19:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:19:49 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 15:19:50 INFO Found: ['5156705']
2026-07-14 15:19:55 INFO [TGCC-IRENE] Submitted job with ID:['5156705']
2026-07-14 15:19:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:19:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 15:19:55 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-14 15:19:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 15:19:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:19:55 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 15:19:55 INFO Queuing job for member 5...
2026-07-14 15:19:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:19:55 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 15:19:56 INFO Found: ['5156706']
2026-07-14 15:20:01 INFO [TGCC-IRENE] Submitted job with ID:['5156706']
2026-07-14 15:20:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 15:20:01 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-14 15:20:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 15:20:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:01 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 15:20:01 INFO Queuing job for member 6...
2026-07-14 15:20:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:01 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 15:20:03 INFO Found: ['5156709']
2026-07-14 15:20:08 INFO [TGCC-IRENE] Submitted job with ID:['5156709']
2026-07-14 15:20:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 15:20:08 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-14 15:20:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 15:20:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:08 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 15:20:08 INFO Queuing job for member 7...
2026-07-14 15:20:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:08 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 15:20:10 INFO Found: ['5156710']
2026-07-14 15:20:15 INFO [TGCC-IRENE] Submitted job with ID:['5156710']
2026-07-14 15:20:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 15:20:15 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-14 15:20:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 15:20:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:15 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 15:20:15 INFO Queuing job for member 8...
2026-07-14 15:20:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:15 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 15:20:16 INFO Found: ['5156712']
2026-07-14 15:20:21 INFO [TGCC-IRENE] Submitted job with ID:['5156712']
2026-07-14 15:20:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 15:20:21 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-14 15:20:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 15:20:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:21 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 15:20:21 INFO Queuing job for member 9...
2026-07-14 15:20:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:21 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 15:20:21 INFO Found: ['5156714']
2026-07-14 15:20:26 INFO [TGCC-IRENE] Submitted job with ID:['5156714']
2026-07-14 15:20:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 15:20:26 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-14 15:20:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 15:20:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:26 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 15:20:26 INFO Queuing job for member 10...
2026-07-14 15:20:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:26 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 15:20:27 INFO Found: ['5156715']
2026-07-14 15:20:32 INFO [TGCC-IRENE] Submitted job with ID:['5156715']
2026-07-14 15:20:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 15:20:32 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-14 15:20:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 15:20:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:32 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 15:20:32 INFO Queuing job for member 11...
2026-07-14 15:20:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:32 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 15:20:33 INFO Found: ['5156716']
2026-07-14 15:20:38 INFO [TGCC-IRENE] Submitted job with ID:['5156716']
2026-07-14 15:20:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 15:20:38 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-14 15:20:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 15:20:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:38 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 15:20:38 INFO Queuing job for member 12...
2026-07-14 15:20:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:38 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 15:20:39 INFO Found: ['5156717']
2026-07-14 15:20:44 INFO [TGCC-IRENE] Submitted job with ID:['5156717']
2026-07-14 15:20:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 15:20:44 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-14 15:20:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 15:20:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:44 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 15:20:44 INFO Queuing job for member 13...
2026-07-14 15:20:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:44 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 15:20:44 INFO Found: ['5156718']
2026-07-14 15:20:49 INFO [TGCC-IRENE] Submitted job with ID:['5156718']
2026-07-14 15:20:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 15:20:49 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-14 15:20:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 15:20:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:49 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 15:20:49 INFO Queuing job for member 14...
2026-07-14 15:20:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:49 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 15:20:50 INFO Found: ['5156719']
2026-07-14 15:20:55 INFO [TGCC-IRENE] Submitted job with ID:['5156719']
2026-07-14 15:20:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:20:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 15:20:55 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-14 15:20:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 15:20:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:20:55 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 15:20:55 INFO Queuing job for member 15...
2026-07-14 15:20:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:20:55 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 15:20:57 INFO Found: ['5156720']
2026-07-14 15:21:02 INFO [TGCC-IRENE] Submitted job with ID:['5156720']
2026-07-14 15:21:02 INFO Checking job status ...
2026-07-14 15:21:02 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:21:02 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:21:02 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:21:17 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:23:14 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:23:14 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:23:29 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:23:29 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:23:29 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:23:44 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:23:44 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:23:45 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:23:45 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:24:00 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:24:00 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:24:00 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:24:15 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:24:15 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:24:15 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:24:30 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:24:30 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:24:31 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:24:31 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:24:46 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:24:46 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:24:46 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:25:01 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:25:01 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:25:01 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:25:16 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:25:16 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:25:17 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:25:17 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:25:32 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156706: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:25:32 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:25:32 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:25:47 INFO None 5156702: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156703: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156704: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156705: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156706: status FINISHED
2026-07-14 15:25:47 INFO None 5156709: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156714: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156715: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156716: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156717: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156718: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156719: status RUNNING/PENDING
2026-07-14 15:25:47 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:25:47 INFO Jobs still running: ['5156702', '5156703', '5156704', '5156705', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720']. Waiting...
2026-07-14 15:28:10 INFO None 5156702: status FINISHED
2026-07-14 15:28:10 INFO None 5156703: status FINISHED
2026-07-14 15:28:10 INFO None 5156704: status FINISHED
2026-07-14 15:28:10 INFO None 5156705: status FINISHED
2026-07-14 15:28:10 INFO None 5156706: status FINISHED
2026-07-14 15:28:10 INFO None 5156709: status FINISHED
2026-07-14 15:28:10 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:28:10 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:28:10 INFO None 5156714: status FINISHED
2026-07-14 15:28:10 INFO None 5156715: status FINISHED
2026-07-14 15:28:10 INFO None 5156716: status FINISHED
2026-07-14 15:28:10 INFO None 5156717: status FINISHED
2026-07-14 15:28:10 INFO None 5156718: status FINISHED
2026-07-14 15:28:10 INFO None 5156719: status FINISHED
2026-07-14 15:28:10 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:28:10 INFO Jobs still running: ['5156710', '5156712', '5156720']. Waiting...
2026-07-14 15:28:25 INFO None 5156702: status FINISHED
2026-07-14 15:28:25 INFO None 5156703: status FINISHED
2026-07-14 15:28:25 INFO None 5156704: status FINISHED
2026-07-14 15:28:25 INFO None 5156705: status FINISHED
2026-07-14 15:28:25 INFO None 5156706: status FINISHED
2026-07-14 15:28:25 INFO None 5156709: status FINISHED
2026-07-14 15:28:25 INFO None 5156710: status RUNNING/PENDING
2026-07-14 15:28:25 INFO None 5156712: status RUNNING/PENDING
2026-07-14 15:28:25 INFO None 5156714: status FINISHED
2026-07-14 15:28:25 INFO None 5156715: status FINISHED
2026-07-14 15:28:25 INFO None 5156716: status FINISHED
2026-07-14 15:28:25 INFO None 5156717: status FINISHED
2026-07-14 15:28:25 INFO None 5156718: status FINISHED
2026-07-14 15:28:25 INFO None 5156719: status FINISHED
2026-07-14 15:28:25 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:28:25 INFO Jobs still running: ['5156710', '5156712', '5156720']. Waiting...
2026-07-14 15:28:40 INFO None 5156702: status FINISHED
2026-07-14 15:28:40 INFO None 5156703: status FINISHED
2026-07-14 15:28:40 INFO None 5156704: status FINISHED
2026-07-14 15:28:40 INFO None 5156705: status FINISHED
2026-07-14 15:28:40 INFO None 5156706: status FINISHED
2026-07-14 15:28:40 INFO None 5156709: status FINISHED
2026-07-14 15:28:40 INFO None 5156710: status FINISHED
2026-07-14 15:28:40 INFO None 5156712: status FINISHED
2026-07-14 15:28:40 INFO None 5156714: status FINISHED
2026-07-14 15:28:40 INFO None 5156715: status FINISHED
2026-07-14 15:28:40 INFO None 5156716: status FINISHED
2026-07-14 15:28:40 INFO None 5156717: status FINISHED
2026-07-14 15:28:40 INFO None 5156718: status FINISHED
2026-07-14 15:28:40 INFO None 5156719: status FINISHED
2026-07-14 15:28:40 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:28:40 INFO Jobs still running: ['5156720']. Waiting...
2026-07-14 15:28:56 INFO None 5156702: status FINISHED
2026-07-14 15:28:56 INFO None 5156703: status FINISHED
2026-07-14 15:28:56 INFO None 5156704: status FINISHED
2026-07-14 15:28:56 INFO None 5156705: status FINISHED
2026-07-14 15:28:56 INFO None 5156706: status FINISHED
2026-07-14 15:28:56 INFO None 5156709: status FINISHED
2026-07-14 15:28:56 INFO None 5156710: status FINISHED
2026-07-14 15:28:56 INFO None 5156712: status FINISHED
2026-07-14 15:28:56 INFO None 5156714: status FINISHED
2026-07-14 15:28:56 INFO None 5156715: status FINISHED
2026-07-14 15:28:56 INFO None 5156716: status FINISHED
2026-07-14 15:28:56 INFO None 5156717: status FINISHED
2026-07-14 15:28:56 INFO None 5156718: status FINISHED
2026-07-14 15:28:56 INFO None 5156719: status FINISHED
2026-07-14 15:28:56 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:28:56 INFO Jobs still running: ['5156720']. Waiting...
2026-07-14 15:29:11 INFO None 5156702: status FINISHED
2026-07-14 15:29:11 INFO None 5156703: status FINISHED
2026-07-14 15:29:11 INFO None 5156704: status FINISHED
2026-07-14 15:29:11 INFO None 5156705: status FINISHED
2026-07-14 15:29:12 INFO None 5156706: status FINISHED
2026-07-14 15:29:12 INFO None 5156709: status FINISHED
2026-07-14 15:29:12 INFO None 5156710: status FINISHED
2026-07-14 15:29:12 INFO None 5156712: status FINISHED
2026-07-14 15:29:12 INFO None 5156714: status FINISHED
2026-07-14 15:29:12 INFO None 5156715: status FINISHED
2026-07-14 15:29:12 INFO None 5156716: status FINISHED
2026-07-14 15:29:12 INFO None 5156717: status FINISHED
2026-07-14 15:29:12 INFO None 5156718: status FINISHED
2026-07-14 15:29:12 INFO None 5156719: status FINISHED
2026-07-14 15:29:12 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:29:12 INFO Jobs still running: ['5156720']. Waiting...
2026-07-14 15:29:27 INFO None 5156702: status FINISHED
2026-07-14 15:29:27 INFO None 5156703: status FINISHED
2026-07-14 15:29:27 INFO None 5156704: status FINISHED
2026-07-14 15:29:27 INFO None 5156705: status FINISHED
2026-07-14 15:29:27 INFO None 5156706: status FINISHED
2026-07-14 15:29:27 INFO None 5156709: status FINISHED
2026-07-14 15:29:27 INFO None 5156710: status FINISHED
2026-07-14 15:29:27 INFO None 5156712: status FINISHED
2026-07-14 15:29:27 INFO None 5156714: status FINISHED
2026-07-14 15:29:27 INFO None 5156715: status FINISHED
2026-07-14 15:29:27 INFO None 5156716: status FINISHED
2026-07-14 15:29:27 INFO None 5156717: status FINISHED
2026-07-14 15:29:27 INFO None 5156718: status FINISHED
2026-07-14 15:29:27 INFO None 5156719: status FINISHED
2026-07-14 15:29:27 INFO None 5156720: status RUNNING/PENDING
2026-07-14 15:29:27 INFO Jobs still running: ['5156720']. Waiting...
2026-07-14 15:29:42 INFO None 5156702: status FINISHED
2026-07-14 15:29:42 INFO None 5156703: status FINISHED
2026-07-14 15:29:42 INFO None 5156704: status FINISHED
2026-07-14 15:29:42 INFO None 5156705: status FINISHED
2026-07-14 15:29:42 INFO None 5156706: status FINISHED
2026-07-14 15:29:42 INFO None 5156709: status FINISHED
2026-07-14 15:29:42 INFO None 5156710: status FINISHED
2026-07-14 15:29:42 INFO None 5156712: status FINISHED
2026-07-14 15:29:42 INFO None 5156714: status FINISHED
2026-07-14 15:29:42 INFO None 5156715: status FINISHED
2026-07-14 15:29:42 INFO None 5156716: status FINISHED
2026-07-14 15:29:42 INFO None 5156717: status FINISHED
2026-07-14 15:29:42 INFO None 5156718: status FINISHED
2026-07-14 15:29:42 INFO None 5156719: status FINISHED
2026-07-14 15:29:42 INFO None 5156720: status FINISHED
2026-07-14 15:29:42 INFO Jobs ['5156702', '5156703', '5156704', '5156705', '5156706', '5156709', '5156710', '5156712', '5156714', '5156715', '5156716', '5156717', '5156718', '5156719', '5156720'] have finished
2026-07-14 15:29:42 INFO Checking restart files were created ...
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc(1002685915 bytes)
2026-07-14 15:29:42 INFO  Run_model() completed successfully.
2026-07-14 15:29:42 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:29:42 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 14:00:00 days=153073 seconds=50400
2026-07-14 15:29:42 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 15:29:42 INFO [TIME] increment current_time 2020-02-07 12:00:00 -> 2020-02-07 14:00:00
2026-07-14 15:29:42 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:29:42 INFO ---------->>> Running process_satellite_data()
2026-07-14 15:29:42 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12019.nc
2026-07-14 15:29:42 INFO ---------->>> Running run_obs_converter()
2026-07-14 15:29:42 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-14 15:29:42 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-14 15:29:42 INFO ---------->>> Running DART
2026-07-14 15:29:42 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-14 15:29:42 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020714_1_out_toDART.nc
2026-07-14 15:29:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020714_1_out_toDART.nc
2026-07-14 15:29:43 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020714_1_out_toDART.nc
2026-07-14 15:29:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020714_1_out_toDART.nc
2026-07-14 15:29:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020714_1_out_toDART.nc
2026-07-14 15:29:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020714_1_out_toDART.nc
2026-07-14 15:29:44 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020714_1_out_toDART.nc
2026-07-14 15:29:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020714_1_out_toDART.nc
2026-07-14 15:29:45 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020714_1_out_toDART.nc
2026-07-14 15:29:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020714_1_out_toDART.nc
2026-07-14 15:29:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020714_1_out_toDART.nc
2026-07-14 15:29:46 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020714_1_out_toDART.nc
2026-07-14 15:29:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020714_1_out_toDART.nc
2026-07-14 15:29:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020714_1_out_toDART.nc
2026-07-14 15:29:47 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020714_1_out_toDART.nc
2026-07-14 15:29:48 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-14 15:29:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-14 15:29:48 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-14 15:29:48 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-14 15:29:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-14 15:29:48 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-14 15:29:57 INFO Found: []
2026-07-14 15:29:57 INFO No job id returned by command ./run_filter.bsh
2026-07-14 15:29:57 INFO No monitoring will be performed
2026-07-14 15:29:57 INFO Moving DART output files to analysis and preassim directories for date 2020020714 if present ...
2026-07-14 15:29:57 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020714'
2026-07-14 15:29:57 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-14 15:29:58 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-14 15:29:58 INFO run_dart() is DONE.
2026-07-14 15:29:58 INFO ---------->>> Running update_pollutant_in_end()
2026-07-14 15:29:58 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:29:58 INFO No previous orbit memory found.
2026-07-14 15:29:58 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:29:58 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020714.nc
2026-07-14 15:29:59 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:29:59 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:29:59 INFO No previous orbit memory found.
2026-07-14 15:29:59 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:29:59 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020714.nc
2026-07-14 15:30:00 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:00 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:00 INFO No previous orbit memory found.
2026-07-14 15:30:00 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:00 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS2_2020020714.nc
2026-07-14 15:30:01 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:01 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:01 INFO No previous orbit memory found.
2026-07-14 15:30:01 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:01 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS2_2020020714.nc
2026-07-14 15:30:02 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:03 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:03 INFO No previous orbit memory found.
2026-07-14 15:30:03 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:03 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS3_2020020714.nc
2026-07-14 15:30:03 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:04 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:04 INFO No previous orbit memory found.
2026-07-14 15:30:04 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:04 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS3_2020020714.nc
2026-07-14 15:30:04 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:05 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:05 INFO No previous orbit memory found.
2026-07-14 15:30:05 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:05 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS4_2020020714.nc
2026-07-14 15:30:06 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:06 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:06 INFO No previous orbit memory found.
2026-07-14 15:30:06 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:06 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS4_2020020714.nc
2026-07-14 15:30:07 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:07 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:07 INFO No previous orbit memory found.
2026-07-14 15:30:07 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:07 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS5_2020020714.nc
2026-07-14 15:30:08 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:08 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:08 INFO No previous orbit memory found.
2026-07-14 15:30:09 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:09 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS5_2020020714.nc
2026-07-14 15:30:09 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:10 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:10 INFO No previous orbit memory found.
2026-07-14 15:30:10 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:10 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS6_2020020714.nc
2026-07-14 15:30:10 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:11 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:11 INFO No previous orbit memory found.
2026-07-14 15:30:11 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:11 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS6_2020020714.nc
2026-07-14 15:30:11 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:12 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:12 INFO No previous orbit memory found.
2026-07-14 15:30:12 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:12 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS7_2020020714.nc
2026-07-14 15:30:12 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:13 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:13 INFO No previous orbit memory found.
2026-07-14 15:30:13 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:13 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS7_2020020714.nc
2026-07-14 15:30:14 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:14 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:14 INFO No previous orbit memory found.
2026-07-14 15:30:14 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:14 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS8_2020020714.nc
2026-07-14 15:30:15 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:15 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:15 INFO No previous orbit memory found.
2026-07-14 15:30:15 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:15 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS8_2020020714.nc
2026-07-14 15:30:16 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:16 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:16 INFO No previous orbit memory found.
2026-07-14 15:30:16 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:16 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS9_2020020714.nc
2026-07-14 15:30:17 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:17 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:17 INFO No previous orbit memory found.
2026-07-14 15:30:17 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:17 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS9_2020020714.nc
2026-07-14 15:30:18 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:18 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:18 INFO No previous orbit memory found.
2026-07-14 15:30:18 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:18 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS10_2020020714.nc
2026-07-14 15:30:19 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:20 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:20 INFO No previous orbit memory found.
2026-07-14 15:30:20 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:20 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS10_2020020714.nc
2026-07-14 15:30:20 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:21 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:21 INFO No previous orbit memory found.
2026-07-14 15:30:21 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:21 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS11_2020020714.nc
2026-07-14 15:30:21 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:22 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:22 INFO No previous orbit memory found.
2026-07-14 15:30:22 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:22 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS11_2020020714.nc
2026-07-14 15:30:22 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:23 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:23 INFO No previous orbit memory found.
2026-07-14 15:30:23 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:23 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS12_2020020714.nc
2026-07-14 15:30:24 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:24 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:24 INFO No previous orbit memory found.
2026-07-14 15:30:24 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:24 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS12_2020020714.nc
2026-07-14 15:30:25 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:26 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:26 INFO No previous orbit memory found.
2026-07-14 15:30:26 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:26 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS13_2020020714.nc
2026-07-14 15:30:26 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:27 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:27 INFO No previous orbit memory found.
2026-07-14 15:30:27 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:27 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS13_2020020714.nc
2026-07-14 15:30:27 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:28 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:28 INFO No previous orbit memory found.
2026-07-14 15:30:28 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:28 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS14_2020020714.nc
2026-07-14 15:30:29 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:29 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:29 INFO No previous orbit memory found.
2026-07-14 15:30:29 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:29 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS14_2020020714.nc
2026-07-14 15:30:30 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:30 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:30 INFO No previous orbit memory found.
2026-07-14 15:30:30 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:30 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS15_2020020714.nc
2026-07-14 15:30:31 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:31 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-14 15:30:31 INFO No previous orbit memory found.
2026-07-14 15:30:31 INFO Emission correction applied with pixel-based damping.
2026-07-14 15:30:31 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS15_2020020714.nc
2026-07-14 15:30:32 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-14 15:30:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-14 15:30:32 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 15:30:32 INFO [TIME] step_end current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:30:32 INFO [TIME] step_start current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 15:30:32 INFO [TIME] window start=2020-02-07 14:00:00 end=2020-02-08 00:00:00 run_hours=10 has_assimilation=False
2026-07-14 15:30:32 INFO Copying EMIS of next day ...
2026-07-14 15:30:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-14 15:30:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:30:34 INFO Hourly dataset computed and listing created
2026-07-14 15:30:52 INFO Hourly dataset computed
2026-07-14 15:30:52 INFO Copying EMIS of next day ...
2026-07-14 15:30:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-14 15:30:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:30:54 INFO Hourly dataset computed and listing created
2026-07-14 15:33:19 INFO Hourly dataset computed
2026-07-14 15:33:19 INFO Copying EMIS of next day ...
2026-07-14 15:33:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-14 15:33:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:33:22 INFO Hourly dataset computed and listing created
2026-07-14 15:33:36 INFO Hourly dataset computed
2026-07-14 15:33:36 INFO Copying EMIS of next day ...
2026-07-14 15:33:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-14 15:33:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:33:38 INFO Hourly dataset computed and listing created
2026-07-14 15:33:48 INFO Hourly dataset computed
2026-07-14 15:33:48 INFO Copying EMIS of next day ...
2026-07-14 15:33:48 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-14 15:33:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:33:50 INFO Hourly dataset computed and listing created
2026-07-14 15:34:03 INFO Hourly dataset computed
2026-07-14 15:34:03 INFO Copying EMIS of next day ...
2026-07-14 15:34:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-14 15:34:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:34:05 INFO Hourly dataset computed and listing created
2026-07-14 15:34:17 INFO Hourly dataset computed
2026-07-14 15:34:17 INFO Copying EMIS of next day ...
2026-07-14 15:34:18 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-14 15:34:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:34:20 INFO Hourly dataset computed and listing created
2026-07-14 15:34:31 INFO Hourly dataset computed
2026-07-14 15:34:31 INFO Copying EMIS of next day ...
2026-07-14 15:34:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-14 15:34:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:34:33 INFO Hourly dataset computed and listing created
2026-07-14 15:34:43 INFO Hourly dataset computed
2026-07-14 15:34:43 INFO Copying EMIS of next day ...
2026-07-14 15:34:43 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-14 15:34:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:34:45 INFO Hourly dataset computed and listing created
2026-07-14 15:34:55 INFO Hourly dataset computed
2026-07-14 15:34:55 INFO Copying EMIS of next day ...
2026-07-14 15:34:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-14 15:34:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:34:57 INFO Hourly dataset computed and listing created
2026-07-14 15:35:06 INFO Hourly dataset computed
2026-07-14 15:35:06 INFO Copying EMIS of next day ...
2026-07-14 15:35:07 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-14 15:35:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:35:09 INFO Hourly dataset computed and listing created
2026-07-14 15:35:18 INFO Hourly dataset computed
2026-07-14 15:35:18 INFO Copying EMIS of next day ...
2026-07-14 15:35:19 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-14 15:35:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:35:21 INFO Hourly dataset computed and listing created
2026-07-14 15:35:32 INFO Hourly dataset computed
2026-07-14 15:35:32 INFO Copying EMIS of next day ...
2026-07-14 15:35:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-14 15:35:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:35:34 INFO Hourly dataset computed and listing created
2026-07-14 15:35:52 INFO Hourly dataset computed
2026-07-14 15:35:52 INFO Copying EMIS of next day ...
2026-07-14 15:35:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-14 15:35:52 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:35:55 INFO Hourly dataset computed and listing created
2026-07-14 15:38:42 INFO Hourly dataset computed
2026-07-14 15:38:44 INFO Copying EMIS of next day ...
2026-07-14 15:38:44 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-14 15:38:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-14 15:38:46 INFO Hourly dataset computed and listing created
2026-07-14 15:39:06 INFO Hourly dataset computed
2026-07-14 15:39:06 INFO ---------->>> Running CHIMERE model from 2020-02-07 14:00:00 to 2020-02-08 00:00:00
2026-07-14 15:39:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-14 15:39:06 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-14 15:39:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-14 15:39:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:06 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-14 15:39:06 INFO Queuing job for member 1...
2026-07-14 15:39:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:06 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-14 15:39:07 INFO Found: ['5156842']
2026-07-14 15:39:12 INFO [TGCC-IRENE] Submitted job with ID:['5156842']
2026-07-14 15:39:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-14 15:39:12 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-14 15:39:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-14 15:39:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:12 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-14 15:39:12 INFO Queuing job for member 2...
2026-07-14 15:39:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:12 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-14 15:39:12 INFO Found: ['5156844']
2026-07-14 15:39:17 INFO [TGCC-IRENE] Submitted job with ID:['5156844']
2026-07-14 15:39:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-14 15:39:17 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-14 15:39:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-14 15:39:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:17 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-14 15:39:17 INFO Queuing job for member 3...
2026-07-14 15:39:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:17 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-14 15:39:18 INFO Found: ['5156846']
2026-07-14 15:39:23 INFO [TGCC-IRENE] Submitted job with ID:['5156846']
2026-07-14 15:39:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-14 15:39:23 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-14 15:39:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-14 15:39:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:23 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-14 15:39:23 INFO Queuing job for member 4...
2026-07-14 15:39:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:23 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-14 15:39:24 INFO Found: ['5156848']
2026-07-14 15:39:29 INFO [TGCC-IRENE] Submitted job with ID:['5156848']
2026-07-14 15:39:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-14 15:39:29 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-14 15:39:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-14 15:39:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:29 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-14 15:39:29 INFO Queuing job for member 5...
2026-07-14 15:39:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:29 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-14 15:39:30 INFO Found: ['5156850']
2026-07-14 15:39:35 INFO [TGCC-IRENE] Submitted job with ID:['5156850']
2026-07-14 15:39:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-14 15:39:35 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-14 15:39:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-14 15:39:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:35 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-14 15:39:35 INFO Queuing job for member 6...
2026-07-14 15:39:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:35 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-14 15:39:35 INFO Found: ['5156851']
2026-07-14 15:39:40 INFO [TGCC-IRENE] Submitted job with ID:['5156851']
2026-07-14 15:39:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-14 15:39:40 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-14 15:39:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-14 15:39:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:40 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-14 15:39:40 INFO Queuing job for member 7...
2026-07-14 15:39:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:40 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-14 15:39:42 INFO Found: ['5156852']
2026-07-14 15:39:47 INFO [TGCC-IRENE] Submitted job with ID:['5156852']
2026-07-14 15:39:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-14 15:39:47 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-14 15:39:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-14 15:39:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:47 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-14 15:39:47 INFO Queuing job for member 8...
2026-07-14 15:39:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:47 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-14 15:39:48 INFO Found: ['5156853']
2026-07-14 15:39:53 INFO [TGCC-IRENE] Submitted job with ID:['5156853']
2026-07-14 15:39:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-14 15:39:53 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-14 15:39:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-14 15:39:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:53 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-14 15:39:53 INFO Queuing job for member 9...
2026-07-14 15:39:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:53 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-14 15:39:54 INFO Found: ['5156854']
2026-07-14 15:39:59 INFO [TGCC-IRENE] Submitted job with ID:['5156854']
2026-07-14 15:39:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:39:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-14 15:39:59 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-14 15:39:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-14 15:39:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:39:59 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-14 15:39:59 INFO Queuing job for member 10...
2026-07-14 15:39:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:39:59 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-14 15:40:00 INFO Found: ['5156855']
2026-07-14 15:40:05 INFO [TGCC-IRENE] Submitted job with ID:['5156855']
2026-07-14 15:40:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:40:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-14 15:40:05 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-14 15:40:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-14 15:40:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:40:05 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-14 15:40:05 INFO Queuing job for member 11...
2026-07-14 15:40:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:40:05 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-14 15:40:06 INFO Found: ['5156857']
2026-07-14 15:40:11 INFO [TGCC-IRENE] Submitted job with ID:['5156857']
2026-07-14 15:40:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:40:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-14 15:40:11 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-14 15:40:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-14 15:40:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:40:11 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-14 15:40:11 INFO Queuing job for member 12...
2026-07-14 15:40:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:40:11 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-14 15:40:11 INFO Found: ['5156858']
2026-07-14 15:40:16 INFO [TGCC-IRENE] Submitted job with ID:['5156858']
2026-07-14 15:40:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:40:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-14 15:40:16 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-14 15:40:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-14 15:40:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:40:16 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-14 15:40:16 INFO Queuing job for member 13...
2026-07-14 15:40:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:40:16 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-14 15:40:17 INFO Found: ['5156859']
2026-07-14 15:40:22 INFO [TGCC-IRENE] Submitted job with ID:['5156859']
2026-07-14 15:40:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:40:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-14 15:40:22 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-14 15:40:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-14 15:40:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:40:22 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-14 15:40:22 INFO Queuing job for member 14...
2026-07-14 15:40:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:40:22 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-14 15:40:23 INFO Found: ['5156860']
2026-07-14 15:40:28 INFO [TGCC-IRENE] Submitted job with ID:['5156860']
2026-07-14 15:40:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-14 15:40:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-14 15:40:28 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-14 15:40:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-14 15:40:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-14 15:40:28 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-14 15:40:28 INFO Queuing job for member 15...
2026-07-14 15:40:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-14 15:40:28 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-14 15:40:30 INFO Found: ['5156861']
2026-07-14 15:40:35 INFO [TGCC-IRENE] Submitted job with ID:['5156861']
2026-07-14 15:40:35 INFO Checking job status ...
2026-07-14 15:40:36 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:40:36 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:40:36 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:40:51 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:40:51 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:40:51 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:43:12 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:43:12 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:43:12 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:43:27 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:43:27 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:43:27 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:43:42 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:43:42 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:43:43 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:43:43 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:43:58 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:43:58 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:43:58 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:44:13 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:44:13 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:44:13 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:44:28 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:44:29 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:44:29 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:44:44 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:44:44 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:44:44 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:44:59 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:44:59 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:44:59 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:45:14 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:45:15 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:45:15 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:45:30 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:45:30 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:45:31 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:45:31 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:45:46 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:45:46 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:45:46 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:48:09 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:48:09 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:48:09 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:48:24 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156850: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:48:25 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:48:25 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:48:40 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156850: status FINISHED
2026-07-14 15:48:40 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:48:40 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:48:40 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:48:55 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156850: status FINISHED
2026-07-14 15:48:55 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:48:55 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:48:55 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:49:10 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156848: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156850: status FINISHED
2026-07-14 15:49:11 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:49:11 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:49:11 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156848', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:49:27 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:49:27 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:49:27 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:49:27 INFO None 5156848: status FINISHED
2026-07-14 15:49:27 INFO None 5156850: status FINISHED
2026-07-14 15:49:27 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:49:28 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:49:28 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:49:43 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156848: status FINISHED
2026-07-14 15:49:43 INFO None 5156850: status FINISHED
2026-07-14 15:49:43 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:49:43 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:49:43 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:49:58 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156848: status FINISHED
2026-07-14 15:49:58 INFO None 5156850: status FINISHED
2026-07-14 15:49:58 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:49:58 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:49:58 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:50:13 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:50:13 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:50:13 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:50:13 INFO None 5156848: status FINISHED
2026-07-14 15:50:13 INFO None 5156850: status FINISHED
2026-07-14 15:50:13 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:50:13 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:50:14 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:50:14 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:50:30 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156848: status FINISHED
2026-07-14 15:50:30 INFO None 5156850: status FINISHED
2026-07-14 15:50:30 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:50:30 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:50:30 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:50:45 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156848: status FINISHED
2026-07-14 15:50:46 INFO None 5156850: status FINISHED
2026-07-14 15:50:46 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:50:46 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:50:46 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:51:01 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156848: status FINISHED
2026-07-14 15:51:01 INFO None 5156850: status FINISHED
2026-07-14 15:51:01 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:51:01 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:51:01 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:51:16 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156848: status FINISHED
2026-07-14 15:51:16 INFO None 5156850: status FINISHED
2026-07-14 15:51:16 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:51:16 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:51:16 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:51:33 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156848: status FINISHED
2026-07-14 15:51:33 INFO None 5156850: status FINISHED
2026-07-14 15:51:33 INFO None 5156851: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156852: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156853: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156854: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156855: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:51:33 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:51:33 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:51:48 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156848: status FINISHED
2026-07-14 15:53:37 INFO None 5156850: status FINISHED
2026-07-14 15:53:37 INFO None 5156851: status FINISHED
2026-07-14 15:53:37 INFO None 5156852: status FINISHED
2026-07-14 15:53:37 INFO None 5156853: status FINISHED
2026-07-14 15:53:37 INFO None 5156854: status FINISHED
2026-07-14 15:53:37 INFO None 5156855: status FINISHED
2026-07-14 15:53:37 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:53:37 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:53:37 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:53:52 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156848: status FINISHED
2026-07-14 15:53:52 INFO None 5156850: status FINISHED
2026-07-14 15:53:52 INFO None 5156851: status FINISHED
2026-07-14 15:53:52 INFO None 5156852: status FINISHED
2026-07-14 15:53:52 INFO None 5156853: status FINISHED
2026-07-14 15:53:52 INFO None 5156854: status FINISHED
2026-07-14 15:53:52 INFO None 5156855: status FINISHED
2026-07-14 15:53:52 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:53:52 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:53:52 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:54:07 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156848: status FINISHED
2026-07-14 15:54:07 INFO None 5156850: status FINISHED
2026-07-14 15:54:07 INFO None 5156851: status FINISHED
2026-07-14 15:54:07 INFO None 5156852: status FINISHED
2026-07-14 15:54:07 INFO None 5156853: status FINISHED
2026-07-14 15:54:07 INFO None 5156854: status FINISHED
2026-07-14 15:54:07 INFO None 5156855: status FINISHED
2026-07-14 15:54:07 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:54:07 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:54:07 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:54:22 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:54:22 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:54:22 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:54:22 INFO None 5156848: status FINISHED
2026-07-14 15:54:23 INFO None 5156850: status FINISHED
2026-07-14 15:54:23 INFO None 5156851: status FINISHED
2026-07-14 15:54:23 INFO None 5156852: status FINISHED
2026-07-14 15:54:23 INFO None 5156853: status FINISHED
2026-07-14 15:54:23 INFO None 5156854: status FINISHED
2026-07-14 15:54:23 INFO None 5156855: status FINISHED
2026-07-14 15:54:23 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:54:23 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:54:23 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:54:23 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:54:23 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:54:23 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:54:38 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156848: status FINISHED
2026-07-14 15:54:38 INFO None 5156850: status FINISHED
2026-07-14 15:54:38 INFO None 5156851: status FINISHED
2026-07-14 15:54:38 INFO None 5156852: status FINISHED
2026-07-14 15:54:38 INFO None 5156853: status FINISHED
2026-07-14 15:54:38 INFO None 5156854: status FINISHED
2026-07-14 15:54:38 INFO None 5156855: status FINISHED
2026-07-14 15:54:38 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:54:38 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:54:38 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:54:53 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:54:53 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:54:53 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:54:53 INFO None 5156848: status FINISHED
2026-07-14 15:54:53 INFO None 5156850: status FINISHED
2026-07-14 15:54:53 INFO None 5156851: status FINISHED
2026-07-14 15:54:53 INFO None 5156852: status FINISHED
2026-07-14 15:54:53 INFO None 5156853: status FINISHED
2026-07-14 15:54:53 INFO None 5156854: status FINISHED
2026-07-14 15:54:53 INFO None 5156855: status FINISHED
2026-07-14 15:54:53 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:54:53 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:54:54 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:54:54 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:54:54 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:54:54 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:55:09 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156848: status FINISHED
2026-07-14 15:55:09 INFO None 5156850: status FINISHED
2026-07-14 15:55:09 INFO None 5156851: status FINISHED
2026-07-14 15:55:09 INFO None 5156852: status FINISHED
2026-07-14 15:55:09 INFO None 5156853: status FINISHED
2026-07-14 15:55:09 INFO None 5156854: status FINISHED
2026-07-14 15:55:09 INFO None 5156855: status FINISHED
2026-07-14 15:55:09 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:55:09 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:55:09 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:55:24 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156846: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156848: status FINISHED
2026-07-14 15:55:24 INFO None 5156850: status FINISHED
2026-07-14 15:55:24 INFO None 5156851: status FINISHED
2026-07-14 15:55:24 INFO None 5156852: status FINISHED
2026-07-14 15:55:24 INFO None 5156853: status FINISHED
2026-07-14 15:55:24 INFO None 5156854: status FINISHED
2026-07-14 15:55:24 INFO None 5156855: status FINISHED
2026-07-14 15:55:24 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:55:24 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:55:24 INFO Jobs still running: ['5156842', '5156844', '5156846', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:55:41 INFO None 5156842: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156846: status FINISHED
2026-07-14 15:55:41 INFO None 5156848: status FINISHED
2026-07-14 15:55:41 INFO None 5156850: status FINISHED
2026-07-14 15:55:41 INFO None 5156851: status FINISHED
2026-07-14 15:55:41 INFO None 5156852: status FINISHED
2026-07-14 15:55:41 INFO None 5156853: status FINISHED
2026-07-14 15:55:41 INFO None 5156854: status FINISHED
2026-07-14 15:55:41 INFO None 5156855: status FINISHED
2026-07-14 15:55:41 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:55:41 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:55:41 INFO Jobs still running: ['5156842', '5156844', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:55:56 INFO None 5156842: status FINISHED
2026-07-14 15:55:56 INFO None 5156844: status RUNNING/PENDING
2026-07-14 15:55:56 INFO None 5156846: status FINISHED
2026-07-14 15:55:56 INFO None 5156848: status FINISHED
2026-07-14 15:55:56 INFO None 5156850: status FINISHED
2026-07-14 15:55:56 INFO None 5156851: status FINISHED
2026-07-14 15:55:56 INFO None 5156852: status FINISHED
2026-07-14 15:55:56 INFO None 5156853: status FINISHED
2026-07-14 15:55:57 INFO None 5156854: status FINISHED
2026-07-14 15:55:57 INFO None 5156855: status FINISHED
2026-07-14 15:55:57 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:55:57 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:55:57 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:55:57 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:55:57 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:55:57 INFO Jobs still running: ['5156844', '5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:56:12 INFO None 5156842: status FINISHED
2026-07-14 15:56:12 INFO None 5156844: status FINISHED
2026-07-14 15:56:12 INFO None 5156846: status FINISHED
2026-07-14 15:56:12 INFO None 5156848: status FINISHED
2026-07-14 15:56:12 INFO None 5156850: status FINISHED
2026-07-14 15:56:12 INFO None 5156851: status FINISHED
2026-07-14 15:56:12 INFO None 5156852: status FINISHED
2026-07-14 15:56:12 INFO None 5156853: status FINISHED
2026-07-14 15:56:12 INFO None 5156854: status FINISHED
2026-07-14 15:56:12 INFO None 5156855: status FINISHED
2026-07-14 15:56:12 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:56:12 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:56:12 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:56:12 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:56:12 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:56:12 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:56:27 INFO None 5156842: status FINISHED
2026-07-14 15:56:27 INFO None 5156844: status FINISHED
2026-07-14 15:56:27 INFO None 5156846: status FINISHED
2026-07-14 15:56:27 INFO None 5156848: status FINISHED
2026-07-14 15:56:27 INFO None 5156850: status FINISHED
2026-07-14 15:56:27 INFO None 5156851: status FINISHED
2026-07-14 15:56:27 INFO None 5156852: status FINISHED
2026-07-14 15:56:27 INFO None 5156853: status FINISHED
2026-07-14 15:56:27 INFO None 5156854: status FINISHED
2026-07-14 15:56:27 INFO None 5156855: status FINISHED
2026-07-14 15:56:27 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:56:27 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:56:27 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:56:27 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:56:27 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:56:27 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:56:44 INFO None 5156842: status FINISHED
2026-07-14 15:56:44 INFO None 5156844: status FINISHED
2026-07-14 15:56:44 INFO None 5156846: status FINISHED
2026-07-14 15:56:44 INFO None 5156848: status FINISHED
2026-07-14 15:56:44 INFO None 5156850: status FINISHED
2026-07-14 15:56:44 INFO None 5156851: status FINISHED
2026-07-14 15:56:44 INFO None 5156852: status FINISHED
2026-07-14 15:56:44 INFO None 5156853: status FINISHED
2026-07-14 15:56:44 INFO None 5156854: status FINISHED
2026-07-14 15:56:44 INFO None 5156855: status FINISHED
2026-07-14 15:56:44 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:56:44 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:56:44 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:56:44 INFO None 5156860: status RUNNING/PENDING
2026-07-14 15:56:44 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:56:44 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156860', '5156861']. Waiting...
2026-07-14 15:56:59 INFO None 5156842: status FINISHED
2026-07-14 15:58:57 INFO None 5156844: status FINISHED
2026-07-14 15:58:57 INFO None 5156846: status FINISHED
2026-07-14 15:58:57 INFO None 5156848: status FINISHED
2026-07-14 15:58:57 INFO None 5156850: status FINISHED
2026-07-14 15:58:57 INFO None 5156851: status FINISHED
2026-07-14 15:58:57 INFO None 5156852: status FINISHED
2026-07-14 15:58:57 INFO None 5156853: status FINISHED
2026-07-14 15:58:57 INFO None 5156854: status FINISHED
2026-07-14 15:58:58 INFO None 5156855: status FINISHED
2026-07-14 15:58:58 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:58:58 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:58:58 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:58:58 INFO None 5156860: status FINISHED
2026-07-14 15:58:58 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:58:58 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156861']. Waiting...
2026-07-14 15:59:13 INFO None 5156842: status FINISHED
2026-07-14 15:59:13 INFO None 5156844: status FINISHED
2026-07-14 15:59:13 INFO None 5156846: status FINISHED
2026-07-14 15:59:13 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 15:59:13 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 15:59:13 INFO None 5156851: status FINISHED
2026-07-14 15:59:13 INFO None 5156852: status FINISHED
2026-07-14 15:59:13 INFO None 5156853: status FINISHED
2026-07-14 15:59:13 INFO None 5156854: status FINISHED
2026-07-14 15:59:13 INFO None 5156855: status FINISHED
2026-07-14 15:59:13 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:59:13 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:59:13 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:59:13 INFO None 5156860: status FINISHED
2026-07-14 15:59:13 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:59:13 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156861']. Waiting...
2026-07-14 15:59:28 INFO None 5156842: status FINISHED
2026-07-14 15:59:28 INFO None 5156844: status FINISHED
2026-07-14 15:59:28 INFO None 5156846: status FINISHED
2026-07-14 15:59:28 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 15:59:28 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 15:59:28 INFO None 5156851: status FINISHED
2026-07-14 15:59:28 INFO None 5156852: status FINISHED
2026-07-14 15:59:28 INFO None 5156853: status FINISHED
2026-07-14 15:59:28 INFO None 5156854: status FINISHED
2026-07-14 15:59:28 INFO None 5156855: status FINISHED
2026-07-14 15:59:28 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:59:28 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:59:28 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:59:28 INFO None 5156860: status FINISHED
2026-07-14 15:59:28 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:59:28 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156861']. Waiting...
2026-07-14 15:59:43 INFO None 5156842: status FINISHED
2026-07-14 15:59:43 INFO None 5156844: status FINISHED
2026-07-14 15:59:43 INFO None 5156846: status FINISHED
2026-07-14 15:59:43 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 15:59:43 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 15:59:43 INFO None 5156851: status FINISHED
2026-07-14 15:59:43 INFO None 5156852: status FINISHED
2026-07-14 15:59:43 INFO None 5156853: status FINISHED
2026-07-14 15:59:43 INFO None 5156854: status FINISHED
2026-07-14 15:59:43 INFO None 5156855: status FINISHED
2026-07-14 15:59:43 INFO None 5156857: status RUNNING/PENDING
2026-07-14 15:59:43 INFO None 5156858: status RUNNING/PENDING
2026-07-14 15:59:43 INFO None 5156859: status RUNNING/PENDING
2026-07-14 15:59:43 INFO None 5156860: status FINISHED
2026-07-14 15:59:44 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:59:44 INFO Jobs still running: ['5156857', '5156858', '5156859', '5156861']. Waiting...
2026-07-14 15:59:59 INFO None 5156842: status FINISHED
2026-07-14 15:59:59 INFO None 5156844: status FINISHED
2026-07-14 15:59:59 INFO None 5156846: status FINISHED
2026-07-14 15:59:59 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 15:59:59 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 15:59:59 INFO None 5156851: status FINISHED
2026-07-14 15:59:59 INFO None 5156852: status FINISHED
2026-07-14 15:59:59 INFO None 5156853: status FINISHED
2026-07-14 15:59:59 INFO None 5156854: status FINISHED
2026-07-14 15:59:59 INFO None 5156855: status FINISHED
2026-07-14 15:59:59 INFO None 5156857: status FINISHED
2026-07-14 15:59:59 INFO None 5156858: status FINISHED
2026-07-14 15:59:59 INFO None 5156859: status FINISHED
2026-07-14 15:59:59 INFO None 5156860: status FINISHED
2026-07-14 15:59:59 INFO None 5156861: status RUNNING/PENDING
2026-07-14 15:59:59 INFO Jobs still running: ['5156861']. Waiting...
2026-07-14 16:00:14 INFO None 5156842: status FINISHED
2026-07-14 16:00:14 INFO None 5156844: status FINISHED
2026-07-14 16:00:14 INFO None 5156846: status FINISHED
2026-07-14 16:00:14 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 16:00:14 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 16:00:14 INFO None 5156851: status FINISHED
2026-07-14 16:00:14 INFO None 5156852: status FINISHED
2026-07-14 16:00:14 INFO None 5156853: status FINISHED
2026-07-14 16:00:14 INFO None 5156854: status FINISHED
2026-07-14 16:00:14 INFO None 5156855: status FINISHED
2026-07-14 16:00:14 INFO None 5156857: status FINISHED
2026-07-14 16:00:14 INFO None 5156858: status FINISHED
2026-07-14 16:00:14 INFO None 5156859: status FINISHED
2026-07-14 16:00:14 INFO None 5156860: status FINISHED
2026-07-14 16:00:14 INFO None 5156861: status RUNNING/PENDING
2026-07-14 16:00:14 INFO Jobs still running: ['5156861']. Waiting...
2026-07-14 16:00:29 INFO None 5156842: status FINISHED
2026-07-14 16:00:29 INFO None 5156844: status FINISHED
2026-07-14 16:00:29 INFO None 5156846: status FINISHED
2026-07-14 16:00:29 INFO None 5156848: status FINISHED (not in squeue)
2026-07-14 16:00:29 INFO None 5156850: status FINISHED (not in squeue)
2026-07-14 16:00:29 INFO None 5156851: status FINISHED
2026-07-14 16:00:29 INFO None 5156852: status FINISHED
2026-07-14 16:00:29 INFO None 5156853: status FINISHED
2026-07-14 16:00:29 INFO None 5156854: status FINISHED
2026-07-14 16:00:29 INFO None 5156855: status FINISHED
2026-07-14 16:00:29 INFO None 5156857: status FINISHED
2026-07-14 16:00:29 INFO None 5156858: status FINISHED
2026-07-14 16:00:29 INFO None 5156859: status FINISHED
2026-07-14 16:00:29 INFO None 5156860: status FINISHED
2026-07-14 16:00:29 INFO None 5156861: status FINISHED
2026-07-14 16:00:29 INFO Jobs ['5156842', '5156844', '5156846', '5156848', '5156850', '5156851', '5156852', '5156853', '5156854', '5156855', '5156857', '5156858', '5156859', '5156860', '5156861'] have finished
2026-07-14 16:00:29 INFO Checking restart files were created ...
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020714_10_ENS1.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020714_10_ENS2.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020714_10_ENS3.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020714_10_ENS4.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020714_10_ENS5.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020714_10_ENS6.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020714_10_ENS7.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020714_10_ENS8.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020714_10_ENS9.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020714_10_ENS10.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020714_10_ENS11.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020714_10_ENS12.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020714_10_ENS13.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020714_10_ENS14.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020714_10_ENS15.nc(3673513755 bytes)
2026-07-14 16:00:29 INFO  Run_model() completed successfully.
2026-07-14 16:00:29 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 14:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:00:29 INFO [TIME] gregorian_conversion simulated_time=2020-02-08 00:00:00 days=153074 seconds=0
2026-07-14 16:00:29 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-14 16:00:29 INFO [TIME] increment current_time 2020-02-07 14:00:00 -> 2020-02-08 00:00:00
2026-07-14 16:00:29 INFO [TIME] after_increment_before_assimilation current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:00:29 INFO ---------->>> Running process_satellite_data()
2026-07-14 16:00:30 INFO [DART] No satellite data found, skipping assimilation
2026-07-14 16:00:30 INFO after_assimilation() skipped
2026-07-14 16:00:30 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-14 16:00:30 INFO [TIME] step_end current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-14 16:00:30 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
