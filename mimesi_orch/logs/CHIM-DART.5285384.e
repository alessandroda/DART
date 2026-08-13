+ /bin/bash -x /tmp/tmp.ZCfCua37gV
+ SCRIPT_PID=2565887
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
+ unset _mlshdbg
+ return 0
+ python -u main.py -c config/config_irene_IM_cp2.yaml
2026-07-27 13:16:47 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-27 13:16:47 INFO [PIPELINE] =======================================
2026-07-27 13:16:47 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-27 13:16:47 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-27 13:16:47 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-27 13:16:47 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260727_131647.log
2026-07-27 13:16:47 INFO [PIPELINE] =======================================
2026-07-27 13:16:47 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-27 13:16:47 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-27 13:16:47 INFO [STEP] ---- TIME LOOP START ----
2026-07-27 13:16:47 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:16:47 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-27 13:16:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:16:56 INFO Hourly dataset computed and listing created
2026-07-27 13:17:01 INFO Hourly dataset computed
2026-07-27 13:17:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:02 INFO Hourly dataset computed and listing created
2026-07-27 13:17:04 INFO Hourly dataset computed
2026-07-27 13:17:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:05 INFO Hourly dataset computed and listing created
2026-07-27 13:17:06 INFO Hourly dataset computed
2026-07-27 13:17:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:07 INFO Hourly dataset computed and listing created
2026-07-27 13:17:08 INFO Hourly dataset computed
2026-07-27 13:17:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:09 INFO Hourly dataset computed and listing created
2026-07-27 13:17:11 INFO Hourly dataset computed
2026-07-27 13:17:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:12 INFO Hourly dataset computed and listing created
2026-07-27 13:17:13 INFO Hourly dataset computed
2026-07-27 13:17:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:14 INFO Hourly dataset computed and listing created
2026-07-27 13:17:16 INFO Hourly dataset computed
2026-07-27 13:17:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:17 INFO Hourly dataset computed and listing created
2026-07-27 13:17:18 INFO Hourly dataset computed
2026-07-27 13:17:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:19 INFO Hourly dataset computed and listing created
2026-07-27 13:17:21 INFO Hourly dataset computed
2026-07-27 13:17:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:22 INFO Hourly dataset computed and listing created
2026-07-27 13:17:23 INFO Hourly dataset computed
2026-07-27 13:17:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:24 INFO Hourly dataset computed and listing created
2026-07-27 13:17:25 INFO Hourly dataset computed
2026-07-27 13:17:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:25 INFO Hourly dataset computed and listing created
2026-07-27 13:17:26 INFO Hourly dataset computed
2026-07-27 13:17:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:27 INFO Hourly dataset computed and listing created
2026-07-27 13:17:27 INFO Hourly dataset computed
2026-07-27 13:17:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:28 INFO Hourly dataset computed and listing created
2026-07-27 13:17:29 INFO Hourly dataset computed
2026-07-27 13:17:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:17:30 INFO Hourly dataset computed and listing created
2026-07-27 13:17:30 INFO Hourly dataset computed
2026-07-27 13:17:30 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-27 13:17:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:17:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 13:17:30 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-27 13:17:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 13:17:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:17:31 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 13:17:31 INFO Queuing job for member 1...
2026-07-27 13:17:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:17:31 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 13:17:31 INFO Found: ['5285391']
2026-07-27 13:17:36 INFO [TGCC-IRENE] Submitted job with ID:['5285391']
2026-07-27 13:17:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:17:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 13:17:36 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-27 13:17:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 13:17:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:17:36 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 13:17:36 INFO Queuing job for member 2...
2026-07-27 13:17:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:17:36 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 13:17:37 INFO Found: ['5285393']
2026-07-27 13:17:42 INFO [TGCC-IRENE] Submitted job with ID:['5285393']
2026-07-27 13:17:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:17:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 13:17:42 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-27 13:17:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 13:17:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:17:42 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 13:17:42 INFO Queuing job for member 3...
2026-07-27 13:17:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:17:42 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 13:17:43 INFO Found: ['5285394']
2026-07-27 13:17:48 INFO [TGCC-IRENE] Submitted job with ID:['5285394']
2026-07-27 13:17:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:17:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 13:17:48 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-27 13:17:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 13:17:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:17:48 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 13:17:48 INFO Queuing job for member 4...
2026-07-27 13:17:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:17:48 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 13:17:49 INFO Found: ['5285396']
2026-07-27 13:17:54 INFO [TGCC-IRENE] Submitted job with ID:['5285396']
2026-07-27 13:17:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:17:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 13:17:54 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-27 13:17:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 13:17:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:17:54 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 13:17:54 INFO Queuing job for member 5...
2026-07-27 13:17:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:17:54 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 13:17:56 INFO Found: ['5285397']
2026-07-27 13:18:01 INFO [TGCC-IRENE] Submitted job with ID:['5285397']
2026-07-27 13:18:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 13:18:01 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-27 13:18:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 13:18:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:01 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 13:18:01 INFO Queuing job for member 6...
2026-07-27 13:18:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:01 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 13:18:03 INFO Found: ['5285398']
2026-07-27 13:18:08 INFO [TGCC-IRENE] Submitted job with ID:['5285398']
2026-07-27 13:18:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 13:18:08 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-27 13:18:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 13:18:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:08 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 13:18:08 INFO Queuing job for member 7...
2026-07-27 13:18:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:08 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 13:18:11 INFO Found: ['5285401']
2026-07-27 13:18:16 INFO [TGCC-IRENE] Submitted job with ID:['5285401']
2026-07-27 13:18:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 13:18:16 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-27 13:18:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 13:18:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:16 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 13:18:16 INFO Queuing job for member 8...
2026-07-27 13:18:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:16 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 13:18:18 INFO Found: ['5285404']
2026-07-27 13:18:23 INFO [TGCC-IRENE] Submitted job with ID:['5285404']
2026-07-27 13:18:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 13:18:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-27 13:18:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 13:18:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 13:18:23 INFO Queuing job for member 9...
2026-07-27 13:18:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 13:18:26 INFO Found: ['5285407']
2026-07-27 13:18:31 INFO [TGCC-IRENE] Submitted job with ID:['5285407']
2026-07-27 13:18:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 13:18:31 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-27 13:18:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 13:18:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:31 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 13:18:31 INFO Queuing job for member 10...
2026-07-27 13:18:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:31 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 13:18:33 INFO Found: ['5285413']
2026-07-27 13:18:38 INFO [TGCC-IRENE] Submitted job with ID:['5285413']
2026-07-27 13:18:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 13:18:38 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-27 13:18:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 13:18:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:38 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 13:18:38 INFO Queuing job for member 11...
2026-07-27 13:18:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:38 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 13:18:39 INFO Found: ['5285416']
2026-07-27 13:18:44 INFO [TGCC-IRENE] Submitted job with ID:['5285416']
2026-07-27 13:18:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 13:18:44 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-27 13:18:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 13:18:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:44 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 13:18:44 INFO Queuing job for member 12...
2026-07-27 13:18:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:44 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 13:18:45 INFO Found: ['5285417']
2026-07-27 13:18:50 INFO [TGCC-IRENE] Submitted job with ID:['5285417']
2026-07-27 13:18:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 13:18:50 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-27 13:18:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 13:18:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:50 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 13:18:50 INFO Queuing job for member 13...
2026-07-27 13:18:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:50 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 13:18:51 INFO Found: ['5285418']
2026-07-27 13:18:56 INFO [TGCC-IRENE] Submitted job with ID:['5285418']
2026-07-27 13:18:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:18:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 13:18:56 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-27 13:18:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 13:18:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:18:56 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 13:18:56 INFO Queuing job for member 14...
2026-07-27 13:18:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:18:56 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 13:18:56 INFO Found: ['5285420']
2026-07-27 13:19:01 INFO [TGCC-IRENE] Submitted job with ID:['5285420']
2026-07-27 13:19:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:19:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 13:19:01 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-27 13:19:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 13:19:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:19:01 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 13:19:01 INFO Queuing job for member 15...
2026-07-27 13:19:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:19:01 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 13:19:02 INFO Found: ['5285422']
2026-07-27 13:19:07 INFO [TGCC-IRENE] Submitted job with ID:['5285422']
2026-07-27 13:19:07 INFO Checking job status ...
2026-07-27 13:19:07 INFO None 5285391: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285393: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285394: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285396: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285398: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285417: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:19:07 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:19:07 INFO Jobs still running: ['5285391', '5285393', '5285394', '5285396', '5285397', '5285398', '5285401', '5285404', '5285407', '5285413', '5285416', '5285417', '5285418', '5285420', '5285422']. Waiting...
2026-07-27 13:19:22 INFO None 5285391: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285393: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285394: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285396: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285398: status RUNNING/PENDING
2026-07-27 13:19:22 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285417: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:19:23 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:19:23 INFO Jobs still running: ['5285391', '5285393', '5285394', '5285396', '5285397', '5285398', '5285401', '5285404', '5285407', '5285413', '5285416', '5285417', '5285418', '5285420', '5285422']. Waiting...
2026-07-27 13:19:39 INFO None 5285391: status FINISHED
2026-07-27 13:19:39 INFO None 5285393: status FINISHED
2026-07-27 13:19:39 INFO None 5285394: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285396: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285398: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285417: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:19:39 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:19:39 INFO Jobs still running: ['5285394', '5285396', '5285397', '5285398', '5285401', '5285404', '5285407', '5285413', '5285416', '5285417', '5285418', '5285420', '5285422']. Waiting...
2026-07-27 13:19:54 INFO None 5285391: status FINISHED
2026-07-27 13:19:54 INFO None 5285393: status FINISHED
2026-07-27 13:19:54 INFO None 5285394: status FINISHED
2026-07-27 13:19:54 INFO None 5285396: status FINISHED
2026-07-27 13:19:54 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285398: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285417: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:19:54 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:19:54 INFO Jobs still running: ['5285397', '5285398', '5285401', '5285404', '5285407', '5285413', '5285416', '5285417', '5285418', '5285420', '5285422']. Waiting...
2026-07-27 13:20:09 INFO None 5285391: status FINISHED
2026-07-27 13:20:09 INFO None 5285393: status FINISHED
2026-07-27 13:20:09 INFO None 5285394: status FINISHED
2026-07-27 13:20:09 INFO None 5285396: status FINISHED
2026-07-27 13:20:09 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285398: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285417: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:20:11 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:20:11 INFO Jobs still running: ['5285397', '5285398', '5285401', '5285404', '5285407', '5285413', '5285416', '5285417', '5285418', '5285420', '5285422']. Waiting...
2026-07-27 13:20:26 INFO None 5285391: status FINISHED
2026-07-27 13:20:26 INFO None 5285393: status FINISHED
2026-07-27 13:20:26 INFO None 5285394: status FINISHED
2026-07-27 13:20:26 INFO None 5285396: status FINISHED
2026-07-27 13:20:27 INFO None 5285397: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285398: status FINISHED
2026-07-27 13:20:27 INFO None 5285401: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285404: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285407: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285413: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285416: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285417: status FINISHED
2026-07-27 13:20:27 INFO None 5285418: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285420: status RUNNING/PENDING
2026-07-27 13:20:27 INFO None 5285422: status RUNNING/PENDING
2026-07-27 13:20:27 INFO Jobs still running: ['5285397', '5285401', '5285404', '5285407', '5285413', '5285416', '5285418', '5285420', '5285422']. Waiting...
[2026-07-27T13:20:42.110] error: *** JOB 5285384 ON irene4113 CANCELLED AT 2026-07-27T13:20:42 DUE to SIGNAL Terminated ***
