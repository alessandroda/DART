+ /bin/bash -x /tmp/tmp.nS2dUY86OR
+ SCRIPT_PID=624134
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
2026-07-27 15:45:37 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-27 15:45:37 INFO [PIPELINE] =======================================
2026-07-27 15:45:37 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-27 15:45:37 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-27 15:45:37 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-27 15:45:37 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260727_154537.log
2026-07-27 15:45:37 INFO [PIPELINE] =======================================
2026-07-27 15:45:37 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-27 15:45:37 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-27 15:45:37 INFO [STEP] ---- TIME LOOP START ----
2026-07-27 15:45:37 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 15:45:37 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-27 15:45:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:45 INFO Hourly dataset computed and listing created
2026-07-27 15:45:51 INFO Hourly dataset computed
2026-07-27 15:45:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:52 INFO Hourly dataset computed and listing created
2026-07-27 15:45:53 INFO Hourly dataset computed
2026-07-27 15:45:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:54 INFO Hourly dataset computed and listing created
2026-07-27 15:45:55 INFO Hourly dataset computed
2026-07-27 15:45:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:55 INFO Hourly dataset computed and listing created
2026-07-27 15:45:56 INFO Hourly dataset computed
2026-07-27 15:45:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:57 INFO Hourly dataset computed and listing created
2026-07-27 15:45:58 INFO Hourly dataset computed
2026-07-27 15:45:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:45:59 INFO Hourly dataset computed and listing created
2026-07-27 15:46:00 INFO Hourly dataset computed
2026-07-27 15:46:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:00 INFO Hourly dataset computed and listing created
2026-07-27 15:46:01 INFO Hourly dataset computed
2026-07-27 15:46:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:02 INFO Hourly dataset computed and listing created
2026-07-27 15:46:03 INFO Hourly dataset computed
2026-07-27 15:46:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:04 INFO Hourly dataset computed and listing created
2026-07-27 15:46:05 INFO Hourly dataset computed
2026-07-27 15:46:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:05 INFO Hourly dataset computed and listing created
2026-07-27 15:46:06 INFO Hourly dataset computed
2026-07-27 15:46:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:07 INFO Hourly dataset computed and listing created
2026-07-27 15:46:08 INFO Hourly dataset computed
2026-07-27 15:46:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:09 INFO Hourly dataset computed and listing created
2026-07-27 15:46:10 INFO Hourly dataset computed
2026-07-27 15:46:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:11 INFO Hourly dataset computed and listing created
2026-07-27 15:46:12 INFO Hourly dataset computed
2026-07-27 15:46:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:13 INFO Hourly dataset computed and listing created
2026-07-27 15:46:13 INFO Hourly dataset computed
2026-07-27 15:46:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:46:14 INFO Hourly dataset computed and listing created
2026-07-27 15:46:15 INFO Hourly dataset computed
2026-07-27 15:46:15 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-27 15:46:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 15:46:15 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-27 15:46:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 15:46:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:15 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 15:46:15 INFO Queuing job for member 1...
2026-07-27 15:46:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:15 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 15:46:16 INFO Found: ['5287539']
2026-07-27 15:46:21 INFO [TGCC-IRENE] Submitted job with ID:['5287539']
2026-07-27 15:46:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 15:46:21 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-27 15:46:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 15:46:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:21 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 15:46:21 INFO Queuing job for member 2...
2026-07-27 15:46:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:21 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 15:46:22 INFO Found: ['5287540']
2026-07-27 15:46:27 INFO [TGCC-IRENE] Submitted job with ID:['5287540']
2026-07-27 15:46:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 15:46:27 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-27 15:46:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 15:46:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:27 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 15:46:27 INFO Queuing job for member 3...
2026-07-27 15:46:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:27 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 15:46:30 INFO Found: ['5287541']
2026-07-27 15:46:35 INFO [TGCC-IRENE] Submitted job with ID:['5287541']
2026-07-27 15:46:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 15:46:35 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-27 15:46:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 15:46:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:35 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 15:46:35 INFO Queuing job for member 4...
2026-07-27 15:46:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:35 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 15:46:37 INFO Found: ['5287543']
2026-07-27 15:46:42 INFO [TGCC-IRENE] Submitted job with ID:['5287543']
2026-07-27 15:46:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 15:46:42 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-27 15:46:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 15:46:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:42 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 15:46:42 INFO Queuing job for member 5...
2026-07-27 15:46:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:42 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 15:46:45 INFO Found: ['5287545']
2026-07-27 15:46:50 INFO [TGCC-IRENE] Submitted job with ID:['5287545']
2026-07-27 15:46:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 15:46:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-27 15:46:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 15:46:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 15:46:50 INFO Queuing job for member 6...
2026-07-27 15:46:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 15:46:52 INFO Found: ['5287547']
2026-07-27 15:46:57 INFO [TGCC-IRENE] Submitted job with ID:['5287547']
2026-07-27 15:46:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:46:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 15:46:57 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-27 15:46:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 15:46:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:46:57 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 15:46:57 INFO Queuing job for member 7...
2026-07-27 15:46:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:46:57 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 15:47:00 INFO Found: ['5287549']
2026-07-27 15:47:05 INFO [TGCC-IRENE] Submitted job with ID:['5287549']
2026-07-27 15:47:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 15:47:05 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-27 15:47:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 15:47:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:05 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 15:47:05 INFO Queuing job for member 8...
2026-07-27 15:47:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:05 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 15:47:07 INFO Found: ['5287552']
2026-07-27 15:47:12 INFO [TGCC-IRENE] Submitted job with ID:['5287552']
2026-07-27 15:47:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 15:47:12 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-27 15:47:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 15:47:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:12 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 15:47:12 INFO Queuing job for member 9...
2026-07-27 15:47:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:12 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 15:47:15 INFO Found: ['5287558']
2026-07-27 15:47:20 INFO [TGCC-IRENE] Submitted job with ID:['5287558']
2026-07-27 15:47:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 15:47:20 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-27 15:47:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 15:47:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:20 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 15:47:20 INFO Queuing job for member 10...
2026-07-27 15:47:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:20 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 15:47:22 INFO Found: ['5287559']
2026-07-27 15:47:27 INFO [TGCC-IRENE] Submitted job with ID:['5287559']
2026-07-27 15:47:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 15:47:27 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-27 15:47:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 15:47:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:27 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 15:47:27 INFO Queuing job for member 11...
2026-07-27 15:47:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:27 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 15:47:30 INFO Found: ['5287560']
2026-07-27 15:47:35 INFO [TGCC-IRENE] Submitted job with ID:['5287560']
2026-07-27 15:47:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 15:47:35 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-27 15:47:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 15:47:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:35 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 15:47:35 INFO Queuing job for member 12...
2026-07-27 15:47:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:35 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 15:47:37 INFO Found: ['5287561']
2026-07-27 15:47:42 INFO [TGCC-IRENE] Submitted job with ID:['5287561']
2026-07-27 15:47:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 15:47:42 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-27 15:47:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 15:47:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:42 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 15:47:42 INFO Queuing job for member 13...
2026-07-27 15:47:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:42 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 15:47:45 INFO Found: ['5287562']
2026-07-27 15:47:50 INFO [TGCC-IRENE] Submitted job with ID:['5287562']
2026-07-27 15:47:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 15:47:50 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-27 15:47:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 15:47:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:50 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 15:47:50 INFO Queuing job for member 14...
2026-07-27 15:47:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:50 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 15:47:51 INFO Found: ['5287563']
2026-07-27 15:47:56 INFO [TGCC-IRENE] Submitted job with ID:['5287563']
2026-07-27 15:47:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:47:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 15:47:56 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-27 15:47:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 15:47:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:47:56 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 15:47:56 INFO Queuing job for member 15...
2026-07-27 15:47:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:47:56 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 15:47:57 INFO Found: ['5287564']
2026-07-27 15:48:02 INFO [TGCC-IRENE] Submitted job with ID:['5287564']
2026-07-27 15:48:02 INFO Checking job status ...
2026-07-27 15:48:02 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:48:02 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:48:02 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:48:17 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:48:17 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:48:17 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:48:32 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:48:32 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:48:32 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:48:48 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:48:48 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:48:48 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:49:03 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:49:03 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:49:03 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:49:18 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:49:18 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:49:18 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:49:18 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:49:20 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:49:20 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:49:35 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:49:35 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:49:35 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:49:35 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:49:35 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:49:36 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:49:36 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:49:51 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:49:51 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:49:51 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:49:51 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:49:51 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:49:51 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:49:53 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:49:53 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:50:08 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:50:08 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:50:08 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:50:23 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:50:23 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:50:24 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:50:24 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:50:39 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:50:39 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:50:39 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:50:55 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:50:55 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:50:55 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:51:10 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:51:10 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:51:12 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:51:12 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:51:12 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:51:12 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:51:27 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:51:27 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:51:27 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:51:43 INFO None 5287539: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:51:43 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:51:43 INFO Jobs still running: ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:51:58 INFO None 5287539: status FINISHED
2026-07-27 15:51:58 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:52:00 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:52:00 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:52:15 INFO None 5287539: status FINISHED
2026-07-27 15:52:15 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:52:15 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:52:15 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:52:30 INFO None 5287539: status FINISHED
2026-07-27 15:52:30 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:52:30 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:52:30 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:52:30 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:52:30 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:52:31 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:52:31 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:52:46 INFO None 5287539: status FINISHED
2026-07-27 15:52:46 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287545: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287547: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287561: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:52:46 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:52:46 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:53:01 INFO None 5287539: status FINISHED
2026-07-27 15:53:01 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287545: status FINISHED
2026-07-27 15:53:01 INFO None 5287547: status FINISHED
2026-07-27 15:53:01 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287560: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287561: status FINISHED
2026-07-27 15:53:01 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287563: status RUNNING/PENDING
2026-07-27 15:53:01 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:53:01 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287549', '5287552', '5287558', '5287559', '5287560', '5287562', '5287563', '5287564']. Waiting...
2026-07-27 15:53:18 INFO None 5287539: status FINISHED
2026-07-27 15:53:18 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287545: status FINISHED
2026-07-27 15:53:18 INFO None 5287547: status FINISHED
2026-07-27 15:53:18 INFO None 5287549: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287552: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287560: status FINISHED
2026-07-27 15:53:18 INFO None 5287561: status FINISHED
2026-07-27 15:53:18 INFO None 5287562: status RUNNING/PENDING
2026-07-27 15:53:18 INFO None 5287563: status FINISHED
2026-07-27 15:53:18 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:53:18 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287549', '5287552', '5287558', '5287559', '5287562', '5287564']. Waiting...
2026-07-27 15:53:33 INFO None 5287539: status FINISHED
2026-07-27 15:53:33 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:53:33 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:53:33 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:53:33 INFO None 5287545: status FINISHED
2026-07-27 15:53:33 INFO None 5287547: status FINISHED
2026-07-27 15:53:33 INFO None 5287549: status FINISHED
2026-07-27 15:53:33 INFO None 5287552: status FINISHED
2026-07-27 15:53:33 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:53:33 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:53:33 INFO None 5287560: status FINISHED
2026-07-27 15:53:33 INFO None 5287561: status FINISHED
2026-07-27 15:53:33 INFO None 5287562: status FINISHED
2026-07-27 15:53:33 INFO None 5287563: status FINISHED
2026-07-27 15:53:33 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:53:33 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287558', '5287559', '5287564']. Waiting...
2026-07-27 15:53:48 INFO None 5287539: status FINISHED
2026-07-27 15:53:48 INFO None 5287540: status RUNNING/PENDING
2026-07-27 15:53:48 INFO None 5287541: status RUNNING/PENDING
2026-07-27 15:53:50 INFO None 5287543: status RUNNING/PENDING
2026-07-27 15:53:50 INFO None 5287545: status FINISHED
2026-07-27 15:53:50 INFO None 5287547: status FINISHED
2026-07-27 15:53:50 INFO None 5287549: status FINISHED
2026-07-27 15:53:50 INFO None 5287552: status FINISHED
2026-07-27 15:53:50 INFO None 5287558: status RUNNING/PENDING
2026-07-27 15:53:50 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:53:50 INFO None 5287560: status FINISHED
2026-07-27 15:53:50 INFO None 5287561: status FINISHED
2026-07-27 15:53:50 INFO None 5287562: status FINISHED
2026-07-27 15:53:50 INFO None 5287563: status FINISHED
2026-07-27 15:53:50 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:53:50 INFO Jobs still running: ['5287540', '5287541', '5287543', '5287558', '5287559', '5287564']. Waiting...
2026-07-27 15:54:05 INFO None 5287539: status FINISHED
2026-07-27 15:54:05 INFO None 5287540: status FINISHED
2026-07-27 15:54:05 INFO None 5287541: status FINISHED
2026-07-27 15:54:05 INFO None 5287543: status FINISHED
2026-07-27 15:54:05 INFO None 5287545: status FINISHED
2026-07-27 15:54:05 INFO None 5287547: status FINISHED
2026-07-27 15:54:06 INFO None 5287549: status FINISHED
2026-07-27 15:54:06 INFO None 5287552: status FINISHED
2026-07-27 15:54:06 INFO None 5287558: status FINISHED
2026-07-27 15:54:06 INFO None 5287559: status RUNNING/PENDING
2026-07-27 15:54:06 INFO None 5287560: status FINISHED
2026-07-27 15:54:06 INFO None 5287561: status FINISHED
2026-07-27 15:54:06 INFO None 5287562: status FINISHED
2026-07-27 15:54:06 INFO None 5287563: status FINISHED
2026-07-27 15:54:06 INFO None 5287564: status RUNNING/PENDING
2026-07-27 15:54:06 INFO Jobs still running: ['5287559', '5287564']. Waiting...
2026-07-27 15:54:21 INFO None 5287539: status FINISHED
2026-07-27 15:54:21 INFO None 5287540: status FINISHED
2026-07-27 15:54:21 INFO None 5287541: status FINISHED
2026-07-27 15:54:21 INFO None 5287543: status FINISHED
2026-07-27 15:54:23 INFO None 5287545: status FINISHED
2026-07-27 15:54:23 INFO None 5287547: status FINISHED
2026-07-27 15:54:23 INFO None 5287549: status FINISHED
2026-07-27 15:54:23 INFO None 5287552: status FINISHED
2026-07-27 15:54:23 INFO None 5287558: status FINISHED
2026-07-27 15:54:23 INFO None 5287559: status FINISHED
2026-07-27 15:54:23 INFO None 5287560: status FINISHED
2026-07-27 15:54:23 INFO None 5287561: status FINISHED
2026-07-27 15:54:23 INFO None 5287562: status FINISHED
2026-07-27 15:54:23 INFO None 5287563: status FINISHED
2026-07-27 15:54:23 INFO None 5287564: status FINISHED
2026-07-27 15:54:23 INFO Jobs ['5287539', '5287540', '5287541', '5287543', '5287545', '5287547', '5287549', '5287552', '5287558', '5287559', '5287560', '5287561', '5287562', '5287563', '5287564'] have finished
2026-07-27 15:54:23 INFO Checking restart files were created ...
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-27 15:54:23 INFO  Run_model() completed successfully.
2026-07-27 15:54:23 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 15:54:23 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-27 15:54:23 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 15:54:23 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-27 15:54:23 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 15:54:23 INFO ---------->>> Running process_satellite_data()
2026-07-27 15:54:23 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-27 15:54:23 INFO ---------->>> Running run_obs_converter()
2026-07-27 15:54:23 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-27 15:54:23 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-27 15:54:23 INFO ---------->>> Running DART
2026-07-27 15:54:23 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-27 15:54:23 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 15:54:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 15:54:23 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 15:54:23 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 15:54:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 15:54:23 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 15:54:52 INFO Found: []
2026-07-27 15:54:52 INFO No job id returned by command ./run_filter.bsh
2026-07-27 15:54:52 INFO No monitoring will be performed
2026-07-27 15:54:52 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-27 15:54:52 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:52 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:52 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 15:54:53 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:53 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 15:54:53 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 15:54:55 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 15:54:55 INFO run_dart() is DONE.
2026-07-27 15:54:55 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 15:54:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-27 15:55:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-27 15:55:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-27 15:55:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-27 15:55:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-27 15:55:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-27 15:55:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-27 15:55:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-27 15:55:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-27 15:55:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-27 15:55:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-27 15:55:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:54 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-27 15:55:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:55:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-27 15:56:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:56:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-27 15:56:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:56:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-27 15:56:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 15:56:15 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 15:56:15 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 15:56:15 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 15:56:15 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-27 15:56:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:17 INFO Hourly dataset computed and listing created
2026-07-27 15:56:22 INFO Hourly dataset computed
2026-07-27 15:56:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:22 INFO Hourly dataset computed and listing created
2026-07-27 15:56:23 INFO Hourly dataset computed
2026-07-27 15:56:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:24 INFO Hourly dataset computed and listing created
2026-07-27 15:56:25 INFO Hourly dataset computed
2026-07-27 15:56:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:26 INFO Hourly dataset computed and listing created
2026-07-27 15:56:27 INFO Hourly dataset computed
2026-07-27 15:56:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:28 INFO Hourly dataset computed and listing created
2026-07-27 15:56:28 INFO Hourly dataset computed
2026-07-27 15:56:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:29 INFO Hourly dataset computed and listing created
2026-07-27 15:56:30 INFO Hourly dataset computed
2026-07-27 15:56:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:31 INFO Hourly dataset computed and listing created
2026-07-27 15:56:32 INFO Hourly dataset computed
2026-07-27 15:56:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:33 INFO Hourly dataset computed and listing created
2026-07-27 15:56:33 INFO Hourly dataset computed
2026-07-27 15:56:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:34 INFO Hourly dataset computed and listing created
2026-07-27 15:56:35 INFO Hourly dataset computed
2026-07-27 15:56:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:36 INFO Hourly dataset computed and listing created
2026-07-27 15:56:37 INFO Hourly dataset computed
2026-07-27 15:56:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:38 INFO Hourly dataset computed and listing created
2026-07-27 15:56:39 INFO Hourly dataset computed
2026-07-27 15:56:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:39 INFO Hourly dataset computed and listing created
2026-07-27 15:56:40 INFO Hourly dataset computed
2026-07-27 15:56:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:41 INFO Hourly dataset computed and listing created
2026-07-27 15:56:42 INFO Hourly dataset computed
2026-07-27 15:56:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:43 INFO Hourly dataset computed and listing created
2026-07-27 15:56:44 INFO Hourly dataset computed
2026-07-27 15:56:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 15:56:45 INFO Hourly dataset computed and listing created
2026-07-27 15:56:46 INFO Hourly dataset computed
2026-07-27 15:56:46 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-27 15:56:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:56:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 15:56:46 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-27 15:56:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 15:56:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:56:46 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 15:56:46 INFO Queuing job for member 1...
2026-07-27 15:56:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:56:46 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 15:56:47 INFO Found: ['5287611']
2026-07-27 15:56:52 INFO [TGCC-IRENE] Submitted job with ID:['5287611']
2026-07-27 15:56:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:56:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 15:56:52 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-27 15:56:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 15:56:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:56:52 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 15:56:53 INFO Queuing job for member 2...
2026-07-27 15:56:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:56:53 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 15:56:53 INFO Found: ['5287613']
2026-07-27 15:56:58 INFO [TGCC-IRENE] Submitted job with ID:['5287613']
2026-07-27 15:56:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:56:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 15:56:58 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-27 15:56:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 15:56:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:56:58 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 15:56:58 INFO Queuing job for member 3...
2026-07-27 15:56:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:56:58 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 15:56:59 INFO Found: ['5287614']
2026-07-27 15:57:04 INFO [TGCC-IRENE] Submitted job with ID:['5287614']
2026-07-27 15:57:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 15:57:04 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-27 15:57:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 15:57:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:04 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 15:57:04 INFO Queuing job for member 4...
2026-07-27 15:57:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:04 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 15:57:05 INFO Found: ['5287616']
2026-07-27 15:57:10 INFO [TGCC-IRENE] Submitted job with ID:['5287616']
2026-07-27 15:57:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 15:57:10 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-27 15:57:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 15:57:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:10 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 15:57:10 INFO Queuing job for member 5...
2026-07-27 15:57:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:10 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 15:57:10 INFO Found: ['5287617']
2026-07-27 15:57:15 INFO [TGCC-IRENE] Submitted job with ID:['5287617']
2026-07-27 15:57:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 15:57:15 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-27 15:57:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 15:57:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:15 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 15:57:15 INFO Queuing job for member 6...
2026-07-27 15:57:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:15 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 15:57:16 INFO Found: ['5287618']
2026-07-27 15:57:21 INFO [TGCC-IRENE] Submitted job with ID:['5287618']
2026-07-27 15:57:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 15:57:21 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-27 15:57:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 15:57:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:21 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 15:57:21 INFO Queuing job for member 7...
2026-07-27 15:57:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:21 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 15:57:22 INFO Found: ['5287619']
2026-07-27 15:57:27 INFO [TGCC-IRENE] Submitted job with ID:['5287619']
2026-07-27 15:57:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 15:57:27 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-27 15:57:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 15:57:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:27 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 15:57:27 INFO Queuing job for member 8...
2026-07-27 15:57:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:27 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 15:57:28 INFO Found: ['5287621']
2026-07-27 15:57:33 INFO [TGCC-IRENE] Submitted job with ID:['5287621']
2026-07-27 15:57:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 15:57:33 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-27 15:57:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 15:57:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:33 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 15:57:33 INFO Queuing job for member 9...
2026-07-27 15:57:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:33 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 15:57:35 INFO Found: ['5287622']
2026-07-27 15:57:40 INFO [TGCC-IRENE] Submitted job with ID:['5287622']
2026-07-27 15:57:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 15:57:40 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-27 15:57:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 15:57:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:40 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 15:57:40 INFO Queuing job for member 10...
2026-07-27 15:57:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:40 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 15:57:43 INFO Found: ['5287623']
2026-07-27 15:57:48 INFO [TGCC-IRENE] Submitted job with ID:['5287623']
2026-07-27 15:57:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 15:57:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-27 15:57:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 15:57:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 15:57:48 INFO Queuing job for member 11...
2026-07-27 15:57:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 15:57:50 INFO Found: ['5287624']
2026-07-27 15:57:55 INFO [TGCC-IRENE] Submitted job with ID:['5287624']
2026-07-27 15:57:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:57:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 15:57:55 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-27 15:57:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 15:57:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:57:55 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 15:57:55 INFO Queuing job for member 12...
2026-07-27 15:57:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:57:55 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 15:57:58 INFO Found: ['5287625']
2026-07-27 15:58:03 INFO [TGCC-IRENE] Submitted job with ID:['5287625']
2026-07-27 15:58:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:58:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 15:58:03 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-27 15:58:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 15:58:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:58:03 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 15:58:03 INFO Queuing job for member 13...
2026-07-27 15:58:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:58:03 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 15:58:05 INFO Found: ['5287628']
2026-07-27 15:58:10 INFO [TGCC-IRENE] Submitted job with ID:['5287628']
2026-07-27 15:58:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:58:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 15:58:10 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-27 15:58:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 15:58:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:58:10 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 15:58:10 INFO Queuing job for member 14...
2026-07-27 15:58:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:58:10 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 15:58:13 INFO Found: ['5287629']
2026-07-27 15:58:18 INFO [TGCC-IRENE] Submitted job with ID:['5287629']
2026-07-27 15:58:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 15:58:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 15:58:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-27 15:58:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 15:58:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 15:58:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 15:58:18 INFO Queuing job for member 15...
2026-07-27 15:58:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 15:58:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 15:58:20 INFO Found: ['5287630']
2026-07-27 15:58:25 INFO [TGCC-IRENE] Submitted job with ID:['5287630']
2026-07-27 15:58:25 INFO Checking job status ...
2026-07-27 15:58:25 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:58:25 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:58:25 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 15:58:41 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:58:41 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:58:43 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:58:43 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 15:58:58 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:58:58 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:58:58 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 15:59:13 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:59:13 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:59:13 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 15:59:28 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:59:28 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:59:28 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:59:28 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:59:28 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:59:29 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:59:29 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 15:59:44 INFO None 5287611: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287613: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287614: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287616: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287617: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287618: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287619: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287621: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287622: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287623: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287624: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287625: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287628: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287629: status RUNNING/PENDING
2026-07-27 15:59:44 INFO None 5287630: status RUNNING/PENDING
2026-07-27 15:59:44 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:00:00 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:00:00 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:00:00 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:00:15 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:00:15 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:00:15 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:00:30 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:00:30 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:00:30 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:00:30 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:00:30 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:00:32 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:00:32 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:00:47 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:00:47 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:00:47 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:00:48 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:00:48 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:01:03 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:01:03 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:01:05 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:01:05 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:01:20 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:01:20 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:01:20 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:01:35 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287617: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:01:35 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:01:36 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:01:36 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:01:51 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287617: status FINISHED
2026-07-27 16:01:51 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287621: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:01:51 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:01:51 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:02:07 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287613: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287614: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287616: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287617: status FINISHED
2026-07-27 16:02:08 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287621: status FINISHED
2026-07-27 16:02:08 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:02:08 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:02:08 INFO Jobs still running: ['5287611', '5287613', '5287614', '5287616', '5287618', '5287619', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:02:23 INFO None 5287611: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287613: status FINISHED
2026-07-27 16:02:23 INFO None 5287614: status FINISHED
2026-07-27 16:02:23 INFO None 5287616: status FINISHED
2026-07-27 16:02:23 INFO None 5287617: status FINISHED
2026-07-27 16:02:23 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287621: status FINISHED
2026-07-27 16:02:23 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:02:23 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:02:23 INFO Jobs still running: ['5287611', '5287618', '5287619', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:02:38 INFO None 5287611: status FINISHED
2026-07-27 16:02:38 INFO None 5287613: status FINISHED
2026-07-27 16:02:38 INFO None 5287614: status FINISHED
2026-07-27 16:02:38 INFO None 5287616: status FINISHED
2026-07-27 16:02:38 INFO None 5287617: status FINISHED
2026-07-27 16:02:40 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287621: status FINISHED
2026-07-27 16:02:40 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:02:40 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:02:40 INFO Jobs still running: ['5287618', '5287619', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:02:55 INFO None 5287611: status FINISHED
2026-07-27 16:02:55 INFO None 5287613: status FINISHED
2026-07-27 16:02:55 INFO None 5287614: status FINISHED
2026-07-27 16:02:55 INFO None 5287616: status FINISHED
2026-07-27 16:02:55 INFO None 5287617: status FINISHED
2026-07-27 16:02:55 INFO None 5287618: status RUNNING/PENDING
2026-07-27 16:02:55 INFO None 5287619: status RUNNING/PENDING
2026-07-27 16:02:55 INFO None 5287621: status FINISHED
2026-07-27 16:02:55 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:02:55 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:02:56 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:02:56 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:02:56 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:02:56 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:02:56 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:02:56 INFO Jobs still running: ['5287618', '5287619', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:03:11 INFO None 5287611: status FINISHED
2026-07-27 16:03:11 INFO None 5287613: status FINISHED
2026-07-27 16:03:11 INFO None 5287614: status FINISHED
2026-07-27 16:03:11 INFO None 5287616: status FINISHED
2026-07-27 16:03:11 INFO None 5287617: status FINISHED
2026-07-27 16:03:11 INFO None 5287618: status FINISHED
2026-07-27 16:03:11 INFO None 5287619: status FINISHED
2026-07-27 16:03:11 INFO None 5287621: status FINISHED
2026-07-27 16:03:13 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:03:13 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:03:13 INFO Jobs still running: ['5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:03:28 INFO None 5287611: status FINISHED
2026-07-27 16:03:28 INFO None 5287613: status FINISHED
2026-07-27 16:03:28 INFO None 5287614: status FINISHED
2026-07-27 16:03:28 INFO None 5287616: status FINISHED
2026-07-27 16:03:28 INFO None 5287617: status FINISHED
2026-07-27 16:03:28 INFO None 5287618: status FINISHED
2026-07-27 16:03:28 INFO None 5287619: status FINISHED
2026-07-27 16:03:28 INFO None 5287621: status FINISHED
2026-07-27 16:03:28 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:03:28 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:03:28 INFO Jobs still running: ['5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:03:43 INFO None 5287611: status FINISHED
2026-07-27 16:03:43 INFO None 5287613: status FINISHED
2026-07-27 16:03:43 INFO None 5287614: status FINISHED
2026-07-27 16:03:43 INFO None 5287616: status FINISHED
2026-07-27 16:03:43 INFO None 5287617: status FINISHED
2026-07-27 16:03:43 INFO None 5287618: status FINISHED
2026-07-27 16:03:43 INFO None 5287619: status FINISHED
2026-07-27 16:03:43 INFO None 5287621: status FINISHED
2026-07-27 16:03:43 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:03:43 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:03:43 INFO Jobs still running: ['5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:03:58 INFO None 5287611: status FINISHED
2026-07-27 16:03:58 INFO None 5287613: status FINISHED
2026-07-27 16:03:58 INFO None 5287614: status FINISHED
2026-07-27 16:03:59 INFO None 5287616: status FINISHED
2026-07-27 16:03:59 INFO None 5287617: status FINISHED
2026-07-27 16:03:59 INFO None 5287618: status FINISHED
2026-07-27 16:03:59 INFO None 5287619: status FINISHED
2026-07-27 16:03:59 INFO None 5287621: status FINISHED
2026-07-27 16:03:59 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287623: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:03:59 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:03:59 INFO Jobs still running: ['5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:04:14 INFO None 5287611: status FINISHED
2026-07-27 16:04:14 INFO None 5287613: status FINISHED
2026-07-27 16:04:14 INFO None 5287614: status FINISHED
2026-07-27 16:04:14 INFO None 5287616: status FINISHED
2026-07-27 16:04:14 INFO None 5287617: status FINISHED
2026-07-27 16:04:14 INFO None 5287618: status FINISHED
2026-07-27 16:04:14 INFO None 5287619: status FINISHED
2026-07-27 16:04:14 INFO None 5287621: status FINISHED
2026-07-27 16:04:14 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:04:14 INFO None 5287623: status FINISHED
2026-07-27 16:04:14 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:04:14 INFO None 5287625: status RUNNING/PENDING
2026-07-27 16:04:14 INFO None 5287628: status RUNNING/PENDING
2026-07-27 16:04:14 INFO None 5287629: status RUNNING/PENDING
2026-07-27 16:04:14 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:04:14 INFO Jobs still running: ['5287622', '5287624', '5287625', '5287628', '5287629', '5287630']. Waiting...
2026-07-27 16:04:30 INFO None 5287611: status FINISHED
2026-07-27 16:04:30 INFO None 5287613: status FINISHED
2026-07-27 16:04:30 INFO None 5287614: status FINISHED
2026-07-27 16:04:30 INFO None 5287616: status FINISHED
2026-07-27 16:04:30 INFO None 5287617: status FINISHED
2026-07-27 16:04:30 INFO None 5287618: status FINISHED
2026-07-27 16:04:30 INFO None 5287619: status FINISHED
2026-07-27 16:04:30 INFO None 5287621: status FINISHED
2026-07-27 16:04:30 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:04:30 INFO None 5287623: status FINISHED
2026-07-27 16:04:30 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:04:30 INFO None 5287625: status FINISHED
2026-07-27 16:04:30 INFO None 5287628: status FINISHED
2026-07-27 16:04:30 INFO None 5287629: status FINISHED
2026-07-27 16:04:30 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:04:30 INFO Jobs still running: ['5287622', '5287624', '5287630']. Waiting...
2026-07-27 16:04:45 INFO None 5287611: status FINISHED
2026-07-27 16:04:45 INFO None 5287613: status FINISHED
2026-07-27 16:04:45 INFO None 5287614: status FINISHED
2026-07-27 16:04:45 INFO None 5287616: status FINISHED
2026-07-27 16:04:45 INFO None 5287617: status FINISHED
2026-07-27 16:04:45 INFO None 5287618: status FINISHED
2026-07-27 16:04:45 INFO None 5287619: status FINISHED
2026-07-27 16:04:45 INFO None 5287621: status FINISHED
2026-07-27 16:04:45 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:04:45 INFO None 5287623: status FINISHED
2026-07-27 16:04:45 INFO None 5287624: status RUNNING/PENDING
2026-07-27 16:04:45 INFO None 5287625: status FINISHED
2026-07-27 16:04:45 INFO None 5287628: status FINISHED
2026-07-27 16:04:45 INFO None 5287629: status FINISHED
2026-07-27 16:04:45 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:04:45 INFO Jobs still running: ['5287622', '5287624', '5287630']. Waiting...
2026-07-27 16:05:00 INFO None 5287611: status FINISHED
2026-07-27 16:05:00 INFO None 5287613: status FINISHED
2026-07-27 16:05:00 INFO None 5287614: status FINISHED
2026-07-27 16:05:00 INFO None 5287616: status FINISHED
2026-07-27 16:05:00 INFO None 5287617: status FINISHED
2026-07-27 16:05:02 INFO None 5287618: status FINISHED
2026-07-27 16:05:02 INFO None 5287619: status FINISHED
2026-07-27 16:05:02 INFO None 5287621: status FINISHED
2026-07-27 16:05:02 INFO None 5287622: status RUNNING/PENDING
2026-07-27 16:05:02 INFO None 5287623: status FINISHED
2026-07-27 16:05:02 INFO None 5287624: status FINISHED
2026-07-27 16:05:02 INFO None 5287625: status FINISHED
2026-07-27 16:05:02 INFO None 5287628: status FINISHED
2026-07-27 16:05:02 INFO None 5287629: status FINISHED
2026-07-27 16:05:02 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:05:02 INFO Jobs still running: ['5287622', '5287630']. Waiting...
2026-07-27 16:05:17 INFO None 5287611: status FINISHED
2026-07-27 16:05:18 INFO None 5287613: status FINISHED
2026-07-27 16:05:18 INFO None 5287614: status FINISHED
2026-07-27 16:05:18 INFO None 5287616: status FINISHED
2026-07-27 16:05:18 INFO None 5287617: status FINISHED
2026-07-27 16:05:18 INFO None 5287618: status FINISHED
2026-07-27 16:05:18 INFO None 5287619: status FINISHED
2026-07-27 16:05:18 INFO None 5287621: status FINISHED
2026-07-27 16:05:18 INFO None 5287622: status FINISHED
2026-07-27 16:05:18 INFO None 5287623: status FINISHED
2026-07-27 16:05:18 INFO None 5287624: status FINISHED
2026-07-27 16:05:18 INFO None 5287625: status FINISHED
2026-07-27 16:05:18 INFO None 5287628: status FINISHED
2026-07-27 16:05:18 INFO None 5287629: status FINISHED
2026-07-27 16:05:18 INFO None 5287630: status RUNNING/PENDING
2026-07-27 16:05:18 INFO Jobs still running: ['5287630']. Waiting...
2026-07-27 16:05:33 INFO None 5287611: status FINISHED
2026-07-27 16:05:33 INFO None 5287613: status FINISHED
2026-07-27 16:05:33 INFO None 5287614: status FINISHED
2026-07-27 16:05:33 INFO None 5287616: status FINISHED
2026-07-27 16:05:33 INFO None 5287617: status FINISHED
2026-07-27 16:05:33 INFO None 5287618: status FINISHED
2026-07-27 16:05:33 INFO None 5287619: status FINISHED
2026-07-27 16:05:35 INFO None 5287621: status FINISHED
2026-07-27 16:05:35 INFO None 5287622: status FINISHED
2026-07-27 16:05:35 INFO None 5287623: status FINISHED
2026-07-27 16:05:35 INFO None 5287624: status FINISHED
2026-07-27 16:05:35 INFO None 5287625: status FINISHED
2026-07-27 16:05:35 INFO None 5287628: status FINISHED
2026-07-27 16:05:35 INFO None 5287629: status FINISHED
2026-07-27 16:05:35 INFO None 5287630: status FINISHED
2026-07-27 16:05:35 INFO Jobs ['5287611', '5287613', '5287614', '5287616', '5287617', '5287618', '5287619', '5287621', '5287622', '5287623', '5287624', '5287625', '5287628', '5287629', '5287630'] have finished
2026-07-27 16:05:35 INFO Checking restart files were created ...
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-27 16:05:35 INFO  Run_model() completed successfully.
2026-07-27 16:05:35 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:05:35 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-27 16:05:35 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 16:05:35 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-27 16:05:35 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:05:35 INFO ---------->>> Running process_satellite_data()
2026-07-27 16:05:35 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-27 16:05:35 INFO ---------->>> Running run_obs_converter()
2026-07-27 16:05:35 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-27 16:05:35 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-27 16:05:35 INFO ---------->>> Running DART
2026-07-27 16:05:35 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-27 16:05:35 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 16:05:35 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 16:05:35 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 16:05:35 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 16:05:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 16:05:35 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 16:06:06 INFO Found: []
2026-07-27 16:06:06 INFO No job id returned by command ./run_filter.bsh
2026-07-27 16:06:06 INFO No monitoring will be performed
2026-07-27 16:06:06 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-27 16:06:06 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:06 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:06 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:06 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:06 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 16:06:07 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 16:06:07 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 16:06:07 INFO run_dart() is DONE.
2026-07-27 16:06:07 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 16:06:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-27 16:06:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-27 16:06:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-27 16:06:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-27 16:06:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:29 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-27 16:06:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-27 16:06:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-27 16:06:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-27 16:06:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-27 16:06:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:06:57 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-27 16:07:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-27 16:07:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-27 16:07:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:12 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-27 16:07:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-27 16:07:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-27 16:07:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:07:28 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 16:07:28 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:07:28 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:07:28 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-27 16:07:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:30 INFO Hourly dataset computed and listing created
2026-07-27 16:07:32 INFO Hourly dataset computed
2026-07-27 16:07:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:33 INFO Hourly dataset computed and listing created
2026-07-27 16:07:33 INFO Hourly dataset computed
2026-07-27 16:07:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:34 INFO Hourly dataset computed and listing created
2026-07-27 16:07:35 INFO Hourly dataset computed
2026-07-27 16:07:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:36 INFO Hourly dataset computed and listing created
2026-07-27 16:07:36 INFO Hourly dataset computed
2026-07-27 16:07:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:37 INFO Hourly dataset computed and listing created
2026-07-27 16:07:38 INFO Hourly dataset computed
2026-07-27 16:07:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:39 INFO Hourly dataset computed and listing created
2026-07-27 16:07:39 INFO Hourly dataset computed
2026-07-27 16:07:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:40 INFO Hourly dataset computed and listing created
2026-07-27 16:07:41 INFO Hourly dataset computed
2026-07-27 16:07:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:42 INFO Hourly dataset computed and listing created
2026-07-27 16:07:42 INFO Hourly dataset computed
2026-07-27 16:07:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:43 INFO Hourly dataset computed and listing created
2026-07-27 16:07:44 INFO Hourly dataset computed
2026-07-27 16:07:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:45 INFO Hourly dataset computed and listing created
2026-07-27 16:07:45 INFO Hourly dataset computed
2026-07-27 16:07:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:46 INFO Hourly dataset computed and listing created
2026-07-27 16:07:47 INFO Hourly dataset computed
2026-07-27 16:07:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:48 INFO Hourly dataset computed and listing created
2026-07-27 16:07:48 INFO Hourly dataset computed
2026-07-27 16:07:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:49 INFO Hourly dataset computed and listing created
2026-07-27 16:07:50 INFO Hourly dataset computed
2026-07-27 16:07:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:51 INFO Hourly dataset computed and listing created
2026-07-27 16:07:51 INFO Hourly dataset computed
2026-07-27 16:07:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:07:52 INFO Hourly dataset computed and listing created
2026-07-27 16:07:53 INFO Hourly dataset computed
2026-07-27 16:07:53 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-27 16:07:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:07:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 16:07:53 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-27 16:07:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 16:07:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:07:53 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 16:07:53 INFO Queuing job for member 1...
2026-07-27 16:07:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:07:53 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 16:07:55 INFO Found: ['5287672']
2026-07-27 16:08:00 INFO [TGCC-IRENE] Submitted job with ID:['5287672']
2026-07-27 16:08:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 16:08:00 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-27 16:08:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 16:08:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:00 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 16:08:00 INFO Queuing job for member 2...
2026-07-27 16:08:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:00 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 16:08:03 INFO Found: ['5287674']
2026-07-27 16:08:08 INFO [TGCC-IRENE] Submitted job with ID:['5287674']
2026-07-27 16:08:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 16:08:08 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-27 16:08:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 16:08:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:08 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 16:08:08 INFO Queuing job for member 3...
2026-07-27 16:08:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:08 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 16:08:08 INFO Found: ['5287675']
2026-07-27 16:08:13 INFO [TGCC-IRENE] Submitted job with ID:['5287675']
2026-07-27 16:08:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 16:08:13 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-27 16:08:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 16:08:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:13 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 16:08:13 INFO Queuing job for member 4...
2026-07-27 16:08:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:13 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 16:08:14 INFO Found: ['5287677']
2026-07-27 16:08:19 INFO [TGCC-IRENE] Submitted job with ID:['5287677']
2026-07-27 16:08:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 16:08:19 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-27 16:08:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 16:08:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:19 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 16:08:19 INFO Queuing job for member 5...
2026-07-27 16:08:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:19 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 16:08:20 INFO Found: ['5287679']
2026-07-27 16:08:25 INFO [TGCC-IRENE] Submitted job with ID:['5287679']
2026-07-27 16:08:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 16:08:25 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-27 16:08:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 16:08:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:25 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 16:08:25 INFO Queuing job for member 6...
2026-07-27 16:08:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:25 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 16:08:26 INFO Found: ['5287681']
2026-07-27 16:08:31 INFO [TGCC-IRENE] Submitted job with ID:['5287681']
2026-07-27 16:08:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 16:08:31 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-27 16:08:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 16:08:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:31 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 16:08:31 INFO Queuing job for member 7...
2026-07-27 16:08:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:31 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 16:08:31 INFO Found: ['5287684']
2026-07-27 16:08:36 INFO [TGCC-IRENE] Submitted job with ID:['5287684']
2026-07-27 16:08:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 16:08:36 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-27 16:08:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 16:08:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:36 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 16:08:36 INFO Queuing job for member 8...
2026-07-27 16:08:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:36 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 16:08:37 INFO Found: ['5287686']
2026-07-27 16:08:42 INFO [TGCC-IRENE] Submitted job with ID:['5287686']
2026-07-27 16:08:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 16:08:42 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-27 16:08:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 16:08:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:42 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 16:08:42 INFO Queuing job for member 9...
2026-07-27 16:08:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:42 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 16:08:43 INFO Found: ['5287687']
2026-07-27 16:08:48 INFO [TGCC-IRENE] Submitted job with ID:['5287687']
2026-07-27 16:08:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 16:08:48 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-27 16:08:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 16:08:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:48 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 16:08:48 INFO Queuing job for member 10...
2026-07-27 16:08:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:48 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 16:08:50 INFO Found: ['5287688']
2026-07-27 16:08:55 INFO [TGCC-IRENE] Submitted job with ID:['5287688']
2026-07-27 16:08:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:08:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 16:08:55 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-27 16:08:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 16:08:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:08:55 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 16:08:55 INFO Queuing job for member 11...
2026-07-27 16:08:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:08:55 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 16:08:57 INFO Found: ['5287690']
2026-07-27 16:09:02 INFO [TGCC-IRENE] Submitted job with ID:['5287690']
2026-07-27 16:09:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:09:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 16:09:02 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-27 16:09:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 16:09:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:09:02 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 16:09:02 INFO Queuing job for member 12...
2026-07-27 16:09:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:09:02 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 16:09:05 INFO Found: ['5287693']
2026-07-27 16:09:10 INFO [TGCC-IRENE] Submitted job with ID:['5287693']
2026-07-27 16:09:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:09:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 16:09:10 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-27 16:09:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 16:09:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:09:10 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 16:09:10 INFO Queuing job for member 13...
2026-07-27 16:09:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:09:10 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 16:09:12 INFO Found: ['5287695']
2026-07-27 16:09:17 INFO [TGCC-IRENE] Submitted job with ID:['5287695']
2026-07-27 16:09:17 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:09:17 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 16:09:17 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-27 16:09:17 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 16:09:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:09:17 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 16:09:17 INFO Queuing job for member 14...
2026-07-27 16:09:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:09:17 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 16:09:20 INFO Found: ['5287696']
2026-07-27 16:09:25 INFO [TGCC-IRENE] Submitted job with ID:['5287696']
2026-07-27 16:09:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:09:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 16:09:25 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-27 16:09:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 16:09:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:09:25 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 16:09:25 INFO Queuing job for member 15...
2026-07-27 16:09:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:09:25 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 16:09:27 INFO Found: ['5287697']
2026-07-27 16:09:32 INFO [TGCC-IRENE] Submitted job with ID:['5287697']
2026-07-27 16:09:32 INFO Checking job status ...
2026-07-27 16:09:32 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:09:32 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:09:32 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:09:32 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:09:32 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:09:33 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:09:33 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:09:48 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:09:48 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:09:48 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:09:48 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:09:48 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:09:50 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:09:50 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:10:05 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:10:05 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:10:05 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:10:20 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:10:20 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:10:21 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:10:21 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:10:21 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:10:21 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:10:21 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:10:21 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:10:36 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:10:36 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:10:36 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:10:51 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:10:51 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:10:51 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:11:07 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:11:07 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:11:07 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:11:07 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:11:08 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:11:08 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:11:23 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:11:23 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:11:23 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:11:38 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:11:38 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:11:38 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:11:38 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:11:40 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:11:40 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:11:55 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287679: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:11:55 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:11:56 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:11:56 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:11:56 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:11:56 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:11:56 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:12:11 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:12:11 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:12:11 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:12:11 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:12:11 INFO None 5287679: status FINISHED
2026-07-27 16:12:11 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287687: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:12:13 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:12:13 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:12:28 INFO None 5287672: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287674: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287675: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287679: status FINISHED
2026-07-27 16:12:28 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287684: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287687: status FINISHED
2026-07-27 16:12:28 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:12:28 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:12:28 INFO Jobs still running: ['5287672', '5287674', '5287675', '5287677', '5287681', '5287684', '5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:12:43 INFO None 5287672: status FINISHED
2026-07-27 16:12:43 INFO None 5287674: status FINISHED
2026-07-27 16:12:43 INFO None 5287675: status FINISHED
2026-07-27 16:12:43 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287679: status FINISHED
2026-07-27 16:12:43 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287684: status FINISHED
2026-07-27 16:12:43 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287687: status FINISHED
2026-07-27 16:12:43 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:12:43 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:12:43 INFO Jobs still running: ['5287677', '5287681', '5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:12:59 INFO None 5287672: status FINISHED
2026-07-27 16:12:59 INFO None 5287674: status FINISHED
2026-07-27 16:12:59 INFO None 5287675: status FINISHED
2026-07-27 16:12:59 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287679: status FINISHED
2026-07-27 16:12:59 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287684: status FINISHED
2026-07-27 16:12:59 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287687: status FINISHED
2026-07-27 16:12:59 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:12:59 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:12:59 INFO Jobs still running: ['5287677', '5287681', '5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:13:14 INFO None 5287672: status FINISHED
2026-07-27 16:13:14 INFO None 5287674: status FINISHED
2026-07-27 16:13:14 INFO None 5287675: status FINISHED
2026-07-27 16:13:14 INFO None 5287677: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287679: status FINISHED
2026-07-27 16:13:14 INFO None 5287681: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287684: status FINISHED
2026-07-27 16:13:14 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287687: status FINISHED
2026-07-27 16:13:14 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:13:14 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:13:14 INFO Jobs still running: ['5287677', '5287681', '5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:13:30 INFO None 5287672: status FINISHED
2026-07-27 16:13:30 INFO None 5287674: status FINISHED
2026-07-27 16:13:30 INFO None 5287675: status FINISHED
2026-07-27 16:13:30 INFO None 5287677: status FINISHED
2026-07-27 16:13:30 INFO None 5287679: status FINISHED
2026-07-27 16:13:30 INFO None 5287681: status FINISHED
2026-07-27 16:13:30 INFO None 5287684: status FINISHED
2026-07-27 16:13:30 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287687: status FINISHED
2026-07-27 16:13:30 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:13:30 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:13:30 INFO Jobs still running: ['5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:13:45 INFO None 5287672: status FINISHED
2026-07-27 16:13:45 INFO None 5287674: status FINISHED
2026-07-27 16:13:45 INFO None 5287675: status FINISHED
2026-07-27 16:13:45 INFO None 5287677: status FINISHED
2026-07-27 16:13:45 INFO None 5287679: status FINISHED
2026-07-27 16:13:45 INFO None 5287681: status FINISHED
2026-07-27 16:13:45 INFO None 5287684: status FINISHED
2026-07-27 16:13:45 INFO None 5287686: status RUNNING/PENDING
2026-07-27 16:13:45 INFO None 5287687: status FINISHED
2026-07-27 16:13:47 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:13:47 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:13:47 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:13:47 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:13:47 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:13:47 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:13:47 INFO Jobs still running: ['5287686', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:14:02 INFO None 5287672: status FINISHED
2026-07-27 16:14:02 INFO None 5287674: status FINISHED
2026-07-27 16:14:02 INFO None 5287675: status FINISHED
2026-07-27 16:14:02 INFO None 5287677: status FINISHED
2026-07-27 16:14:02 INFO None 5287679: status FINISHED
2026-07-27 16:14:02 INFO None 5287681: status FINISHED
2026-07-27 16:14:02 INFO None 5287684: status FINISHED
2026-07-27 16:14:02 INFO None 5287686: status FINISHED
2026-07-27 16:14:02 INFO None 5287687: status FINISHED
2026-07-27 16:14:03 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:14:03 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:14:03 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:14:03 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:14:03 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:14:03 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:14:03 INFO Jobs still running: ['5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:14:18 INFO None 5287672: status FINISHED
2026-07-27 16:14:18 INFO None 5287674: status FINISHED
2026-07-27 16:14:18 INFO None 5287675: status FINISHED
2026-07-27 16:14:18 INFO None 5287677: status FINISHED
2026-07-27 16:14:18 INFO None 5287679: status FINISHED
2026-07-27 16:14:18 INFO None 5287681: status FINISHED
2026-07-27 16:14:18 INFO None 5287684: status FINISHED
2026-07-27 16:14:18 INFO None 5287686: status FINISHED
2026-07-27 16:14:18 INFO None 5287687: status FINISHED
2026-07-27 16:14:18 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:14:20 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:14:20 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:14:20 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:14:20 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:14:20 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:14:20 INFO Jobs still running: ['5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:14:35 INFO None 5287672: status FINISHED
2026-07-27 16:14:35 INFO None 5287674: status FINISHED
2026-07-27 16:14:35 INFO None 5287675: status FINISHED
2026-07-27 16:14:35 INFO None 5287677: status FINISHED
2026-07-27 16:14:35 INFO None 5287679: status FINISHED
2026-07-27 16:14:35 INFO None 5287681: status FINISHED
2026-07-27 16:14:35 INFO None 5287684: status FINISHED
2026-07-27 16:14:35 INFO None 5287686: status FINISHED
2026-07-27 16:14:35 INFO None 5287687: status FINISHED
2026-07-27 16:14:35 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:14:35 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:14:35 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:14:35 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:14:35 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:14:35 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:14:35 INFO Jobs still running: ['5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:14:50 INFO None 5287672: status FINISHED
2026-07-27 16:14:50 INFO None 5287674: status FINISHED
2026-07-27 16:14:50 INFO None 5287675: status FINISHED
2026-07-27 16:14:50 INFO None 5287677: status FINISHED
2026-07-27 16:14:50 INFO None 5287679: status FINISHED
2026-07-27 16:14:50 INFO None 5287681: status FINISHED
2026-07-27 16:14:50 INFO None 5287684: status FINISHED
2026-07-27 16:14:50 INFO None 5287686: status FINISHED
2026-07-27 16:14:50 INFO None 5287687: status FINISHED
2026-07-27 16:14:50 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:14:50 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:14:50 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:14:50 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:14:50 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:14:51 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:14:51 INFO Jobs still running: ['5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:15:06 INFO None 5287672: status FINISHED
2026-07-27 16:15:06 INFO None 5287674: status FINISHED
2026-07-27 16:15:06 INFO None 5287675: status FINISHED
2026-07-27 16:15:06 INFO None 5287677: status FINISHED
2026-07-27 16:15:06 INFO None 5287679: status FINISHED
2026-07-27 16:15:06 INFO None 5287681: status FINISHED
2026-07-27 16:15:06 INFO None 5287684: status FINISHED
2026-07-27 16:15:06 INFO None 5287686: status FINISHED
2026-07-27 16:15:06 INFO None 5287687: status FINISHED
2026-07-27 16:15:06 INFO None 5287688: status RUNNING/PENDING
2026-07-27 16:15:06 INFO None 5287690: status RUNNING/PENDING
2026-07-27 16:15:06 INFO None 5287693: status RUNNING/PENDING
2026-07-27 16:15:06 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:15:06 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:15:06 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:15:06 INFO Jobs still running: ['5287688', '5287690', '5287693', '5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:15:21 INFO None 5287672: status FINISHED
2026-07-27 16:15:21 INFO None 5287674: status FINISHED
2026-07-27 16:15:21 INFO None 5287675: status FINISHED
2026-07-27 16:15:21 INFO None 5287677: status FINISHED
2026-07-27 16:15:21 INFO None 5287679: status FINISHED
2026-07-27 16:15:21 INFO None 5287681: status FINISHED
2026-07-27 16:15:21 INFO None 5287684: status FINISHED
2026-07-27 16:15:21 INFO None 5287686: status FINISHED
2026-07-27 16:15:21 INFO None 5287687: status FINISHED
2026-07-27 16:15:21 INFO None 5287688: status FINISHED
2026-07-27 16:15:21 INFO None 5287690: status FINISHED
2026-07-27 16:15:21 INFO None 5287693: status FINISHED
2026-07-27 16:15:21 INFO None 5287695: status RUNNING/PENDING
2026-07-27 16:15:21 INFO None 5287696: status RUNNING/PENDING
2026-07-27 16:15:21 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:15:21 INFO Jobs still running: ['5287695', '5287696', '5287697']. Waiting...
2026-07-27 16:15:38 INFO None 5287672: status FINISHED
2026-07-27 16:15:38 INFO None 5287674: status FINISHED
2026-07-27 16:15:38 INFO None 5287675: status FINISHED
2026-07-27 16:15:38 INFO None 5287677: status FINISHED
2026-07-27 16:15:38 INFO None 5287679: status FINISHED
2026-07-27 16:15:38 INFO None 5287681: status FINISHED
2026-07-27 16:15:38 INFO None 5287684: status FINISHED
2026-07-27 16:15:38 INFO None 5287686: status FINISHED
2026-07-27 16:15:38 INFO None 5287687: status FINISHED
2026-07-27 16:15:38 INFO None 5287688: status FINISHED
2026-07-27 16:15:38 INFO None 5287690: status FINISHED
2026-07-27 16:15:38 INFO None 5287693: status FINISHED
2026-07-27 16:15:38 INFO None 5287695: status FINISHED
2026-07-27 16:15:38 INFO None 5287696: status FINISHED
2026-07-27 16:15:38 INFO None 5287697: status RUNNING/PENDING
2026-07-27 16:15:38 INFO Jobs still running: ['5287697']. Waiting...
2026-07-27 16:15:53 INFO None 5287672: status FINISHED
2026-07-27 16:15:53 INFO None 5287674: status FINISHED
2026-07-27 16:15:53 INFO None 5287675: status FINISHED
2026-07-27 16:15:53 INFO None 5287677: status FINISHED
2026-07-27 16:15:53 INFO None 5287679: status FINISHED
2026-07-27 16:15:53 INFO None 5287681: status FINISHED
2026-07-27 16:15:53 INFO None 5287684: status FINISHED
2026-07-27 16:15:53 INFO None 5287686: status FINISHED
2026-07-27 16:15:53 INFO None 5287687: status FINISHED
2026-07-27 16:15:53 INFO None 5287688: status FINISHED
2026-07-27 16:15:53 INFO None 5287690: status FINISHED
2026-07-27 16:15:53 INFO None 5287693: status FINISHED
2026-07-27 16:15:53 INFO None 5287695: status FINISHED
2026-07-27 16:15:53 INFO None 5287696: status FINISHED
2026-07-27 16:15:53 INFO None 5287697: status FINISHED
2026-07-27 16:15:53 INFO Jobs ['5287672', '5287674', '5287675', '5287677', '5287679', '5287681', '5287684', '5287686', '5287687', '5287688', '5287690', '5287693', '5287695', '5287696', '5287697'] have finished
2026-07-27 16:15:53 INFO Checking restart files were created ...
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-27 16:15:53 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-27 16:15:53 INFO  Run_model() completed successfully.
2026-07-27 16:15:53 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:15:53 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-27 16:15:53 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 16:15:53 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-27 16:15:53 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:15:53 INFO ---------->>> Running process_satellite_data()
2026-07-27 16:15:53 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-27 16:15:53 INFO ---------->>> Running run_obs_converter()
2026-07-27 16:15:53 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-27 16:15:53 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-27 16:15:53 INFO ---------->>> Running DART
2026-07-27 16:15:53 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-27 16:15:53 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 16:15:53 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 16:15:53 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 16:15:53 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 16:15:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 16:15:53 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 16:16:08 INFO Found: []
2026-07-27 16:16:08 INFO No job id returned by command ./run_filter.bsh
2026-07-27 16:16:08 INFO No monitoring will be performed
2026-07-27 16:16:08 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-27 16:16:08 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:08 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:08 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 16:16:09 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 16:16:09 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 16:16:09 INFO run_dart() is DONE.
2026-07-27 16:16:09 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 16:16:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-27 16:16:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-27 16:16:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-27 16:16:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-27 16:16:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-27 16:16:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-27 16:16:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:31 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-27 16:16:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-27 16:16:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-27 16:16:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-27 16:16:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-27 16:16:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-27 16:16:54 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-27 16:16:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:16:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-27 16:17:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:17:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-27 16:17:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 16:17:06 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 16:17:06 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:17:06 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:17:06 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-27 16:17:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:17:07 INFO Hourly dataset computed and listing created
2026-07-27 16:17:28 INFO Hourly dataset computed
2026-07-27 16:17:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:17:30 INFO Hourly dataset computed and listing created
2026-07-27 16:17:49 INFO Hourly dataset computed
2026-07-27 16:17:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:17:50 INFO Hourly dataset computed and listing created
2026-07-27 16:18:02 INFO Hourly dataset computed
2026-07-27 16:18:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:03 INFO Hourly dataset computed and listing created
2026-07-27 16:18:06 INFO Hourly dataset computed
2026-07-27 16:18:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:07 INFO Hourly dataset computed and listing created
2026-07-27 16:18:09 INFO Hourly dataset computed
2026-07-27 16:18:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:10 INFO Hourly dataset computed and listing created
2026-07-27 16:18:13 INFO Hourly dataset computed
2026-07-27 16:18:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:14 INFO Hourly dataset computed and listing created
2026-07-27 16:18:17 INFO Hourly dataset computed
2026-07-27 16:18:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:18 INFO Hourly dataset computed and listing created
2026-07-27 16:18:21 INFO Hourly dataset computed
2026-07-27 16:18:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:22 INFO Hourly dataset computed and listing created
2026-07-27 16:18:24 INFO Hourly dataset computed
2026-07-27 16:18:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:25 INFO Hourly dataset computed and listing created
2026-07-27 16:18:28 INFO Hourly dataset computed
2026-07-27 16:18:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:29 INFO Hourly dataset computed and listing created
2026-07-27 16:18:32 INFO Hourly dataset computed
2026-07-27 16:18:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:33 INFO Hourly dataset computed and listing created
2026-07-27 16:18:36 INFO Hourly dataset computed
2026-07-27 16:18:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:37 INFO Hourly dataset computed and listing created
2026-07-27 16:18:39 INFO Hourly dataset computed
2026-07-27 16:18:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:41 INFO Hourly dataset computed and listing created
2026-07-27 16:18:43 INFO Hourly dataset computed
2026-07-27 16:18:43 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:18:45 INFO Hourly dataset computed and listing created
2026-07-27 16:18:47 INFO Hourly dataset computed
2026-07-27 16:18:47 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-27 16:18:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:18:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 16:18:47 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-27 16:18:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 16:18:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:18:47 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 16:18:47 INFO Queuing job for member 1...
2026-07-27 16:18:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:18:47 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 16:18:50 INFO Found: ['5287777']
2026-07-27 16:18:55 INFO [TGCC-IRENE] Submitted job with ID:['5287777']
2026-07-27 16:18:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:18:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 16:18:55 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-27 16:18:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 16:18:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:18:55 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 16:18:55 INFO Queuing job for member 2...
2026-07-27 16:18:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:18:55 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 16:18:57 INFO Found: ['5287778']
2026-07-27 16:19:03 INFO [TGCC-IRENE] Submitted job with ID:['5287778']
2026-07-27 16:19:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 16:19:03 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-27 16:19:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 16:19:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:03 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 16:19:03 INFO Queuing job for member 3...
2026-07-27 16:19:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:03 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 16:19:05 INFO Found: ['5287782']
2026-07-27 16:19:10 INFO [TGCC-IRENE] Submitted job with ID:['5287782']
2026-07-27 16:19:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 16:19:10 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-27 16:19:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 16:19:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:10 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 16:19:10 INFO Queuing job for member 4...
2026-07-27 16:19:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:10 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 16:19:13 INFO Found: ['5287784']
2026-07-27 16:19:18 INFO [TGCC-IRENE] Submitted job with ID:['5287784']
2026-07-27 16:19:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 16:19:18 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-27 16:19:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 16:19:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:18 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 16:19:18 INFO Queuing job for member 5...
2026-07-27 16:19:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:18 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 16:19:18 INFO Found: ['5287785']
2026-07-27 16:19:23 INFO [TGCC-IRENE] Submitted job with ID:['5287785']
2026-07-27 16:19:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 16:19:23 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-27 16:19:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 16:19:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:23 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 16:19:23 INFO Queuing job for member 6...
2026-07-27 16:19:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:23 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 16:19:24 INFO Found: ['5287786']
2026-07-27 16:19:29 INFO [TGCC-IRENE] Submitted job with ID:['5287786']
2026-07-27 16:19:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 16:19:29 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-27 16:19:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 16:19:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:29 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 16:19:29 INFO Queuing job for member 7...
2026-07-27 16:19:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:29 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 16:19:30 INFO Found: ['5287789']
2026-07-27 16:19:35 INFO [TGCC-IRENE] Submitted job with ID:['5287789']
2026-07-27 16:19:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 16:19:35 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-27 16:19:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 16:19:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:35 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 16:19:35 INFO Queuing job for member 8...
2026-07-27 16:19:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:35 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 16:19:36 INFO Found: ['5287790']
2026-07-27 16:19:41 INFO [TGCC-IRENE] Submitted job with ID:['5287790']
2026-07-27 16:19:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 16:19:41 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-27 16:19:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 16:19:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:41 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 16:19:41 INFO Queuing job for member 9...
2026-07-27 16:19:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:41 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 16:19:41 INFO Found: ['5287791']
2026-07-27 16:19:46 INFO [TGCC-IRENE] Submitted job with ID:['5287791']
2026-07-27 16:19:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 16:19:46 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-27 16:19:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 16:19:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:46 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 16:19:46 INFO Queuing job for member 10...
2026-07-27 16:19:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:46 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 16:19:47 INFO Found: ['5287793']
2026-07-27 16:19:52 INFO [TGCC-IRENE] Submitted job with ID:['5287793']
2026-07-27 16:19:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 16:19:52 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-27 16:19:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 16:19:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:52 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 16:19:52 INFO Queuing job for member 11...
2026-07-27 16:19:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:52 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 16:19:53 INFO Found: ['5287795']
2026-07-27 16:19:58 INFO [TGCC-IRENE] Submitted job with ID:['5287795']
2026-07-27 16:19:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:19:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 16:19:58 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-27 16:19:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 16:19:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:19:58 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 16:19:58 INFO Queuing job for member 12...
2026-07-27 16:19:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:19:58 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 16:19:58 INFO Found: ['5287796']
2026-07-27 16:20:03 INFO [TGCC-IRENE] Submitted job with ID:['5287796']
2026-07-27 16:20:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:20:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 16:20:03 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-27 16:20:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 16:20:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:20:03 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 16:20:04 INFO Queuing job for member 13...
2026-07-27 16:20:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:20:04 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 16:20:05 INFO Found: ['5287798']
2026-07-27 16:20:10 INFO [TGCC-IRENE] Submitted job with ID:['5287798']
2026-07-27 16:20:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:20:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 16:20:10 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-27 16:20:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 16:20:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:20:10 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 16:20:10 INFO Queuing job for member 14...
2026-07-27 16:20:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:20:10 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 16:20:13 INFO Found: ['5287799']
2026-07-27 16:20:18 INFO [TGCC-IRENE] Submitted job with ID:['5287799']
2026-07-27 16:20:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:20:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 16:20:18 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-27 16:20:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 16:20:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:20:18 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 16:20:18 INFO Queuing job for member 15...
2026-07-27 16:20:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:20:18 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 16:20:20 INFO Found: ['5287800']
2026-07-27 16:20:25 INFO [TGCC-IRENE] Submitted job with ID:['5287800']
2026-07-27 16:20:25 INFO Checking job status ...
2026-07-27 16:20:25 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:20:25 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:20:25 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:20:40 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:20:41 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:20:43 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:20:43 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:20:43 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:20:43 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:20:43 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:20:43 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:20:58 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:20:58 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:20:58 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:21:13 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:21:13 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:21:15 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:21:15 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:21:15 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:21:15 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:21:15 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:21:16 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:21:16 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:21:16 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:21:31 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:21:31 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:21:31 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:21:46 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:21:46 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:21:46 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:22:01 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:22:01 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:22:01 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:22:16 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:22:16 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:22:16 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:22:16 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:22:16 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:22:17 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:22:17 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:22:32 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:22:32 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:22:32 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:22:47 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:22:47 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:22:48 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:22:48 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:23:03 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:23:03 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:23:03 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:23:05 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:23:05 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:23:20 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:23:20 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:23:22 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:23:22 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:23:37 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:23:37 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:23:38 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:23:38 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:23:53 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:23:53 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:23:53 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:24:08 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:24:08 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:24:08 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:24:23 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:24:23 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:24:24 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:24:24 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:24:40 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:24:40 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:24:40 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:24:55 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:24:55 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:24:57 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:24:58 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:24:58 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:24:58 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:25:13 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:25:13 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:25:13 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:25:28 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:25:28 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:25:28 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:25:28 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:25:28 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:25:30 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:25:30 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:25:45 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:25:45 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:25:45 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:25:46 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:25:48 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:25:48 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:26:03 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:26:03 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:26:03 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:26:18 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:26:18 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:26:18 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:26:33 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:26:33 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:26:34 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:26:34 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:26:34 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:26:34 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:26:34 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:26:34 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:26:50 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:26:50 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:26:50 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:27:05 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:27:05 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:27:07 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:27:07 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:27:07 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:27:07 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:27:22 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:27:22 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:27:23 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:27:23 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:27:38 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:27:38 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:27:40 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:27:40 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:27:55 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:27:55 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:27:55 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:28:10 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:28:10 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:28:10 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:28:10 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:28:10 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:28:13 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:28:13 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:28:28 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:28:28 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:28:28 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:28:43 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:28:43 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:28:43 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:28:58 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:28:58 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:28:59 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:28:59 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:29:15 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:29:15 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:29:15 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:29:30 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:29:30 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:29:30 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:29:30 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:29:31 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:29:33 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:29:33 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:29:33 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:29:33 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:29:33 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:29:33 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:29:48 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:29:48 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:29:48 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:30:03 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:30:03 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:30:05 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:30:05 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:30:05 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:30:05 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:30:05 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:30:05 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:30:20 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:30:20 INFO None 5287778: status RUNNING/PENDING
2026-07-27 16:30:20 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:30:20 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:30:20 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:30:21 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:30:21 INFO Jobs still running: ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:30:36 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287778: status FINISHED
2026-07-27 16:30:36 INFO None 5287782: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287785: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:30:36 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:30:36 INFO Jobs still running: ['5287777', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:30:51 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287778: status FINISHED
2026-07-27 16:30:51 INFO None 5287782: status FINISHED
2026-07-27 16:30:51 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287785: status FINISHED
2026-07-27 16:30:51 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287789: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:30:51 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:30:51 INFO Jobs still running: ['5287777', '5287784', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:31:06 INFO None 5287777: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287778: status FINISHED
2026-07-27 16:31:06 INFO None 5287782: status FINISHED
2026-07-27 16:31:06 INFO None 5287784: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287785: status FINISHED
2026-07-27 16:31:06 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287789: status FINISHED
2026-07-27 16:31:06 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:31:06 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:31:07 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:31:07 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:31:07 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:31:07 INFO Jobs still running: ['5287777', '5287784', '5287786', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:31:22 INFO None 5287777: status FINISHED
2026-07-27 16:31:22 INFO None 5287778: status FINISHED
2026-07-27 16:31:22 INFO None 5287782: status FINISHED
2026-07-27 16:31:22 INFO None 5287784: status FINISHED
2026-07-27 16:31:22 INFO None 5287785: status FINISHED
2026-07-27 16:31:22 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287789: status FINISHED
2026-07-27 16:31:22 INFO None 5287790: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287791: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287793: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287795: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:31:22 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:31:22 INFO Jobs still running: ['5287786', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:31:37 INFO None 5287777: status FINISHED
2026-07-27 16:31:37 INFO None 5287778: status FINISHED
2026-07-27 16:31:37 INFO None 5287782: status FINISHED
2026-07-27 16:31:37 INFO None 5287784: status FINISHED
2026-07-27 16:31:37 INFO None 5287785: status FINISHED
2026-07-27 16:31:37 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:31:37 INFO None 5287789: status FINISHED
2026-07-27 16:31:37 INFO None 5287790: status FINISHED
2026-07-27 16:31:37 INFO None 5287791: status FINISHED
2026-07-27 16:31:37 INFO None 5287793: status FINISHED
2026-07-27 16:31:38 INFO None 5287795: status FINISHED
2026-07-27 16:31:38 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:31:38 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:31:38 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:31:38 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:31:38 INFO Jobs still running: ['5287786', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:31:55 INFO None 5287777: status FINISHED
2026-07-27 16:31:55 INFO None 5287778: status FINISHED
2026-07-27 16:31:55 INFO None 5287782: status FINISHED
2026-07-27 16:31:55 INFO None 5287784: status FINISHED
2026-07-27 16:31:55 INFO None 5287785: status FINISHED
2026-07-27 16:31:55 INFO None 5287786: status RUNNING/PENDING
2026-07-27 16:31:55 INFO None 5287789: status FINISHED
2026-07-27 16:31:55 INFO None 5287790: status FINISHED
2026-07-27 16:31:55 INFO None 5287791: status FINISHED
2026-07-27 16:31:55 INFO None 5287793: status FINISHED
2026-07-27 16:31:55 INFO None 5287795: status FINISHED
2026-07-27 16:31:55 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:31:55 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:31:55 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:31:55 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:31:55 INFO Jobs still running: ['5287786', '5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:32:10 INFO None 5287777: status FINISHED
2026-07-27 16:32:10 INFO None 5287778: status FINISHED
2026-07-27 16:32:10 INFO None 5287782: status FINISHED
2026-07-27 16:32:10 INFO None 5287784: status FINISHED
2026-07-27 16:32:10 INFO None 5287785: status FINISHED
2026-07-27 16:32:10 INFO None 5287786: status FINISHED
2026-07-27 16:32:10 INFO None 5287789: status FINISHED
2026-07-27 16:32:10 INFO None 5287790: status FINISHED
2026-07-27 16:32:10 INFO None 5287791: status FINISHED
2026-07-27 16:32:10 INFO None 5287793: status FINISHED
2026-07-27 16:32:10 INFO None 5287795: status FINISHED
2026-07-27 16:32:10 INFO None 5287796: status RUNNING/PENDING
2026-07-27 16:32:12 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:32:12 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:32:12 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:32:12 INFO Jobs still running: ['5287796', '5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:32:27 INFO None 5287777: status FINISHED
2026-07-27 16:32:27 INFO None 5287778: status FINISHED
2026-07-27 16:32:27 INFO None 5287782: status FINISHED
2026-07-27 16:32:27 INFO None 5287784: status FINISHED
2026-07-27 16:32:27 INFO None 5287785: status FINISHED
2026-07-27 16:32:27 INFO None 5287786: status FINISHED
2026-07-27 16:32:27 INFO None 5287789: status FINISHED
2026-07-27 16:32:27 INFO None 5287790: status FINISHED
2026-07-27 16:32:28 INFO None 5287791: status FINISHED
2026-07-27 16:32:28 INFO None 5287793: status FINISHED
2026-07-27 16:32:28 INFO None 5287795: status FINISHED
2026-07-27 16:32:28 INFO None 5287796: status FINISHED
2026-07-27 16:32:28 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:32:28 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:32:28 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:32:28 INFO Jobs still running: ['5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:32:43 INFO None 5287777: status FINISHED
2026-07-27 16:32:43 INFO None 5287778: status FINISHED
2026-07-27 16:32:43 INFO None 5287782: status FINISHED
2026-07-27 16:32:43 INFO None 5287784: status FINISHED
2026-07-27 16:32:43 INFO None 5287785: status FINISHED
2026-07-27 16:32:43 INFO None 5287786: status FINISHED
2026-07-27 16:32:43 INFO None 5287789: status FINISHED
2026-07-27 16:32:43 INFO None 5287790: status FINISHED
2026-07-27 16:32:43 INFO None 5287791: status FINISHED
2026-07-27 16:32:43 INFO None 5287793: status FINISHED
2026-07-27 16:32:43 INFO None 5287795: status FINISHED
2026-07-27 16:32:43 INFO None 5287796: status FINISHED
2026-07-27 16:32:43 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:32:45 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:32:45 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:32:45 INFO Jobs still running: ['5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:33:00 INFO None 5287777: status FINISHED
2026-07-27 16:33:00 INFO None 5287778: status FINISHED
2026-07-27 16:33:00 INFO None 5287782: status FINISHED
2026-07-27 16:33:00 INFO None 5287784: status FINISHED
2026-07-27 16:33:00 INFO None 5287785: status FINISHED
2026-07-27 16:33:00 INFO None 5287786: status FINISHED
2026-07-27 16:33:00 INFO None 5287789: status FINISHED
2026-07-27 16:33:00 INFO None 5287790: status FINISHED
2026-07-27 16:33:00 INFO None 5287791: status FINISHED
2026-07-27 16:33:00 INFO None 5287793: status FINISHED
2026-07-27 16:33:00 INFO None 5287795: status FINISHED
2026-07-27 16:33:00 INFO None 5287796: status FINISHED
2026-07-27 16:33:00 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:33:00 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:33:00 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:33:00 INFO Jobs still running: ['5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:33:15 INFO None 5287777: status FINISHED
2026-07-27 16:33:15 INFO None 5287778: status FINISHED
2026-07-27 16:33:15 INFO None 5287782: status FINISHED
2026-07-27 16:33:15 INFO None 5287784: status FINISHED
2026-07-27 16:33:15 INFO None 5287785: status FINISHED
2026-07-27 16:33:15 INFO None 5287786: status FINISHED
2026-07-27 16:33:15 INFO None 5287789: status FINISHED
2026-07-27 16:33:15 INFO None 5287790: status FINISHED
2026-07-27 16:33:15 INFO None 5287791: status FINISHED
2026-07-27 16:33:15 INFO None 5287793: status FINISHED
2026-07-27 16:33:15 INFO None 5287795: status FINISHED
2026-07-27 16:33:15 INFO None 5287796: status FINISHED
2026-07-27 16:33:16 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:33:16 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:33:16 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:33:16 INFO Jobs still running: ['5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:33:31 INFO None 5287777: status FINISHED
2026-07-27 16:33:31 INFO None 5287778: status FINISHED
2026-07-27 16:33:31 INFO None 5287782: status FINISHED
2026-07-27 16:33:31 INFO None 5287784: status FINISHED
2026-07-27 16:33:31 INFO None 5287785: status FINISHED
2026-07-27 16:33:31 INFO None 5287786: status FINISHED
2026-07-27 16:33:31 INFO None 5287789: status FINISHED
2026-07-27 16:33:31 INFO None 5287790: status FINISHED
2026-07-27 16:33:31 INFO None 5287791: status FINISHED
2026-07-27 16:33:31 INFO None 5287793: status FINISHED
2026-07-27 16:33:31 INFO None 5287795: status FINISHED
2026-07-27 16:33:31 INFO None 5287796: status FINISHED
2026-07-27 16:33:31 INFO None 5287798: status RUNNING/PENDING
2026-07-27 16:33:31 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:33:31 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:33:31 INFO Jobs still running: ['5287798', '5287799', '5287800']. Waiting...
2026-07-27 16:33:47 INFO None 5287777: status FINISHED
2026-07-27 16:33:48 INFO None 5287778: status FINISHED
2026-07-27 16:33:48 INFO None 5287782: status FINISHED
2026-07-27 16:33:48 INFO None 5287784: status FINISHED
2026-07-27 16:33:48 INFO None 5287785: status FINISHED
2026-07-27 16:33:48 INFO None 5287786: status FINISHED
2026-07-27 16:33:48 INFO None 5287789: status FINISHED
2026-07-27 16:33:48 INFO None 5287790: status FINISHED
2026-07-27 16:33:48 INFO None 5287791: status FINISHED
2026-07-27 16:33:48 INFO None 5287793: status FINISHED
2026-07-27 16:33:48 INFO None 5287795: status FINISHED
2026-07-27 16:33:48 INFO None 5287796: status FINISHED
2026-07-27 16:33:48 INFO None 5287798: status FINISHED
2026-07-27 16:33:48 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:33:48 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:33:48 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:34:03 INFO None 5287777: status FINISHED
2026-07-27 16:34:03 INFO None 5287778: status FINISHED
2026-07-27 16:34:03 INFO None 5287782: status FINISHED
2026-07-27 16:34:03 INFO None 5287784: status FINISHED
2026-07-27 16:34:03 INFO None 5287785: status FINISHED
2026-07-27 16:34:03 INFO None 5287786: status FINISHED
2026-07-27 16:34:03 INFO None 5287789: status FINISHED
2026-07-27 16:34:03 INFO None 5287790: status FINISHED
2026-07-27 16:34:03 INFO None 5287791: status FINISHED
2026-07-27 16:34:03 INFO None 5287793: status FINISHED
2026-07-27 16:34:03 INFO None 5287795: status FINISHED
2026-07-27 16:34:03 INFO None 5287796: status FINISHED
2026-07-27 16:34:03 INFO None 5287798: status FINISHED
2026-07-27 16:34:03 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:34:03 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:34:03 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:34:18 INFO None 5287777: status FINISHED
2026-07-27 16:34:18 INFO None 5287778: status FINISHED
2026-07-27 16:34:18 INFO None 5287782: status FINISHED
2026-07-27 16:34:20 INFO None 5287784: status FINISHED
2026-07-27 16:34:20 INFO None 5287785: status FINISHED
2026-07-27 16:34:20 INFO None 5287786: status FINISHED
2026-07-27 16:34:20 INFO None 5287789: status FINISHED
2026-07-27 16:34:20 INFO None 5287790: status FINISHED
2026-07-27 16:34:20 INFO None 5287791: status FINISHED
2026-07-27 16:34:20 INFO None 5287793: status FINISHED
2026-07-27 16:34:20 INFO None 5287795: status FINISHED
2026-07-27 16:34:20 INFO None 5287796: status FINISHED
2026-07-27 16:34:20 INFO None 5287798: status FINISHED
2026-07-27 16:34:20 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:34:20 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:34:20 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:34:35 INFO None 5287777: status FINISHED
2026-07-27 16:34:35 INFO None 5287778: status FINISHED
2026-07-27 16:34:35 INFO None 5287782: status FINISHED
2026-07-27 16:34:35 INFO None 5287784: status FINISHED
2026-07-27 16:34:35 INFO None 5287785: status FINISHED
2026-07-27 16:34:35 INFO None 5287786: status FINISHED
2026-07-27 16:34:35 INFO None 5287789: status FINISHED
2026-07-27 16:34:35 INFO None 5287790: status FINISHED
2026-07-27 16:34:36 INFO None 5287791: status FINISHED
2026-07-27 16:34:36 INFO None 5287793: status FINISHED
2026-07-27 16:34:36 INFO None 5287795: status FINISHED
2026-07-27 16:34:36 INFO None 5287796: status FINISHED
2026-07-27 16:34:36 INFO None 5287798: status FINISHED
2026-07-27 16:34:36 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:34:36 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:34:36 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:34:51 INFO None 5287777: status FINISHED
2026-07-27 16:34:51 INFO None 5287778: status FINISHED
2026-07-27 16:34:51 INFO None 5287782: status FINISHED
2026-07-27 16:34:51 INFO None 5287784: status FINISHED
2026-07-27 16:34:51 INFO None 5287785: status FINISHED
2026-07-27 16:34:51 INFO None 5287786: status FINISHED
2026-07-27 16:34:51 INFO None 5287789: status FINISHED
2026-07-27 16:34:51 INFO None 5287790: status FINISHED
2026-07-27 16:34:53 INFO None 5287791: status FINISHED
2026-07-27 16:34:53 INFO None 5287793: status FINISHED
2026-07-27 16:34:53 INFO None 5287795: status FINISHED
2026-07-27 16:34:53 INFO None 5287796: status FINISHED
2026-07-27 16:34:53 INFO None 5287798: status FINISHED
2026-07-27 16:34:53 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:34:53 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:34:53 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:35:08 INFO None 5287777: status FINISHED
2026-07-27 16:35:08 INFO None 5287778: status FINISHED
2026-07-27 16:35:08 INFO None 5287782: status FINISHED
2026-07-27 16:35:08 INFO None 5287784: status FINISHED
2026-07-27 16:35:08 INFO None 5287785: status FINISHED
2026-07-27 16:35:08 INFO None 5287786: status FINISHED
2026-07-27 16:35:08 INFO None 5287789: status FINISHED
2026-07-27 16:35:08 INFO None 5287790: status FINISHED
2026-07-27 16:35:08 INFO None 5287791: status FINISHED
2026-07-27 16:35:08 INFO None 5287793: status FINISHED
2026-07-27 16:35:08 INFO None 5287795: status FINISHED
2026-07-27 16:35:08 INFO None 5287796: status FINISHED
2026-07-27 16:35:08 INFO None 5287798: status FINISHED
2026-07-27 16:35:08 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:35:08 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:35:08 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:35:23 INFO None 5287777: status FINISHED
2026-07-27 16:35:23 INFO None 5287778: status FINISHED
2026-07-27 16:35:23 INFO None 5287782: status FINISHED
2026-07-27 16:35:23 INFO None 5287784: status FINISHED
2026-07-27 16:35:23 INFO None 5287785: status FINISHED
2026-07-27 16:35:23 INFO None 5287786: status FINISHED
2026-07-27 16:35:23 INFO None 5287789: status FINISHED
2026-07-27 16:35:23 INFO None 5287790: status FINISHED
2026-07-27 16:35:23 INFO None 5287791: status FINISHED
2026-07-27 16:35:23 INFO None 5287793: status FINISHED
2026-07-27 16:35:23 INFO None 5287795: status FINISHED
2026-07-27 16:35:23 INFO None 5287796: status FINISHED
2026-07-27 16:35:23 INFO None 5287798: status FINISHED
2026-07-27 16:35:23 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:35:23 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:35:23 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:35:39 INFO None 5287777: status FINISHED
2026-07-27 16:35:39 INFO None 5287778: status FINISHED
2026-07-27 16:35:39 INFO None 5287782: status FINISHED
2026-07-27 16:35:39 INFO None 5287784: status FINISHED
2026-07-27 16:35:39 INFO None 5287785: status FINISHED
2026-07-27 16:35:39 INFO None 5287786: status FINISHED
2026-07-27 16:35:39 INFO None 5287789: status FINISHED
2026-07-27 16:35:39 INFO None 5287790: status FINISHED
2026-07-27 16:35:39 INFO None 5287791: status FINISHED
2026-07-27 16:35:39 INFO None 5287793: status FINISHED
2026-07-27 16:35:39 INFO None 5287795: status FINISHED
2026-07-27 16:35:39 INFO None 5287796: status FINISHED
2026-07-27 16:35:39 INFO None 5287798: status FINISHED
2026-07-27 16:35:39 INFO None 5287799: status RUNNING/PENDING
2026-07-27 16:35:39 INFO None 5287800: status RUNNING/PENDING
2026-07-27 16:35:39 INFO Jobs still running: ['5287799', '5287800']. Waiting...
2026-07-27 16:35:55 INFO None 5287777: status FINISHED
2026-07-27 16:35:55 INFO None 5287778: status FINISHED
2026-07-27 16:35:55 INFO None 5287782: status FINISHED
2026-07-27 16:35:55 INFO None 5287784: status FINISHED
2026-07-27 16:35:55 INFO None 5287785: status FINISHED
2026-07-27 16:35:55 INFO None 5287786: status FINISHED
2026-07-27 16:35:55 INFO None 5287789: status FINISHED
2026-07-27 16:35:55 INFO None 5287790: status FINISHED
2026-07-27 16:35:55 INFO None 5287791: status FINISHED
2026-07-27 16:35:55 INFO None 5287793: status FINISHED
2026-07-27 16:35:55 INFO None 5287795: status FINISHED
2026-07-27 16:35:55 INFO None 5287796: status FINISHED
2026-07-27 16:35:55 INFO None 5287798: status FINISHED
2026-07-27 16:35:55 INFO None 5287799: status FINISHED
2026-07-27 16:35:55 INFO None 5287800: status FINISHED
2026-07-27 16:35:55 INFO Jobs ['5287777', '5287778', '5287782', '5287784', '5287785', '5287786', '5287789', '5287790', '5287791', '5287793', '5287795', '5287796', '5287798', '5287799', '5287800'] have finished
2026-07-27 16:35:55 INFO Checking restart files were created ...
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-27 16:35:55 INFO  Run_model() completed successfully.
2026-07-27 16:35:55 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:35:55 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-27 16:35:55 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 16:35:55 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-27 16:35:55 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:35:55 INFO ---------->>> Running process_satellite_data()
2026-07-27 16:35:55 INFO [DART] No satellite data found, skipping assimilation
2026-07-27 16:35:55 INFO after_assimilation() skipped
2026-07-27 16:35:55 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 16:35:55 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:35:55 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:35:55 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-27 16:35:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:35:56 INFO Hourly dataset computed and listing created
2026-07-27 16:36:01 INFO Hourly dataset computed
2026-07-27 16:36:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:02 INFO Hourly dataset computed and listing created
2026-07-27 16:36:03 INFO Hourly dataset computed
2026-07-27 16:36:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:04 INFO Hourly dataset computed and listing created
2026-07-27 16:36:04 INFO Hourly dataset computed
2026-07-27 16:36:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:05 INFO Hourly dataset computed and listing created
2026-07-27 16:36:06 INFO Hourly dataset computed
2026-07-27 16:36:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:07 INFO Hourly dataset computed and listing created
2026-07-27 16:36:07 INFO Hourly dataset computed
2026-07-27 16:36:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:08 INFO Hourly dataset computed and listing created
2026-07-27 16:36:09 INFO Hourly dataset computed
2026-07-27 16:36:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:10 INFO Hourly dataset computed and listing created
2026-07-27 16:36:10 INFO Hourly dataset computed
2026-07-27 16:36:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:11 INFO Hourly dataset computed and listing created
2026-07-27 16:36:12 INFO Hourly dataset computed
2026-07-27 16:36:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:13 INFO Hourly dataset computed and listing created
2026-07-27 16:36:13 INFO Hourly dataset computed
2026-07-27 16:36:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:14 INFO Hourly dataset computed and listing created
2026-07-27 16:36:15 INFO Hourly dataset computed
2026-07-27 16:36:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:16 INFO Hourly dataset computed and listing created
2026-07-27 16:36:16 INFO Hourly dataset computed
2026-07-27 16:36:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:17 INFO Hourly dataset computed and listing created
2026-07-27 16:36:18 INFO Hourly dataset computed
2026-07-27 16:36:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:19 INFO Hourly dataset computed and listing created
2026-07-27 16:36:19 INFO Hourly dataset computed
2026-07-27 16:36:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:20 INFO Hourly dataset computed and listing created
2026-07-27 16:36:21 INFO Hourly dataset computed
2026-07-27 16:36:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:36:22 INFO Hourly dataset computed and listing created
2026-07-27 16:36:23 INFO Hourly dataset computed
2026-07-27 16:36:23 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-27 16:36:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:36:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 16:36:23 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-27 16:36:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 16:36:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:36:23 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 16:36:23 INFO Queuing job for member 1...
2026-07-27 16:36:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:36:23 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 16:36:25 INFO Found: ['5288006']
2026-07-27 16:36:30 INFO [TGCC-IRENE] Submitted job with ID:['5288006']
2026-07-27 16:36:30 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:36:30 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 16:36:30 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-27 16:36:30 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 16:36:30 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:36:30 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 16:36:30 INFO Queuing job for member 2...
2026-07-27 16:36:30 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:36:30 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 16:36:32 INFO Found: ['5288007']
2026-07-27 16:36:37 INFO [TGCC-IRENE] Submitted job with ID:['5288007']
2026-07-27 16:36:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:36:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 16:36:37 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-27 16:36:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 16:36:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:36:37 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 16:36:37 INFO Queuing job for member 3...
2026-07-27 16:36:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:36:37 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 16:36:40 INFO Found: ['5288008']
2026-07-27 16:36:45 INFO [TGCC-IRENE] Submitted job with ID:['5288008']
2026-07-27 16:36:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:36:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 16:36:45 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-27 16:36:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 16:36:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:36:45 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 16:36:45 INFO Queuing job for member 4...
2026-07-27 16:36:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:36:45 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 16:36:47 INFO Found: ['5288009']
2026-07-27 16:36:52 INFO [TGCC-IRENE] Submitted job with ID:['5288009']
2026-07-27 16:36:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:36:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 16:36:52 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-27 16:36:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 16:36:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:36:52 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 16:36:52 INFO Queuing job for member 5...
2026-07-27 16:36:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:36:52 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 16:36:55 INFO Found: ['5288010']
2026-07-27 16:37:00 INFO [TGCC-IRENE] Submitted job with ID:['5288010']
2026-07-27 16:37:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 16:37:00 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-27 16:37:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 16:37:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:00 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 16:37:00 INFO Queuing job for member 6...
2026-07-27 16:37:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:00 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 16:37:02 INFO Found: ['5288012']
2026-07-27 16:37:07 INFO [TGCC-IRENE] Submitted job with ID:['5288012']
2026-07-27 16:37:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 16:37:07 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-27 16:37:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 16:37:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:07 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 16:37:07 INFO Queuing job for member 7...
2026-07-27 16:37:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:07 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 16:37:10 INFO Found: ['5288018']
2026-07-27 16:37:15 INFO [TGCC-IRENE] Submitted job with ID:['5288018']
2026-07-27 16:37:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 16:37:15 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-27 16:37:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 16:37:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:15 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 16:37:15 INFO Queuing job for member 8...
2026-07-27 16:37:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:15 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 16:37:18 INFO Found: ['5288020']
2026-07-27 16:37:23 INFO [TGCC-IRENE] Submitted job with ID:['5288020']
2026-07-27 16:37:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 16:37:23 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-27 16:37:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 16:37:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:23 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 16:37:23 INFO Queuing job for member 9...
2026-07-27 16:37:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:23 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 16:37:23 INFO Found: ['5288021']
2026-07-27 16:37:28 INFO [TGCC-IRENE] Submitted job with ID:['5288021']
2026-07-27 16:37:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 16:37:28 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-27 16:37:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 16:37:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:28 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 16:37:28 INFO Queuing job for member 10...
2026-07-27 16:37:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:28 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 16:37:29 INFO Found: ['5288022']
2026-07-27 16:37:34 INFO [TGCC-IRENE] Submitted job with ID:['5288022']
2026-07-27 16:37:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 16:37:34 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-27 16:37:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 16:37:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:34 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 16:37:34 INFO Queuing job for member 11...
2026-07-27 16:37:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:34 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 16:37:35 INFO Found: ['5288023']
2026-07-27 16:37:40 INFO [TGCC-IRENE] Submitted job with ID:['5288023']
2026-07-27 16:37:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 16:37:40 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-27 16:37:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 16:37:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:40 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 16:37:40 INFO Queuing job for member 12...
2026-07-27 16:37:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:40 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 16:37:40 INFO Found: ['5288024']
2026-07-27 16:37:45 INFO [TGCC-IRENE] Submitted job with ID:['5288024']
2026-07-27 16:37:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 16:37:45 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-27 16:37:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 16:37:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:45 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 16:37:45 INFO Queuing job for member 13...
2026-07-27 16:37:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:45 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 16:37:46 INFO Found: ['5288025']
2026-07-27 16:37:51 INFO [TGCC-IRENE] Submitted job with ID:['5288025']
2026-07-27 16:37:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 16:37:51 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-27 16:37:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 16:37:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:51 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 16:37:51 INFO Queuing job for member 14...
2026-07-27 16:37:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:51 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 16:37:52 INFO Found: ['5288027']
2026-07-27 16:37:57 INFO [TGCC-IRENE] Submitted job with ID:['5288027']
2026-07-27 16:37:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:37:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 16:37:57 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-27 16:37:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 16:37:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:37:57 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 16:37:57 INFO Queuing job for member 15...
2026-07-27 16:37:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:37:57 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 16:37:58 INFO Found: ['5288028']
2026-07-27 16:38:03 INFO [TGCC-IRENE] Submitted job with ID:['5288028']
2026-07-27 16:38:03 INFO Checking job status ...
2026-07-27 16:38:03 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:38:03 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:38:03 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:38:18 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:38:18 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:38:18 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:38:18 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:38:20 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:38:20 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:38:35 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:38:35 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:38:36 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:38:36 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:38:36 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:38:36 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:38:36 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:38:36 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:38:51 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:38:51 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:38:51 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:38:53 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:38:53 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:39:08 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:39:08 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:39:10 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:39:10 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:39:10 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:39:10 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:39:10 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:39:25 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:39:25 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:39:25 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:39:25 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:39:25 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:39:25 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:39:26 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:39:26 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:39:41 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:39:41 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:39:41 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:39:56 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:39:56 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:39:56 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:40:11 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:40:11 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:40:12 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:40:12 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:40:12 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:40:12 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:40:12 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:40:12 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:40:27 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:40:27 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:40:27 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:40:43 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:40:43 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:40:45 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:40:45 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:41:00 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:41:00 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:41:02 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:41:02 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:41:17 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:41:17 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:41:18 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:41:18 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:41:18 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:41:18 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:41:18 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:41:18 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:41:33 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:41:33 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:41:35 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:41:35 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:41:50 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:41:50 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:41:50 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:42:05 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:42:05 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:42:05 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:42:21 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:42:21 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:42:21 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:42:37 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:42:37 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:42:38 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:42:38 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:42:53 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:42:53 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:42:53 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:43:08 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:43:08 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:43:08 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:43:10 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:43:10 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:43:25 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:43:25 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:43:26 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:43:26 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:43:41 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:43:41 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:43:41 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:43:43 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:43:43 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:43:58 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:43:58 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:44:00 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:44:00 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:44:15 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:44:15 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:44:16 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:44:16 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:44:31 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:44:31 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:44:31 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:44:46 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:44:46 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:44:46 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:45:02 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:45:02 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:45:02 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:45:17 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:45:17 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:45:17 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:45:17 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:45:18 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:45:20 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:45:20 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:45:20 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:45:35 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:45:35 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:45:35 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:45:35 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:45:37 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:45:38 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:45:38 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:45:38 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:45:38 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:45:38 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:45:53 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:45:53 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:45:55 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:45:55 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:45:55 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:45:55 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:45:55 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:46:10 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:46:10 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:46:10 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:46:25 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:46:25 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:46:25 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:46:40 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:46:41 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:46:41 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:46:56 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:46:56 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:46:56 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:47:13 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:47:13 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:47:13 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:47:28 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:47:28 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:47:28 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:47:28 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:47:28 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:47:30 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:47:30 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:47:45 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:47:45 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:47:46 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:47:48 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:47:48 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:47:48 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:47:48 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:47:48 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:47:48 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:48:03 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288023: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:48:03 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:48:05 INFO None 5288028: status RUNNING/PENDING
2026-07-27 16:48:05 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028']. Waiting...
2026-07-27 16:48:20 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288007: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288008: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:48:20 INFO None 5288020: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288021: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288022: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288023: status FINISHED
2026-07-27 16:48:21 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:48:21 INFO None 5288028: status FINISHED
2026-07-27 16:48:21 INFO Jobs still running: ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288024', '5288025', '5288027']. Waiting...
2026-07-27 16:48:36 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288007: status FINISHED
2026-07-27 16:48:36 INFO None 5288008: status FINISHED
2026-07-27 16:48:36 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288020: status FINISHED
2026-07-27 16:48:36 INFO None 5288021: status FINISHED
2026-07-27 16:48:36 INFO None 5288022: status FINISHED
2026-07-27 16:48:36 INFO None 5288023: status FINISHED
2026-07-27 16:48:36 INFO None 5288024: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:48:36 INFO None 5288028: status FINISHED
2026-07-27 16:48:36 INFO Jobs still running: ['5288006', '5288009', '5288010', '5288012', '5288018', '5288024', '5288025', '5288027']. Waiting...
2026-07-27 16:48:51 INFO None 5288006: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288007: status FINISHED
2026-07-27 16:48:51 INFO None 5288008: status FINISHED
2026-07-27 16:48:51 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288012: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288018: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288020: status FINISHED
2026-07-27 16:48:51 INFO None 5288021: status FINISHED
2026-07-27 16:48:51 INFO None 5288022: status FINISHED
2026-07-27 16:48:51 INFO None 5288023: status FINISHED
2026-07-27 16:48:51 INFO None 5288024: status FINISHED
2026-07-27 16:48:51 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:48:51 INFO None 5288028: status FINISHED
2026-07-27 16:48:51 INFO Jobs still running: ['5288006', '5288009', '5288010', '5288012', '5288018', '5288025', '5288027']. Waiting...
2026-07-27 16:49:06 INFO None 5288006: status FINISHED
2026-07-27 16:49:06 INFO None 5288007: status FINISHED
2026-07-27 16:49:06 INFO None 5288008: status FINISHED
2026-07-27 16:49:06 INFO None 5288009: status RUNNING/PENDING
2026-07-27 16:49:06 INFO None 5288010: status RUNNING/PENDING
2026-07-27 16:49:06 INFO None 5288012: status FINISHED
2026-07-27 16:49:06 INFO None 5288018: status FINISHED
2026-07-27 16:49:06 INFO None 5288020: status FINISHED
2026-07-27 16:49:06 INFO None 5288021: status FINISHED
2026-07-27 16:49:06 INFO None 5288022: status FINISHED
2026-07-27 16:49:06 INFO None 5288023: status FINISHED
2026-07-27 16:49:06 INFO None 5288024: status FINISHED
2026-07-27 16:49:06 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:49:06 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:49:06 INFO None 5288028: status FINISHED
2026-07-27 16:49:06 INFO Jobs still running: ['5288009', '5288010', '5288025', '5288027']. Waiting...
2026-07-27 16:49:22 INFO None 5288006: status FINISHED
2026-07-27 16:49:22 INFO None 5288007: status FINISHED
2026-07-27 16:49:22 INFO None 5288008: status FINISHED
2026-07-27 16:49:22 INFO None 5288009: status FINISHED
2026-07-27 16:49:22 INFO None 5288010: status FINISHED
2026-07-27 16:49:22 INFO None 5288012: status FINISHED
2026-07-27 16:49:22 INFO None 5288018: status FINISHED
2026-07-27 16:49:22 INFO None 5288020: status FINISHED
2026-07-27 16:49:22 INFO None 5288021: status FINISHED
2026-07-27 16:49:22 INFO None 5288022: status FINISHED
2026-07-27 16:49:22 INFO None 5288023: status FINISHED
2026-07-27 16:49:22 INFO None 5288024: status FINISHED
2026-07-27 16:49:22 INFO None 5288025: status RUNNING/PENDING
2026-07-27 16:49:22 INFO None 5288027: status RUNNING/PENDING
2026-07-27 16:49:22 INFO None 5288028: status FINISHED
2026-07-27 16:49:22 INFO Jobs still running: ['5288025', '5288027']. Waiting...
2026-07-27 16:49:37 INFO None 5288006: status FINISHED
2026-07-27 16:49:37 INFO None 5288007: status FINISHED
2026-07-27 16:49:37 INFO None 5288008: status FINISHED
2026-07-27 16:49:37 INFO None 5288009: status FINISHED
2026-07-27 16:49:37 INFO None 5288010: status FINISHED
2026-07-27 16:49:37 INFO None 5288012: status FINISHED
2026-07-27 16:49:37 INFO None 5288018: status FINISHED
2026-07-27 16:49:37 INFO None 5288020: status FINISHED
2026-07-27 16:49:37 INFO None 5288021: status FINISHED
2026-07-27 16:49:38 INFO None 5288022: status FINISHED
2026-07-27 16:49:38 INFO None 5288023: status FINISHED
2026-07-27 16:49:38 INFO None 5288024: status FINISHED
2026-07-27 16:49:38 INFO None 5288025: status FINISHED
2026-07-27 16:49:38 INFO None 5288027: status FINISHED
2026-07-27 16:49:38 INFO None 5288028: status FINISHED
2026-07-27 16:49:38 INFO Jobs ['5288006', '5288007', '5288008', '5288009', '5288010', '5288012', '5288018', '5288020', '5288021', '5288022', '5288023', '5288024', '5288025', '5288027', '5288028'] have finished
2026-07-27 16:49:38 INFO Checking restart files were created ...
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc(668832435 bytes)
2026-07-27 16:49:38 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc(668832435 bytes)
2026-07-27 16:49:38 INFO  Run_model() completed successfully.
2026-07-27 16:49:38 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:49:38 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 01:00:00 days=153073 seconds=3600
2026-07-27 16:49:38 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 16:49:38 INFO [TIME] increment current_time 2020-02-07 00:00:00 -> 2020-02-07 01:00:00
2026-07-27 16:49:38 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:49:38 INFO ---------->>> Running process_satellite_data()
2026-07-27 16:49:38 INFO [DART] No satellite data found, skipping assimilation
2026-07-27 16:49:38 INFO after_assimilation() skipped
2026-07-27 16:49:38 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 16:49:38 INFO [TIME] step_end current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:49:38 INFO [TIME] step_start current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 16:49:38 INFO [TIME] window start=2020-02-07 01:00:00 end=2020-02-07 09:00:00 run_hours=8 has_assimilation=True
2026-07-27 16:49:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:49:39 INFO Hourly dataset computed and listing created
2026-07-27 16:49:54 INFO Hourly dataset computed
2026-07-27 16:49:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:49:56 INFO Hourly dataset computed and listing created
2026-07-27 16:49:58 INFO Hourly dataset computed
2026-07-27 16:49:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:49:59 INFO Hourly dataset computed and listing created
2026-07-27 16:50:01 INFO Hourly dataset computed
2026-07-27 16:50:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:02 INFO Hourly dataset computed and listing created
2026-07-27 16:50:04 INFO Hourly dataset computed
2026-07-27 16:50:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:06 INFO Hourly dataset computed and listing created
2026-07-27 16:50:08 INFO Hourly dataset computed
2026-07-27 16:50:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:09 INFO Hourly dataset computed and listing created
2026-07-27 16:50:11 INFO Hourly dataset computed
2026-07-27 16:50:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:12 INFO Hourly dataset computed and listing created
2026-07-27 16:50:14 INFO Hourly dataset computed
2026-07-27 16:50:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:16 INFO Hourly dataset computed and listing created
2026-07-27 16:50:18 INFO Hourly dataset computed
2026-07-27 16:50:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:19 INFO Hourly dataset computed and listing created
2026-07-27 16:50:21 INFO Hourly dataset computed
2026-07-27 16:50:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:22 INFO Hourly dataset computed and listing created
2026-07-27 16:50:24 INFO Hourly dataset computed
2026-07-27 16:50:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:25 INFO Hourly dataset computed and listing created
2026-07-27 16:50:27 INFO Hourly dataset computed
2026-07-27 16:50:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:28 INFO Hourly dataset computed and listing created
2026-07-27 16:50:30 INFO Hourly dataset computed
2026-07-27 16:50:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:31 INFO Hourly dataset computed and listing created
2026-07-27 16:50:34 INFO Hourly dataset computed
2026-07-27 16:50:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:35 INFO Hourly dataset computed and listing created
2026-07-27 16:50:37 INFO Hourly dataset computed
2026-07-27 16:50:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 16:50:38 INFO Hourly dataset computed and listing created
2026-07-27 16:50:40 INFO Hourly dataset computed
2026-07-27 16:50:40 INFO ---------->>> Running CHIMERE model from 2020-02-07 01:00:00 to 2020-02-07 09:00:00
2026-07-27 16:50:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:50:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 16:50:40 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020700_1_ENS1.nc
2026-07-27 16:50:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 16:50:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:50:40 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 16:50:40 INFO Queuing job for member 1...
2026-07-27 16:50:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:50:40 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 16:50:42 INFO Found: ['5288182']
2026-07-27 16:50:47 INFO [TGCC-IRENE] Submitted job with ID:['5288182']
2026-07-27 16:50:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:50:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 16:50:47 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020700_1_ENS2.nc
2026-07-27 16:50:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 16:50:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:50:48 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 16:50:48 INFO Queuing job for member 2...
2026-07-27 16:50:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:50:48 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 16:50:48 INFO Found: ['5288183']
2026-07-27 16:50:53 INFO [TGCC-IRENE] Submitted job with ID:['5288183']
2026-07-27 16:50:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:50:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 16:50:53 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020700_1_ENS3.nc
2026-07-27 16:50:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 16:50:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:50:53 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 16:50:53 INFO Queuing job for member 3...
2026-07-27 16:50:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:50:53 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 16:50:54 INFO Found: ['5288184']
2026-07-27 16:50:59 INFO [TGCC-IRENE] Submitted job with ID:['5288184']
2026-07-27 16:50:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:50:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 16:50:59 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020700_1_ENS4.nc
2026-07-27 16:50:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 16:50:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:50:59 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 16:50:59 INFO Queuing job for member 4...
2026-07-27 16:50:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:50:59 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 16:51:00 INFO Found: ['5288185']
2026-07-27 16:51:05 INFO [TGCC-IRENE] Submitted job with ID:['5288185']
2026-07-27 16:51:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 16:51:05 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020700_1_ENS5.nc
2026-07-27 16:51:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 16:51:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:05 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 16:51:05 INFO Queuing job for member 5...
2026-07-27 16:51:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:05 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 16:51:05 INFO Found: ['5288188']
2026-07-27 16:51:10 INFO [TGCC-IRENE] Submitted job with ID:['5288188']
2026-07-27 16:51:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 16:51:10 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020700_1_ENS6.nc
2026-07-27 16:51:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 16:51:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:10 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 16:51:10 INFO Queuing job for member 6...
2026-07-27 16:51:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:10 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 16:51:11 INFO Found: ['5288189']
2026-07-27 16:51:16 INFO [TGCC-IRENE] Submitted job with ID:['5288189']
2026-07-27 16:51:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 16:51:16 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020700_1_ENS7.nc
2026-07-27 16:51:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 16:51:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:16 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 16:51:16 INFO Queuing job for member 7...
2026-07-27 16:51:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:16 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 16:51:17 INFO Found: ['5288190']
2026-07-27 16:51:22 INFO [TGCC-IRENE] Submitted job with ID:['5288190']
2026-07-27 16:51:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 16:51:22 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020700_1_ENS8.nc
2026-07-27 16:51:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 16:51:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:22 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 16:51:22 INFO Queuing job for member 8...
2026-07-27 16:51:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:22 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 16:51:23 INFO Found: ['5288191']
2026-07-27 16:51:28 INFO [TGCC-IRENE] Submitted job with ID:['5288191']
2026-07-27 16:51:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 16:51:28 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020700_1_ENS9.nc
2026-07-27 16:51:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 16:51:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:28 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 16:51:28 INFO Queuing job for member 9...
2026-07-27 16:51:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:28 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 16:51:28 INFO Found: ['5288192']
2026-07-27 16:51:33 INFO [TGCC-IRENE] Submitted job with ID:['5288192']
2026-07-27 16:51:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 16:51:33 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020700_1_ENS10.nc
2026-07-27 16:51:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 16:51:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:33 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 16:51:33 INFO Queuing job for member 10...
2026-07-27 16:51:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:33 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 16:51:35 INFO Found: ['5288200']
2026-07-27 16:51:40 INFO [TGCC-IRENE] Submitted job with ID:['5288200']
2026-07-27 16:51:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 16:51:40 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020700_1_ENS11.nc
2026-07-27 16:51:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 16:51:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:40 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 16:51:40 INFO Queuing job for member 11...
2026-07-27 16:51:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:40 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 16:51:43 INFO Found: ['5288205']
2026-07-27 16:51:48 INFO [TGCC-IRENE] Submitted job with ID:['5288205']
2026-07-27 16:51:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 16:51:48 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020700_1_ENS12.nc
2026-07-27 16:51:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 16:51:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:48 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 16:51:48 INFO Queuing job for member 12...
2026-07-27 16:51:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:48 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 16:51:50 INFO Found: ['5288209']
2026-07-27 16:51:55 INFO [TGCC-IRENE] Submitted job with ID:['5288209']
2026-07-27 16:51:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:51:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 16:51:55 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020700_1_ENS13.nc
2026-07-27 16:51:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 16:51:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:51:55 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 16:51:55 INFO Queuing job for member 13...
2026-07-27 16:51:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:51:55 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 16:51:58 INFO Found: ['5288214']
2026-07-27 16:52:03 INFO [TGCC-IRENE] Submitted job with ID:['5288214']
2026-07-27 16:52:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:52:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 16:52:03 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020700_1_ENS14.nc
2026-07-27 16:52:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 16:52:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:52:03 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 16:52:03 INFO Queuing job for member 14...
2026-07-27 16:52:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:52:03 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 16:52:05 INFO Found: ['5288219']
2026-07-27 16:52:10 INFO [TGCC-IRENE] Submitted job with ID:['5288219']
2026-07-27 16:52:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 16:52:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 16:52:10 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020700_1_ENS15.nc
2026-07-27 16:52:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 16:52:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 16:52:10 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 16:52:10 INFO Queuing job for member 15...
2026-07-27 16:52:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 16:52:10 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 16:52:13 INFO Found: ['5288225']
2026-07-27 16:52:18 INFO [TGCC-IRENE] Submitted job with ID:['5288225']
2026-07-27 16:52:18 INFO Checking job status ...
2026-07-27 16:52:18 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:52:18 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:52:20 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:52:20 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:52:20 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:52:20 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:52:35 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:52:35 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:52:36 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:52:36 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:52:36 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:52:36 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:52:36 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:52:51 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:52:51 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:52:53 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:52:53 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:52:53 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:52:53 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:52:53 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:53:08 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:53:08 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:53:08 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:53:23 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:53:23 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:53:24 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:53:24 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:53:39 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:53:39 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:53:39 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:53:55 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:53:55 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:53:55 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:54:10 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:54:10 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:54:10 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:54:10 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:54:12 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:54:12 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:54:27 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:54:27 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:54:27 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:54:27 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:54:27 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:54:27 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:54:28 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:54:28 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:54:43 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:54:43 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:54:43 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:54:43 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:54:45 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:54:45 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:55:00 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:55:00 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:55:02 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:55:02 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:55:02 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:55:17 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:55:17 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:55:18 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:55:18 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:55:33 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:55:33 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:55:33 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:55:48 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:55:48 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:55:48 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:56:05 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:56:05 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:56:05 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:56:20 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:56:20 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:56:21 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:56:21 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:56:21 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:56:23 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:56:23 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:56:23 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:56:38 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:56:38 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:56:38 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:56:53 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:56:53 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:56:55 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:56:55 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:56:55 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:56:55 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:56:55 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:57:10 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:57:10 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:57:10 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:57:11 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:57:11 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:57:26 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:57:26 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:57:28 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:57:28 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:57:43 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:57:43 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:57:43 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:57:58 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:57:58 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:57:58 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:57:58 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:57:58 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:57:58 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:57:59 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:57:59 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:58:14 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:58:14 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:58:14 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:58:30 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:58:30 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:58:30 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:58:45 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:58:45 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:58:45 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:59:00 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:59:00 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:59:00 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:59:03 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:59:03 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:59:18 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:59:18 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:59:18 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:59:18 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:59:18 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:59:20 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:59:20 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:59:35 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:59:35 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:59:36 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:59:36 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:59:36 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 16:59:51 INFO None 5288182: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288183: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288184: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288185: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288188: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288189: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288190: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288191: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288192: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288200: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288205: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288209: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288214: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288219: status RUNNING/PENDING
2026-07-27 16:59:51 INFO None 5288225: status RUNNING/PENDING
2026-07-27 16:59:51 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:00:06 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:00:06 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:00:06 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:00:21 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:00:21 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:00:21 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:00:37 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:00:37 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:00:37 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:00:52 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:00:52 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:00:53 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:00:53 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:01:08 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:01:08 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:01:08 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:01:10 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:01:10 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:01:25 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:01:25 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:01:25 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:01:40 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:01:40 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:01:40 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:01:40 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:01:40 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:01:42 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:01:42 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:01:42 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:01:42 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:01:42 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:01:43 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:01:43 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:01:43 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:01:43 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:01:43 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:01:43 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:01:58 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:01:58 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:01:58 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:02:13 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:02:13 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:02:13 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:02:28 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:02:28 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:02:28 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:02:43 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:02:43 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:02:44 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:02:44 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:03:00 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:03:00 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:03:00 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:03:15 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:03:15 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:03:16 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:03:16 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:03:18 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:03:18 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:03:18 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:03:33 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:03:33 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:03:33 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:03:48 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:03:48 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:03:48 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:04:03 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:04:03 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:04:03 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:04:05 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:04:05 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:04:20 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:04:20 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:04:21 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:04:21 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:04:36 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:04:36 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:04:36 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:04:51 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:04:51 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:04:51 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:05:07 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288209: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:05:07 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:05:07 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:05:22 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:05:22 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288209: status FINISHED
2026-07-27 17:05:23 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:05:23 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:05:23 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:05:38 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:05:38 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288190: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288191: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288200: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288205: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288209: status FINISHED
2026-07-27 17:05:40 INFO None 5288214: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288219: status RUNNING/PENDING
2026-07-27 17:05:40 INFO None 5288225: status RUNNING/PENDING
2026-07-27 17:05:40 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288214', '5288219', '5288225']. Waiting...
2026-07-27 17:05:55 INFO None 5288182: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288189: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288190: status FINISHED
2026-07-27 17:05:55 INFO None 5288191: status FINISHED
2026-07-27 17:05:55 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:05:55 INFO None 5288200: status FINISHED
2026-07-27 17:05:55 INFO None 5288205: status FINISHED
2026-07-27 17:05:55 INFO None 5288209: status FINISHED
2026-07-27 17:05:55 INFO None 5288214: status FINISHED
2026-07-27 17:05:57 INFO None 5288219: status FINISHED
2026-07-27 17:05:57 INFO None 5288225: status FINISHED
2026-07-27 17:05:57 INFO Jobs still running: ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288192']. Waiting...
2026-07-27 17:06:12 INFO None 5288182: status FINISHED
2026-07-27 17:06:12 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:06:12 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:06:12 INFO None 5288185: status RUNNING/PENDING
2026-07-27 17:06:12 INFO None 5288188: status RUNNING/PENDING
2026-07-27 17:06:13 INFO None 5288189: status FINISHED
2026-07-27 17:06:13 INFO None 5288190: status FINISHED
2026-07-27 17:06:13 INFO None 5288191: status FINISHED
2026-07-27 17:06:13 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:06:13 INFO None 5288200: status FINISHED
2026-07-27 17:06:13 INFO None 5288205: status FINISHED
2026-07-27 17:06:13 INFO None 5288209: status FINISHED
2026-07-27 17:06:13 INFO None 5288214: status FINISHED
2026-07-27 17:06:13 INFO None 5288219: status FINISHED
2026-07-27 17:06:13 INFO None 5288225: status FINISHED
2026-07-27 17:06:13 INFO Jobs still running: ['5288183', '5288184', '5288185', '5288188', '5288192']. Waiting...
2026-07-27 17:06:28 INFO None 5288182: status FINISHED
2026-07-27 17:06:28 INFO None 5288183: status RUNNING/PENDING
2026-07-27 17:06:28 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:06:28 INFO None 5288185: status FINISHED
2026-07-27 17:06:28 INFO None 5288188: status FINISHED
2026-07-27 17:06:28 INFO None 5288189: status FINISHED
2026-07-27 17:06:28 INFO None 5288190: status FINISHED
2026-07-27 17:06:28 INFO None 5288191: status FINISHED
2026-07-27 17:06:28 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:06:28 INFO None 5288200: status FINISHED
2026-07-27 17:06:28 INFO None 5288205: status FINISHED
2026-07-27 17:06:28 INFO None 5288209: status FINISHED
2026-07-27 17:06:28 INFO None 5288214: status FINISHED
2026-07-27 17:06:28 INFO None 5288219: status FINISHED
2026-07-27 17:06:30 INFO None 5288225: status FINISHED
2026-07-27 17:06:30 INFO Jobs still running: ['5288183', '5288184', '5288192']. Waiting...
2026-07-27 17:06:45 INFO None 5288182: status FINISHED
2026-07-27 17:06:45 INFO None 5288183: status FINISHED
2026-07-27 17:06:45 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:06:45 INFO None 5288185: status FINISHED
2026-07-27 17:06:45 INFO None 5288188: status FINISHED
2026-07-27 17:06:45 INFO None 5288189: status FINISHED
2026-07-27 17:06:45 INFO None 5288190: status FINISHED
2026-07-27 17:06:45 INFO None 5288191: status FINISHED
2026-07-27 17:06:45 INFO None 5288192: status RUNNING/PENDING
2026-07-27 17:06:45 INFO None 5288200: status FINISHED
2026-07-27 17:06:45 INFO None 5288205: status FINISHED
2026-07-27 17:06:45 INFO None 5288209: status FINISHED
2026-07-27 17:06:45 INFO None 5288214: status FINISHED
2026-07-27 17:06:45 INFO None 5288219: status FINISHED
2026-07-27 17:06:45 INFO None 5288225: status FINISHED
2026-07-27 17:06:45 INFO Jobs still running: ['5288184', '5288192']. Waiting...
2026-07-27 17:07:00 INFO None 5288182: status FINISHED
2026-07-27 17:07:00 INFO None 5288183: status FINISHED
2026-07-27 17:07:00 INFO None 5288184: status RUNNING/PENDING
2026-07-27 17:07:00 INFO None 5288185: status FINISHED
2026-07-27 17:07:00 INFO None 5288188: status FINISHED
2026-07-27 17:07:00 INFO None 5288189: status FINISHED
2026-07-27 17:07:00 INFO None 5288190: status FINISHED
2026-07-27 17:07:00 INFO None 5288191: status FINISHED
2026-07-27 17:07:00 INFO None 5288192: status FINISHED
2026-07-27 17:07:00 INFO None 5288200: status FINISHED
2026-07-27 17:07:00 INFO None 5288205: status FINISHED
2026-07-27 17:07:01 INFO None 5288209: status FINISHED
2026-07-27 17:07:01 INFO None 5288214: status FINISHED
2026-07-27 17:07:01 INFO None 5288219: status FINISHED
2026-07-27 17:07:01 INFO None 5288225: status FINISHED
2026-07-27 17:07:01 INFO Jobs still running: ['5288184']. Waiting...
2026-07-27 17:07:16 INFO None 5288182: status FINISHED
2026-07-27 17:07:16 INFO None 5288183: status FINISHED
2026-07-27 17:07:16 INFO None 5288184: status FINISHED
2026-07-27 17:07:16 INFO None 5288185: status FINISHED
2026-07-27 17:07:16 INFO None 5288188: status FINISHED
2026-07-27 17:07:16 INFO None 5288189: status FINISHED
2026-07-27 17:07:16 INFO None 5288190: status FINISHED
2026-07-27 17:07:16 INFO None 5288191: status FINISHED
2026-07-27 17:07:16 INFO None 5288192: status FINISHED
2026-07-27 17:07:16 INFO None 5288200: status FINISHED
2026-07-27 17:07:16 INFO None 5288205: status FINISHED
2026-07-27 17:07:16 INFO None 5288209: status FINISHED
2026-07-27 17:07:16 INFO None 5288214: status FINISHED
2026-07-27 17:07:16 INFO None 5288219: status FINISHED
2026-07-27 17:07:16 INFO None 5288225: status FINISHED
2026-07-27 17:07:16 INFO Jobs ['5288182', '5288183', '5288184', '5288185', '5288188', '5288189', '5288190', '5288191', '5288192', '5288200', '5288205', '5288209', '5288214', '5288219', '5288225'] have finished
2026-07-27 17:07:16 INFO Checking restart files were created ...
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc(3005806795 bytes)
2026-07-27 17:07:16 INFO  Run_model() completed successfully.
2026-07-27 17:07:16 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 01:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:07:16 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 09:00:00 days=153073 seconds=32400
2026-07-27 17:07:16 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 17:07:16 INFO [TIME] increment current_time 2020-02-07 01:00:00 -> 2020-02-07 09:00:00
2026-07-27 17:07:16 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:07:16 INFO ---------->>> Running process_satellite_data()
2026-07-27 17:07:16 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12016.nc
2026-07-27 17:07:16 INFO ---------->>> Running run_obs_converter()
2026-07-27 17:07:16 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-27 17:07:16 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_32090_153073.out
2026-07-27 17:07:16 INFO ---------->>> Running DART
2026-07-27 17:07:16 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 17:07:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_1_out_toDART.nc
2026-07-27 17:07:16 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_1_out_toDART.nc
2026-07-27 17:07:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_1_out_toDART.nc
2026-07-27 17:07:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_1_out_toDART.nc
2026-07-27 17:07:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_1_out_toDART.nc
2026-07-27 17:07:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_1_out_toDART.nc
2026-07-27 17:07:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_1_out_toDART.nc
2026-07-27 17:07:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_1_out_toDART.nc
2026-07-27 17:07:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_1_out_toDART.nc
2026-07-27 17:07:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_1_out_toDART.nc
2026-07-27 17:07:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_1_out_toDART.nc
2026-07-27 17:07:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_1_out_toDART.nc
2026-07-27 17:07:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_1_out_toDART.nc
2026-07-27 17:07:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_1_out_toDART.nc
2026-07-27 17:07:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020701_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_1_out_toDART.nc
2026-07-27 17:07:22 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 17:07:22 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 17:07:22 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 17:07:22 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 17:07:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 17:07:22 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 17:07:31 INFO Found: []
2026-07-27 17:07:31 INFO No job id returned by command ./run_filter.bsh
2026-07-27 17:07:31 INFO No monitoring will be performed
2026-07-27 17:07:31 INFO Moving DART output files to analysis and preassim directories for date 2020020709 if present ...
2026-07-27 17:07:31 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020709'
2026-07-27 17:07:31 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 17:07:31 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 17:07:31 INFO run_dart() is DONE.
2026-07-27 17:07:31 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 17:07:32 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-07-27 17:07:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-27 17:07:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:07:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-07-27 17:07:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-27 17:07:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:07:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-07-27 17:07:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-27 17:08:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:08:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-07-27 17:08:14 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-27 17:08:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:08:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-07-27 17:08:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-27 17:08:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:08:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-07-27 17:08:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-27 17:08:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:08:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-07-27 17:08:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-27 17:09:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:09:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-07-27 17:09:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-27 17:09:23 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:09:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-07-27 17:09:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-27 17:09:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:09:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-07-27 17:09:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-27 17:09:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:09:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-07-27 17:09:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-27 17:10:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:10:13 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-07-27 17:10:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-27 17:10:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:10:27 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-07-27 17:10:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-27 17:10:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:10:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-07-27 17:10:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-27 17:10:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:10:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-07-27 17:10:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-27 17:11:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:11:13 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 17:11:13 INFO [TIME] step_end current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:11:13 INFO [TIME] step_start current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:11:13 INFO [TIME] window start=2020-02-07 09:00:00 end=2020-02-07 11:00:00 run_hours=2 has_assimilation=True
2026-07-27 17:11:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:15 INFO Hourly dataset computed and listing created
2026-07-27 17:11:18 INFO Hourly dataset computed
2026-07-27 17:11:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:19 INFO Hourly dataset computed and listing created
2026-07-27 17:11:20 INFO Hourly dataset computed
2026-07-27 17:11:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:21 INFO Hourly dataset computed and listing created
2026-07-27 17:11:22 INFO Hourly dataset computed
2026-07-27 17:11:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:23 INFO Hourly dataset computed and listing created
2026-07-27 17:11:23 INFO Hourly dataset computed
2026-07-27 17:11:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:24 INFO Hourly dataset computed and listing created
2026-07-27 17:11:25 INFO Hourly dataset computed
2026-07-27 17:11:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:26 INFO Hourly dataset computed and listing created
2026-07-27 17:11:27 INFO Hourly dataset computed
2026-07-27 17:11:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:28 INFO Hourly dataset computed and listing created
2026-07-27 17:11:29 INFO Hourly dataset computed
2026-07-27 17:11:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:30 INFO Hourly dataset computed and listing created
2026-07-27 17:11:31 INFO Hourly dataset computed
2026-07-27 17:11:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:32 INFO Hourly dataset computed and listing created
2026-07-27 17:11:32 INFO Hourly dataset computed
2026-07-27 17:11:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:33 INFO Hourly dataset computed and listing created
2026-07-27 17:11:34 INFO Hourly dataset computed
2026-07-27 17:11:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:35 INFO Hourly dataset computed and listing created
2026-07-27 17:11:36 INFO Hourly dataset computed
2026-07-27 17:11:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:37 INFO Hourly dataset computed and listing created
2026-07-27 17:11:38 INFO Hourly dataset computed
2026-07-27 17:11:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:39 INFO Hourly dataset computed and listing created
2026-07-27 17:11:40 INFO Hourly dataset computed
2026-07-27 17:11:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:40 INFO Hourly dataset computed and listing created
2026-07-27 17:11:41 INFO Hourly dataset computed
2026-07-27 17:11:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:11:42 INFO Hourly dataset computed and listing created
2026-07-27 17:11:43 INFO Hourly dataset computed
2026-07-27 17:11:43 INFO ---------->>> Running CHIMERE model from 2020-02-07 09:00:00 to 2020-02-07 11:00:00
2026-07-27 17:11:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:11:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 17:11:43 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020701_8_ENS1.nc
2026-07-27 17:11:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 17:11:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:11:43 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 17:11:43 INFO Queuing job for member 1...
2026-07-27 17:11:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:11:43 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 17:11:44 INFO Found: ['5288415']
2026-07-27 17:11:49 INFO [TGCC-IRENE] Submitted job with ID:['5288415']
2026-07-27 17:11:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:11:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 17:11:49 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020701_8_ENS2.nc
2026-07-27 17:11:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 17:11:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:11:49 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 17:11:49 INFO Queuing job for member 2...
2026-07-27 17:11:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:11:49 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 17:11:50 INFO Found: ['5288416']
2026-07-27 17:11:55 INFO [TGCC-IRENE] Submitted job with ID:['5288416']
2026-07-27 17:11:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:11:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 17:11:55 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020701_8_ENS3.nc
2026-07-27 17:11:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 17:11:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:11:55 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 17:11:55 INFO Queuing job for member 3...
2026-07-27 17:11:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:11:55 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 17:11:58 INFO Found: ['5288417']
2026-07-27 17:12:03 INFO [TGCC-IRENE] Submitted job with ID:['5288417']
2026-07-27 17:12:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 17:12:03 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020701_8_ENS4.nc
2026-07-27 17:12:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 17:12:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:03 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 17:12:03 INFO Queuing job for member 4...
2026-07-27 17:12:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:03 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 17:12:05 INFO Found: ['5288419']
2026-07-27 17:12:10 INFO [TGCC-IRENE] Submitted job with ID:['5288419']
2026-07-27 17:12:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 17:12:10 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020701_8_ENS5.nc
2026-07-27 17:12:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 17:12:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:10 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 17:12:10 INFO Queuing job for member 5...
2026-07-27 17:12:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:10 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 17:12:13 INFO Found: ['5288420']
2026-07-27 17:12:18 INFO [TGCC-IRENE] Submitted job with ID:['5288420']
2026-07-27 17:12:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 17:12:18 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020701_8_ENS6.nc
2026-07-27 17:12:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 17:12:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:18 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 17:12:18 INFO Queuing job for member 6...
2026-07-27 17:12:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:18 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 17:12:20 INFO Found: ['5288421']
2026-07-27 17:12:25 INFO [TGCC-IRENE] Submitted job with ID:['5288421']
2026-07-27 17:12:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 17:12:25 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020701_8_ENS7.nc
2026-07-27 17:12:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 17:12:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:25 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 17:12:25 INFO Queuing job for member 7...
2026-07-27 17:12:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:25 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 17:12:28 INFO Found: ['5288422']
2026-07-27 17:12:33 INFO [TGCC-IRENE] Submitted job with ID:['5288422']
2026-07-27 17:12:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 17:12:33 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020701_8_ENS8.nc
2026-07-27 17:12:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 17:12:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:33 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 17:12:33 INFO Queuing job for member 8...
2026-07-27 17:12:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:33 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 17:12:35 INFO Found: ['5288423']
2026-07-27 17:12:40 INFO [TGCC-IRENE] Submitted job with ID:['5288423']
2026-07-27 17:12:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 17:12:40 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020701_8_ENS9.nc
2026-07-27 17:12:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 17:12:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:40 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 17:12:40 INFO Queuing job for member 9...
2026-07-27 17:12:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:40 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 17:12:43 INFO Found: ['5288425']
2026-07-27 17:12:48 INFO [TGCC-IRENE] Submitted job with ID:['5288425']
2026-07-27 17:12:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 17:12:48 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020701_8_ENS10.nc
2026-07-27 17:12:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 17:12:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:48 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 17:12:48 INFO Queuing job for member 10...
2026-07-27 17:12:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:48 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 17:12:50 INFO Found: ['5288428']
2026-07-27 17:12:55 INFO [TGCC-IRENE] Submitted job with ID:['5288428']
2026-07-27 17:12:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:12:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 17:12:55 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020701_8_ENS11.nc
2026-07-27 17:12:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 17:12:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:12:55 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 17:12:55 INFO Queuing job for member 11...
2026-07-27 17:12:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:12:55 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 17:12:58 INFO Found: ['5288430']
2026-07-27 17:13:03 INFO [TGCC-IRENE] Submitted job with ID:['5288430']
2026-07-27 17:13:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:13:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 17:13:03 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020701_8_ENS12.nc
2026-07-27 17:13:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 17:13:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:13:03 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 17:13:03 INFO Queuing job for member 12...
2026-07-27 17:13:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:13:03 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 17:13:05 INFO Found: ['5288433']
2026-07-27 17:13:10 INFO [TGCC-IRENE] Submitted job with ID:['5288433']
2026-07-27 17:13:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:13:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 17:13:10 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020701_8_ENS13.nc
2026-07-27 17:13:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 17:13:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:13:10 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 17:13:10 INFO Queuing job for member 13...
2026-07-27 17:13:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:13:10 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 17:13:13 INFO Found: ['5288437']
2026-07-27 17:13:18 INFO [TGCC-IRENE] Submitted job with ID:['5288437']
2026-07-27 17:13:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:13:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 17:13:18 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020701_8_ENS14.nc
2026-07-27 17:13:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 17:13:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:13:18 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 17:13:18 INFO Queuing job for member 14...
2026-07-27 17:13:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:13:18 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 17:13:19 INFO Found: ['5288443']
2026-07-27 17:13:24 INFO [TGCC-IRENE] Submitted job with ID:['5288443']
2026-07-27 17:13:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:13:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 17:13:24 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020701_8_ENS15.nc
2026-07-27 17:13:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 17:13:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:13:24 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 17:13:24 INFO Queuing job for member 15...
2026-07-27 17:13:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:13:24 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 17:13:24 INFO Found: ['5288450']
2026-07-27 17:13:29 INFO [TGCC-IRENE] Submitted job with ID:['5288450']
2026-07-27 17:13:29 INFO Checking job status ...
2026-07-27 17:13:29 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:13:29 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:13:29 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:13:29 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:13:30 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:13:30 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:13:45 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:13:45 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:13:45 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:14:00 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:14:00 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:14:00 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:14:17 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:14:17 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:14:17 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:14:32 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:14:33 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:14:35 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:14:35 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:14:35 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:14:35 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:14:35 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:14:50 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:14:50 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:14:50 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:15:05 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:15:05 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:15:05 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:15:05 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:15:05 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:15:05 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:15:07 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:15:08 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:15:08 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:15:08 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:15:08 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:15:23 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:15:23 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:15:23 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:15:38 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:15:38 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:15:38 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:15:53 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:15:53 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:15:53 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:16:08 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:16:08 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:16:08 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:16:08 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:16:08 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:16:09 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:16:09 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:16:25 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:16:25 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:16:26 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:16:26 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:16:41 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:16:41 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:16:41 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:16:41 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:16:43 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:16:43 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:16:58 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:16:58 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:17:00 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:17:00 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:17:15 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:17:15 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:17:15 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:17:16 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:17:16 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:17:31 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:17:31 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:17:33 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:17:33 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:17:48 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:17:48 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:17:48 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:18:03 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:18:03 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:18:04 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:18:04 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:18:19 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:18:19 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:18:19 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:18:34 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:18:34 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:18:34 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:18:49 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:18:49 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:18:50 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:18:52 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:18:52 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:19:07 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:19:07 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:19:07 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:19:22 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:19:22 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:19:24 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:19:24 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:19:24 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:19:24 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:19:24 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:19:39 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:19:39 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:19:40 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:19:40 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:19:40 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:19:40 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:19:55 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:19:55 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:19:55 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:20:10 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:20:10 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:20:10 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:20:25 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:20:25 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:20:25 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:20:41 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:20:41 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:20:41 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:20:56 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:20:56 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:20:56 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:20:57 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:20:57 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:20:57 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:20:57 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:20:57 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:20:59 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:20:59 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:21:14 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:21:14 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:21:14 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:21:29 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:21:29 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:21:29 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:21:29 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:21:29 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:21:29 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:21:31 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:21:32 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:21:32 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:21:47 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:21:47 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:21:47 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:22:02 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:22:02 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:22:02 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:22:04 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:22:04 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:22:19 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:22:19 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:22:19 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:22:35 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:22:35 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:22:35 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:22:50 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:22:50 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:22:50 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:23:06 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:23:06 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:23:07 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:23:07 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:23:22 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:23:22 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:23:22 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:23:37 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:23:37 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:23:37 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:23:37 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:23:39 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:23:40 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:23:40 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:23:40 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:23:40 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:23:40 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:23:55 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:23:55 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:23:57 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:23:57 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:24:12 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:24:12 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:24:12 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:24:27 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:24:27 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:24:28 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:24:28 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:24:28 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:24:43 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288443: status RUNNING/PENDING
2026-07-27 17:24:43 INFO None 5288450: status RUNNING/PENDING
2026-07-27 17:24:43 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450']. Waiting...
2026-07-27 17:24:58 INFO None 5288415: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:24:58 INFO None 5288443: status FINISHED
2026-07-27 17:24:58 INFO None 5288450: status FINISHED
2026-07-27 17:24:58 INFO Jobs still running: ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:25:13 INFO None 5288415: status FINISHED
2026-07-27 17:25:13 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:25:13 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288419: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288420: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288421: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288422: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288423: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288425: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:25:14 INFO None 5288443: status FINISHED
2026-07-27 17:25:16 INFO None 5288450: status FINISHED
2026-07-27 17:25:16 INFO Jobs still running: ['5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:25:31 INFO None 5288415: status FINISHED
2026-07-27 17:25:31 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288417: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288419: status FINISHED
2026-07-27 17:25:31 INFO None 5288420: status FINISHED
2026-07-27 17:25:31 INFO None 5288421: status FINISHED
2026-07-27 17:25:31 INFO None 5288422: status FINISHED
2026-07-27 17:25:31 INFO None 5288423: status FINISHED
2026-07-27 17:25:31 INFO None 5288425: status FINISHED
2026-07-27 17:25:31 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:25:31 INFO None 5288443: status FINISHED
2026-07-27 17:25:31 INFO None 5288450: status FINISHED
2026-07-27 17:25:31 INFO Jobs still running: ['5288416', '5288417', '5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:25:46 INFO None 5288415: status FINISHED
2026-07-27 17:25:46 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:25:46 INFO None 5288417: status FINISHED
2026-07-27 17:25:46 INFO None 5288419: status FINISHED
2026-07-27 17:25:46 INFO None 5288420: status FINISHED
2026-07-27 17:25:46 INFO None 5288421: status FINISHED
2026-07-27 17:25:46 INFO None 5288422: status FINISHED
2026-07-27 17:25:46 INFO None 5288423: status FINISHED
2026-07-27 17:25:46 INFO None 5288425: status FINISHED
2026-07-27 17:25:46 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:25:46 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:25:47 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:25:47 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:25:47 INFO None 5288443: status FINISHED
2026-07-27 17:25:47 INFO None 5288450: status FINISHED
2026-07-27 17:25:47 INFO Jobs still running: ['5288416', '5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:26:02 INFO None 5288415: status FINISHED
2026-07-27 17:26:02 INFO None 5288416: status RUNNING/PENDING
2026-07-27 17:26:02 INFO None 5288417: status FINISHED
2026-07-27 17:26:04 INFO None 5288419: status FINISHED
2026-07-27 17:26:04 INFO None 5288420: status FINISHED
2026-07-27 17:26:04 INFO None 5288421: status FINISHED
2026-07-27 17:26:04 INFO None 5288422: status FINISHED
2026-07-27 17:26:04 INFO None 5288423: status FINISHED
2026-07-27 17:26:04 INFO None 5288425: status FINISHED
2026-07-27 17:26:04 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:26:04 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:26:04 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:26:04 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:26:04 INFO None 5288443: status FINISHED
2026-07-27 17:26:04 INFO None 5288450: status FINISHED
2026-07-27 17:26:04 INFO Jobs still running: ['5288416', '5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:26:19 INFO None 5288415: status FINISHED
2026-07-27 17:26:19 INFO None 5288416: status FINISHED
2026-07-27 17:26:19 INFO None 5288417: status FINISHED
2026-07-27 17:26:19 INFO None 5288419: status FINISHED
2026-07-27 17:26:19 INFO None 5288420: status FINISHED
2026-07-27 17:26:19 INFO None 5288421: status FINISHED
2026-07-27 17:26:19 INFO None 5288422: status FINISHED
2026-07-27 17:26:19 INFO None 5288423: status FINISHED
2026-07-27 17:26:19 INFO None 5288425: status FINISHED
2026-07-27 17:26:19 INFO None 5288428: status RUNNING/PENDING
2026-07-27 17:26:19 INFO None 5288430: status RUNNING/PENDING
2026-07-27 17:26:21 INFO None 5288433: status RUNNING/PENDING
2026-07-27 17:26:21 INFO None 5288437: status RUNNING/PENDING
2026-07-27 17:26:21 INFO None 5288443: status FINISHED
2026-07-27 17:26:21 INFO None 5288450: status FINISHED
2026-07-27 17:26:21 INFO Jobs still running: ['5288428', '5288430', '5288433', '5288437']. Waiting...
2026-07-27 17:26:36 INFO None 5288415: status FINISHED
2026-07-27 17:26:37 INFO None 5288416: status FINISHED
2026-07-27 17:26:37 INFO None 5288417: status FINISHED
2026-07-27 17:26:37 INFO None 5288419: status FINISHED
2026-07-27 17:26:37 INFO None 5288420: status FINISHED
2026-07-27 17:26:37 INFO None 5288421: status FINISHED
2026-07-27 17:26:37 INFO None 5288422: status FINISHED
2026-07-27 17:26:37 INFO None 5288423: status FINISHED
2026-07-27 17:26:37 INFO None 5288425: status FINISHED
2026-07-27 17:26:37 INFO None 5288428: status FINISHED
2026-07-27 17:26:37 INFO None 5288430: status FINISHED
2026-07-27 17:26:37 INFO None 5288433: status FINISHED
2026-07-27 17:26:37 INFO None 5288437: status FINISHED
2026-07-27 17:26:37 INFO None 5288443: status FINISHED
2026-07-27 17:26:37 INFO None 5288450: status FINISHED
2026-07-27 17:26:37 INFO Jobs ['5288415', '5288416', '5288417', '5288419', '5288420', '5288421', '5288422', '5288423', '5288425', '5288428', '5288430', '5288433', '5288437', '5288443', '5288450'] have finished
2026-07-27 17:26:37 INFO Checking restart files were created ...
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc(1002685915 bytes)
2026-07-27 17:26:37 INFO  Run_model() completed successfully.
2026-07-27 17:26:37 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 09:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:26:37 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 11:00:00 days=153073 seconds=39600
2026-07-27 17:26:37 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 17:26:37 INFO [TIME] increment current_time 2020-02-07 09:00:00 -> 2020-02-07 11:00:00
2026-07-27 17:26:37 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:26:37 INFO ---------->>> Running process_satellite_data()
2026-07-27 17:26:37 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12017.nc
2026-07-27 17:26:37 INFO ---------->>> Running run_obs_converter()
2026-07-27 17:26:37 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-27 17:26:37 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_37907_153073.out
2026-07-27 17:26:37 INFO ---------->>> Running DART
2026-07-27 17:26:37 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 17:26:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out_toDART.nc
2026-07-27 17:26:37 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out_toDART.nc
2026-07-27 17:26:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out_toDART.nc
2026-07-27 17:26:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out_toDART.nc
2026-07-27 17:26:38 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out_toDART.nc
2026-07-27 17:26:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out_toDART.nc
2026-07-27 17:26:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out_toDART.nc
2026-07-27 17:26:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out_toDART.nc
2026-07-27 17:26:39 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out_toDART.nc
2026-07-27 17:26:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out_toDART.nc
2026-07-27 17:26:40 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out_toDART.nc
2026-07-27 17:26:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out_toDART.nc
2026-07-27 17:26:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out_toDART.nc
2026-07-27 17:26:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out_toDART.nc
2026-07-27 17:26:41 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020709_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out_toDART.nc
2026-07-27 17:26:42 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 17:26:42 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 17:26:42 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 17:26:42 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 17:26:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 17:26:42 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 17:27:04 INFO Found: []
2026-07-27 17:27:04 INFO No job id returned by command ./run_filter.bsh
2026-07-27 17:27:04 INFO No monitoring will be performed
2026-07-27 17:27:04 INFO Moving DART output files to analysis and preassim directories for date 2020020711 if present ...
2026-07-27 17:27:04 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:04 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020711'
2026-07-27 17:27:05 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 17:27:05 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 17:27:05 INFO run_dart() is DONE.
2026-07-27 17:27:05 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 17:27:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-27 17:27:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-27 17:27:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-27 17:27:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-27 17:27:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-27 17:27:31 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-27 17:27:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-27 17:27:42 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-27 17:27:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-27 17:27:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-27 17:27:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:27:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-27 17:28:05 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:28:05 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-27 17:28:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:28:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-27 17:28:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:28:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-27 17:28:21 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:28:21 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-27 17:28:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:28:26 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 17:28:26 INFO [TIME] step_end current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:28:26 INFO [TIME] step_start current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:28:26 INFO [TIME] window start=2020-02-07 11:00:00 end=2020-02-07 12:00:00 run_hours=1 has_assimilation=True
2026-07-27 17:28:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:27 INFO Hourly dataset computed and listing created
2026-07-27 17:28:30 INFO Hourly dataset computed
2026-07-27 17:28:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:31 INFO Hourly dataset computed and listing created
2026-07-27 17:28:32 INFO Hourly dataset computed
2026-07-27 17:28:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:33 INFO Hourly dataset computed and listing created
2026-07-27 17:28:33 INFO Hourly dataset computed
2026-07-27 17:28:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:34 INFO Hourly dataset computed and listing created
2026-07-27 17:28:35 INFO Hourly dataset computed
2026-07-27 17:28:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:36 INFO Hourly dataset computed and listing created
2026-07-27 17:28:36 INFO Hourly dataset computed
2026-07-27 17:28:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:37 INFO Hourly dataset computed and listing created
2026-07-27 17:28:38 INFO Hourly dataset computed
2026-07-27 17:28:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:39 INFO Hourly dataset computed and listing created
2026-07-27 17:28:39 INFO Hourly dataset computed
2026-07-27 17:28:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:40 INFO Hourly dataset computed and listing created
2026-07-27 17:28:41 INFO Hourly dataset computed
2026-07-27 17:28:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:42 INFO Hourly dataset computed and listing created
2026-07-27 17:28:42 INFO Hourly dataset computed
2026-07-27 17:28:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:43 INFO Hourly dataset computed and listing created
2026-07-27 17:28:44 INFO Hourly dataset computed
2026-07-27 17:28:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:45 INFO Hourly dataset computed and listing created
2026-07-27 17:28:45 INFO Hourly dataset computed
2026-07-27 17:28:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:46 INFO Hourly dataset computed and listing created
2026-07-27 17:28:47 INFO Hourly dataset computed
2026-07-27 17:28:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:48 INFO Hourly dataset computed and listing created
2026-07-27 17:28:48 INFO Hourly dataset computed
2026-07-27 17:28:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:49 INFO Hourly dataset computed and listing created
2026-07-27 17:28:50 INFO Hourly dataset computed
2026-07-27 17:28:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:28:51 INFO Hourly dataset computed and listing created
2026-07-27 17:28:51 INFO Hourly dataset computed
2026-07-27 17:28:51 INFO ---------->>> Running CHIMERE model from 2020-02-07 11:00:00 to 2020-02-07 12:00:00
2026-07-27 17:28:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:28:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 17:28:51 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020709_2_ENS1.nc
2026-07-27 17:28:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 17:28:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:28:51 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 17:28:51 INFO Queuing job for member 1...
2026-07-27 17:28:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:28:51 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 17:28:52 INFO Found: ['5288829']
2026-07-27 17:28:57 INFO [TGCC-IRENE] Submitted job with ID:['5288829']
2026-07-27 17:28:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:28:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 17:28:57 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020709_2_ENS2.nc
2026-07-27 17:28:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 17:28:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:28:57 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 17:28:57 INFO Queuing job for member 2...
2026-07-27 17:28:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:28:57 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 17:28:58 INFO Found: ['5288830']
2026-07-27 17:29:03 INFO [TGCC-IRENE] Submitted job with ID:['5288830']
2026-07-27 17:29:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 17:29:03 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020709_2_ENS3.nc
2026-07-27 17:29:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 17:29:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:03 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 17:29:03 INFO Queuing job for member 3...
2026-07-27 17:29:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:03 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 17:29:04 INFO Found: ['5288839']
2026-07-27 17:29:09 INFO [TGCC-IRENE] Submitted job with ID:['5288839']
2026-07-27 17:29:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 17:29:09 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020709_2_ENS4.nc
2026-07-27 17:29:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 17:29:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:09 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 17:29:09 INFO Queuing job for member 4...
2026-07-27 17:29:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:09 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 17:29:09 INFO Found: ['5288846']
2026-07-27 17:29:14 INFO [TGCC-IRENE] Submitted job with ID:['5288846']
2026-07-27 17:29:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 17:29:14 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020709_2_ENS5.nc
2026-07-27 17:29:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 17:29:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:14 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 17:29:14 INFO Queuing job for member 5...
2026-07-27 17:29:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:14 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 17:29:15 INFO Found: ['5288853']
2026-07-27 17:29:20 INFO [TGCC-IRENE] Submitted job with ID:['5288853']
2026-07-27 17:29:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 17:29:20 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020709_2_ENS6.nc
2026-07-27 17:29:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 17:29:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:20 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 17:29:20 INFO Queuing job for member 6...
2026-07-27 17:29:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:20 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 17:29:21 INFO Found: ['5288857']
2026-07-27 17:29:26 INFO [TGCC-IRENE] Submitted job with ID:['5288857']
2026-07-27 17:29:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 17:29:26 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020709_2_ENS7.nc
2026-07-27 17:29:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 17:29:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:26 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 17:29:26 INFO Queuing job for member 7...
2026-07-27 17:29:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:26 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 17:29:27 INFO Found: ['5288858']
2026-07-27 17:29:32 INFO [TGCC-IRENE] Submitted job with ID:['5288858']
2026-07-27 17:29:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 17:29:32 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020709_2_ENS8.nc
2026-07-27 17:29:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 17:29:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:32 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 17:29:32 INFO Queuing job for member 8...
2026-07-27 17:29:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:32 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 17:29:32 INFO Found: ['5288859']
2026-07-27 17:29:37 INFO [TGCC-IRENE] Submitted job with ID:['5288859']
2026-07-27 17:29:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 17:29:37 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020709_2_ENS9.nc
2026-07-27 17:29:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 17:29:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:37 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 17:29:37 INFO Queuing job for member 9...
2026-07-27 17:29:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:37 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 17:29:38 INFO Found: ['5288860']
2026-07-27 17:29:43 INFO [TGCC-IRENE] Submitted job with ID:['5288860']
2026-07-27 17:29:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 17:29:43 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020709_2_ENS10.nc
2026-07-27 17:29:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 17:29:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:43 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 17:29:43 INFO Queuing job for member 10...
2026-07-27 17:29:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:43 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 17:29:46 INFO Found: ['5288864']
2026-07-27 17:29:51 INFO [TGCC-IRENE] Submitted job with ID:['5288864']
2026-07-27 17:29:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 17:29:51 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020709_2_ENS11.nc
2026-07-27 17:29:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 17:29:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:51 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 17:29:51 INFO Queuing job for member 11...
2026-07-27 17:29:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:51 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 17:29:53 INFO Found: ['5288867']
2026-07-27 17:29:58 INFO [TGCC-IRENE] Submitted job with ID:['5288867']
2026-07-27 17:29:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:29:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 17:29:58 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020709_2_ENS12.nc
2026-07-27 17:29:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 17:29:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:29:58 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 17:29:58 INFO Queuing job for member 12...
2026-07-27 17:29:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:29:58 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 17:30:01 INFO Found: ['5288870']
2026-07-27 17:30:06 INFO [TGCC-IRENE] Submitted job with ID:['5288870']
2026-07-27 17:30:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:30:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 17:30:06 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020709_2_ENS13.nc
2026-07-27 17:30:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 17:30:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:30:06 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 17:30:06 INFO Queuing job for member 13...
2026-07-27 17:30:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:30:06 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 17:30:08 INFO Found: ['5288876']
2026-07-27 17:30:13 INFO [TGCC-IRENE] Submitted job with ID:['5288876']
2026-07-27 17:30:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:30:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 17:30:13 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020709_2_ENS14.nc
2026-07-27 17:30:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 17:30:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:30:13 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 17:30:13 INFO Queuing job for member 14...
2026-07-27 17:30:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:30:13 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 17:30:16 INFO Found: ['5288880']
2026-07-27 17:30:21 INFO [TGCC-IRENE] Submitted job with ID:['5288880']
2026-07-27 17:30:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:30:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 17:30:21 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020709_2_ENS15.nc
2026-07-27 17:30:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 17:30:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:30:21 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 17:30:21 INFO Queuing job for member 15...
2026-07-27 17:30:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:30:21 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 17:30:23 INFO Found: ['5288883']
2026-07-27 17:30:28 INFO [TGCC-IRENE] Submitted job with ID:['5288883']
2026-07-27 17:30:28 INFO Checking job status ...
2026-07-27 17:30:28 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:30:28 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:30:29 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:30:29 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:30:29 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:30:29 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:30:44 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:30:44 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:30:46 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:30:46 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:31:01 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:31:01 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:31:01 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:31:16 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:31:16 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:31:16 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:31:31 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:31:31 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:31:31 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:31:32 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:31:32 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:31:47 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:31:47 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:31:47 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:32:03 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:32:03 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:32:03 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:32:18 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:32:18 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:32:20 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:32:20 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:32:20 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:32:35 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:32:35 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:32:36 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:32:36 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:32:36 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:32:36 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:32:51 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:32:51 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:32:53 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:32:53 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:32:53 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:32:53 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:32:53 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:33:08 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:33:08 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:33:08 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:33:23 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288857: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288858: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288859: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:33:23 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:33:24 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:33:24 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:33:24 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:33:24 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:33:24 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:33:24 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:33:39 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288839: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288846: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288853: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288857: status FINISHED
2026-07-27 17:33:39 INFO None 5288858: status FINISHED
2026-07-27 17:33:39 INFO None 5288859: status FINISHED
2026-07-27 17:33:39 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:33:39 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:33:39 INFO Jobs still running: ['5288829', '5288830', '5288839', '5288846', '5288853', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:33:54 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288839: status FINISHED
2026-07-27 17:33:54 INFO None 5288846: status FINISHED
2026-07-27 17:33:54 INFO None 5288853: status FINISHED
2026-07-27 17:33:54 INFO None 5288857: status FINISHED
2026-07-27 17:33:54 INFO None 5288858: status FINISHED
2026-07-27 17:33:54 INFO None 5288859: status FINISHED
2026-07-27 17:33:54 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:33:54 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:33:54 INFO Jobs still running: ['5288829', '5288830', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:34:10 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288839: status FINISHED
2026-07-27 17:34:11 INFO None 5288846: status FINISHED
2026-07-27 17:34:11 INFO None 5288853: status FINISHED
2026-07-27 17:34:11 INFO None 5288857: status FINISHED
2026-07-27 17:34:11 INFO None 5288858: status FINISHED
2026-07-27 17:34:11 INFO None 5288859: status FINISHED
2026-07-27 17:34:11 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:34:11 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:34:11 INFO Jobs still running: ['5288829', '5288830', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:34:26 INFO None 5288829: status RUNNING/PENDING
2026-07-27 17:34:26 INFO None 5288830: status RUNNING/PENDING
2026-07-27 17:34:26 INFO None 5288839: status FINISHED
2026-07-27 17:34:26 INFO None 5288846: status FINISHED
2026-07-27 17:34:26 INFO None 5288853: status FINISHED
2026-07-27 17:34:26 INFO None 5288857: status FINISHED
2026-07-27 17:34:26 INFO None 5288858: status FINISHED
2026-07-27 17:34:26 INFO None 5288859: status FINISHED
2026-07-27 17:34:28 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:34:28 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:34:28 INFO Jobs still running: ['5288829', '5288830', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:34:43 INFO None 5288829: status FINISHED
2026-07-27 17:34:43 INFO None 5288830: status FINISHED
2026-07-27 17:34:43 INFO None 5288839: status FINISHED
2026-07-27 17:34:43 INFO None 5288846: status FINISHED
2026-07-27 17:34:43 INFO None 5288853: status FINISHED
2026-07-27 17:34:43 INFO None 5288857: status FINISHED
2026-07-27 17:34:43 INFO None 5288858: status FINISHED
2026-07-27 17:34:43 INFO None 5288859: status FINISHED
2026-07-27 17:34:43 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:34:43 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:34:43 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:34:58 INFO None 5288829: status FINISHED
2026-07-27 17:34:58 INFO None 5288830: status FINISHED
2026-07-27 17:34:58 INFO None 5288839: status FINISHED
2026-07-27 17:34:58 INFO None 5288846: status FINISHED
2026-07-27 17:34:59 INFO None 5288853: status FINISHED
2026-07-27 17:34:59 INFO None 5288857: status FINISHED
2026-07-27 17:34:59 INFO None 5288858: status FINISHED
2026-07-27 17:34:59 INFO None 5288859: status FINISHED
2026-07-27 17:34:59 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:34:59 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:34:59 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:34:59 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:35:01 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:35:01 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:35:01 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:35:01 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:35:16 INFO None 5288829: status FINISHED
2026-07-27 17:35:16 INFO None 5288830: status FINISHED
2026-07-27 17:35:16 INFO None 5288839: status FINISHED
2026-07-27 17:35:16 INFO None 5288846: status FINISHED
2026-07-27 17:35:16 INFO None 5288853: status FINISHED
2026-07-27 17:35:16 INFO None 5288857: status FINISHED
2026-07-27 17:35:16 INFO None 5288858: status FINISHED
2026-07-27 17:35:16 INFO None 5288859: status FINISHED
2026-07-27 17:35:16 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:35:16 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:35:16 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:35:31 INFO None 5288829: status FINISHED
2026-07-27 17:35:31 INFO None 5288830: status FINISHED
2026-07-27 17:35:31 INFO None 5288839: status FINISHED
2026-07-27 17:35:31 INFO None 5288846: status FINISHED
2026-07-27 17:35:31 INFO None 5288853: status FINISHED
2026-07-27 17:35:31 INFO None 5288857: status FINISHED
2026-07-27 17:35:31 INFO None 5288858: status FINISHED
2026-07-27 17:35:31 INFO None 5288859: status FINISHED
2026-07-27 17:35:31 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:35:31 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:35:33 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:35:33 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:35:33 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:35:33 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:35:33 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:35:33 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:35:48 INFO None 5288829: status FINISHED
2026-07-27 17:35:48 INFO None 5288830: status FINISHED
2026-07-27 17:35:48 INFO None 5288839: status FINISHED
2026-07-27 17:35:48 INFO None 5288846: status FINISHED
2026-07-27 17:35:49 INFO None 5288853: status FINISHED
2026-07-27 17:35:49 INFO None 5288857: status FINISHED
2026-07-27 17:35:49 INFO None 5288858: status FINISHED
2026-07-27 17:35:49 INFO None 5288859: status FINISHED
2026-07-27 17:35:49 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:35:49 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:35:49 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:36:04 INFO None 5288829: status FINISHED
2026-07-27 17:36:04 INFO None 5288830: status FINISHED
2026-07-27 17:36:04 INFO None 5288839: status FINISHED
2026-07-27 17:36:04 INFO None 5288846: status FINISHED
2026-07-27 17:36:04 INFO None 5288853: status FINISHED
2026-07-27 17:36:04 INFO None 5288857: status FINISHED
2026-07-27 17:36:04 INFO None 5288858: status FINISHED
2026-07-27 17:36:04 INFO None 5288859: status FINISHED
2026-07-27 17:36:04 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:36:04 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:36:04 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:36:19 INFO None 5288829: status FINISHED
2026-07-27 17:36:19 INFO None 5288830: status FINISHED
2026-07-27 17:36:19 INFO None 5288839: status FINISHED
2026-07-27 17:36:19 INFO None 5288846: status FINISHED
2026-07-27 17:36:19 INFO None 5288853: status FINISHED
2026-07-27 17:36:19 INFO None 5288857: status FINISHED
2026-07-27 17:36:19 INFO None 5288858: status FINISHED
2026-07-27 17:36:19 INFO None 5288859: status FINISHED
2026-07-27 17:36:19 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:36:19 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:36:19 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:36:35 INFO None 5288829: status FINISHED
2026-07-27 17:36:35 INFO None 5288830: status FINISHED
2026-07-27 17:36:35 INFO None 5288839: status FINISHED
2026-07-27 17:36:35 INFO None 5288846: status FINISHED
2026-07-27 17:36:35 INFO None 5288853: status FINISHED
2026-07-27 17:36:35 INFO None 5288857: status FINISHED
2026-07-27 17:36:35 INFO None 5288858: status FINISHED
2026-07-27 17:36:35 INFO None 5288859: status FINISHED
2026-07-27 17:36:35 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:36:35 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:36:35 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:36:50 INFO None 5288829: status FINISHED
2026-07-27 17:36:50 INFO None 5288830: status FINISHED
2026-07-27 17:36:50 INFO None 5288839: status FINISHED
2026-07-27 17:36:50 INFO None 5288846: status FINISHED
2026-07-27 17:36:50 INFO None 5288853: status FINISHED
2026-07-27 17:36:50 INFO None 5288857: status FINISHED
2026-07-27 17:36:50 INFO None 5288858: status FINISHED
2026-07-27 17:36:50 INFO None 5288859: status FINISHED
2026-07-27 17:36:50 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:36:51 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:36:51 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:37:06 INFO None 5288829: status FINISHED
2026-07-27 17:37:06 INFO None 5288830: status FINISHED
2026-07-27 17:37:06 INFO None 5288839: status FINISHED
2026-07-27 17:37:06 INFO None 5288846: status FINISHED
2026-07-27 17:37:08 INFO None 5288853: status FINISHED
2026-07-27 17:37:08 INFO None 5288857: status FINISHED
2026-07-27 17:37:08 INFO None 5288858: status FINISHED
2026-07-27 17:37:08 INFO None 5288859: status FINISHED
2026-07-27 17:37:08 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:37:08 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:37:08 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:37:23 INFO None 5288829: status FINISHED
2026-07-27 17:37:23 INFO None 5288830: status FINISHED
2026-07-27 17:37:23 INFO None 5288839: status FINISHED
2026-07-27 17:37:23 INFO None 5288846: status FINISHED
2026-07-27 17:37:23 INFO None 5288853: status FINISHED
2026-07-27 17:37:23 INFO None 5288857: status FINISHED
2026-07-27 17:37:23 INFO None 5288858: status FINISHED
2026-07-27 17:37:23 INFO None 5288859: status FINISHED
2026-07-27 17:37:23 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:37:23 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:37:23 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:37:38 INFO None 5288829: status FINISHED
2026-07-27 17:37:38 INFO None 5288830: status FINISHED
2026-07-27 17:37:38 INFO None 5288839: status FINISHED
2026-07-27 17:37:38 INFO None 5288846: status FINISHED
2026-07-27 17:37:40 INFO None 5288853: status FINISHED
2026-07-27 17:37:40 INFO None 5288857: status FINISHED
2026-07-27 17:37:40 INFO None 5288858: status FINISHED
2026-07-27 17:37:40 INFO None 5288859: status FINISHED
2026-07-27 17:37:40 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:37:40 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:37:40 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:37:41 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:37:41 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:37:41 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:37:41 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:37:41 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:37:56 INFO None 5288829: status FINISHED
2026-07-27 17:37:56 INFO None 5288830: status FINISHED
2026-07-27 17:37:56 INFO None 5288839: status FINISHED
2026-07-27 17:37:56 INFO None 5288846: status FINISHED
2026-07-27 17:37:56 INFO None 5288853: status FINISHED
2026-07-27 17:37:56 INFO None 5288857: status FINISHED
2026-07-27 17:37:56 INFO None 5288858: status FINISHED
2026-07-27 17:37:56 INFO None 5288859: status FINISHED
2026-07-27 17:37:56 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:37:56 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:37:56 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:38:11 INFO None 5288829: status FINISHED
2026-07-27 17:38:11 INFO None 5288830: status FINISHED
2026-07-27 17:38:11 INFO None 5288839: status FINISHED
2026-07-27 17:38:11 INFO None 5288846: status FINISHED
2026-07-27 17:38:11 INFO None 5288853: status FINISHED
2026-07-27 17:38:11 INFO None 5288857: status FINISHED
2026-07-27 17:38:11 INFO None 5288858: status FINISHED
2026-07-27 17:38:11 INFO None 5288859: status FINISHED
2026-07-27 17:38:11 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:38:11 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:38:11 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:38:27 INFO None 5288829: status FINISHED
2026-07-27 17:38:27 INFO None 5288830: status FINISHED
2026-07-27 17:38:27 INFO None 5288839: status FINISHED
2026-07-27 17:38:27 INFO None 5288846: status FINISHED
2026-07-27 17:38:27 INFO None 5288853: status FINISHED
2026-07-27 17:38:27 INFO None 5288857: status FINISHED
2026-07-27 17:38:27 INFO None 5288858: status FINISHED
2026-07-27 17:38:27 INFO None 5288859: status FINISHED
2026-07-27 17:38:27 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:38:27 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:38:27 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:38:43 INFO None 5288829: status FINISHED
2026-07-27 17:38:43 INFO None 5288830: status FINISHED
2026-07-27 17:38:43 INFO None 5288839: status FINISHED
2026-07-27 17:38:43 INFO None 5288846: status FINISHED
2026-07-27 17:38:43 INFO None 5288853: status FINISHED
2026-07-27 17:38:43 INFO None 5288857: status FINISHED
2026-07-27 17:38:43 INFO None 5288858: status FINISHED
2026-07-27 17:38:43 INFO None 5288859: status FINISHED
2026-07-27 17:38:43 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:38:43 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:38:43 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:38:58 INFO None 5288829: status FINISHED
2026-07-27 17:38:58 INFO None 5288830: status FINISHED
2026-07-27 17:38:58 INFO None 5288839: status FINISHED
2026-07-27 17:38:58 INFO None 5288846: status FINISHED
2026-07-27 17:38:58 INFO None 5288853: status FINISHED
2026-07-27 17:38:58 INFO None 5288857: status FINISHED
2026-07-27 17:38:58 INFO None 5288858: status FINISHED
2026-07-27 17:38:58 INFO None 5288859: status FINISHED
2026-07-27 17:38:58 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:38:58 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:38:58 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:38:58 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:38:58 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:38:59 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:39:01 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:39:01 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:39:16 INFO None 5288829: status FINISHED
2026-07-27 17:39:16 INFO None 5288830: status FINISHED
2026-07-27 17:39:16 INFO None 5288839: status FINISHED
2026-07-27 17:39:16 INFO None 5288846: status FINISHED
2026-07-27 17:39:16 INFO None 5288853: status FINISHED
2026-07-27 17:39:16 INFO None 5288857: status FINISHED
2026-07-27 17:39:16 INFO None 5288858: status FINISHED
2026-07-27 17:39:16 INFO None 5288859: status FINISHED
2026-07-27 17:39:16 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:39:16 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:39:16 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:39:31 INFO None 5288829: status FINISHED
2026-07-27 17:39:31 INFO None 5288830: status FINISHED
2026-07-27 17:39:31 INFO None 5288839: status FINISHED
2026-07-27 17:39:31 INFO None 5288846: status FINISHED
2026-07-27 17:39:31 INFO None 5288853: status FINISHED
2026-07-27 17:39:31 INFO None 5288857: status FINISHED
2026-07-27 17:39:31 INFO None 5288858: status FINISHED
2026-07-27 17:39:31 INFO None 5288859: status FINISHED
2026-07-27 17:39:31 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:39:31 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:39:31 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:39:33 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:39:33 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:39:33 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:39:33 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:39:33 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:39:48 INFO None 5288829: status FINISHED
2026-07-27 17:39:48 INFO None 5288830: status FINISHED
2026-07-27 17:39:48 INFO None 5288839: status FINISHED
2026-07-27 17:39:48 INFO None 5288846: status FINISHED
2026-07-27 17:39:48 INFO None 5288853: status FINISHED
2026-07-27 17:39:48 INFO None 5288857: status FINISHED
2026-07-27 17:39:48 INFO None 5288858: status FINISHED
2026-07-27 17:39:48 INFO None 5288859: status FINISHED
2026-07-27 17:39:48 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:39:48 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:39:49 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:39:49 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:39:49 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:39:49 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:39:49 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:39:49 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:40:04 INFO None 5288829: status FINISHED
2026-07-27 17:40:04 INFO None 5288830: status FINISHED
2026-07-27 17:40:04 INFO None 5288839: status FINISHED
2026-07-27 17:40:04 INFO None 5288846: status FINISHED
2026-07-27 17:40:04 INFO None 5288853: status FINISHED
2026-07-27 17:40:04 INFO None 5288857: status FINISHED
2026-07-27 17:40:04 INFO None 5288858: status FINISHED
2026-07-27 17:40:04 INFO None 5288859: status FINISHED
2026-07-27 17:40:06 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:40:06 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:40:06 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:40:21 INFO None 5288829: status FINISHED
2026-07-27 17:40:21 INFO None 5288830: status FINISHED
2026-07-27 17:40:21 INFO None 5288839: status FINISHED
2026-07-27 17:40:21 INFO None 5288846: status FINISHED
2026-07-27 17:40:21 INFO None 5288853: status FINISHED
2026-07-27 17:40:21 INFO None 5288857: status FINISHED
2026-07-27 17:40:21 INFO None 5288858: status FINISHED
2026-07-27 17:40:21 INFO None 5288859: status FINISHED
2026-07-27 17:40:21 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:40:21 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:40:21 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:40:36 INFO None 5288829: status FINISHED
2026-07-27 17:40:36 INFO None 5288830: status FINISHED
2026-07-27 17:40:36 INFO None 5288839: status FINISHED
2026-07-27 17:40:36 INFO None 5288846: status FINISHED
2026-07-27 17:40:36 INFO None 5288853: status FINISHED
2026-07-27 17:40:36 INFO None 5288857: status FINISHED
2026-07-27 17:40:36 INFO None 5288858: status FINISHED
2026-07-27 17:40:36 INFO None 5288859: status FINISHED
2026-07-27 17:40:36 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:40:36 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:40:36 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:40:53 INFO None 5288829: status FINISHED
2026-07-27 17:40:53 INFO None 5288830: status FINISHED
2026-07-27 17:40:53 INFO None 5288839: status FINISHED
2026-07-27 17:40:53 INFO None 5288846: status FINISHED
2026-07-27 17:40:53 INFO None 5288853: status FINISHED
2026-07-27 17:40:53 INFO None 5288857: status FINISHED
2026-07-27 17:40:53 INFO None 5288858: status FINISHED
2026-07-27 17:40:53 INFO None 5288859: status FINISHED
2026-07-27 17:40:53 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:40:53 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:40:53 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:41:08 INFO None 5288829: status FINISHED
2026-07-27 17:41:08 INFO None 5288830: status FINISHED
2026-07-27 17:41:08 INFO None 5288839: status FINISHED
2026-07-27 17:41:08 INFO None 5288846: status FINISHED
2026-07-27 17:41:08 INFO None 5288853: status FINISHED
2026-07-27 17:41:08 INFO None 5288857: status FINISHED
2026-07-27 17:41:08 INFO None 5288858: status FINISHED
2026-07-27 17:41:08 INFO None 5288859: status FINISHED
2026-07-27 17:41:08 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:41:08 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:41:08 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:41:08 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:41:08 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:41:10 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:41:10 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:41:10 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:41:25 INFO None 5288829: status FINISHED
2026-07-27 17:41:25 INFO None 5288830: status FINISHED
2026-07-27 17:41:25 INFO None 5288839: status FINISHED
2026-07-27 17:41:25 INFO None 5288846: status FINISHED
2026-07-27 17:41:25 INFO None 5288853: status FINISHED
2026-07-27 17:41:25 INFO None 5288857: status FINISHED
2026-07-27 17:41:25 INFO None 5288858: status FINISHED
2026-07-27 17:41:25 INFO None 5288859: status FINISHED
2026-07-27 17:41:25 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:41:25 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:41:25 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:41:40 INFO None 5288829: status FINISHED
2026-07-27 17:41:40 INFO None 5288830: status FINISHED
2026-07-27 17:41:41 INFO None 5288839: status FINISHED
2026-07-27 17:41:41 INFO None 5288846: status FINISHED
2026-07-27 17:41:41 INFO None 5288853: status FINISHED
2026-07-27 17:41:41 INFO None 5288857: status FINISHED
2026-07-27 17:41:41 INFO None 5288858: status FINISHED
2026-07-27 17:41:41 INFO None 5288859: status FINISHED
2026-07-27 17:41:41 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:41:41 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:41:41 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:41:41 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:41:41 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:41:41 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:41:43 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:41:43 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:41:58 INFO None 5288829: status FINISHED
2026-07-27 17:41:58 INFO None 5288830: status FINISHED
2026-07-27 17:41:58 INFO None 5288839: status FINISHED
2026-07-27 17:41:58 INFO None 5288846: status FINISHED
2026-07-27 17:41:58 INFO None 5288853: status FINISHED
2026-07-27 17:41:58 INFO None 5288857: status FINISHED
2026-07-27 17:41:58 INFO None 5288858: status FINISHED
2026-07-27 17:41:58 INFO None 5288859: status FINISHED
2026-07-27 17:41:58 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:41:58 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:41:58 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:41:58 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:41:58 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:42:00 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:42:00 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:42:00 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:42:15 INFO None 5288829: status FINISHED
2026-07-27 17:42:15 INFO None 5288830: status FINISHED
2026-07-27 17:42:16 INFO None 5288839: status FINISHED
2026-07-27 17:42:16 INFO None 5288846: status FINISHED
2026-07-27 17:42:16 INFO None 5288853: status FINISHED
2026-07-27 17:42:16 INFO None 5288857: status FINISHED
2026-07-27 17:42:16 INFO None 5288858: status FINISHED
2026-07-27 17:42:16 INFO None 5288859: status FINISHED
2026-07-27 17:42:16 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:42:16 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:42:16 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:42:31 INFO None 5288829: status FINISHED
2026-07-27 17:42:31 INFO None 5288830: status FINISHED
2026-07-27 17:42:31 INFO None 5288839: status FINISHED
2026-07-27 17:42:31 INFO None 5288846: status FINISHED
2026-07-27 17:42:31 INFO None 5288853: status FINISHED
2026-07-27 17:42:31 INFO None 5288857: status FINISHED
2026-07-27 17:42:31 INFO None 5288858: status FINISHED
2026-07-27 17:42:31 INFO None 5288859: status FINISHED
2026-07-27 17:42:31 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:42:31 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:42:31 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:42:46 INFO None 5288829: status FINISHED
2026-07-27 17:42:46 INFO None 5288830: status FINISHED
2026-07-27 17:42:46 INFO None 5288839: status FINISHED
2026-07-27 17:42:46 INFO None 5288846: status FINISHED
2026-07-27 17:42:46 INFO None 5288853: status FINISHED
2026-07-27 17:42:46 INFO None 5288857: status FINISHED
2026-07-27 17:42:46 INFO None 5288858: status FINISHED
2026-07-27 17:42:46 INFO None 5288859: status FINISHED
2026-07-27 17:42:46 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:42:46 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:42:46 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:43:01 INFO None 5288829: status FINISHED
2026-07-27 17:43:01 INFO None 5288830: status FINISHED
2026-07-27 17:43:01 INFO None 5288839: status FINISHED
2026-07-27 17:43:01 INFO None 5288846: status FINISHED
2026-07-27 17:43:01 INFO None 5288853: status FINISHED
2026-07-27 17:43:01 INFO None 5288857: status FINISHED
2026-07-27 17:43:01 INFO None 5288858: status FINISHED
2026-07-27 17:43:01 INFO None 5288859: status FINISHED
2026-07-27 17:43:02 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:43:02 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:43:02 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:43:18 INFO None 5288829: status FINISHED
2026-07-27 17:43:18 INFO None 5288830: status FINISHED
2026-07-27 17:43:18 INFO None 5288839: status FINISHED
2026-07-27 17:43:18 INFO None 5288846: status FINISHED
2026-07-27 17:43:18 INFO None 5288853: status FINISHED
2026-07-27 17:43:18 INFO None 5288857: status FINISHED
2026-07-27 17:43:18 INFO None 5288858: status FINISHED
2026-07-27 17:43:18 INFO None 5288859: status FINISHED
2026-07-27 17:43:18 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:43:18 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:43:18 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:43:34 INFO None 5288829: status FINISHED
2026-07-27 17:43:34 INFO None 5288830: status FINISHED
2026-07-27 17:43:36 INFO None 5288839: status FINISHED
2026-07-27 17:43:36 INFO None 5288846: status FINISHED
2026-07-27 17:43:36 INFO None 5288853: status FINISHED
2026-07-27 17:43:36 INFO None 5288857: status FINISHED
2026-07-27 17:43:36 INFO None 5288858: status FINISHED
2026-07-27 17:43:36 INFO None 5288859: status FINISHED
2026-07-27 17:43:36 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:43:36 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:43:36 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:43:51 INFO None 5288829: status FINISHED
2026-07-27 17:43:51 INFO None 5288830: status FINISHED
2026-07-27 17:43:51 INFO None 5288839: status FINISHED
2026-07-27 17:43:51 INFO None 5288846: status FINISHED
2026-07-27 17:43:51 INFO None 5288853: status FINISHED
2026-07-27 17:43:51 INFO None 5288857: status FINISHED
2026-07-27 17:43:51 INFO None 5288858: status FINISHED
2026-07-27 17:43:51 INFO None 5288859: status FINISHED
2026-07-27 17:43:51 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:43:51 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:43:51 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:43:51 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:43:51 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:43:51 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:43:53 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:43:53 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:44:08 INFO None 5288829: status FINISHED
2026-07-27 17:44:08 INFO None 5288830: status FINISHED
2026-07-27 17:44:08 INFO None 5288839: status FINISHED
2026-07-27 17:44:08 INFO None 5288846: status FINISHED
2026-07-27 17:44:08 INFO None 5288853: status FINISHED
2026-07-27 17:44:09 INFO None 5288857: status FINISHED
2026-07-27 17:44:09 INFO None 5288858: status FINISHED
2026-07-27 17:44:09 INFO None 5288859: status FINISHED
2026-07-27 17:44:09 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:44:09 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:44:09 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:44:09 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:44:09 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:44:09 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:44:11 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:44:11 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:44:26 INFO None 5288829: status FINISHED
2026-07-27 17:44:26 INFO None 5288830: status FINISHED
2026-07-27 17:44:26 INFO None 5288839: status FINISHED
2026-07-27 17:44:26 INFO None 5288846: status FINISHED
2026-07-27 17:44:26 INFO None 5288853: status FINISHED
2026-07-27 17:44:26 INFO None 5288857: status FINISHED
2026-07-27 17:44:26 INFO None 5288858: status FINISHED
2026-07-27 17:44:26 INFO None 5288859: status FINISHED
2026-07-27 17:44:26 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:44:26 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:44:26 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:44:41 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:44:41 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:44:41 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:44:42 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:44:42 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:44:42 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:44:42 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:44:42 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:44:42 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:44:57 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:44:57 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:44:57 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:44:57 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:45:12 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:45:12 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:45:12 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:45:12 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:45:28 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:45:28 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:45:28 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:45:28 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:45:43 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:45:43 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:45:43 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:45:43 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:45:58 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:45:58 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:45:58 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:45:58 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:45:58 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:45:58 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:46:00 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:46:00 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:46:00 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:46:00 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:46:00 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:46:15 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:46:16 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:46:16 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:46:16 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:46:31 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:46:31 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:46:31 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:46:31 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:46:31 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:46:31 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:46:33 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:46:33 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:46:33 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288867: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288870: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:46:33 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:46:33 INFO Jobs still running: ['5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:46:48 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:46:48 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:46:48 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:46:48 INFO None 5288867: status FINISHED
2026-07-27 17:46:48 INFO None 5288870: status FINISHED
2026-07-27 17:46:48 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:46:50 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:46:50 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:46:50 INFO Jobs still running: ['5288860', '5288864', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:47:06 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:47:06 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:47:06 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:47:06 INFO None 5288867: status FINISHED
2026-07-27 17:47:06 INFO None 5288870: status FINISHED
2026-07-27 17:47:06 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:47:06 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:47:06 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:47:06 INFO Jobs still running: ['5288860', '5288864', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:47:21 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:47:21 INFO None 5288860: status RUNNING/PENDING
2026-07-27 17:47:21 INFO None 5288864: status RUNNING/PENDING
2026-07-27 17:47:21 INFO None 5288867: status FINISHED
2026-07-27 17:47:21 INFO None 5288870: status FINISHED
2026-07-27 17:47:21 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:47:21 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:47:21 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:47:21 INFO Jobs still running: ['5288860', '5288864', '5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:47:38 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:47:38 INFO None 5288860: status FINISHED
2026-07-27 17:47:38 INFO None 5288864: status FINISHED
2026-07-27 17:47:38 INFO None 5288867: status FINISHED
2026-07-27 17:47:38 INFO None 5288870: status FINISHED
2026-07-27 17:47:38 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:47:38 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:47:38 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:47:38 INFO Jobs still running: ['5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:47:53 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:47:53 INFO None 5288860: status FINISHED
2026-07-27 17:47:53 INFO None 5288864: status FINISHED
2026-07-27 17:47:53 INFO None 5288867: status FINISHED
2026-07-27 17:47:53 INFO None 5288870: status FINISHED
2026-07-27 17:47:54 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:47:54 INFO None 5288880: status RUNNING/PENDING
2026-07-27 17:47:54 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:47:54 INFO Jobs still running: ['5288876', '5288880', '5288883']. Waiting...
2026-07-27 17:48:09 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:48:09 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:48:09 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:48:09 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:48:11 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:48:11 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:48:11 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:48:11 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:48:11 INFO None 5288860: status FINISHED
2026-07-27 17:48:11 INFO None 5288864: status FINISHED
2026-07-27 17:48:11 INFO None 5288867: status FINISHED
2026-07-27 17:48:11 INFO None 5288870: status FINISHED
2026-07-27 17:48:11 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:48:11 INFO None 5288880: status FINISHED
2026-07-27 17:48:11 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:48:11 INFO Jobs still running: ['5288876', '5288883']. Waiting...
2026-07-27 17:48:26 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:48:26 INFO None 5288860: status FINISHED
2026-07-27 17:48:26 INFO None 5288864: status FINISHED
2026-07-27 17:48:26 INFO None 5288867: status FINISHED
2026-07-27 17:48:26 INFO None 5288870: status FINISHED
2026-07-27 17:48:26 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:48:26 INFO None 5288880: status FINISHED
2026-07-27 17:48:26 INFO None 5288883: status RUNNING/PENDING
2026-07-27 17:48:26 INFO Jobs still running: ['5288876', '5288883']. Waiting...
2026-07-27 17:48:41 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:48:41 INFO None 5288860: status FINISHED
2026-07-27 17:48:41 INFO None 5288864: status FINISHED
2026-07-27 17:48:44 INFO None 5288867: status FINISHED
2026-07-27 17:48:44 INFO None 5288870: status FINISHED
2026-07-27 17:48:44 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:48:44 INFO None 5288880: status FINISHED
2026-07-27 17:48:44 INFO None 5288883: status FINISHED
2026-07-27 17:48:44 INFO Jobs still running: ['5288876']. Waiting...
2026-07-27 17:48:59 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:48:59 INFO None 5288860: status FINISHED
2026-07-27 17:48:59 INFO None 5288864: status FINISHED
2026-07-27 17:48:59 INFO None 5288867: status FINISHED
2026-07-27 17:48:59 INFO None 5288870: status FINISHED
2026-07-27 17:48:59 INFO None 5288876: status RUNNING/PENDING
2026-07-27 17:49:01 INFO None 5288880: status FINISHED
2026-07-27 17:49:01 INFO None 5288883: status FINISHED
2026-07-27 17:49:01 INFO Jobs still running: ['5288876']. Waiting...
2026-07-27 17:49:16 INFO None 5288829: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288830: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288839: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288846: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288853: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288857: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288858: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288859: status FINISHED (not in squeue)
2026-07-27 17:49:16 INFO None 5288860: status FINISHED
2026-07-27 17:49:16 INFO None 5288864: status FINISHED
2026-07-27 17:49:16 INFO None 5288867: status FINISHED
2026-07-27 17:49:16 INFO None 5288870: status FINISHED
2026-07-27 17:49:17 INFO None 5288876: status FINISHED
2026-07-27 17:49:17 INFO None 5288880: status FINISHED
2026-07-27 17:49:17 INFO None 5288883: status FINISHED
2026-07-27 17:49:17 INFO Jobs ['5288829', '5288830', '5288839', '5288846', '5288853', '5288857', '5288858', '5288859', '5288860', '5288864', '5288867', '5288870', '5288876', '5288880', '5288883'] have finished
2026-07-27 17:49:17 INFO Checking restart files were created ...
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc(668832435 bytes)
2026-07-27 17:49:17 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc(668832435 bytes)
2026-07-27 17:49:17 INFO  Run_model() completed successfully.
2026-07-27 17:49:17 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 11:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:49:17 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 12:00:00 days=153073 seconds=43200
2026-07-27 17:49:17 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 17:49:17 INFO [TIME] increment current_time 2020-02-07 11:00:00 -> 2020-02-07 12:00:00
2026-07-27 17:49:17 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:49:17 INFO ---------->>> Running process_satellite_data()
2026-07-27 17:49:17 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12018.nc
2026-07-27 17:49:17 INFO ---------->>> Running run_obs_converter()
2026-07-27 17:49:17 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-27 17:49:17 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_43993_153073.out
2026-07-27 17:49:17 INFO ---------->>> Running DART
2026-07-27 17:49:17 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 17:49:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_1_out_toDART.nc
2026-07-27 17:49:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_1_out_toDART.nc
2026-07-27 17:49:17 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_1_out_toDART.nc
2026-07-27 17:49:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_1_out_toDART.nc
2026-07-27 17:49:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_1_out_toDART.nc
2026-07-27 17:49:18 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_1_out_toDART.nc
2026-07-27 17:49:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_1_out_toDART.nc
2026-07-27 17:49:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_1_out_toDART.nc
2026-07-27 17:49:19 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_1_out_toDART.nc
2026-07-27 17:49:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_1_out_toDART.nc
2026-07-27 17:49:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_1_out_toDART.nc
2026-07-27 17:49:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_1_out_toDART.nc
2026-07-27 17:49:20 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_1_out_toDART.nc
2026-07-27 17:49:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_1_out_toDART.nc
2026-07-27 17:49:21 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020711_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_1_out_toDART.nc
2026-07-27 17:49:21 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 17:49:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 17:49:21 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 17:49:21 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 17:49:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 17:49:21 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 17:49:54 INFO Found: []
2026-07-27 17:49:54 INFO No job id returned by command ./run_filter.bsh
2026-07-27 17:49:54 INFO No monitoring will be performed
2026-07-27 17:49:54 INFO Moving DART output files to analysis and preassim directories for date 2020020712 if present ...
2026-07-27 17:49:54 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:54 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:55 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020712'
2026-07-27 17:49:55 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020712'
2026-07-27 17:49:55 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 17:49:55 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 17:49:55 INFO run_dart() is DONE.
2026-07-27 17:49:55 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 17:49:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-27 17:49:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:49:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-27 17:50:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-27 17:50:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-27 17:50:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-27 17:50:14 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:14 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-27 17:50:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:18 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-27 17:50:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-27 17:50:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-27 17:50:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-27 17:50:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-27 17:50:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:37 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-27 17:50:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-27 17:50:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-27 17:50:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-27 17:50:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 17:50:53 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 17:50:53 INFO [TIME] step_end current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:50:53 INFO [TIME] step_start current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 12:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 17:50:53 INFO [TIME] window start=2020-02-07 12:00:00 end=2020-02-07 14:00:00 run_hours=2 has_assimilation=True
2026-07-27 17:50:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:50:55 INFO Hourly dataset computed and listing created
2026-07-27 17:51:06 INFO Hourly dataset computed
2026-07-27 17:51:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:07 INFO Hourly dataset computed and listing created
2026-07-27 17:51:08 INFO Hourly dataset computed
2026-07-27 17:51:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:09 INFO Hourly dataset computed and listing created
2026-07-27 17:51:10 INFO Hourly dataset computed
2026-07-27 17:51:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:10 INFO Hourly dataset computed and listing created
2026-07-27 17:51:11 INFO Hourly dataset computed
2026-07-27 17:51:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:12 INFO Hourly dataset computed and listing created
2026-07-27 17:51:13 INFO Hourly dataset computed
2026-07-27 17:51:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:14 INFO Hourly dataset computed and listing created
2026-07-27 17:51:16 INFO Hourly dataset computed
2026-07-27 17:51:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:17 INFO Hourly dataset computed and listing created
2026-07-27 17:51:18 INFO Hourly dataset computed
2026-07-27 17:51:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:19 INFO Hourly dataset computed and listing created
2026-07-27 17:51:19 INFO Hourly dataset computed
2026-07-27 17:51:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:20 INFO Hourly dataset computed and listing created
2026-07-27 17:51:21 INFO Hourly dataset computed
2026-07-27 17:51:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:22 INFO Hourly dataset computed and listing created
2026-07-27 17:51:23 INFO Hourly dataset computed
2026-07-27 17:51:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:24 INFO Hourly dataset computed and listing created
2026-07-27 17:51:25 INFO Hourly dataset computed
2026-07-27 17:51:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:26 INFO Hourly dataset computed and listing created
2026-07-27 17:51:27 INFO Hourly dataset computed
2026-07-27 17:51:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:27 INFO Hourly dataset computed and listing created
2026-07-27 17:51:28 INFO Hourly dataset computed
2026-07-27 17:51:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:29 INFO Hourly dataset computed and listing created
2026-07-27 17:51:30 INFO Hourly dataset computed
2026-07-27 17:51:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 17:51:31 INFO Hourly dataset computed and listing created
2026-07-27 17:51:32 INFO Hourly dataset computed
2026-07-27 17:51:32 INFO ---------->>> Running CHIMERE model from 2020-02-07 12:00:00 to 2020-02-07 14:00:00
2026-07-27 17:51:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:51:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 17:51:32 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020711_1_ENS1.nc
2026-07-27 17:51:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 17:51:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:51:32 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 17:51:32 INFO Queuing job for member 1...
2026-07-27 17:51:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:51:32 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 17:51:33 INFO Found: ['5289275']
2026-07-27 17:51:38 INFO [TGCC-IRENE] Submitted job with ID:['5289275']
2026-07-27 17:51:38 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:51:38 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 17:51:38 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020711_1_ENS2.nc
2026-07-27 17:51:38 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 17:51:38 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:51:38 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 17:51:38 INFO Queuing job for member 2...
2026-07-27 17:51:38 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:51:38 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 17:51:38 INFO Found: ['5289276']
2026-07-27 17:51:43 INFO [TGCC-IRENE] Submitted job with ID:['5289276']
2026-07-27 17:51:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:51:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 17:51:43 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020711_1_ENS3.nc
2026-07-27 17:51:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 17:51:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:51:43 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 17:51:43 INFO Queuing job for member 3...
2026-07-27 17:51:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:51:43 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 17:51:44 INFO Found: ['5289278']
2026-07-27 17:51:49 INFO [TGCC-IRENE] Submitted job with ID:['5289278']
2026-07-27 17:51:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:51:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 17:51:49 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020711_1_ENS4.nc
2026-07-27 17:51:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 17:51:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:51:49 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 17:51:49 INFO Queuing job for member 4...
2026-07-27 17:51:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:51:49 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 17:51:50 INFO Found: ['5289279']
2026-07-27 17:51:55 INFO [TGCC-IRENE] Submitted job with ID:['5289279']
2026-07-27 17:51:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:51:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 17:51:55 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020711_1_ENS5.nc
2026-07-27 17:51:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 17:51:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:51:55 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 17:51:55 INFO Queuing job for member 5...
2026-07-27 17:51:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:51:55 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 17:51:56 INFO Found: ['5289280']
2026-07-27 17:52:01 INFO [TGCC-IRENE] Submitted job with ID:['5289280']
2026-07-27 17:52:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 17:52:01 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020711_1_ENS6.nc
2026-07-27 17:52:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 17:52:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:01 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 17:52:01 INFO Queuing job for member 6...
2026-07-27 17:52:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:01 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 17:52:01 INFO Found: ['5289282']
2026-07-27 17:52:06 INFO [TGCC-IRENE] Submitted job with ID:['5289282']
2026-07-27 17:52:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 17:52:06 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020711_1_ENS7.nc
2026-07-27 17:52:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 17:52:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:06 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 17:52:06 INFO Queuing job for member 7...
2026-07-27 17:52:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:06 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 17:52:08 INFO Found: ['5289286']
2026-07-27 17:52:13 INFO [TGCC-IRENE] Submitted job with ID:['5289286']
2026-07-27 17:52:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 17:52:13 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020711_1_ENS8.nc
2026-07-27 17:52:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 17:52:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:13 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 17:52:13 INFO Queuing job for member 8...
2026-07-27 17:52:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:13 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 17:52:16 INFO Found: ['5289288']
2026-07-27 17:52:21 INFO [TGCC-IRENE] Submitted job with ID:['5289288']
2026-07-27 17:52:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 17:52:21 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020711_1_ENS9.nc
2026-07-27 17:52:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 17:52:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:21 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 17:52:21 INFO Queuing job for member 9...
2026-07-27 17:52:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:21 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 17:52:23 INFO Found: ['5289290']
2026-07-27 17:52:28 INFO [TGCC-IRENE] Submitted job with ID:['5289290']
2026-07-27 17:52:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 17:52:28 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020711_1_ENS10.nc
2026-07-27 17:52:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 17:52:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:28 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 17:52:28 INFO Queuing job for member 10...
2026-07-27 17:52:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:28 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 17:52:31 INFO Found: ['5289292']
2026-07-27 17:52:36 INFO [TGCC-IRENE] Submitted job with ID:['5289292']
2026-07-27 17:52:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 17:52:36 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020711_1_ENS11.nc
2026-07-27 17:52:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 17:52:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:36 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 17:52:36 INFO Queuing job for member 11...
2026-07-27 17:52:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:36 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 17:52:38 INFO Found: ['5289293']
2026-07-27 17:52:43 INFO [TGCC-IRENE] Submitted job with ID:['5289293']
2026-07-27 17:52:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 17:52:43 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020711_1_ENS12.nc
2026-07-27 17:52:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 17:52:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:43 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 17:52:43 INFO Queuing job for member 12...
2026-07-27 17:52:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:43 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 17:52:46 INFO Found: ['5289294']
2026-07-27 17:52:51 INFO [TGCC-IRENE] Submitted job with ID:['5289294']
2026-07-27 17:52:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 17:52:51 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020711_1_ENS13.nc
2026-07-27 17:52:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 17:52:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:51 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 17:52:51 INFO Queuing job for member 13...
2026-07-27 17:52:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:51 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 17:52:53 INFO Found: ['5289295']
2026-07-27 17:52:58 INFO [TGCC-IRENE] Submitted job with ID:['5289295']
2026-07-27 17:52:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:52:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 17:52:58 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020711_1_ENS14.nc
2026-07-27 17:52:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 17:52:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:52:58 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 17:52:58 INFO Queuing job for member 14...
2026-07-27 17:52:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:52:58 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 17:53:01 INFO Found: ['5289296']
2026-07-27 17:53:06 INFO [TGCC-IRENE] Submitted job with ID:['5289296']
2026-07-27 17:53:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 17:53:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 17:53:06 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020711_1_ENS15.nc
2026-07-27 17:53:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 17:53:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 17:53:06 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 17:53:06 INFO Queuing job for member 15...
2026-07-27 17:53:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 17:53:06 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 17:53:08 INFO Found: ['5289299']
2026-07-27 17:53:13 INFO [TGCC-IRENE] Submitted job with ID:['5289299']
2026-07-27 17:53:13 INFO Checking job status ...
2026-07-27 17:53:13 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:53:13 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:53:13 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:53:14 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:53:14 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:53:29 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:53:29 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:53:31 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:53:31 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:53:46 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:53:46 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:53:46 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:54:01 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:54:01 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:54:02 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:54:02 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:54:17 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:54:17 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:54:17 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:54:33 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:54:33 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:54:33 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:54:48 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:54:48 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:54:50 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:54:50 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:55:05 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:55:05 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:55:06 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:55:06 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:55:21 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:55:21 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:55:23 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:55:23 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:55:38 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:55:38 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:55:40 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:55:40 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:55:55 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:55:55 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:55:56 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:55:56 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:56:11 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289276: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:56:11 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:56:11 INFO Jobs still running: ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:56:26 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289276: status FINISHED
2026-07-27 17:56:26 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:56:26 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:56:26 INFO Jobs still running: ['5289275', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:56:43 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289276: status FINISHED
2026-07-27 17:56:43 INFO None 5289278: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289279: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:56:43 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:56:43 INFO Jobs still running: ['5289275', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:56:58 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:56:58 INFO None 5289276: status FINISHED
2026-07-27 17:56:58 INFO None 5289278: status FINISHED
2026-07-27 17:56:58 INFO None 5289279: status FINISHED
2026-07-27 17:56:59 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:56:59 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:57:01 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:57:01 INFO Jobs still running: ['5289275', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:57:16 INFO None 5289275: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289276: status FINISHED
2026-07-27 17:57:16 INFO None 5289278: status FINISHED
2026-07-27 17:57:16 INFO None 5289279: status FINISHED
2026-07-27 17:57:16 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:57:16 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:57:16 INFO Jobs still running: ['5289275', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:57:31 INFO None 5289275: status FINISHED
2026-07-27 17:57:31 INFO None 5289276: status FINISHED
2026-07-27 17:57:31 INFO None 5289278: status FINISHED
2026-07-27 17:57:31 INFO None 5289279: status FINISHED
2026-07-27 17:57:31 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:57:31 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:57:31 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:57:31 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:57:31 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:57:33 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:57:33 INFO Jobs still running: ['5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:57:48 INFO None 5289275: status FINISHED
2026-07-27 17:57:48 INFO None 5289276: status FINISHED
2026-07-27 17:57:48 INFO None 5289278: status FINISHED
2026-07-27 17:57:48 INFO None 5289279: status FINISHED
2026-07-27 17:57:48 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:57:48 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:57:48 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:57:48 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:57:49 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:57:49 INFO Jobs still running: ['5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:58:04 INFO None 5289275: status FINISHED
2026-07-27 17:58:04 INFO None 5289276: status FINISHED
2026-07-27 17:58:04 INFO None 5289278: status FINISHED
2026-07-27 17:58:04 INFO None 5289279: status FINISHED
2026-07-27 17:58:04 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:58:04 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:58:04 INFO None 5289286: status RUNNING/PENDING
2026-07-27 17:58:04 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:58:04 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:58:04 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:58:06 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:58:06 INFO None 5289294: status RUNNING/PENDING
2026-07-27 17:58:06 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:58:06 INFO None 5289296: status RUNNING/PENDING
2026-07-27 17:58:06 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:58:06 INFO Jobs still running: ['5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299']. Waiting...
2026-07-27 17:58:21 INFO None 5289275: status FINISHED
2026-07-27 17:58:21 INFO None 5289276: status FINISHED
2026-07-27 17:58:21 INFO None 5289278: status FINISHED
2026-07-27 17:58:21 INFO None 5289279: status FINISHED
2026-07-27 17:58:21 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289286: status FINISHED
2026-07-27 17:58:21 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289292: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289293: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289294: status FINISHED
2026-07-27 17:58:21 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:58:21 INFO None 5289296: status FINISHED
2026-07-27 17:58:21 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:58:21 INFO Jobs still running: ['5289280', '5289282', '5289288', '5289290', '5289292', '5289293', '5289295', '5289299']. Waiting...
2026-07-27 17:58:36 INFO None 5289275: status FINISHED
2026-07-27 17:58:36 INFO None 5289276: status FINISHED
2026-07-27 17:58:36 INFO None 5289278: status FINISHED
2026-07-27 17:58:36 INFO None 5289279: status FINISHED
2026-07-27 17:58:36 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:58:36 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:58:36 INFO None 5289286: status FINISHED
2026-07-27 17:58:36 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:58:36 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:58:36 INFO None 5289292: status FINISHED
2026-07-27 17:58:36 INFO None 5289293: status FINISHED
2026-07-27 17:58:36 INFO None 5289294: status FINISHED
2026-07-27 17:58:36 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:58:36 INFO None 5289296: status FINISHED
2026-07-27 17:58:36 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:58:36 INFO Jobs still running: ['5289280', '5289282', '5289288', '5289290', '5289295', '5289299']. Waiting...
2026-07-27 17:58:53 INFO None 5289275: status FINISHED
2026-07-27 17:58:53 INFO None 5289276: status FINISHED
2026-07-27 17:58:53 INFO None 5289278: status FINISHED
2026-07-27 17:58:53 INFO None 5289279: status FINISHED
2026-07-27 17:58:53 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:58:53 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:58:53 INFO None 5289286: status FINISHED
2026-07-27 17:58:53 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:58:53 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:58:53 INFO None 5289292: status FINISHED
2026-07-27 17:58:53 INFO None 5289293: status FINISHED
2026-07-27 17:58:53 INFO None 5289294: status FINISHED
2026-07-27 17:58:53 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:58:53 INFO None 5289296: status FINISHED
2026-07-27 17:58:53 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:58:53 INFO Jobs still running: ['5289280', '5289282', '5289288', '5289290', '5289295', '5289299']. Waiting...
2026-07-27 17:59:08 INFO None 5289275: status FINISHED
2026-07-27 17:59:08 INFO None 5289276: status FINISHED
2026-07-27 17:59:08 INFO None 5289278: status FINISHED
2026-07-27 17:59:08 INFO None 5289279: status FINISHED
2026-07-27 17:59:08 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:59:08 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:59:08 INFO None 5289286: status FINISHED
2026-07-27 17:59:08 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:59:08 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:59:08 INFO None 5289292: status FINISHED
2026-07-27 17:59:08 INFO None 5289293: status FINISHED
2026-07-27 17:59:08 INFO None 5289294: status FINISHED
2026-07-27 17:59:08 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:59:08 INFO None 5289296: status FINISHED
2026-07-27 17:59:08 INFO None 5289299: status RUNNING/PENDING
2026-07-27 17:59:08 INFO Jobs still running: ['5289280', '5289282', '5289288', '5289290', '5289295', '5289299']. Waiting...
2026-07-27 17:59:23 INFO None 5289275: status FINISHED
2026-07-27 17:59:23 INFO None 5289276: status FINISHED
2026-07-27 17:59:23 INFO None 5289278: status FINISHED
2026-07-27 17:59:25 INFO None 5289279: status FINISHED
2026-07-27 17:59:25 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:59:25 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:59:25 INFO None 5289286: status FINISHED
2026-07-27 17:59:25 INFO None 5289288: status RUNNING/PENDING
2026-07-27 17:59:25 INFO None 5289290: status RUNNING/PENDING
2026-07-27 17:59:25 INFO None 5289292: status FINISHED
2026-07-27 17:59:25 INFO None 5289293: status FINISHED
2026-07-27 17:59:25 INFO None 5289294: status FINISHED
2026-07-27 17:59:25 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:59:25 INFO None 5289296: status FINISHED
2026-07-27 17:59:25 INFO None 5289299: status FINISHED
2026-07-27 17:59:25 INFO Jobs still running: ['5289280', '5289282', '5289288', '5289290', '5289295']. Waiting...
2026-07-27 17:59:40 INFO None 5289275: status FINISHED
2026-07-27 17:59:40 INFO None 5289276: status FINISHED
2026-07-27 17:59:40 INFO None 5289278: status FINISHED
2026-07-27 17:59:41 INFO None 5289279: status FINISHED
2026-07-27 17:59:41 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:59:41 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:59:41 INFO None 5289286: status FINISHED
2026-07-27 17:59:41 INFO None 5289288: status FINISHED
2026-07-27 17:59:41 INFO None 5289290: status FINISHED
2026-07-27 17:59:41 INFO None 5289292: status FINISHED
2026-07-27 17:59:41 INFO None 5289293: status FINISHED
2026-07-27 17:59:41 INFO None 5289294: status FINISHED
2026-07-27 17:59:41 INFO None 5289295: status RUNNING/PENDING
2026-07-27 17:59:41 INFO None 5289296: status FINISHED
2026-07-27 17:59:43 INFO None 5289299: status FINISHED
2026-07-27 17:59:43 INFO Jobs still running: ['5289280', '5289282', '5289295']. Waiting...
2026-07-27 17:59:58 INFO None 5289275: status FINISHED
2026-07-27 17:59:58 INFO None 5289276: status FINISHED
2026-07-27 17:59:58 INFO None 5289278: status FINISHED
2026-07-27 17:59:58 INFO None 5289279: status FINISHED
2026-07-27 17:59:58 INFO None 5289280: status RUNNING/PENDING
2026-07-27 17:59:58 INFO None 5289282: status RUNNING/PENDING
2026-07-27 17:59:58 INFO None 5289286: status FINISHED
2026-07-27 17:59:58 INFO None 5289288: status FINISHED
2026-07-27 17:59:58 INFO None 5289290: status FINISHED
2026-07-27 17:59:58 INFO None 5289292: status FINISHED
2026-07-27 17:59:58 INFO None 5289293: status FINISHED
2026-07-27 17:59:58 INFO None 5289294: status FINISHED
2026-07-27 17:59:58 INFO None 5289295: status FINISHED
2026-07-27 17:59:58 INFO None 5289296: status FINISHED
2026-07-27 17:59:58 INFO None 5289299: status FINISHED
2026-07-27 17:59:58 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:00:13 INFO None 5289275: status FINISHED
2026-07-27 18:00:13 INFO None 5289276: status FINISHED
2026-07-27 18:00:13 INFO None 5289278: status FINISHED
2026-07-27 18:00:13 INFO None 5289279: status FINISHED
2026-07-27 18:00:13 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:00:13 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:00:13 INFO None 5289286: status FINISHED
2026-07-27 18:00:13 INFO None 5289288: status FINISHED
2026-07-27 18:00:13 INFO None 5289290: status FINISHED
2026-07-27 18:00:13 INFO None 5289292: status FINISHED
2026-07-27 18:00:13 INFO None 5289293: status FINISHED
2026-07-27 18:00:15 INFO None 5289294: status FINISHED
2026-07-27 18:00:15 INFO None 5289295: status FINISHED
2026-07-27 18:00:15 INFO None 5289296: status FINISHED
2026-07-27 18:00:16 INFO None 5289299: status FINISHED
2026-07-27 18:00:16 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:00:31 INFO None 5289275: status FINISHED
2026-07-27 18:00:31 INFO None 5289276: status FINISHED
2026-07-27 18:00:31 INFO None 5289278: status FINISHED
2026-07-27 18:00:31 INFO None 5289279: status FINISHED
2026-07-27 18:00:31 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:00:31 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:00:31 INFO None 5289286: status FINISHED
2026-07-27 18:00:31 INFO None 5289288: status FINISHED
2026-07-27 18:00:31 INFO None 5289290: status FINISHED
2026-07-27 18:00:31 INFO None 5289292: status FINISHED
2026-07-27 18:00:31 INFO None 5289293: status FINISHED
2026-07-27 18:00:31 INFO None 5289294: status FINISHED
2026-07-27 18:00:31 INFO None 5289295: status FINISHED
2026-07-27 18:00:31 INFO None 5289296: status FINISHED
2026-07-27 18:00:31 INFO None 5289299: status FINISHED
2026-07-27 18:00:31 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:00:46 INFO None 5289275: status FINISHED
2026-07-27 18:00:46 INFO None 5289276: status FINISHED
2026-07-27 18:00:46 INFO None 5289278: status FINISHED
2026-07-27 18:00:46 INFO None 5289279: status FINISHED
2026-07-27 18:00:46 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:00:46 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:00:46 INFO None 5289286: status FINISHED
2026-07-27 18:00:46 INFO None 5289288: status FINISHED
2026-07-27 18:00:46 INFO None 5289290: status FINISHED
2026-07-27 18:00:46 INFO None 5289292: status FINISHED
2026-07-27 18:00:46 INFO None 5289293: status FINISHED
2026-07-27 18:00:46 INFO None 5289294: status FINISHED
2026-07-27 18:00:46 INFO None 5289295: status FINISHED
2026-07-27 18:00:46 INFO None 5289296: status FINISHED
2026-07-27 18:00:46 INFO None 5289299: status FINISHED
2026-07-27 18:00:46 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:01:01 INFO None 5289275: status FINISHED
2026-07-27 18:01:01 INFO None 5289276: status FINISHED
2026-07-27 18:01:01 INFO None 5289278: status FINISHED
2026-07-27 18:01:01 INFO None 5289279: status FINISHED
2026-07-27 18:01:01 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:01:01 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:01:01 INFO None 5289286: status FINISHED
2026-07-27 18:01:01 INFO None 5289288: status FINISHED
2026-07-27 18:01:01 INFO None 5289290: status FINISHED
2026-07-27 18:01:01 INFO None 5289292: status FINISHED
2026-07-27 18:01:01 INFO None 5289293: status FINISHED
2026-07-27 18:01:01 INFO None 5289294: status FINISHED
2026-07-27 18:01:01 INFO None 5289295: status FINISHED
2026-07-27 18:01:01 INFO None 5289296: status FINISHED
2026-07-27 18:01:01 INFO None 5289299: status FINISHED
2026-07-27 18:01:01 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:01:18 INFO None 5289275: status FINISHED
2026-07-27 18:01:18 INFO None 5289276: status FINISHED
2026-07-27 18:01:18 INFO None 5289278: status FINISHED
2026-07-27 18:01:18 INFO None 5289279: status FINISHED
2026-07-27 18:01:18 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:01:18 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:01:18 INFO None 5289286: status FINISHED
2026-07-27 18:01:18 INFO None 5289288: status FINISHED
2026-07-27 18:01:18 INFO None 5289290: status FINISHED
2026-07-27 18:01:18 INFO None 5289292: status FINISHED
2026-07-27 18:01:18 INFO None 5289293: status FINISHED
2026-07-27 18:01:18 INFO None 5289294: status FINISHED
2026-07-27 18:01:18 INFO None 5289295: status FINISHED
2026-07-27 18:01:18 INFO None 5289296: status FINISHED
2026-07-27 18:01:18 INFO None 5289299: status FINISHED
2026-07-27 18:01:18 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:01:33 INFO None 5289275: status FINISHED
2026-07-27 18:01:33 INFO None 5289276: status FINISHED
2026-07-27 18:01:33 INFO None 5289278: status FINISHED
2026-07-27 18:01:33 INFO None 5289279: status FINISHED
2026-07-27 18:01:33 INFO None 5289280: status RUNNING/PENDING
2026-07-27 18:01:33 INFO None 5289282: status RUNNING/PENDING
2026-07-27 18:01:33 INFO None 5289286: status FINISHED
2026-07-27 18:01:36 INFO None 5289288: status FINISHED
2026-07-27 18:01:36 INFO None 5289290: status FINISHED
2026-07-27 18:01:36 INFO None 5289292: status FINISHED
2026-07-27 18:01:36 INFO None 5289293: status FINISHED
2026-07-27 18:01:36 INFO None 5289294: status FINISHED
2026-07-27 18:01:36 INFO None 5289295: status FINISHED
2026-07-27 18:01:36 INFO None 5289296: status FINISHED
2026-07-27 18:01:36 INFO None 5289299: status FINISHED
2026-07-27 18:01:36 INFO Jobs still running: ['5289280', '5289282']. Waiting...
2026-07-27 18:01:51 INFO None 5289275: status FINISHED
2026-07-27 18:01:51 INFO None 5289276: status FINISHED
2026-07-27 18:01:51 INFO None 5289278: status FINISHED
2026-07-27 18:01:51 INFO None 5289279: status FINISHED
2026-07-27 18:01:51 INFO None 5289280: status FINISHED
2026-07-27 18:01:51 INFO None 5289282: status FINISHED
2026-07-27 18:01:51 INFO None 5289286: status FINISHED
2026-07-27 18:01:51 INFO None 5289288: status FINISHED
2026-07-27 18:01:51 INFO None 5289290: status FINISHED
2026-07-27 18:01:51 INFO None 5289292: status FINISHED
2026-07-27 18:01:51 INFO None 5289293: status FINISHED
2026-07-27 18:01:51 INFO None 5289294: status FINISHED
2026-07-27 18:01:51 INFO None 5289295: status FINISHED
2026-07-27 18:01:51 INFO None 5289296: status FINISHED
2026-07-27 18:01:51 INFO None 5289299: status FINISHED
2026-07-27 18:01:51 INFO Jobs ['5289275', '5289276', '5289278', '5289279', '5289280', '5289282', '5289286', '5289288', '5289290', '5289292', '5289293', '5289294', '5289295', '5289296', '5289299'] have finished
2026-07-27 18:01:51 INFO Checking restart files were created ...
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc(1002685915 bytes)
2026-07-27 18:01:51 INFO  Run_model() completed successfully.
2026-07-27 18:01:51 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 12:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:01:51 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 14:00:00 days=153073 seconds=50400
2026-07-27 18:01:51 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 18:01:51 INFO [TIME] increment current_time 2020-02-07 12:00:00 -> 2020-02-07 14:00:00
2026-07-27 18:01:51 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:01:51 INFO ---------->>> Running process_satellite_data()
2026-07-27 18:01:51 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12019.nc
2026-07-27 18:01:51 INFO ---------->>> Running run_obs_converter()
2026-07-27 18:01:51 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-27 18:01:51 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_50083_153073.out
2026-07-27 18:01:51 INFO ---------->>> Running DART
2026-07-27 18:01:51 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 18:01:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020714_1_out_toDART.nc
2026-07-27 18:01:51 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020714_1_out_toDART.nc
2026-07-27 18:01:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020714_1_out_toDART.nc
2026-07-27 18:01:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020714_1_out_toDART.nc
2026-07-27 18:01:52 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020714_1_out_toDART.nc
2026-07-27 18:01:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020714_1_out_toDART.nc
2026-07-27 18:01:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020714_1_out_toDART.nc
2026-07-27 18:01:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020714_1_out_toDART.nc
2026-07-27 18:01:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020714_1_out_toDART.nc
2026-07-27 18:01:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020714_1_out_toDART.nc
2026-07-27 18:01:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020714_1_out_toDART.nc
2026-07-27 18:01:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020714_1_out_toDART.nc
2026-07-27 18:01:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020714_1_out_toDART.nc
2026-07-27 18:01:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020714_1_out_toDART.nc
2026-07-27 18:01:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020712_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020714_1_out_toDART.nc
2026-07-27 18:01:56 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 18:01:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 18:01:56 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 18:01:56 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 18:01:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 18:01:56 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 18:02:13 INFO Found: []
2026-07-27 18:02:13 INFO No job id returned by command ./run_filter.bsh
2026-07-27 18:02:13 INFO No monitoring will be performed
2026-07-27 18:02:13 INFO Moving DART output files to analysis and preassim directories for date 2020020714 if present ...
2026-07-27 18:02:13 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_postinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:13 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_postinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020714'
2026-07-27 18:02:14 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 18:02:14 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 18:02:14 INFO run_dart() is DONE.
2026-07-27 18:02:14 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 18:02:14 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-27 18:02:19 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-27 18:02:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:25 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-27 18:02:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-27 18:02:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-27 18:02:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-27 18:02:47 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:47 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-27 18:02:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-27 18:02:59 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:02:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-27 18:03:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-27 18:03:10 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-27 18:03:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:16 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-27 18:03:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-27 18:03:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-27 18:03:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-27 18:03:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 18:03:38 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 18:03:38 INFO [TIME] step_end current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:03:38 INFO [TIME] step_start current_time=2020-02-07 14:00:00 simulated_time=2020-02-07 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:03:39 INFO [TIME] window start=2020-02-07 14:00:00 end=2020-02-08 00:00:00 run_hours=10 has_assimilation=False
2026-07-27 18:03:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:03:40 INFO Hourly dataset computed and listing created
2026-07-27 18:03:57 INFO Hourly dataset computed
2026-07-27 18:03:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:03:58 INFO Hourly dataset computed and listing created
2026-07-27 18:04:01 INFO Hourly dataset computed
2026-07-27 18:04:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:02 INFO Hourly dataset computed and listing created
2026-07-27 18:04:05 INFO Hourly dataset computed
2026-07-27 18:04:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:06 INFO Hourly dataset computed and listing created
2026-07-27 18:04:09 INFO Hourly dataset computed
2026-07-27 18:04:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:10 INFO Hourly dataset computed and listing created
2026-07-27 18:04:12 INFO Hourly dataset computed
2026-07-27 18:04:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:14 INFO Hourly dataset computed and listing created
2026-07-27 18:04:17 INFO Hourly dataset computed
2026-07-27 18:04:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:18 INFO Hourly dataset computed and listing created
2026-07-27 18:04:21 INFO Hourly dataset computed
2026-07-27 18:04:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:22 INFO Hourly dataset computed and listing created
2026-07-27 18:04:25 INFO Hourly dataset computed
2026-07-27 18:04:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:26 INFO Hourly dataset computed and listing created
2026-07-27 18:04:28 INFO Hourly dataset computed
2026-07-27 18:04:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:30 INFO Hourly dataset computed and listing created
2026-07-27 18:04:32 INFO Hourly dataset computed
2026-07-27 18:04:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:33 INFO Hourly dataset computed and listing created
2026-07-27 18:04:36 INFO Hourly dataset computed
2026-07-27 18:04:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:37 INFO Hourly dataset computed and listing created
2026-07-27 18:04:40 INFO Hourly dataset computed
2026-07-27 18:04:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:41 INFO Hourly dataset computed and listing created
2026-07-27 18:04:44 INFO Hourly dataset computed
2026-07-27 18:04:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:45 INFO Hourly dataset computed and listing created
2026-07-27 18:04:47 INFO Hourly dataset computed
2026-07-27 18:04:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 18:04:49 INFO Hourly dataset computed and listing created
2026-07-27 18:04:51 INFO Hourly dataset computed
2026-07-27 18:04:51 INFO ---------->>> Running CHIMERE model from 2020-02-07 14:00:00 to 2020-02-08 00:00:00
2026-07-27 18:04:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:04:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 18:04:51 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020712_2_ENS1.nc
2026-07-27 18:04:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 18:04:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:04:51 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 18:04:51 INFO Queuing job for member 1...
2026-07-27 18:04:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:04:51 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 18:04:52 INFO Found: ['5289368']
2026-07-27 18:04:57 INFO [TGCC-IRENE] Submitted job with ID:['5289368']
2026-07-27 18:04:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:04:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 18:04:57 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020712_2_ENS2.nc
2026-07-27 18:04:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 18:04:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:04:57 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 18:04:57 INFO Queuing job for member 2...
2026-07-27 18:04:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:04:57 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 18:04:58 INFO Found: ['5289369']
2026-07-27 18:05:03 INFO [TGCC-IRENE] Submitted job with ID:['5289369']
2026-07-27 18:05:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 18:05:03 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020712_2_ENS3.nc
2026-07-27 18:05:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 18:05:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:03 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 18:05:03 INFO Queuing job for member 3...
2026-07-27 18:05:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:03 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 18:05:04 INFO Found: ['5289376']
2026-07-27 18:05:09 INFO [TGCC-IRENE] Submitted job with ID:['5289376']
2026-07-27 18:05:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 18:05:09 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020712_2_ENS4.nc
2026-07-27 18:05:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 18:05:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:09 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 18:05:09 INFO Queuing job for member 4...
2026-07-27 18:05:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:09 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 18:05:09 INFO Found: ['5289380']
2026-07-27 18:05:14 INFO [TGCC-IRENE] Submitted job with ID:['5289380']
2026-07-27 18:05:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 18:05:14 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020712_2_ENS5.nc
2026-07-27 18:05:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 18:05:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:14 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 18:05:14 INFO Queuing job for member 5...
2026-07-27 18:05:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:14 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 18:05:15 INFO Found: ['5289381']
2026-07-27 18:05:20 INFO [TGCC-IRENE] Submitted job with ID:['5289381']
2026-07-27 18:05:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 18:05:20 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020712_2_ENS6.nc
2026-07-27 18:05:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 18:05:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:20 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 18:05:20 INFO Queuing job for member 6...
2026-07-27 18:05:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:20 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 18:05:21 INFO Found: ['5289382']
2026-07-27 18:05:26 INFO [TGCC-IRENE] Submitted job with ID:['5289382']
2026-07-27 18:05:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 18:05:26 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020712_2_ENS7.nc
2026-07-27 18:05:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 18:05:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:26 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 18:05:26 INFO Queuing job for member 7...
2026-07-27 18:05:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:26 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 18:05:26 INFO Found: ['5289383']
2026-07-27 18:05:31 INFO [TGCC-IRENE] Submitted job with ID:['5289383']
2026-07-27 18:05:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 18:05:31 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020712_2_ENS8.nc
2026-07-27 18:05:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 18:05:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:31 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 18:05:31 INFO Queuing job for member 8...
2026-07-27 18:05:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:31 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 18:05:32 INFO Found: ['5289384']
2026-07-27 18:05:37 INFO [TGCC-IRENE] Submitted job with ID:['5289384']
2026-07-27 18:05:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 18:05:37 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020712_2_ENS9.nc
2026-07-27 18:05:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 18:05:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:37 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 18:05:37 INFO Queuing job for member 9...
2026-07-27 18:05:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:37 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 18:05:38 INFO Found: ['5289385']
2026-07-27 18:05:43 INFO [TGCC-IRENE] Submitted job with ID:['5289385']
2026-07-27 18:05:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 18:05:43 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020712_2_ENS10.nc
2026-07-27 18:05:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 18:05:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:43 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 18:05:43 INFO Queuing job for member 10...
2026-07-27 18:05:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:43 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 18:05:45 INFO Found: ['5289386']
2026-07-27 18:05:50 INFO [TGCC-IRENE] Submitted job with ID:['5289386']
2026-07-27 18:05:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 18:05:50 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020712_2_ENS11.nc
2026-07-27 18:05:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 18:05:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:50 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 18:05:50 INFO Queuing job for member 11...
2026-07-27 18:05:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:50 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 18:05:53 INFO Found: ['5289388']
2026-07-27 18:05:58 INFO [TGCC-IRENE] Submitted job with ID:['5289388']
2026-07-27 18:05:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:05:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 18:05:58 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020712_2_ENS12.nc
2026-07-27 18:05:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 18:05:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:05:58 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 18:05:58 INFO Queuing job for member 12...
2026-07-27 18:05:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:05:58 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 18:06:00 INFO Found: ['5289390']
2026-07-27 18:06:05 INFO [TGCC-IRENE] Submitted job with ID:['5289390']
2026-07-27 18:06:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:06:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 18:06:05 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020712_2_ENS13.nc
2026-07-27 18:06:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 18:06:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:06:05 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 18:06:05 INFO Queuing job for member 13...
2026-07-27 18:06:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:06:05 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 18:06:08 INFO Found: ['5289392']
2026-07-27 18:06:13 INFO [TGCC-IRENE] Submitted job with ID:['5289392']
2026-07-27 18:06:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:06:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 18:06:13 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020712_2_ENS14.nc
2026-07-27 18:06:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 18:06:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:06:13 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 18:06:13 INFO Queuing job for member 14...
2026-07-27 18:06:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:06:13 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 18:06:15 INFO Found: ['5289394']
2026-07-27 18:06:20 INFO [TGCC-IRENE] Submitted job with ID:['5289394']
2026-07-27 18:06:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 18:06:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 18:06:20 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020712_2_ENS15.nc
2026-07-27 18:06:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 18:06:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 18:06:20 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 18:06:20 INFO Queuing job for member 15...
2026-07-27 18:06:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 18:06:20 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 18:06:23 INFO Found: ['5289396']
2026-07-27 18:06:28 INFO [TGCC-IRENE] Submitted job with ID:['5289396']
2026-07-27 18:06:28 INFO Checking job status ...
2026-07-27 18:06:28 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:06:28 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:06:28 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:06:43 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:06:43 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:06:45 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:06:46 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:06:46 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:06:46 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:07:01 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:07:01 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:07:01 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:07:16 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:07:16 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:07:16 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:07:31 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:07:31 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:07:31 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:07:46 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:07:46 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:07:46 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:07:46 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:07:46 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:07:47 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:07:47 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:08:03 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:08:03 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:08:03 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:08:18 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:08:18 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:08:21 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:08:21 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:08:36 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:08:36 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:08:36 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:08:51 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:08:51 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:08:51 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:09:06 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:09:06 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:09:06 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:09:08 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:09:08 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:09:23 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:09:23 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:09:24 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:09:24 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:09:39 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:09:39 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:09:39 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:09:54 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:09:54 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:09:54 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:10:10 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:10:10 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:10:10 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:10:25 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:10:25 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:10:26 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:10:26 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:10:26 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:10:26 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:10:26 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:10:28 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:10:28 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:10:43 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:10:43 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:10:43 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:10:58 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:10:58 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:11:00 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:11:00 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:11:15 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:11:15 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:11:15 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:11:15 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:11:15 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:11:16 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:11:16 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:11:31 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:11:31 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:11:33 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:11:33 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:11:48 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:11:48 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:11:48 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:12:03 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:12:03 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:12:04 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:12:04 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:12:04 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:12:04 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:12:04 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:12:04 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:12:19 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:12:19 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:12:19 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:12:35 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:12:36 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:12:36 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:12:51 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:12:51 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:12:53 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:12:53 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:13:08 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:13:08 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:13:09 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:13:09 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:13:09 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:13:24 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:13:24 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:13:24 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:13:24 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:13:24 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:13:26 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:13:26 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:13:41 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:13:41 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:13:41 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:13:56 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:13:56 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:13:56 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:14:11 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:14:11 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:14:12 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:14:12 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:14:27 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:14:27 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:14:27 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:14:43 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:14:43 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:14:43 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:14:58 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:14:58 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:15:00 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:15:00 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:15:00 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:15:15 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:15:15 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:15:16 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:15:16 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:15:16 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:15:31 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:15:31 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:15:33 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:15:33 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:15:33 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:15:33 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:15:33 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:15:33 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:15:48 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:15:48 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:15:48 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:16:03 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:16:03 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:16:05 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:16:05 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:16:05 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:16:05 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:16:06 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:16:06 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:16:21 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:16:21 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:16:21 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:16:36 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:16:36 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:16:36 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:16:53 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289369: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289381: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289382: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:16:53 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:16:53 INFO Jobs still running: ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:17:08 INFO None 5289368: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289369: status FINISHED
2026-07-27 18:17:08 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289381: status FINISHED
2026-07-27 18:17:08 INFO None 5289382: status FINISHED
2026-07-27 18:17:08 INFO None 5289383: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289384: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:17:08 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:17:10 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:17:11 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:17:11 INFO Jobs still running: ['5289368', '5289376', '5289380', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:17:26 INFO None 5289368: status FINISHED
2026-07-27 18:17:26 INFO None 5289369: status FINISHED
2026-07-27 18:17:26 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289380: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289381: status FINISHED
2026-07-27 18:17:26 INFO None 5289382: status FINISHED
2026-07-27 18:17:26 INFO None 5289383: status FINISHED
2026-07-27 18:17:26 INFO None 5289384: status FINISHED
2026-07-27 18:17:26 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:17:26 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:17:26 INFO Jobs still running: ['5289376', '5289380', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:17:41 INFO None 5289368: status FINISHED
2026-07-27 18:17:41 INFO None 5289369: status FINISHED
2026-07-27 18:17:41 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:17:41 INFO None 5289380: status FINISHED
2026-07-27 18:17:41 INFO None 5289381: status FINISHED
2026-07-27 18:17:41 INFO None 5289382: status FINISHED
2026-07-27 18:17:41 INFO None 5289383: status FINISHED
2026-07-27 18:17:41 INFO None 5289384: status FINISHED
2026-07-27 18:17:41 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:17:41 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:17:41 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:17:41 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:17:43 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:17:43 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:17:43 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:17:43 INFO Jobs still running: ['5289376', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:17:58 INFO None 5289368: status FINISHED
2026-07-27 18:17:58 INFO None 5289369: status FINISHED
2026-07-27 18:17:58 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289380: status FINISHED
2026-07-27 18:17:58 INFO None 5289381: status FINISHED
2026-07-27 18:17:58 INFO None 5289382: status FINISHED
2026-07-27 18:17:58 INFO None 5289383: status FINISHED
2026-07-27 18:17:58 INFO None 5289384: status FINISHED
2026-07-27 18:17:58 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:17:58 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:17:58 INFO Jobs still running: ['5289376', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:18:13 INFO None 5289368: status FINISHED
2026-07-27 18:18:13 INFO None 5289369: status FINISHED
2026-07-27 18:18:14 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:18:14 INFO None 5289380: status FINISHED
2026-07-27 18:18:14 INFO None 5289381: status FINISHED
2026-07-27 18:18:14 INFO None 5289382: status FINISHED
2026-07-27 18:18:14 INFO None 5289383: status FINISHED
2026-07-27 18:18:14 INFO None 5289384: status FINISHED
2026-07-27 18:18:14 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:18:14 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:18:16 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:18:16 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:18:16 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:18:16 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:18:16 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:18:16 INFO Jobs still running: ['5289376', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:18:31 INFO None 5289368: status FINISHED
2026-07-27 18:18:31 INFO None 5289369: status FINISHED
2026-07-27 18:18:31 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289380: status FINISHED
2026-07-27 18:18:31 INFO None 5289381: status FINISHED
2026-07-27 18:18:31 INFO None 5289382: status FINISHED
2026-07-27 18:18:31 INFO None 5289383: status FINISHED
2026-07-27 18:18:31 INFO None 5289384: status FINISHED
2026-07-27 18:18:31 INFO None 5289385: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289386: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289390: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289392: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289394: status RUNNING/PENDING
2026-07-27 18:18:31 INFO None 5289396: status RUNNING/PENDING
2026-07-27 18:18:31 INFO Jobs still running: ['5289376', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396']. Waiting...
2026-07-27 18:18:46 INFO None 5289368: status FINISHED
2026-07-27 18:18:46 INFO None 5289369: status FINISHED
2026-07-27 18:18:46 INFO None 5289376: status RUNNING/PENDING
2026-07-27 18:18:46 INFO None 5289380: status FINISHED
2026-07-27 18:18:46 INFO None 5289381: status FINISHED
2026-07-27 18:18:46 INFO None 5289382: status FINISHED
2026-07-27 18:18:46 INFO None 5289383: status FINISHED
2026-07-27 18:18:46 INFO None 5289384: status FINISHED
2026-07-27 18:18:46 INFO None 5289385: status FINISHED
2026-07-27 18:18:46 INFO None 5289386: status FINISHED
2026-07-27 18:18:47 INFO None 5289388: status RUNNING/PENDING
2026-07-27 18:18:47 INFO None 5289390: status FINISHED
2026-07-27 18:18:47 INFO None 5289392: status FINISHED
2026-07-27 18:18:47 INFO None 5289394: status FINISHED
2026-07-27 18:18:47 INFO None 5289396: status FINISHED
2026-07-27 18:18:47 INFO Jobs still running: ['5289376', '5289388']. Waiting...
2026-07-27 18:19:02 INFO None 5289368: status FINISHED
2026-07-27 18:19:02 INFO None 5289369: status FINISHED
2026-07-27 18:19:02 INFO None 5289376: status FINISHED
2026-07-27 18:19:02 INFO None 5289380: status FINISHED
2026-07-27 18:19:02 INFO None 5289381: status FINISHED
2026-07-27 18:19:02 INFO None 5289382: status FINISHED
2026-07-27 18:19:02 INFO None 5289383: status FINISHED
2026-07-27 18:19:02 INFO None 5289384: status FINISHED
2026-07-27 18:19:02 INFO None 5289385: status FINISHED
2026-07-27 18:19:02 INFO None 5289386: status FINISHED
2026-07-27 18:19:02 INFO None 5289388: status FINISHED
2026-07-27 18:19:02 INFO None 5289390: status FINISHED
2026-07-27 18:19:02 INFO None 5289392: status FINISHED
2026-07-27 18:19:02 INFO None 5289394: status FINISHED
2026-07-27 18:19:02 INFO None 5289396: status FINISHED
2026-07-27 18:19:02 INFO Jobs ['5289368', '5289369', '5289376', '5289380', '5289381', '5289382', '5289383', '5289384', '5289385', '5289386', '5289388', '5289390', '5289392', '5289394', '5289396'] have finished
2026-07-27 18:19:02 INFO Checking restart files were created ...
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020714_10_ENS1.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020714_10_ENS2.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020714_10_ENS3.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020714_10_ENS4.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020714_10_ENS5.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020714_10_ENS6.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020714_10_ENS7.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020714_10_ENS8.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020714_10_ENS9.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020714_10_ENS10.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020714_10_ENS11.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020714_10_ENS12.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020714_10_ENS13.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020714_10_ENS14.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020714_10_ENS15.nc(3673513755 bytes)
2026-07-27 18:19:02 INFO  Run_model() completed successfully.
2026-07-27 18:19:02 INFO [TIME] after_model_set_simulated_time current_time=2020-02-07 14:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:19:02 INFO [TIME] gregorian_conversion simulated_time=2020-02-08 00:00:00 days=153074 seconds=0
2026-07-27 18:19:02 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 18:19:02 INFO [TIME] increment current_time 2020-02-07 14:00:00 -> 2020-02-08 00:00:00
2026-07-27 18:19:02 INFO [TIME] after_increment_before_assimilation current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:19:02 INFO ---------->>> Running process_satellite_data()
2026-07-27 18:19:02 INFO [DART] No satellite data found, skipping assimilation
2026-07-27 18:19:02 INFO after_assimilation() skipped
2026-07-27 18:19:02 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 18:19:02 INFO [TIME] step_end current_time=2020-02-08 00:00:00 simulated_time=2020-02-08 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 18:19:02 INFO [PIPELINE] ---- TIME LOOP END ----
+ exit 0
