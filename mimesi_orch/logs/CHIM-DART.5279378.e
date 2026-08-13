+ /bin/bash -x /tmp/tmp.2JBdxuUynW
+ SCRIPT_PID=2012210
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
2026-07-26 11:41:09 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-26 11:41:09 INFO [PIPELINE] =======================================
2026-07-26 11:41:09 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-26 11:41:09 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-26 11:41:09 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-26 11:41:09 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260726_114109.log
2026-07-26 11:41:09 INFO [PIPELINE] =======================================
2026-07-26 11:41:09 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-26 11:41:09 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-26 11:41:09 INFO [STEP] ---- TIME LOOP START ----
2026-07-26 11:41:09 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 11:41:09 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-26 11:41:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:17 INFO Hourly dataset computed and listing created
2026-07-26 11:41:23 INFO Hourly dataset computed
2026-07-26 11:41:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:24 INFO Hourly dataset computed and listing created
2026-07-26 11:41:25 INFO Hourly dataset computed
2026-07-26 11:41:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:26 INFO Hourly dataset computed and listing created
2026-07-26 11:41:26 INFO Hourly dataset computed
2026-07-26 11:41:27 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:27 INFO Hourly dataset computed and listing created
2026-07-26 11:41:28 INFO Hourly dataset computed
2026-07-26 11:41:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:29 INFO Hourly dataset computed and listing created
2026-07-26 11:41:30 INFO Hourly dataset computed
2026-07-26 11:41:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:31 INFO Hourly dataset computed and listing created
2026-07-26 11:41:32 INFO Hourly dataset computed
2026-07-26 11:41:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:33 INFO Hourly dataset computed and listing created
2026-07-26 11:41:33 INFO Hourly dataset computed
2026-07-26 11:41:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:34 INFO Hourly dataset computed and listing created
2026-07-26 11:41:35 INFO Hourly dataset computed
2026-07-26 11:41:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:36 INFO Hourly dataset computed and listing created
2026-07-26 11:41:37 INFO Hourly dataset computed
2026-07-26 11:41:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:38 INFO Hourly dataset computed and listing created
2026-07-26 11:41:39 INFO Hourly dataset computed
2026-07-26 11:41:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:39 INFO Hourly dataset computed and listing created
2026-07-26 11:41:40 INFO Hourly dataset computed
2026-07-26 11:41:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:41 INFO Hourly dataset computed and listing created
2026-07-26 11:41:42 INFO Hourly dataset computed
2026-07-26 11:41:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:43 INFO Hourly dataset computed and listing created
2026-07-26 11:41:44 INFO Hourly dataset computed
2026-07-26 11:41:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:45 INFO Hourly dataset computed and listing created
2026-07-26 11:41:46 INFO Hourly dataset computed
2026-07-26 11:41:46 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:41:46 INFO Hourly dataset computed and listing created
2026-07-26 11:41:47 INFO Hourly dataset computed
2026-07-26 11:41:47 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-26 11:41:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:41:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 11:41:47 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-26 11:41:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 11:41:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:41:47 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 11:41:47 INFO Queuing job for member 1...
2026-07-26 11:41:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:41:47 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 11:41:48 INFO Found: ['5279386']
2026-07-26 11:41:53 INFO [TGCC-IRENE] Submitted job with ID:['5279386']
2026-07-26 11:41:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:41:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 11:41:53 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-26 11:41:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 11:41:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:41:53 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 11:41:53 INFO Queuing job for member 2...
2026-07-26 11:41:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:41:53 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 11:41:54 INFO Found: ['5279387']
2026-07-26 11:41:59 INFO [TGCC-IRENE] Submitted job with ID:['5279387']
2026-07-26 11:41:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:41:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 11:41:59 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-26 11:41:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 11:41:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:41:59 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 11:41:59 INFO Queuing job for member 3...
2026-07-26 11:41:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:41:59 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 11:42:00 INFO Found: ['5279388']
2026-07-26 11:42:05 INFO [TGCC-IRENE] Submitted job with ID:['5279388']
2026-07-26 11:42:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 11:42:05 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-26 11:42:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 11:42:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:05 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 11:42:05 INFO Queuing job for member 4...
2026-07-26 11:42:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:05 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 11:42:05 INFO Found: ['5279390']
2026-07-26 11:42:10 INFO [TGCC-IRENE] Submitted job with ID:['5279390']
2026-07-26 11:42:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 11:42:10 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-26 11:42:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 11:42:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:10 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 11:42:10 INFO Queuing job for member 5...
2026-07-26 11:42:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:10 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 11:42:11 INFO Found: ['5279391']
2026-07-26 11:42:16 INFO [TGCC-IRENE] Submitted job with ID:['5279391']
2026-07-26 11:42:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 11:42:16 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-26 11:42:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 11:42:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:16 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 11:42:16 INFO Queuing job for member 6...
2026-07-26 11:42:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:16 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 11:42:17 INFO Found: ['5279392']
2026-07-26 11:42:22 INFO [TGCC-IRENE] Submitted job with ID:['5279392']
2026-07-26 11:42:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 11:42:22 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-26 11:42:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 11:42:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:22 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 11:42:22 INFO Queuing job for member 7...
2026-07-26 11:42:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:22 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 11:42:23 INFO Found: ['5279393']
2026-07-26 11:42:28 INFO [TGCC-IRENE] Submitted job with ID:['5279393']
2026-07-26 11:42:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 11:42:28 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-26 11:42:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 11:42:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:28 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 11:42:28 INFO Queuing job for member 8...
2026-07-26 11:42:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:28 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 11:42:28 INFO Found: ['5279394']
2026-07-26 11:42:33 INFO [TGCC-IRENE] Submitted job with ID:['5279394']
2026-07-26 11:42:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 11:42:33 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-26 11:42:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 11:42:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:33 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 11:42:33 INFO Queuing job for member 9...
2026-07-26 11:42:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:33 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 11:42:36 INFO Found: ['5279396']
2026-07-26 11:42:41 INFO [TGCC-IRENE] Submitted job with ID:['5279396']
2026-07-26 11:42:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 11:42:41 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-26 11:42:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 11:42:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:41 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 11:42:41 INFO Queuing job for member 10...
2026-07-26 11:42:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:41 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 11:42:43 INFO Found: ['5279397']
2026-07-26 11:42:48 INFO [TGCC-IRENE] Submitted job with ID:['5279397']
2026-07-26 11:42:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 11:42:48 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-26 11:42:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 11:42:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:48 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 11:42:48 INFO Queuing job for member 11...
2026-07-26 11:42:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:48 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 11:42:51 INFO Found: ['5279398']
2026-07-26 11:42:56 INFO [TGCC-IRENE] Submitted job with ID:['5279398']
2026-07-26 11:42:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:42:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 11:42:56 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-26 11:42:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 11:42:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:42:56 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 11:42:56 INFO Queuing job for member 12...
2026-07-26 11:42:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:42:56 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 11:42:58 INFO Found: ['5279399']
2026-07-26 11:43:03 INFO [TGCC-IRENE] Submitted job with ID:['5279399']
2026-07-26 11:43:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:43:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 11:43:03 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-26 11:43:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 11:43:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:43:03 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 11:43:03 INFO Queuing job for member 13...
2026-07-26 11:43:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:43:03 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 11:43:06 INFO Found: ['5279401']
2026-07-26 11:43:11 INFO [TGCC-IRENE] Submitted job with ID:['5279401']
2026-07-26 11:43:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:43:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 11:43:11 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-26 11:43:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 11:43:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:43:11 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 11:43:11 INFO Queuing job for member 14...
2026-07-26 11:43:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:43:11 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 11:43:14 INFO Found: ['5279402']
2026-07-26 11:43:19 INFO [TGCC-IRENE] Submitted job with ID:['5279402']
2026-07-26 11:43:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:43:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 11:43:19 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-26 11:43:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 11:43:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:43:19 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 11:43:19 INFO Queuing job for member 15...
2026-07-26 11:43:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:43:19 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 11:43:21 INFO Found: ['5279403']
2026-07-26 11:43:26 INFO [TGCC-IRENE] Submitted job with ID:['5279403']
2026-07-26 11:43:26 INFO Checking job status ...
2026-07-26 11:43:26 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:43:26 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:43:26 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:43:41 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:43:41 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:43:41 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:43:41 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:43:41 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:43:41 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:43:42 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:43:42 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:43:44 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:43:44 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:43:59 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:43:59 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:43:59 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:44:14 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:44:14 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:44:14 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:44:29 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:44:29 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:44:30 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:44:30 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:44:30 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:44:45 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:44:45 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:44:47 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:44:47 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:45:02 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:45:02 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:45:02 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:45:17 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:45:17 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:45:19 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:45:19 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:45:19 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:45:19 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:45:34 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:45:34 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:45:35 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:45:35 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:45:50 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:45:50 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:45:52 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:45:52 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:45:52 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:45:52 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:46:07 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:46:07 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:46:07 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:46:22 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:46:22 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:46:23 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:46:23 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:46:38 INFO None 5279386: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279387: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279388: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279390: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279391: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279392: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:46:38 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:46:38 INFO Jobs still running: ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:46:53 INFO None 5279386: status FINISHED
2026-07-26 11:46:55 INFO None 5279387: status FINISHED
2026-07-26 11:46:55 INFO None 5279388: status FINISHED
2026-07-26 11:46:55 INFO None 5279390: status FINISHED
2026-07-26 11:46:55 INFO None 5279391: status FINISHED
2026-07-26 11:46:55 INFO None 5279392: status FINISHED
2026-07-26 11:46:55 INFO None 5279393: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:46:55 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:46:55 INFO Jobs still running: ['5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:47:10 INFO None 5279386: status FINISHED
2026-07-26 11:47:10 INFO None 5279387: status FINISHED
2026-07-26 11:47:10 INFO None 5279388: status FINISHED
2026-07-26 11:47:10 INFO None 5279390: status FINISHED
2026-07-26 11:47:10 INFO None 5279391: status FINISHED
2026-07-26 11:47:10 INFO None 5279392: status FINISHED
2026-07-26 11:47:10 INFO None 5279393: status FINISHED
2026-07-26 11:47:10 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:47:10 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:47:11 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:47:11 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:47:26 INFO None 5279386: status FINISHED
2026-07-26 11:47:26 INFO None 5279387: status FINISHED
2026-07-26 11:47:26 INFO None 5279388: status FINISHED
2026-07-26 11:47:26 INFO None 5279390: status FINISHED
2026-07-26 11:47:28 INFO None 5279391: status FINISHED
2026-07-26 11:47:28 INFO None 5279392: status FINISHED
2026-07-26 11:47:28 INFO None 5279393: status FINISHED
2026-07-26 11:47:28 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:47:28 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:47:28 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:47:43 INFO None 5279386: status FINISHED
2026-07-26 11:47:43 INFO None 5279387: status FINISHED
2026-07-26 11:47:43 INFO None 5279388: status FINISHED
2026-07-26 11:47:43 INFO None 5279390: status FINISHED
2026-07-26 11:47:43 INFO None 5279391: status FINISHED
2026-07-26 11:47:43 INFO None 5279392: status FINISHED
2026-07-26 11:47:43 INFO None 5279393: status FINISHED
2026-07-26 11:47:43 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:47:43 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:47:43 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:47:58 INFO None 5279386: status FINISHED
2026-07-26 11:47:58 INFO None 5279387: status FINISHED
2026-07-26 11:47:58 INFO None 5279388: status FINISHED
2026-07-26 11:47:58 INFO None 5279390: status FINISHED
2026-07-26 11:47:58 INFO None 5279391: status FINISHED
2026-07-26 11:47:58 INFO None 5279392: status FINISHED
2026-07-26 11:47:58 INFO None 5279393: status FINISHED
2026-07-26 11:48:00 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:48:00 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:48:00 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:48:15 INFO None 5279386: status FINISHED
2026-07-26 11:48:15 INFO None 5279387: status FINISHED
2026-07-26 11:48:15 INFO None 5279388: status FINISHED
2026-07-26 11:48:15 INFO None 5279390: status FINISHED
2026-07-26 11:48:15 INFO None 5279391: status FINISHED
2026-07-26 11:48:16 INFO None 5279392: status FINISHED
2026-07-26 11:48:16 INFO None 5279393: status FINISHED
2026-07-26 11:48:16 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:48:16 INFO None 5279403: status RUNNING/PENDING
2026-07-26 11:48:16 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403']. Waiting...
2026-07-26 11:48:31 INFO None 5279386: status FINISHED
2026-07-26 11:48:31 INFO None 5279387: status FINISHED
2026-07-26 11:48:31 INFO None 5279388: status FINISHED
2026-07-26 11:48:31 INFO None 5279390: status FINISHED
2026-07-26 11:48:31 INFO None 5279391: status FINISHED
2026-07-26 11:48:31 INFO None 5279392: status FINISHED
2026-07-26 11:48:31 INFO None 5279393: status FINISHED
2026-07-26 11:48:31 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:48:31 INFO None 5279403: status FINISHED
2026-07-26 11:48:31 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:48:46 INFO None 5279386: status FINISHED
2026-07-26 11:48:46 INFO None 5279387: status FINISHED
2026-07-26 11:48:46 INFO None 5279388: status FINISHED
2026-07-26 11:48:46 INFO None 5279390: status FINISHED
2026-07-26 11:48:46 INFO None 5279391: status FINISHED
2026-07-26 11:48:46 INFO None 5279392: status FINISHED
2026-07-26 11:48:46 INFO None 5279393: status FINISHED
2026-07-26 11:48:46 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:48:46 INFO None 5279403: status FINISHED
2026-07-26 11:48:46 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:49:01 INFO None 5279386: status FINISHED
2026-07-26 11:49:01 INFO None 5279387: status FINISHED
2026-07-26 11:49:01 INFO None 5279388: status FINISHED
2026-07-26 11:49:01 INFO None 5279390: status FINISHED
2026-07-26 11:49:01 INFO None 5279391: status FINISHED
2026-07-26 11:49:01 INFO None 5279392: status FINISHED
2026-07-26 11:49:01 INFO None 5279393: status FINISHED
2026-07-26 11:49:01 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:49:01 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:49:01 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:49:01 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:49:01 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:49:02 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:49:02 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:49:02 INFO None 5279403: status FINISHED
2026-07-26 11:49:02 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:49:17 INFO None 5279386: status FINISHED
2026-07-26 11:49:17 INFO None 5279387: status FINISHED
2026-07-26 11:49:17 INFO None 5279388: status FINISHED
2026-07-26 11:49:17 INFO None 5279390: status FINISHED
2026-07-26 11:49:17 INFO None 5279391: status FINISHED
2026-07-26 11:49:17 INFO None 5279392: status FINISHED
2026-07-26 11:49:17 INFO None 5279393: status FINISHED
2026-07-26 11:49:17 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279398: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:49:17 INFO None 5279403: status FINISHED
2026-07-26 11:49:17 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:49:32 INFO None 5279386: status FINISHED
2026-07-26 11:49:32 INFO None 5279387: status FINISHED
2026-07-26 11:49:32 INFO None 5279388: status FINISHED
2026-07-26 11:49:32 INFO None 5279390: status FINISHED
2026-07-26 11:49:32 INFO None 5279391: status FINISHED
2026-07-26 11:49:32 INFO None 5279392: status FINISHED
2026-07-26 11:49:32 INFO None 5279393: status FINISHED
2026-07-26 11:49:32 INFO None 5279394: status RUNNING/PENDING
2026-07-26 11:49:32 INFO None 5279396: status RUNNING/PENDING
2026-07-26 11:49:33 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:49:33 INFO None 5279398: status FINISHED
2026-07-26 11:49:35 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:49:35 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:49:35 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:49:35 INFO None 5279403: status FINISHED
2026-07-26 11:49:35 INFO Jobs still running: ['5279394', '5279396', '5279397', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:49:50 INFO None 5279386: status FINISHED
2026-07-26 11:49:50 INFO None 5279387: status FINISHED
2026-07-26 11:49:50 INFO None 5279388: status FINISHED
2026-07-26 11:49:50 INFO None 5279390: status FINISHED
2026-07-26 11:49:50 INFO None 5279391: status FINISHED
2026-07-26 11:49:50 INFO None 5279392: status FINISHED
2026-07-26 11:49:50 INFO None 5279393: status FINISHED
2026-07-26 11:49:50 INFO None 5279394: status FINISHED
2026-07-26 11:49:50 INFO None 5279396: status FINISHED
2026-07-26 11:49:50 INFO None 5279397: status RUNNING/PENDING
2026-07-26 11:49:50 INFO None 5279398: status FINISHED
2026-07-26 11:49:50 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:49:50 INFO None 5279401: status RUNNING/PENDING
2026-07-26 11:49:50 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:49:50 INFO None 5279403: status FINISHED
2026-07-26 11:49:50 INFO Jobs still running: ['5279397', '5279399', '5279401', '5279402']. Waiting...
2026-07-26 11:50:05 INFO None 5279386: status FINISHED
2026-07-26 11:50:05 INFO None 5279387: status FINISHED
2026-07-26 11:50:05 INFO None 5279388: status FINISHED
2026-07-26 11:50:05 INFO None 5279390: status FINISHED
2026-07-26 11:50:05 INFO None 5279391: status FINISHED
2026-07-26 11:50:05 INFO None 5279392: status FINISHED
2026-07-26 11:50:05 INFO None 5279393: status FINISHED
2026-07-26 11:50:05 INFO None 5279394: status FINISHED
2026-07-26 11:50:05 INFO None 5279396: status FINISHED
2026-07-26 11:50:05 INFO None 5279397: status FINISHED
2026-07-26 11:50:05 INFO None 5279398: status FINISHED
2026-07-26 11:50:05 INFO None 5279399: status RUNNING/PENDING
2026-07-26 11:50:05 INFO None 5279401: status FINISHED
2026-07-26 11:50:05 INFO None 5279402: status RUNNING/PENDING
2026-07-26 11:50:07 INFO None 5279403: status FINISHED
2026-07-26 11:50:07 INFO Jobs still running: ['5279399', '5279402']. Waiting...
2026-07-26 11:50:22 INFO None 5279386: status FINISHED
2026-07-26 11:50:22 INFO None 5279387: status FINISHED
2026-07-26 11:50:22 INFO None 5279388: status FINISHED
2026-07-26 11:50:22 INFO None 5279390: status FINISHED
2026-07-26 11:50:22 INFO None 5279391: status FINISHED
2026-07-26 11:50:22 INFO None 5279392: status FINISHED
2026-07-26 11:50:22 INFO None 5279393: status FINISHED
2026-07-26 11:50:22 INFO None 5279394: status FINISHED
2026-07-26 11:50:22 INFO None 5279396: status FINISHED
2026-07-26 11:50:22 INFO None 5279397: status FINISHED
2026-07-26 11:50:22 INFO None 5279398: status FINISHED
2026-07-26 11:50:22 INFO None 5279399: status FINISHED
2026-07-26 11:50:22 INFO None 5279401: status FINISHED
2026-07-26 11:50:22 INFO None 5279402: status FINISHED
2026-07-26 11:50:22 INFO None 5279403: status FINISHED
2026-07-26 11:50:22 INFO Jobs ['5279386', '5279387', '5279388', '5279390', '5279391', '5279392', '5279393', '5279394', '5279396', '5279397', '5279398', '5279399', '5279401', '5279402', '5279403'] have finished
2026-07-26 11:50:22 INFO Checking restart files were created ...
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-26 11:50:22 INFO  Run_model() completed successfully.
2026-07-26 11:50:22 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 11:50:22 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-26 11:50:22 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 11:50:22 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-26 11:50:22 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 11:50:22 INFO ---------->>> Running process_satellite_data()
2026-07-26 11:50:23 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-26 11:50:23 INFO ---------->>> Running run_obs_converter()
2026-07-26 11:50:23 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-26 11:50:23 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-26 11:50:23 INFO ---------->>> Running DART
2026-07-26 11:50:23 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 11:50:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-26 11:50:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-26 11:50:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-26 11:50:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-26 11:50:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-26 11:50:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-26 11:50:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-26 11:50:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-26 11:50:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-26 11:50:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-26 11:50:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-26 11:50:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-26 11:50:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-26 11:50:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-26 11:50:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-26 11:50:27 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 11:50:27 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 11:50:27 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 11:50:27 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 11:50:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 11:50:27 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 11:50:42 INFO Found: []
2026-07-26 11:50:42 INFO No job id returned by command ./run_filter.bsh
2026-07-26 11:50:42 INFO No monitoring will be performed
2026-07-26 11:50:42 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-26 11:50:42 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020611'
2026-07-26 11:50:42 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 11:50:45 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 11:50:45 INFO run_dart() is DONE.
2026-07-26 11:50:45 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 11:50:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-26 11:50:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:50:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-26 11:50:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:50:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-26 11:51:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-26 11:51:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-26 11:51:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-26 11:51:17 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-26 11:51:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-26 11:51:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-26 11:51:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-26 11:51:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-26 11:51:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-26 11:51:50 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-26 11:51:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:51:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-26 11:52:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:52:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-26 11:52:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 11:52:06 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 11:52:06 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 11:52:06 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 11:52:06 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-26 11:52:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:07 INFO Hourly dataset computed and listing created
2026-07-26 11:52:12 INFO Hourly dataset computed
2026-07-26 11:52:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:13 INFO Hourly dataset computed and listing created
2026-07-26 11:52:14 INFO Hourly dataset computed
2026-07-26 11:52:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:15 INFO Hourly dataset computed and listing created
2026-07-26 11:52:15 INFO Hourly dataset computed
2026-07-26 11:52:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:16 INFO Hourly dataset computed and listing created
2026-07-26 11:52:17 INFO Hourly dataset computed
2026-07-26 11:52:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:18 INFO Hourly dataset computed and listing created
2026-07-26 11:52:19 INFO Hourly dataset computed
2026-07-26 11:52:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:20 INFO Hourly dataset computed and listing created
2026-07-26 11:52:21 INFO Hourly dataset computed
2026-07-26 11:52:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:22 INFO Hourly dataset computed and listing created
2026-07-26 11:52:22 INFO Hourly dataset computed
2026-07-26 11:52:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:23 INFO Hourly dataset computed and listing created
2026-07-26 11:52:24 INFO Hourly dataset computed
2026-07-26 11:52:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:25 INFO Hourly dataset computed and listing created
2026-07-26 11:52:26 INFO Hourly dataset computed
2026-07-26 11:52:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:27 INFO Hourly dataset computed and listing created
2026-07-26 11:52:28 INFO Hourly dataset computed
2026-07-26 11:52:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:28 INFO Hourly dataset computed and listing created
2026-07-26 11:52:29 INFO Hourly dataset computed
2026-07-26 11:52:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:30 INFO Hourly dataset computed and listing created
2026-07-26 11:52:31 INFO Hourly dataset computed
2026-07-26 11:52:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:32 INFO Hourly dataset computed and listing created
2026-07-26 11:52:33 INFO Hourly dataset computed
2026-07-26 11:52:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:34 INFO Hourly dataset computed and listing created
2026-07-26 11:52:34 INFO Hourly dataset computed
2026-07-26 11:52:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 11:52:35 INFO Hourly dataset computed and listing created
2026-07-26 11:52:36 INFO Hourly dataset computed
2026-07-26 11:52:36 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-26 11:52:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:52:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 11:52:36 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-26 11:52:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 11:52:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:52:36 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 11:52:36 INFO Queuing job for member 1...
2026-07-26 11:52:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:52:36 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 11:52:37 INFO Found: ['5279427']
2026-07-26 11:52:42 INFO [TGCC-IRENE] Submitted job with ID:['5279427']
2026-07-26 11:52:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:52:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 11:52:42 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-26 11:52:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 11:52:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:52:42 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 11:52:42 INFO Queuing job for member 2...
2026-07-26 11:52:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:52:42 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 11:52:43 INFO Found: ['5279428']
2026-07-26 11:52:48 INFO [TGCC-IRENE] Submitted job with ID:['5279428']
2026-07-26 11:52:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:52:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 11:52:48 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-26 11:52:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 11:52:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:52:48 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 11:52:48 INFO Queuing job for member 3...
2026-07-26 11:52:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:52:48 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 11:52:48 INFO Found: ['5279430']
2026-07-26 11:52:53 INFO [TGCC-IRENE] Submitted job with ID:['5279430']
2026-07-26 11:52:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:52:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 11:52:53 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-26 11:52:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 11:52:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:52:54 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 11:52:54 INFO Queuing job for member 4...
2026-07-26 11:52:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:52:54 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 11:52:54 INFO Found: ['5279431']
2026-07-26 11:52:59 INFO [TGCC-IRENE] Submitted job with ID:['5279431']
2026-07-26 11:52:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:52:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 11:52:59 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-26 11:52:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 11:52:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:52:59 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 11:52:59 INFO Queuing job for member 5...
2026-07-26 11:52:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:52:59 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 11:53:00 INFO Found: ['5279432']
2026-07-26 11:53:05 INFO [TGCC-IRENE] Submitted job with ID:['5279432']
2026-07-26 11:53:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 11:53:05 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-26 11:53:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 11:53:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:05 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 11:53:05 INFO Queuing job for member 6...
2026-07-26 11:53:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:05 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 11:53:06 INFO Found: ['5279434']
2026-07-26 11:53:11 INFO [TGCC-IRENE] Submitted job with ID:['5279434']
2026-07-26 11:53:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 11:53:11 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-26 11:53:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 11:53:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:11 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 11:53:11 INFO Queuing job for member 7...
2026-07-26 11:53:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:11 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 11:53:11 INFO Found: ['5279435']
2026-07-26 11:53:16 INFO [TGCC-IRENE] Submitted job with ID:['5279435']
2026-07-26 11:53:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 11:53:16 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-26 11:53:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 11:53:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:16 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 11:53:16 INFO Queuing job for member 8...
2026-07-26 11:53:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:16 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 11:53:17 INFO Found: ['5279436']
2026-07-26 11:53:22 INFO [TGCC-IRENE] Submitted job with ID:['5279436']
2026-07-26 11:53:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 11:53:22 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-26 11:53:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 11:53:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:22 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 11:53:22 INFO Queuing job for member 9...
2026-07-26 11:53:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:22 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 11:53:23 INFO Found: ['5279438']
2026-07-26 11:53:28 INFO [TGCC-IRENE] Submitted job with ID:['5279438']
2026-07-26 11:53:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 11:53:28 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-26 11:53:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 11:53:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:28 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 11:53:28 INFO Queuing job for member 10...
2026-07-26 11:53:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:28 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 11:53:30 INFO Found: ['5279439']
2026-07-26 11:53:35 INFO [TGCC-IRENE] Submitted job with ID:['5279439']
2026-07-26 11:53:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 11:53:35 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-26 11:53:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 11:53:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:35 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 11:53:35 INFO Queuing job for member 11...
2026-07-26 11:53:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:35 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 11:53:38 INFO Found: ['5279441']
2026-07-26 11:53:43 INFO [TGCC-IRENE] Submitted job with ID:['5279441']
2026-07-26 11:53:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 11:53:43 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-26 11:53:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 11:53:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:43 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 11:53:43 INFO Queuing job for member 12...
2026-07-26 11:53:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:43 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 11:53:45 INFO Found: ['5279442']
2026-07-26 11:53:50 INFO [TGCC-IRENE] Submitted job with ID:['5279442']
2026-07-26 11:53:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 11:53:50 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-26 11:53:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 11:53:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:50 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 11:53:51 INFO Queuing job for member 13...
2026-07-26 11:53:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:51 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 11:53:53 INFO Found: ['5279443']
2026-07-26 11:53:58 INFO [TGCC-IRENE] Submitted job with ID:['5279443']
2026-07-26 11:53:58 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:53:58 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 11:53:58 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-26 11:53:58 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 11:53:58 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:53:58 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 11:53:58 INFO Queuing job for member 14...
2026-07-26 11:53:58 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:53:58 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 11:54:00 INFO Found: ['5279444']
2026-07-26 11:54:05 INFO [TGCC-IRENE] Submitted job with ID:['5279444']
2026-07-26 11:54:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 11:54:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 11:54:05 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-26 11:54:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 11:54:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 11:54:05 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 11:54:05 INFO Queuing job for member 15...
2026-07-26 11:54:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 11:54:05 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 11:54:08 INFO Found: ['5279446']
2026-07-26 11:54:13 INFO [TGCC-IRENE] Submitted job with ID:['5279446']
2026-07-26 11:54:13 INFO Checking job status ...
2026-07-26 11:54:13 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:54:13 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:54:13 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:54:28 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:54:28 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:54:30 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:54:30 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:54:31 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:54:31 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:54:31 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:54:31 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:54:31 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:54:46 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:54:46 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:54:46 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:55:01 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:55:01 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:55:01 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:55:16 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:55:16 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:55:16 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:55:31 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:55:31 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:55:31 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:55:32 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:55:32 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:55:47 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:55:47 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:55:47 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:56:02 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:56:02 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:56:03 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:56:05 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:56:05 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:56:20 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:56:20 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:56:20 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:56:35 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:56:35 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:56:37 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:56:37 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:56:37 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:56:37 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:56:37 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:56:37 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:56:52 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:56:52 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:56:52 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:56:52 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279436: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:56:53 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:56:53 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:57:08 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279436: status FINISHED
2026-07-26 11:57:08 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:57:08 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:57:08 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:57:23 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279430: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279431: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279432: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279434: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279436: status FINISHED
2026-07-26 11:57:23 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:57:23 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:57:23 INFO Jobs still running: ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:57:38 INFO None 5279427: status RUNNING/PENDING
2026-07-26 11:57:38 INFO None 5279428: status RUNNING/PENDING
2026-07-26 11:57:38 INFO None 5279430: status FINISHED
2026-07-26 11:57:38 INFO None 5279431: status FINISHED
2026-07-26 11:57:38 INFO None 5279432: status FINISHED
2026-07-26 11:57:38 INFO None 5279434: status FINISHED
2026-07-26 11:57:38 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:57:38 INFO None 5279436: status FINISHED
2026-07-26 11:57:38 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:57:38 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:57:39 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:57:39 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:57:39 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:57:39 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:57:39 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:57:39 INFO Jobs still running: ['5279427', '5279428', '5279435', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:57:54 INFO None 5279427: status FINISHED
2026-07-26 11:57:54 INFO None 5279428: status FINISHED
2026-07-26 11:57:54 INFO None 5279430: status FINISHED
2026-07-26 11:57:54 INFO None 5279431: status FINISHED
2026-07-26 11:57:54 INFO None 5279432: status FINISHED
2026-07-26 11:57:54 INFO None 5279434: status FINISHED
2026-07-26 11:57:54 INFO None 5279435: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279436: status FINISHED
2026-07-26 11:57:54 INFO None 5279438: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:57:54 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:57:54 INFO Jobs still running: ['5279435', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:58:10 INFO None 5279427: status FINISHED
2026-07-26 11:58:10 INFO None 5279428: status FINISHED
2026-07-26 11:58:10 INFO None 5279430: status FINISHED
2026-07-26 11:58:10 INFO None 5279431: status FINISHED
2026-07-26 11:58:10 INFO None 5279432: status FINISHED
2026-07-26 11:58:10 INFO None 5279434: status FINISHED
2026-07-26 11:58:10 INFO None 5279435: status FINISHED
2026-07-26 11:58:10 INFO None 5279436: status FINISHED
2026-07-26 11:58:10 INFO None 5279438: status FINISHED
2026-07-26 11:58:10 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:58:10 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:58:10 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:58:10 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:58:10 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:58:10 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:58:10 INFO Jobs still running: ['5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:58:25 INFO None 5279427: status FINISHED
2026-07-26 11:58:25 INFO None 5279428: status FINISHED
2026-07-26 11:58:27 INFO None 5279430: status FINISHED
2026-07-26 11:58:27 INFO None 5279431: status FINISHED
2026-07-26 11:58:27 INFO None 5279432: status FINISHED
2026-07-26 11:58:27 INFO None 5279434: status FINISHED
2026-07-26 11:58:27 INFO None 5279435: status FINISHED
2026-07-26 11:58:27 INFO None 5279436: status FINISHED
2026-07-26 11:58:27 INFO None 5279438: status FINISHED
2026-07-26 11:58:27 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:58:27 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:58:27 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:58:27 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:58:27 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:58:27 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:58:27 INFO Jobs still running: ['5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:58:42 INFO None 5279427: status FINISHED
2026-07-26 11:58:42 INFO None 5279428: status FINISHED
2026-07-26 11:58:42 INFO None 5279430: status FINISHED
2026-07-26 11:58:42 INFO None 5279431: status FINISHED
2026-07-26 11:58:42 INFO None 5279432: status FINISHED
2026-07-26 11:58:42 INFO None 5279434: status FINISHED
2026-07-26 11:58:42 INFO None 5279435: status FINISHED
2026-07-26 11:58:42 INFO None 5279436: status FINISHED
2026-07-26 11:58:42 INFO None 5279438: status FINISHED
2026-07-26 11:58:42 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:58:42 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:58:42 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:58:42 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:58:42 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:58:42 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:58:42 INFO Jobs still running: ['5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:58:57 INFO None 5279427: status FINISHED
2026-07-26 11:58:57 INFO None 5279428: status FINISHED
2026-07-26 11:58:57 INFO None 5279430: status FINISHED
2026-07-26 11:58:57 INFO None 5279431: status FINISHED
2026-07-26 11:58:57 INFO None 5279432: status FINISHED
2026-07-26 11:58:59 INFO None 5279434: status FINISHED
2026-07-26 11:59:00 INFO None 5279435: status FINISHED
2026-07-26 11:59:00 INFO None 5279436: status FINISHED
2026-07-26 11:59:00 INFO None 5279438: status FINISHED
2026-07-26 11:59:00 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:59:00 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:59:00 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:59:00 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:59:00 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:59:00 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:59:00 INFO Jobs still running: ['5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:59:15 INFO None 5279427: status FINISHED
2026-07-26 11:59:15 INFO None 5279428: status FINISHED
2026-07-26 11:59:15 INFO None 5279430: status FINISHED
2026-07-26 11:59:15 INFO None 5279431: status FINISHED
2026-07-26 11:59:15 INFO None 5279432: status FINISHED
2026-07-26 11:59:15 INFO None 5279434: status FINISHED
2026-07-26 11:59:15 INFO None 5279435: status FINISHED
2026-07-26 11:59:15 INFO None 5279436: status FINISHED
2026-07-26 11:59:15 INFO None 5279438: status FINISHED
2026-07-26 11:59:15 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:59:15 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:59:15 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:59:15 INFO None 5279443: status RUNNING/PENDING
2026-07-26 11:59:15 INFO None 5279444: status RUNNING/PENDING
2026-07-26 11:59:15 INFO None 5279446: status RUNNING/PENDING
2026-07-26 11:59:15 INFO Jobs still running: ['5279439', '5279441', '5279442', '5279443', '5279444', '5279446']. Waiting...
2026-07-26 11:59:30 INFO None 5279427: status FINISHED
2026-07-26 11:59:30 INFO None 5279428: status FINISHED
2026-07-26 11:59:30 INFO None 5279430: status FINISHED
2026-07-26 11:59:30 INFO None 5279431: status FINISHED
2026-07-26 11:59:30 INFO None 5279432: status FINISHED
2026-07-26 11:59:30 INFO None 5279434: status FINISHED
2026-07-26 11:59:30 INFO None 5279435: status FINISHED
2026-07-26 11:59:30 INFO None 5279436: status FINISHED
2026-07-26 11:59:30 INFO None 5279438: status FINISHED
2026-07-26 11:59:30 INFO None 5279439: status RUNNING/PENDING
2026-07-26 11:59:30 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:59:30 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:59:30 INFO None 5279443: status FINISHED
2026-07-26 11:59:30 INFO None 5279444: status FINISHED
2026-07-26 11:59:30 INFO None 5279446: status FINISHED
2026-07-26 11:59:30 INFO Jobs still running: ['5279439', '5279441', '5279442']. Waiting...
2026-07-26 11:59:45 INFO None 5279427: status FINISHED
2026-07-26 11:59:45 INFO None 5279428: status FINISHED
2026-07-26 11:59:45 INFO None 5279430: status FINISHED
2026-07-26 11:59:45 INFO None 5279431: status FINISHED
2026-07-26 11:59:45 INFO None 5279432: status FINISHED
2026-07-26 11:59:45 INFO None 5279434: status FINISHED
2026-07-26 11:59:45 INFO None 5279435: status FINISHED
2026-07-26 11:59:45 INFO None 5279436: status FINISHED
2026-07-26 11:59:45 INFO None 5279438: status FINISHED
2026-07-26 11:59:45 INFO None 5279439: status FINISHED
2026-07-26 11:59:45 INFO None 5279441: status RUNNING/PENDING
2026-07-26 11:59:45 INFO None 5279442: status RUNNING/PENDING
2026-07-26 11:59:45 INFO None 5279443: status FINISHED
2026-07-26 11:59:45 INFO None 5279444: status FINISHED
2026-07-26 11:59:46 INFO None 5279446: status FINISHED
2026-07-26 11:59:46 INFO Jobs still running: ['5279441', '5279442']. Waiting...
2026-07-26 12:00:01 INFO None 5279427: status FINISHED
2026-07-26 12:00:01 INFO None 5279428: status FINISHED
2026-07-26 12:00:01 INFO None 5279430: status FINISHED
2026-07-26 12:00:01 INFO None 5279431: status FINISHED
2026-07-26 12:00:01 INFO None 5279432: status FINISHED
2026-07-26 12:00:01 INFO None 5279434: status FINISHED
2026-07-26 12:00:01 INFO None 5279435: status FINISHED
2026-07-26 12:00:01 INFO None 5279436: status FINISHED
2026-07-26 12:00:01 INFO None 5279438: status FINISHED
2026-07-26 12:00:01 INFO None 5279439: status FINISHED
2026-07-26 12:00:01 INFO None 5279441: status RUNNING/PENDING
2026-07-26 12:00:01 INFO None 5279442: status RUNNING/PENDING
2026-07-26 12:00:01 INFO None 5279443: status FINISHED
2026-07-26 12:00:01 INFO None 5279444: status FINISHED
2026-07-26 12:00:01 INFO None 5279446: status FINISHED
2026-07-26 12:00:01 INFO Jobs still running: ['5279441', '5279442']. Waiting...
2026-07-26 12:00:16 INFO None 5279427: status FINISHED
2026-07-26 12:00:16 INFO None 5279428: status FINISHED
2026-07-26 12:00:16 INFO None 5279430: status FINISHED
2026-07-26 12:00:16 INFO None 5279431: status FINISHED
2026-07-26 12:00:16 INFO None 5279432: status FINISHED
2026-07-26 12:00:16 INFO None 5279434: status FINISHED
2026-07-26 12:00:16 INFO None 5279435: status FINISHED
2026-07-26 12:00:16 INFO None 5279436: status FINISHED
2026-07-26 12:00:16 INFO None 5279438: status FINISHED
2026-07-26 12:00:16 INFO None 5279439: status FINISHED
2026-07-26 12:00:16 INFO None 5279441: status RUNNING/PENDING
2026-07-26 12:00:16 INFO None 5279442: status RUNNING/PENDING
2026-07-26 12:00:16 INFO None 5279443: status FINISHED
2026-07-26 12:00:16 INFO None 5279444: status FINISHED
2026-07-26 12:00:16 INFO None 5279446: status FINISHED
2026-07-26 12:00:16 INFO Jobs still running: ['5279441', '5279442']. Waiting...
2026-07-26 12:00:31 INFO None 5279427: status FINISHED
2026-07-26 12:00:31 INFO None 5279428: status FINISHED
2026-07-26 12:00:31 INFO None 5279430: status FINISHED
2026-07-26 12:00:31 INFO None 5279431: status FINISHED
2026-07-26 12:00:31 INFO None 5279432: status FINISHED
2026-07-26 12:00:31 INFO None 5279434: status FINISHED
2026-07-26 12:00:32 INFO None 5279435: status FINISHED
2026-07-26 12:00:32 INFO None 5279436: status FINISHED
2026-07-26 12:00:32 INFO None 5279438: status FINISHED
2026-07-26 12:00:32 INFO None 5279439: status FINISHED
2026-07-26 12:00:32 INFO None 5279441: status RUNNING/PENDING
2026-07-26 12:00:32 INFO None 5279442: status RUNNING/PENDING
2026-07-26 12:00:32 INFO None 5279443: status FINISHED
2026-07-26 12:00:34 INFO None 5279444: status FINISHED
2026-07-26 12:00:34 INFO None 5279446: status FINISHED
2026-07-26 12:00:34 INFO Jobs still running: ['5279441', '5279442']. Waiting...
2026-07-26 12:00:49 INFO None 5279427: status FINISHED
2026-07-26 12:00:49 INFO None 5279428: status FINISHED
2026-07-26 12:00:49 INFO None 5279430: status FINISHED
2026-07-26 12:00:49 INFO None 5279431: status FINISHED
2026-07-26 12:00:49 INFO None 5279432: status FINISHED
2026-07-26 12:00:49 INFO None 5279434: status FINISHED
2026-07-26 12:00:49 INFO None 5279435: status FINISHED
2026-07-26 12:00:49 INFO None 5279436: status FINISHED
2026-07-26 12:00:49 INFO None 5279438: status FINISHED
2026-07-26 12:00:49 INFO None 5279439: status FINISHED
2026-07-26 12:00:49 INFO None 5279441: status RUNNING/PENDING
2026-07-26 12:00:49 INFO None 5279442: status RUNNING/PENDING
2026-07-26 12:00:49 INFO None 5279443: status FINISHED
2026-07-26 12:00:49 INFO None 5279444: status FINISHED
2026-07-26 12:00:49 INFO None 5279446: status FINISHED
2026-07-26 12:00:49 INFO Jobs still running: ['5279441', '5279442']. Waiting...
2026-07-26 12:01:04 INFO None 5279427: status FINISHED
2026-07-26 12:01:04 INFO None 5279428: status FINISHED
2026-07-26 12:01:04 INFO None 5279430: status FINISHED
2026-07-26 12:01:04 INFO None 5279431: status FINISHED
2026-07-26 12:01:04 INFO None 5279432: status FINISHED
2026-07-26 12:01:04 INFO None 5279434: status FINISHED
2026-07-26 12:01:04 INFO None 5279435: status FINISHED
2026-07-26 12:01:04 INFO None 5279436: status FINISHED
2026-07-26 12:01:04 INFO None 5279438: status FINISHED
2026-07-26 12:01:04 INFO None 5279439: status FINISHED
2026-07-26 12:01:04 INFO None 5279441: status RUNNING/PENDING
2026-07-26 12:01:04 INFO None 5279442: status FINISHED
2026-07-26 12:01:04 INFO None 5279443: status FINISHED
2026-07-26 12:01:04 INFO None 5279444: status FINISHED
2026-07-26 12:01:04 INFO None 5279446: status FINISHED
2026-07-26 12:01:04 INFO Jobs still running: ['5279441']. Waiting...
2026-07-26 12:01:19 INFO None 5279427: status FINISHED
2026-07-26 12:01:19 INFO None 5279428: status FINISHED
2026-07-26 12:01:19 INFO None 5279430: status FINISHED
2026-07-26 12:01:21 INFO None 5279431: status FINISHED
2026-07-26 12:01:21 INFO None 5279432: status FINISHED
2026-07-26 12:01:21 INFO None 5279434: status FINISHED
2026-07-26 12:01:21 INFO None 5279435: status FINISHED
2026-07-26 12:01:21 INFO None 5279436: status FINISHED
2026-07-26 12:01:21 INFO None 5279438: status FINISHED
2026-07-26 12:01:21 INFO None 5279439: status FINISHED
2026-07-26 12:01:21 INFO None 5279441: status FINISHED
2026-07-26 12:01:21 INFO None 5279442: status FINISHED
2026-07-26 12:01:21 INFO None 5279443: status FINISHED
2026-07-26 12:01:21 INFO None 5279444: status FINISHED
2026-07-26 12:01:22 INFO None 5279446: status FINISHED
2026-07-26 12:01:22 INFO Jobs ['5279427', '5279428', '5279430', '5279431', '5279432', '5279434', '5279435', '5279436', '5279438', '5279439', '5279441', '5279442', '5279443', '5279444', '5279446'] have finished
2026-07-26 12:01:22 INFO Checking restart files were created ...
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-26 12:01:22 INFO  Run_model() completed successfully.
2026-07-26 12:01:22 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:01:22 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-26 12:01:22 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-26 12:01:22 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-26 12:01:22 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:01:22 INFO ---------->>> Running process_satellite_data()
2026-07-26 12:01:22 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-26 12:01:22 INFO ---------->>> Running run_obs_converter()
2026-07-26 12:01:22 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-26 12:01:22 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-26 12:01:22 INFO ---------->>> Running DART
2026-07-26 12:01:22 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-26 12:01:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-26 12:01:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-26 12:01:22 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-26 12:01:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-26 12:01:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-26 12:01:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-26 12:01:23 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-26 12:01:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-26 12:01:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-26 12:01:24 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-26 12:01:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-26 12:01:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-26 12:01:25 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-26 12:01:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-26 12:01:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-26 12:01:26 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-26 12:01:26 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-26 12:01:26 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-26 12:01:26 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-26 12:01:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-26 12:01:26 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-26 12:01:41 INFO Found: []
2026-07-26 12:01:41 INFO No job id returned by command ./run_filter.bsh
2026-07-26 12:01:41 INFO No monitoring will be performed
2026-07-26 12:01:41 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-26 12:01:41 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/analysis/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl_0607_15m_low_v2/preassim/2020020613'
2026-07-26 12:01:41 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-26 12:01:41 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-26 12:01:41 INFO run_dart() is DONE.
2026-07-26 12:01:41 INFO ---------->>> Running update_pollutant_in_end()
2026-07-26 12:01:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-26 12:01:46 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:01:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-26 12:01:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:01:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-26 12:01:57 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:01:57 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-26 12:02:02 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:02 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-26 12:02:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:08 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-26 12:02:13 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-26 12:02:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-26 12:02:24 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:24 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-26 12:02:29 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-26 12:02:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:35 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-26 12:02:40 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:40 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-26 12:02:45 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:46 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-26 12:02:51 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:51 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-26 12:02:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:02:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-26 12:03:01 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-26 12:03:01 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-26 12:03:01 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:03:01 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-26 12:03:02 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-26 12:03:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:03 INFO Hourly dataset computed and listing created
2026-07-26 12:03:05 INFO Hourly dataset computed
2026-07-26 12:03:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:06 INFO Hourly dataset computed and listing created
2026-07-26 12:03:07 INFO Hourly dataset computed
2026-07-26 12:03:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:07 INFO Hourly dataset computed and listing created
2026-07-26 12:03:08 INFO Hourly dataset computed
2026-07-26 12:03:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:09 INFO Hourly dataset computed and listing created
2026-07-26 12:03:09 INFO Hourly dataset computed
2026-07-26 12:03:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:10 INFO Hourly dataset computed and listing created
2026-07-26 12:03:11 INFO Hourly dataset computed
2026-07-26 12:03:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:12 INFO Hourly dataset computed and listing created
2026-07-26 12:03:12 INFO Hourly dataset computed
2026-07-26 12:03:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:13 INFO Hourly dataset computed and listing created
2026-07-26 12:03:14 INFO Hourly dataset computed
2026-07-26 12:03:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:15 INFO Hourly dataset computed and listing created
2026-07-26 12:03:15 INFO Hourly dataset computed
2026-07-26 12:03:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:16 INFO Hourly dataset computed and listing created
2026-07-26 12:03:17 INFO Hourly dataset computed
2026-07-26 12:03:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:18 INFO Hourly dataset computed and listing created
2026-07-26 12:03:18 INFO Hourly dataset computed
2026-07-26 12:03:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:19 INFO Hourly dataset computed and listing created
2026-07-26 12:03:20 INFO Hourly dataset computed
2026-07-26 12:03:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:21 INFO Hourly dataset computed and listing created
2026-07-26 12:03:21 INFO Hourly dataset computed
2026-07-26 12:03:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:22 INFO Hourly dataset computed and listing created
2026-07-26 12:03:23 INFO Hourly dataset computed
2026-07-26 12:03:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:24 INFO Hourly dataset computed and listing created
2026-07-26 12:03:24 INFO Hourly dataset computed
2026-07-26 12:03:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-26 12:03:25 INFO Hourly dataset computed and listing created
2026-07-26 12:03:26 INFO Hourly dataset computed
2026-07-26 12:03:26 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-26 12:03:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1
2026-07-26 12:03:26 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-26 12:03:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-26 12:03:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:26 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-26 12:03:26 INFO Queuing job for member 1...
2026-07-26 12:03:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:26 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-26 12:03:27 INFO Found: ['5279472']
2026-07-26 12:03:32 INFO [TGCC-IRENE] Submitted job with ID:['5279472']
2026-07-26 12:03:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2
2026-07-26 12:03:32 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-26 12:03:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-26 12:03:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:32 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-26 12:03:32 INFO Queuing job for member 2...
2026-07-26 12:03:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:32 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-26 12:03:35 INFO Found: ['5279474']
2026-07-26 12:03:40 INFO [TGCC-IRENE] Submitted job with ID:['5279474']
2026-07-26 12:03:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3
2026-07-26 12:03:40 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-26 12:03:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-26 12:03:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:40 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-26 12:03:40 INFO Queuing job for member 3...
2026-07-26 12:03:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:40 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-26 12:03:42 INFO Found: ['5279475']
2026-07-26 12:03:47 INFO [TGCC-IRENE] Submitted job with ID:['5279475']
2026-07-26 12:03:47 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:47 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4
2026-07-26 12:03:47 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-26 12:03:47 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-26 12:03:47 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:47 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-26 12:03:47 INFO Queuing job for member 4...
2026-07-26 12:03:47 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:47 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-26 12:03:48 INFO Found: ['5279477']
2026-07-26 12:03:53 INFO [TGCC-IRENE] Submitted job with ID:['5279477']
2026-07-26 12:03:53 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:53 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5
2026-07-26 12:03:53 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-26 12:03:53 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-26 12:03:53 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:53 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-26 12:03:53 INFO Queuing job for member 5...
2026-07-26 12:03:53 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:53 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-26 12:03:54 INFO Found: ['5279478']
2026-07-26 12:03:59 INFO [TGCC-IRENE] Submitted job with ID:['5279478']
2026-07-26 12:03:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:03:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6
2026-07-26 12:03:59 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-26 12:03:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-26 12:03:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:03:59 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-26 12:03:59 INFO Queuing job for member 6...
2026-07-26 12:03:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:03:59 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-26 12:04:00 INFO Found: ['5279480']
2026-07-26 12:04:05 INFO [TGCC-IRENE] Submitted job with ID:['5279480']
2026-07-26 12:04:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7
2026-07-26 12:04:05 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-26 12:04:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-26 12:04:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:05 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-26 12:04:05 INFO Queuing job for member 7...
2026-07-26 12:04:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:05 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-26 12:04:05 INFO Found: ['5279482']
2026-07-26 12:04:10 INFO [TGCC-IRENE] Submitted job with ID:['5279482']
2026-07-26 12:04:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8
2026-07-26 12:04:10 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-26 12:04:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-26 12:04:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:10 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-26 12:04:10 INFO Queuing job for member 8...
2026-07-26 12:04:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:10 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-26 12:04:11 INFO Found: ['5279483']
2026-07-26 12:04:16 INFO [TGCC-IRENE] Submitted job with ID:['5279483']
2026-07-26 12:04:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9
2026-07-26 12:04:16 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-26 12:04:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-26 12:04:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:16 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-26 12:04:16 INFO Queuing job for member 9...
2026-07-26 12:04:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:16 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-26 12:04:17 INFO Found: ['5279485']
2026-07-26 12:04:22 INFO [TGCC-IRENE] Submitted job with ID:['5279485']
2026-07-26 12:04:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10
2026-07-26 12:04:22 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-26 12:04:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-26 12:04:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:22 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-26 12:04:22 INFO Queuing job for member 10...
2026-07-26 12:04:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:22 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-26 12:04:23 INFO Found: ['5279486']
2026-07-26 12:04:28 INFO [TGCC-IRENE] Submitted job with ID:['5279486']
2026-07-26 12:04:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11
2026-07-26 12:04:28 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-26 12:04:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-26 12:04:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:28 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-26 12:04:28 INFO Queuing job for member 11...
2026-07-26 12:04:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:28 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-26 12:04:28 INFO Found: ['5279487']
2026-07-26 12:04:33 INFO [TGCC-IRENE] Submitted job with ID:['5279487']
2026-07-26 12:04:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12
2026-07-26 12:04:33 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-26 12:04:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-26 12:04:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:33 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-26 12:04:33 INFO Queuing job for member 12...
2026-07-26 12:04:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:33 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-26 12:04:34 INFO Found: ['5279489']
2026-07-26 12:04:39 INFO [TGCC-IRENE] Submitted job with ID:['5279489']
2026-07-26 12:04:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13
2026-07-26 12:04:39 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-26 12:04:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-26 12:04:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:39 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-26 12:04:39 INFO Queuing job for member 13...
2026-07-26 12:04:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:39 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-26 12:04:41 INFO Found: ['5279491']
2026-07-26 12:04:46 INFO [TGCC-IRENE] Submitted job with ID:['5279491']
2026-07-26 12:04:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14
2026-07-26 12:04:46 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-26 12:04:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-26 12:04:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:46 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-26 12:04:46 INFO Queuing job for member 14...
2026-07-26 12:04:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:46 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-26 12:04:49 INFO Found: ['5279493']
2026-07-26 12:04:54 INFO [TGCC-IRENE] Submitted job with ID:['5279493']
2026-07-26 12:04:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-26 12:04:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15
2026-07-26 12:04:54 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-26 12:04:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-26 12:04:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-26 12:04:54 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-26 12:04:54 INFO Queuing job for member 15...
2026-07-26 12:04:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-26 12:04:54 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-26 12:04:56 INFO Found: ['5279494']
2026-07-26 12:05:01 INFO [TGCC-IRENE] Submitted job with ID:['5279494']
2026-07-26 12:05:01 INFO Checking job status ...
2026-07-26 12:05:01 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:05:01 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:05:02 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:05:02 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:05:02 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:05:02 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:05:02 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:05:02 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:05:17 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:05:17 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:05:19 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:05:19 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:05:34 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:05:34 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:05:34 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:05:49 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:05:49 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:05:51 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:05:52 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:05:52 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:05:52 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:06:07 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:06:07 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:06:07 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:06:22 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:06:22 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:06:22 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:06:37 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:06:37 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:06:38 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:06:38 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:06:38 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:06:54 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:06:54 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:06:54 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:07:09 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:07:09 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:07:10 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:07:10 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:07:10 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:07:10 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:07:10 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:07:12 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:07:12 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:07:27 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:07:27 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:07:27 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:07:42 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279475: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:07:42 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:07:44 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:07:44 INFO Jobs still running: ['5279472', '5279474', '5279475', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:07:59 INFO None 5279472: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279474: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279475: status FINISHED
2026-07-26 12:07:59 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:07:59 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:08:00 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:08:00 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:08:00 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:08:00 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:08:00 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:08:00 INFO Jobs still running: ['5279472', '5279474', '5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:08:15 INFO None 5279472: status FINISHED
2026-07-26 12:08:15 INFO None 5279474: status FINISHED
2026-07-26 12:08:15 INFO None 5279475: status FINISHED
2026-07-26 12:08:15 INFO None 5279477: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279478: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279480: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279482: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279483: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279485: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:08:15 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:08:15 INFO Jobs still running: ['5279477', '5279478', '5279480', '5279482', '5279483', '5279485', '5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:08:30 INFO None 5279472: status FINISHED
2026-07-26 12:08:30 INFO None 5279474: status FINISHED
2026-07-26 12:08:30 INFO None 5279475: status FINISHED
2026-07-26 12:08:30 INFO None 5279477: status FINISHED
2026-07-26 12:08:30 INFO None 5279478: status FINISHED
2026-07-26 12:08:30 INFO None 5279480: status FINISHED
2026-07-26 12:08:30 INFO None 5279482: status FINISHED
2026-07-26 12:08:30 INFO None 5279483: status FINISHED
2026-07-26 12:08:30 INFO None 5279485: status FINISHED
2026-07-26 12:08:30 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:08:30 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:08:30 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:08:30 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:08:30 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:08:30 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:08:30 INFO Jobs still running: ['5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:08:45 INFO None 5279472: status FINISHED
2026-07-26 12:08:45 INFO None 5279474: status FINISHED
2026-07-26 12:08:45 INFO None 5279475: status FINISHED
2026-07-26 12:08:45 INFO None 5279477: status FINISHED
2026-07-26 12:08:45 INFO None 5279478: status FINISHED
2026-07-26 12:08:45 INFO None 5279480: status FINISHED
2026-07-26 12:08:45 INFO None 5279482: status FINISHED
2026-07-26 12:08:45 INFO None 5279483: status FINISHED
2026-07-26 12:08:45 INFO None 5279485: status FINISHED
2026-07-26 12:08:45 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:08:45 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:08:45 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:08:45 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:08:45 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:08:45 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:08:45 INFO Jobs still running: ['5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:09:01 INFO None 5279472: status FINISHED
2026-07-26 12:09:01 INFO None 5279474: status FINISHED
2026-07-26 12:09:01 INFO None 5279475: status FINISHED
2026-07-26 12:09:01 INFO None 5279477: status FINISHED
2026-07-26 12:09:01 INFO None 5279478: status FINISHED
2026-07-26 12:09:01 INFO None 5279480: status FINISHED
2026-07-26 12:09:01 INFO None 5279482: status FINISHED
2026-07-26 12:09:01 INFO None 5279483: status FINISHED
2026-07-26 12:09:01 INFO None 5279485: status FINISHED
2026-07-26 12:09:01 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:09:01 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:09:01 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:09:01 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:09:01 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:09:01 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:09:01 INFO Jobs still running: ['5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
2026-07-26 12:09:16 INFO None 5279472: status FINISHED
2026-07-26 12:09:16 INFO None 5279474: status FINISHED
2026-07-26 12:09:16 INFO None 5279475: status FINISHED
2026-07-26 12:09:16 INFO None 5279477: status FINISHED
2026-07-26 12:09:16 INFO None 5279478: status FINISHED
2026-07-26 12:09:16 INFO None 5279480: status FINISHED
2026-07-26 12:09:16 INFO None 5279482: status FINISHED
2026-07-26 12:09:16 INFO None 5279483: status FINISHED
2026-07-26 12:09:16 INFO None 5279485: status FINISHED
2026-07-26 12:09:16 INFO None 5279486: status RUNNING/PENDING
2026-07-26 12:09:17 INFO None 5279487: status RUNNING/PENDING
2026-07-26 12:09:17 INFO None 5279489: status RUNNING/PENDING
2026-07-26 12:09:17 INFO None 5279491: status RUNNING/PENDING
2026-07-26 12:09:17 INFO None 5279493: status RUNNING/PENDING
2026-07-26 12:09:17 INFO None 5279494: status RUNNING/PENDING
2026-07-26 12:09:17 INFO Jobs still running: ['5279486', '5279487', '5279489', '5279491', '5279493', '5279494']. Waiting...
[2026-07-26T12:09:29.378] error: *** JOB 5279378 ON irene4135 CANCELLED AT 2026-07-26T12:09:29 DUE to SIGNAL Terminated ***
