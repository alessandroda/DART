+ SCRIPT_PID=542090
+ /bin/bash -x /tmp/tmp.A2JwCCnYg5
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
2026-07-27 13:21:34 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-27 13:21:34 INFO [PIPELINE] =======================================
2026-07-27 13:21:34 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-27 13:21:34 INFO [PIPELINE] Config file: config/config_irene_IM_cp2.yaml
2026-07-27 13:21:34 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-27 13:21:34 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260727_132134.log
2026-07-27 13:21:34 INFO [PIPELINE] =======================================
2026-07-27 13:21:34 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-27 13:21:34 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-27 13:21:34 INFO [STEP] ---- TIME LOOP START ----
2026-07-27 13:21:34 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:21:34 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-27 13:21:34 INFO Asked to restart from control run ...
2026-07-27 13:21:34 INFO Copying EMIS ...
2026-07-27 13:21:34 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-27 13:21:34 INFO Linking first END ...
2026-07-27 13:21:34 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:34 INFO >> Checking links...
2026-07-27 13:21:34 INFO >> All links are good for ENS1  ...
2026-07-27 13:21:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:21:42 INFO Hourly dataset computed and listing created
2026-07-27 13:21:47 INFO Hourly dataset computed
2026-07-27 13:21:47 INFO Asked to restart from control run ...
2026-07-27 13:21:47 INFO Copying EMIS ...
2026-07-27 13:21:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-27 13:21:47 INFO Linking first END ...
2026-07-27 13:21:47 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:47 INFO >> Checking links...
2026-07-27 13:21:47 INFO >> All links are good for ENS2  ...
2026-07-27 13:21:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:21:48 INFO Hourly dataset computed and listing created
2026-07-27 13:21:50 INFO Hourly dataset computed
2026-07-27 13:21:50 INFO Asked to restart from control run ...
2026-07-27 13:21:50 INFO Copying EMIS ...
2026-07-27 13:21:50 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-27 13:21:50 INFO Linking first END ...
2026-07-27 13:21:50 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:50 INFO >> Checking links...
2026-07-27 13:21:50 INFO >> All links are good for ENS3  ...
2026-07-27 13:21:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:21:51 INFO Hourly dataset computed and listing created
2026-07-27 13:21:52 INFO Hourly dataset computed
2026-07-27 13:21:53 INFO Asked to restart from control run ...
2026-07-27 13:21:53 INFO Copying EMIS ...
2026-07-27 13:21:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-27 13:21:53 INFO Linking first END ...
2026-07-27 13:21:53 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:53 INFO >> Checking links...
2026-07-27 13:21:53 INFO >> All links are good for ENS4  ...
2026-07-27 13:21:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:21:54 INFO Hourly dataset computed and listing created
2026-07-27 13:21:56 INFO Hourly dataset computed
2026-07-27 13:21:56 INFO Asked to restart from control run ...
2026-07-27 13:21:56 INFO Copying EMIS ...
2026-07-27 13:21:56 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-27 13:21:56 INFO Linking first END ...
2026-07-27 13:21:56 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:56 INFO >> Checking links...
2026-07-27 13:21:56 INFO >> All links are good for ENS5  ...
2026-07-27 13:21:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:21:57 INFO Hourly dataset computed and listing created
2026-07-27 13:21:59 INFO Hourly dataset computed
2026-07-27 13:21:59 INFO Asked to restart from control run ...
2026-07-27 13:21:59 INFO Copying EMIS ...
2026-07-27 13:21:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-27 13:21:59 INFO Linking first END ...
2026-07-27 13:21:59 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:21:59 INFO >> Checking links...
2026-07-27 13:22:00 INFO >> All links are good for ENS6  ...
2026-07-27 13:22:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:01 INFO Hourly dataset computed and listing created
2026-07-27 13:22:02 INFO Hourly dataset computed
2026-07-27 13:22:02 INFO Asked to restart from control run ...
2026-07-27 13:22:02 INFO Copying EMIS ...
2026-07-27 13:22:02 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-27 13:22:02 INFO Linking first END ...
2026-07-27 13:22:02 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:02 INFO >> Checking links...
2026-07-27 13:22:02 INFO >> All links are good for ENS7  ...
2026-07-27 13:22:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:03 INFO Hourly dataset computed and listing created
2026-07-27 13:22:05 INFO Hourly dataset computed
2026-07-27 13:22:05 INFO Asked to restart from control run ...
2026-07-27 13:22:05 INFO Copying EMIS ...
2026-07-27 13:22:05 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-27 13:22:05 INFO Linking first END ...
2026-07-27 13:22:05 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:05 INFO >> Checking links...
2026-07-27 13:22:05 INFO >> All links are good for ENS8  ...
2026-07-27 13:22:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:06 INFO Hourly dataset computed and listing created
2026-07-27 13:22:08 INFO Hourly dataset computed
2026-07-27 13:22:08 INFO Asked to restart from control run ...
2026-07-27 13:22:08 INFO Copying EMIS ...
2026-07-27 13:22:08 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-27 13:22:08 INFO Linking first END ...
2026-07-27 13:22:08 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:08 INFO >> Checking links...
2026-07-27 13:22:08 INFO >> All links are good for ENS9  ...
2026-07-27 13:22:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:09 INFO Hourly dataset computed and listing created
2026-07-27 13:22:11 INFO Hourly dataset computed
2026-07-27 13:22:11 INFO Asked to restart from control run ...
2026-07-27 13:22:11 INFO Copying EMIS ...
2026-07-27 13:22:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-27 13:22:11 INFO Linking first END ...
2026-07-27 13:22:11 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:11 INFO >> Checking links...
2026-07-27 13:22:11 INFO >> All links are good for ENS10  ...
2026-07-27 13:22:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:12 INFO Hourly dataset computed and listing created
2026-07-27 13:22:13 INFO Hourly dataset computed
2026-07-27 13:22:14 INFO Asked to restart from control run ...
2026-07-27 13:22:14 INFO Copying EMIS ...
2026-07-27 13:22:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-27 13:22:14 INFO Linking first END ...
2026-07-27 13:22:14 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:14 INFO >> Checking links...
2026-07-27 13:22:14 INFO >> All links are good for ENS11  ...
2026-07-27 13:22:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:15 INFO Hourly dataset computed and listing created
2026-07-27 13:22:16 INFO Hourly dataset computed
2026-07-27 13:22:16 INFO Asked to restart from control run ...
2026-07-27 13:22:16 INFO Copying EMIS ...
2026-07-27 13:22:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-27 13:22:17 INFO Linking first END ...
2026-07-27 13:22:17 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:17 INFO >> Checking links...
2026-07-27 13:22:17 INFO >> All links are good for ENS12  ...
2026-07-27 13:22:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:18 INFO Hourly dataset computed and listing created
2026-07-27 13:22:19 INFO Hourly dataset computed
2026-07-27 13:22:19 INFO Asked to restart from control run ...
2026-07-27 13:22:19 INFO Copying EMIS ...
2026-07-27 13:22:20 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-27 13:22:20 INFO Linking first END ...
2026-07-27 13:22:20 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:20 INFO >> Checking links...
2026-07-27 13:22:20 INFO >> All links are good for ENS13  ...
2026-07-27 13:22:20 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:21 INFO Hourly dataset computed and listing created
2026-07-27 13:22:22 INFO Hourly dataset computed
2026-07-27 13:22:22 INFO Asked to restart from control run ...
2026-07-27 13:22:22 INFO Copying EMIS ...
2026-07-27 13:22:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-27 13:22:23 INFO Linking first END ...
2026-07-27 13:22:23 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:23 INFO >> Checking links...
2026-07-27 13:22:23 INFO >> All links are good for ENS14  ...
2026-07-27 13:22:23 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:24 INFO Hourly dataset computed and listing created
2026-07-27 13:22:25 INFO Hourly dataset computed
2026-07-27 13:22:25 INFO Asked to restart from control run ...
2026-07-27 13:22:25 INFO Copying EMIS ...
2026-07-27 13:22:25 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-27 13:22:25 INFO Linking first END ...
2026-07-27 13:22:25 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-27 13:22:25 INFO >> Checking links...
2026-07-27 13:22:25 INFO >> All links are good for ENS15  ...
2026-07-27 13:22:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:22:26 INFO Hourly dataset computed and listing created
2026-07-27 13:22:28 INFO Hourly dataset computed
2026-07-27 13:22:28 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-27 13:22:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:22:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 13:22:28 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-27 13:22:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 13:22:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:22:28 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 13:22:28 INFO Queuing job for member 1...
2026-07-27 13:22:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:22:28 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 13:22:29 INFO Found: ['5285440']
2026-07-27 13:22:34 INFO [TGCC-IRENE] Submitted job with ID:['5285440']
2026-07-27 13:22:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:22:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 13:22:34 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-27 13:22:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 13:22:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:22:34 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 13:22:34 INFO Queuing job for member 2...
2026-07-27 13:22:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:22:34 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 13:22:35 INFO Found: ['5285441']
2026-07-27 13:22:40 INFO [TGCC-IRENE] Submitted job with ID:['5285441']
2026-07-27 13:22:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:22:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 13:22:40 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-27 13:22:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 13:22:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:22:40 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 13:22:40 INFO Queuing job for member 3...
2026-07-27 13:22:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:22:40 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 13:22:43 INFO Found: ['5285442']
2026-07-27 13:22:48 INFO [TGCC-IRENE] Submitted job with ID:['5285442']
2026-07-27 13:22:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:22:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 13:22:48 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-27 13:22:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 13:22:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:22:48 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 13:22:48 INFO Queuing job for member 4...
2026-07-27 13:22:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:22:48 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 13:22:50 INFO Found: ['5285443']
2026-07-27 13:22:55 INFO [TGCC-IRENE] Submitted job with ID:['5285443']
2026-07-27 13:22:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:22:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 13:22:55 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-27 13:22:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 13:22:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:22:55 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 13:22:55 INFO Queuing job for member 5...
2026-07-27 13:22:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:22:55 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 13:22:58 INFO Found: ['5285444']
2026-07-27 13:23:03 INFO [TGCC-IRENE] Submitted job with ID:['5285444']
2026-07-27 13:23:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 13:23:03 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-27 13:23:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 13:23:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:03 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 13:23:03 INFO Queuing job for member 6...
2026-07-27 13:23:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:03 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 13:23:05 INFO Found: ['5285446']
2026-07-27 13:23:10 INFO [TGCC-IRENE] Submitted job with ID:['5285446']
2026-07-27 13:23:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 13:23:10 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-27 13:23:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 13:23:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:10 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 13:23:11 INFO Queuing job for member 7...
2026-07-27 13:23:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:11 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 13:23:13 INFO Found: ['5285447']
2026-07-27 13:23:18 INFO [TGCC-IRENE] Submitted job with ID:['5285447']
2026-07-27 13:23:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 13:23:18 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-27 13:23:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 13:23:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:18 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 13:23:18 INFO Queuing job for member 8...
2026-07-27 13:23:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:18 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 13:23:21 INFO Found: ['5285448']
2026-07-27 13:23:26 INFO [TGCC-IRENE] Submitted job with ID:['5285448']
2026-07-27 13:23:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 13:23:26 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-27 13:23:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 13:23:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:26 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 13:23:26 INFO Queuing job for member 9...
2026-07-27 13:23:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:26 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 13:23:26 INFO Found: ['5285449']
2026-07-27 13:23:31 INFO [TGCC-IRENE] Submitted job with ID:['5285449']
2026-07-27 13:23:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 13:23:31 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-27 13:23:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 13:23:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:31 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 13:23:31 INFO Queuing job for member 10...
2026-07-27 13:23:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:31 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 13:23:32 INFO Found: ['5285450']
2026-07-27 13:23:37 INFO [TGCC-IRENE] Submitted job with ID:['5285450']
2026-07-27 13:23:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 13:23:37 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-27 13:23:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 13:23:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:37 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 13:23:37 INFO Queuing job for member 11...
2026-07-27 13:23:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:37 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 13:23:38 INFO Found: ['5285452']
2026-07-27 13:23:43 INFO [TGCC-IRENE] Submitted job with ID:['5285452']
2026-07-27 13:23:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 13:23:43 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-27 13:23:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 13:23:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:43 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 13:23:43 INFO Queuing job for member 12...
2026-07-27 13:23:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:43 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 13:23:43 INFO Found: ['5285453']
2026-07-27 13:23:48 INFO [TGCC-IRENE] Submitted job with ID:['5285453']
2026-07-27 13:23:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 13:23:48 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-27 13:23:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 13:23:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:48 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 13:23:48 INFO Queuing job for member 13...
2026-07-27 13:23:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:48 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 13:23:49 INFO Found: ['5285454']
2026-07-27 13:23:54 INFO [TGCC-IRENE] Submitted job with ID:['5285454']
2026-07-27 13:23:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:23:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 13:23:54 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-27 13:23:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 13:23:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:23:54 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 13:23:54 INFO Queuing job for member 14...
2026-07-27 13:23:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:23:54 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 13:23:55 INFO Found: ['5285455']
2026-07-27 13:24:00 INFO [TGCC-IRENE] Submitted job with ID:['5285455']
2026-07-27 13:24:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:24:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 13:24:00 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-27 13:24:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 13:24:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:24:00 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 13:24:00 INFO Queuing job for member 15...
2026-07-27 13:24:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:24:00 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 13:24:01 INFO Found: ['5285456']
2026-07-27 13:24:06 INFO [TGCC-IRENE] Submitted job with ID:['5285456']
2026-07-27 13:24:06 INFO Checking job status ...
2026-07-27 13:24:06 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:24:06 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:24:06 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:24:22 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:24:22 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:24:23 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:24:23 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:24:23 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:24:23 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:24:23 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:24:23 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:24:38 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:24:38 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:24:40 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:24:40 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:24:40 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:24:55 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:24:55 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:24:55 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:25:10 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:25:10 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:25:11 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:25:11 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:25:11 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:25:11 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:25:26 INFO None 5285440: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:25:26 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:25:26 INFO Jobs still running: ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:25:41 INFO None 5285440: status FINISHED
2026-07-27 13:25:41 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:25:41 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:25:41 INFO Jobs still running: ['5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:25:57 INFO None 5285440: status FINISHED
2026-07-27 13:25:57 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:25:57 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:25:57 INFO Jobs still running: ['5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:26:12 INFO None 5285440: status FINISHED
2026-07-27 13:26:12 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:26:12 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:26:14 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:26:14 INFO Jobs still running: ['5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:26:30 INFO None 5285440: status FINISHED
2026-07-27 13:26:30 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285442: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:26:30 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:26:30 INFO Jobs still running: ['5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:26:45 INFO None 5285440: status FINISHED
2026-07-27 13:26:45 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285442: status FINISHED
2026-07-27 13:26:45 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285449: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285452: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285453: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:26:45 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:26:45 INFO Jobs still running: ['5285441', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:27:00 INFO None 5285440: status FINISHED
2026-07-27 13:27:00 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285442: status FINISHED
2026-07-27 13:27:00 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285449: status FINISHED
2026-07-27 13:27:00 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285452: status FINISHED
2026-07-27 13:27:00 INFO None 5285453: status FINISHED
2026-07-27 13:27:00 INFO None 5285454: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:27:00 INFO None 5285456: status RUNNING/PENDING
2026-07-27 13:27:00 INFO Jobs still running: ['5285441', '5285443', '5285444', '5285446', '5285447', '5285448', '5285450', '5285454', '5285455', '5285456']. Waiting...
2026-07-27 13:27:15 INFO None 5285440: status FINISHED
2026-07-27 13:27:15 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285442: status FINISHED
2026-07-27 13:27:16 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285447: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285449: status FINISHED
2026-07-27 13:27:16 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285452: status FINISHED
2026-07-27 13:27:16 INFO None 5285453: status FINISHED
2026-07-27 13:27:16 INFO None 5285454: status FINISHED
2026-07-27 13:27:16 INFO None 5285455: status RUNNING/PENDING
2026-07-27 13:27:16 INFO None 5285456: status FINISHED
2026-07-27 13:27:16 INFO Jobs still running: ['5285441', '5285443', '5285444', '5285446', '5285447', '5285448', '5285450', '5285455']. Waiting...
2026-07-27 13:27:31 INFO None 5285440: status FINISHED
2026-07-27 13:27:31 INFO None 5285441: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285442: status FINISHED
2026-07-27 13:27:31 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285447: status FINISHED
2026-07-27 13:27:31 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285449: status FINISHED
2026-07-27 13:27:31 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:27:31 INFO None 5285452: status FINISHED
2026-07-27 13:27:31 INFO None 5285453: status FINISHED
2026-07-27 13:27:31 INFO None 5285454: status FINISHED
2026-07-27 13:27:31 INFO None 5285455: status FINISHED
2026-07-27 13:27:31 INFO None 5285456: status FINISHED
2026-07-27 13:27:31 INFO Jobs still running: ['5285441', '5285443', '5285444', '5285446', '5285448', '5285450']. Waiting...
2026-07-27 13:27:46 INFO None 5285440: status FINISHED
2026-07-27 13:27:46 INFO None 5285441: status FINISHED
2026-07-27 13:27:46 INFO None 5285442: status FINISHED
2026-07-27 13:27:46 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:27:46 INFO None 5285444: status RUNNING/PENDING
2026-07-27 13:27:46 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:27:46 INFO None 5285447: status FINISHED
2026-07-27 13:27:46 INFO None 5285448: status RUNNING/PENDING
2026-07-27 13:27:46 INFO None 5285449: status FINISHED
2026-07-27 13:27:46 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:27:46 INFO None 5285452: status FINISHED
2026-07-27 13:27:46 INFO None 5285453: status FINISHED
2026-07-27 13:27:46 INFO None 5285454: status FINISHED
2026-07-27 13:27:46 INFO None 5285455: status FINISHED
2026-07-27 13:27:46 INFO None 5285456: status FINISHED
2026-07-27 13:27:46 INFO Jobs still running: ['5285443', '5285444', '5285446', '5285448', '5285450']. Waiting...
2026-07-27 13:28:01 INFO None 5285440: status FINISHED
2026-07-27 13:28:01 INFO None 5285441: status FINISHED
2026-07-27 13:28:01 INFO None 5285442: status FINISHED
2026-07-27 13:28:01 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:28:04 INFO None 5285444: status FINISHED
2026-07-27 13:28:04 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:28:04 INFO None 5285447: status FINISHED
2026-07-27 13:28:04 INFO None 5285448: status FINISHED
2026-07-27 13:28:04 INFO None 5285449: status FINISHED
2026-07-27 13:28:04 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:28:04 INFO None 5285452: status FINISHED
2026-07-27 13:28:04 INFO None 5285453: status FINISHED
2026-07-27 13:28:04 INFO None 5285454: status FINISHED
2026-07-27 13:28:04 INFO None 5285455: status FINISHED
2026-07-27 13:28:04 INFO None 5285456: status FINISHED
2026-07-27 13:28:04 INFO Jobs still running: ['5285443', '5285446', '5285450']. Waiting...
2026-07-27 13:28:19 INFO None 5285440: status FINISHED
2026-07-27 13:28:19 INFO None 5285441: status FINISHED
2026-07-27 13:28:19 INFO None 5285442: status FINISHED
2026-07-27 13:28:19 INFO None 5285443: status RUNNING/PENDING
2026-07-27 13:28:19 INFO None 5285444: status FINISHED
2026-07-27 13:28:19 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:28:19 INFO None 5285447: status FINISHED
2026-07-27 13:28:19 INFO None 5285448: status FINISHED
2026-07-27 13:28:19 INFO None 5285449: status FINISHED
2026-07-27 13:28:19 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:28:19 INFO None 5285452: status FINISHED
2026-07-27 13:28:19 INFO None 5285453: status FINISHED
2026-07-27 13:28:19 INFO None 5285454: status FINISHED
2026-07-27 13:28:19 INFO None 5285455: status FINISHED
2026-07-27 13:28:19 INFO None 5285456: status FINISHED
2026-07-27 13:28:19 INFO Jobs still running: ['5285443', '5285446', '5285450']. Waiting...
2026-07-27 13:28:34 INFO None 5285440: status FINISHED
2026-07-27 13:28:34 INFO None 5285441: status FINISHED
2026-07-27 13:28:34 INFO None 5285442: status FINISHED
2026-07-27 13:28:34 INFO None 5285443: status FINISHED
2026-07-27 13:28:34 INFO None 5285444: status FINISHED
2026-07-27 13:28:34 INFO None 5285446: status RUNNING/PENDING
2026-07-27 13:28:34 INFO None 5285447: status FINISHED
2026-07-27 13:28:34 INFO None 5285448: status FINISHED
2026-07-27 13:28:34 INFO None 5285449: status FINISHED
2026-07-27 13:28:34 INFO None 5285450: status RUNNING/PENDING
2026-07-27 13:28:34 INFO None 5285452: status FINISHED
2026-07-27 13:28:34 INFO None 5285453: status FINISHED
2026-07-27 13:28:34 INFO None 5285454: status FINISHED
2026-07-27 13:28:34 INFO None 5285455: status FINISHED
2026-07-27 13:28:34 INFO None 5285456: status FINISHED
2026-07-27 13:28:34 INFO Jobs still running: ['5285446', '5285450']. Waiting...
2026-07-27 13:28:49 INFO None 5285440: status FINISHED
2026-07-27 13:28:49 INFO None 5285441: status FINISHED
2026-07-27 13:28:49 INFO None 5285442: status FINISHED
2026-07-27 13:28:49 INFO None 5285443: status FINISHED
2026-07-27 13:28:49 INFO None 5285444: status FINISHED
2026-07-27 13:28:49 INFO None 5285446: status FINISHED
2026-07-27 13:28:49 INFO None 5285447: status FINISHED
2026-07-27 13:28:49 INFO None 5285448: status FINISHED
2026-07-27 13:28:49 INFO None 5285449: status FINISHED
2026-07-27 13:28:49 INFO None 5285450: status FINISHED
2026-07-27 13:28:50 INFO None 5285452: status FINISHED
2026-07-27 13:28:50 INFO None 5285453: status FINISHED
2026-07-27 13:28:50 INFO None 5285454: status FINISHED
2026-07-27 13:28:50 INFO None 5285455: status FINISHED
2026-07-27 13:28:50 INFO None 5285456: status FINISHED
2026-07-27 13:28:50 INFO Jobs ['5285440', '5285441', '5285442', '5285443', '5285444', '5285446', '5285447', '5285448', '5285449', '5285450', '5285452', '5285453', '5285454', '5285455', '5285456'] have finished
2026-07-27 13:28:50 INFO Checking restart files were created ...
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-27 13:28:50 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-27 13:28:50 INFO  Run_model() completed successfully.
2026-07-27 13:28:50 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:28:50 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-27 13:28:50 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 13:28:50 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-27 13:28:50 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:28:50 INFO ---------->>> Running process_satellite_data()
2026-07-27 13:28:50 INFO [DART] No satellite data found, skipping assimilation
2026-07-27 13:28:50 INFO after_assimilation() skipped
2026-07-27 13:28:50 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 13:28:50 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:28:50 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:28:50 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-27 13:28:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:28:51 INFO Hourly dataset computed and listing created
2026-07-27 13:29:09 INFO Hourly dataset computed
2026-07-27 13:29:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:29:10 INFO Hourly dataset computed and listing created
2026-07-27 13:29:22 INFO Hourly dataset computed
2026-07-27 13:29:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:29:24 INFO Hourly dataset computed and listing created
2026-07-27 13:29:37 INFO Hourly dataset computed
2026-07-27 13:29:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:29:39 INFO Hourly dataset computed and listing created
2026-07-27 13:29:51 INFO Hourly dataset computed
2026-07-27 13:29:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:29:53 INFO Hourly dataset computed and listing created
2026-07-27 13:30:04 INFO Hourly dataset computed
2026-07-27 13:30:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:30:06 INFO Hourly dataset computed and listing created
2026-07-27 13:30:21 INFO Hourly dataset computed
2026-07-27 13:30:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:30:22 INFO Hourly dataset computed and listing created
2026-07-27 13:30:36 INFO Hourly dataset computed
2026-07-27 13:30:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:30:37 INFO Hourly dataset computed and listing created
2026-07-27 13:30:50 INFO Hourly dataset computed
2026-07-27 13:30:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:30:51 INFO Hourly dataset computed and listing created
2026-07-27 13:31:02 INFO Hourly dataset computed
2026-07-27 13:31:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:31:04 INFO Hourly dataset computed and listing created
2026-07-27 13:31:18 INFO Hourly dataset computed
2026-07-27 13:31:18 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:31:19 INFO Hourly dataset computed and listing created
2026-07-27 13:31:32 INFO Hourly dataset computed
2026-07-27 13:31:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:31:34 INFO Hourly dataset computed and listing created
2026-07-27 13:31:45 INFO Hourly dataset computed
2026-07-27 13:31:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:31:47 INFO Hourly dataset computed and listing created
2026-07-27 13:32:02 INFO Hourly dataset computed
2026-07-27 13:32:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:32:03 INFO Hourly dataset computed and listing created
2026-07-27 13:32:15 INFO Hourly dataset computed
2026-07-27 13:32:15 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:32:17 INFO Hourly dataset computed and listing created
2026-07-27 13:32:28 INFO Hourly dataset computed
2026-07-27 13:32:28 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-27 13:32:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:32:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 13:32:28 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-27 13:32:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 13:32:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:32:28 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 13:32:28 INFO Queuing job for member 1...
2026-07-27 13:32:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:32:28 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 13:32:30 INFO Found: ['5285657']
2026-07-27 13:32:35 INFO [TGCC-IRENE] Submitted job with ID:['5285657']
2026-07-27 13:32:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:32:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 13:32:35 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-27 13:32:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 13:32:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:32:35 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 13:32:35 INFO Queuing job for member 2...
2026-07-27 13:32:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:32:35 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 13:32:37 INFO Found: ['5285658']
2026-07-27 13:32:42 INFO [TGCC-IRENE] Submitted job with ID:['5285658']
2026-07-27 13:32:42 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:32:42 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 13:32:42 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-27 13:32:42 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 13:32:42 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:32:42 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 13:32:42 INFO Queuing job for member 3...
2026-07-27 13:32:42 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:32:42 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 13:32:45 INFO Found: ['5285660']
2026-07-27 13:32:50 INFO [TGCC-IRENE] Submitted job with ID:['5285660']
2026-07-27 13:32:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:32:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 13:32:50 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-27 13:32:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 13:32:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:32:50 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 13:32:50 INFO Queuing job for member 4...
2026-07-27 13:32:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:32:50 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 13:32:52 INFO Found: ['5285661']
2026-07-27 13:32:57 INFO [TGCC-IRENE] Submitted job with ID:['5285661']
2026-07-27 13:32:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:32:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 13:32:57 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-27 13:32:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 13:32:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:32:57 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 13:32:57 INFO Queuing job for member 5...
2026-07-27 13:32:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:32:57 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 13:33:00 INFO Found: ['5285662']
2026-07-27 13:33:05 INFO [TGCC-IRENE] Submitted job with ID:['5285662']
2026-07-27 13:33:05 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:05 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 13:33:05 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-27 13:33:05 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 13:33:05 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:05 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 13:33:05 INFO Queuing job for member 6...
2026-07-27 13:33:05 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:05 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 13:33:07 INFO Found: ['5285664']
2026-07-27 13:33:12 INFO [TGCC-IRENE] Submitted job with ID:['5285664']
2026-07-27 13:33:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 13:33:12 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-27 13:33:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 13:33:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:12 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 13:33:12 INFO Queuing job for member 7...
2026-07-27 13:33:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:12 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 13:33:15 INFO Found: ['5285665']
2026-07-27 13:33:20 INFO [TGCC-IRENE] Submitted job with ID:['5285665']
2026-07-27 13:33:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 13:33:20 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-27 13:33:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 13:33:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:20 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 13:33:20 INFO Queuing job for member 8...
2026-07-27 13:33:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:20 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 13:33:20 INFO Found: ['5285668']
2026-07-27 13:33:25 INFO [TGCC-IRENE] Submitted job with ID:['5285668']
2026-07-27 13:33:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 13:33:25 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-27 13:33:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 13:33:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:25 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 13:33:25 INFO Queuing job for member 9...
2026-07-27 13:33:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:25 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 13:33:26 INFO Found: ['5285671']
2026-07-27 13:33:31 INFO [TGCC-IRENE] Submitted job with ID:['5285671']
2026-07-27 13:33:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 13:33:31 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-27 13:33:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 13:33:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:31 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 13:33:31 INFO Queuing job for member 10...
2026-07-27 13:33:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:31 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 13:33:32 INFO Found: ['5285672']
2026-07-27 13:33:37 INFO [TGCC-IRENE] Submitted job with ID:['5285672']
2026-07-27 13:33:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 13:33:37 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-27 13:33:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 13:33:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:37 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 13:33:37 INFO Queuing job for member 11...
2026-07-27 13:33:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:37 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 13:33:38 INFO Found: ['5285673']
2026-07-27 13:33:43 INFO [TGCC-IRENE] Submitted job with ID:['5285673']
2026-07-27 13:33:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 13:33:43 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-27 13:33:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 13:33:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:43 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 13:33:43 INFO Queuing job for member 12...
2026-07-27 13:33:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:43 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 13:33:43 INFO Found: ['5285675']
2026-07-27 13:33:48 INFO [TGCC-IRENE] Submitted job with ID:['5285675']
2026-07-27 13:33:48 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:48 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 13:33:48 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-27 13:33:48 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 13:33:48 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:48 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 13:33:48 INFO Queuing job for member 13...
2026-07-27 13:33:48 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:48 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 13:33:49 INFO Found: ['5285678']
2026-07-27 13:33:54 INFO [TGCC-IRENE] Submitted job with ID:['5285678']
2026-07-27 13:33:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:33:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 13:33:54 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-27 13:33:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 13:33:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:33:54 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 13:33:54 INFO Queuing job for member 14...
2026-07-27 13:33:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:33:54 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 13:33:55 INFO Found: ['5285682']
2026-07-27 13:34:00 INFO [TGCC-IRENE] Submitted job with ID:['5285682']
2026-07-27 13:34:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:34:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 13:34:00 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-27 13:34:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 13:34:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:34:00 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 13:34:00 INFO Queuing job for member 15...
2026-07-27 13:34:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:34:00 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 13:34:01 INFO Found: ['5285683']
2026-07-27 13:34:06 INFO [TGCC-IRENE] Submitted job with ID:['5285683']
2026-07-27 13:34:06 INFO Checking job status ...
2026-07-27 13:34:06 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:34:06 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:34:06 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:34:21 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:34:21 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:34:23 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:34:23 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:34:38 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:34:38 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:34:38 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:34:53 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:34:53 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:34:53 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:34:56 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:34:56 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:35:11 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:35:11 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:35:11 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:35:26 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:35:26 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:35:26 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:35:41 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:35:41 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:35:42 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:35:42 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:35:57 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:35:57 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:35:57 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:36:12 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:36:12 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:36:14 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:36:14 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:36:29 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:36:29 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:36:30 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:36:32 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:36:32 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:36:32 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:36:32 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:36:32 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:36:32 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:36:47 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:36:47 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:36:47 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:37:02 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:37:02 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:37:02 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:37:17 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:37:17 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:37:18 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:37:18 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:37:34 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:37:34 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:37:35 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:37:35 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:37:35 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:37:35 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:37:35 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:37:35 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:37:50 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:37:50 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:37:52 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:37:52 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:37:52 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:37:52 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:37:52 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:37:52 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:38:07 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:38:07 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:38:07 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:38:22 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:38:22 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:38:22 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:38:38 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:38:38 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:38:38 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:38:53 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:38:53 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:38:53 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:39:09 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:39:09 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:39:09 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:39:24 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:39:24 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:39:24 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:39:39 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:39:39 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:39:39 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:39:39 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:39:41 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:39:41 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:39:56 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:39:56 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:39:56 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:39:56 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:39:56 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:39:56 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:39:57 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:39:57 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:40:12 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:40:12 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:40:12 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:40:27 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:40:27 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:40:27 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:40:42 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:40:42 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:40:43 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:40:43 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:40:43 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:40:43 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:40:58 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:40:58 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:40:58 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:41:00 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:41:00 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:41:15 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:41:15 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:41:15 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:41:30 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:41:30 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:41:30 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:41:30 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:41:30 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:41:32 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:41:33 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:41:33 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:41:48 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:41:48 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:41:48 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:42:03 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:42:03 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:42:03 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:42:18 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:42:18 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:42:18 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:42:34 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:42:34 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:42:34 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:42:49 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:42:49 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:42:50 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:42:52 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:42:52 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:43:07 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:43:07 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:43:07 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:43:22 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:43:22 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:43:24 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:43:24 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:43:24 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:43:24 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:43:39 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:43:39 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:43:40 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:43:40 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:43:40 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:43:40 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:43:40 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:43:40 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:43:55 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:43:55 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:43:55 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:44:10 INFO None 5285657: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:44:10 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:44:10 INFO Jobs still running: ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:44:27 INFO None 5285657: status FINISHED
2026-07-27 13:44:27 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:44:27 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:44:27 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:44:42 INFO None 5285657: status FINISHED
2026-07-27 13:44:42 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:44:42 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:44:45 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:44:45 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:45:00 INFO None 5285657: status FINISHED
2026-07-27 13:45:00 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:45:00 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:45:00 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:45:15 INFO None 5285657: status FINISHED
2026-07-27 13:45:15 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:45:15 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:45:15 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:45:30 INFO None 5285657: status FINISHED
2026-07-27 13:45:30 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:45:30 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:45:31 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:45:31 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:45:46 INFO None 5285657: status FINISHED
2026-07-27 13:45:46 INFO None 5285658: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:45:46 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:45:46 INFO Jobs still running: ['5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:46:02 INFO None 5285657: status FINISHED
2026-07-27 13:46:02 INFO None 5285658: status FINISHED
2026-07-27 13:46:02 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:46:02 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:46:02 INFO Jobs still running: ['5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:46:17 INFO None 5285657: status FINISHED
2026-07-27 13:46:17 INFO None 5285658: status FINISHED
2026-07-27 13:46:17 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:46:17 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:46:19 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:46:19 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:46:19 INFO Jobs still running: ['5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:46:34 INFO None 5285657: status FINISHED
2026-07-27 13:46:34 INFO None 5285658: status FINISHED
2026-07-27 13:46:35 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:46:35 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:46:35 INFO Jobs still running: ['5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:46:50 INFO None 5285657: status FINISHED
2026-07-27 13:46:50 INFO None 5285658: status FINISHED
2026-07-27 13:46:50 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:46:50 INFO None 5285661: status RUNNING/PENDING
2026-07-27 13:46:50 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:46:51 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:46:51 INFO Jobs still running: ['5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:47:06 INFO None 5285657: status FINISHED
2026-07-27 13:47:06 INFO None 5285658: status FINISHED
2026-07-27 13:47:06 INFO None 5285660: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285661: status FINISHED
2026-07-27 13:47:06 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285668: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:47:06 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:47:06 INFO Jobs still running: ['5285660', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:47:21 INFO None 5285657: status FINISHED
2026-07-27 13:47:21 INFO None 5285658: status FINISHED
2026-07-27 13:47:21 INFO None 5285660: status FINISHED
2026-07-27 13:47:21 INFO None 5285661: status FINISHED
2026-07-27 13:47:21 INFO None 5285662: status RUNNING/PENDING
2026-07-27 13:47:21 INFO None 5285664: status RUNNING/PENDING
2026-07-27 13:47:21 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:47:21 INFO None 5285668: status FINISHED
2026-07-27 13:47:21 INFO None 5285671: status RUNNING/PENDING
2026-07-27 13:47:21 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:47:21 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:47:22 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:47:22 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:47:22 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:47:22 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:47:22 INFO Jobs still running: ['5285662', '5285664', '5285665', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:47:38 INFO None 5285657: status FINISHED
2026-07-27 13:47:38 INFO None 5285658: status FINISHED
2026-07-27 13:47:38 INFO None 5285660: status FINISHED
2026-07-27 13:47:38 INFO None 5285661: status FINISHED
2026-07-27 13:47:38 INFO None 5285662: status FINISHED
2026-07-27 13:47:38 INFO None 5285664: status FINISHED
2026-07-27 13:47:38 INFO None 5285665: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285668: status FINISHED
2026-07-27 13:47:38 INFO None 5285671: status FINISHED
2026-07-27 13:47:38 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:47:38 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:47:38 INFO Jobs still running: ['5285665', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:47:53 INFO None 5285657: status FINISHED
2026-07-27 13:47:53 INFO None 5285658: status FINISHED
2026-07-27 13:47:53 INFO None 5285660: status FINISHED
2026-07-27 13:47:53 INFO None 5285661: status FINISHED
2026-07-27 13:47:53 INFO None 5285662: status FINISHED
2026-07-27 13:47:53 INFO None 5285664: status FINISHED
2026-07-27 13:47:53 INFO None 5285665: status FINISHED
2026-07-27 13:47:53 INFO None 5285668: status FINISHED
2026-07-27 13:47:53 INFO None 5285671: status FINISHED
2026-07-27 13:47:53 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:47:53 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:47:53 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:47:53 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:47:53 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:47:53 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:47:53 INFO Jobs still running: ['5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:48:08 INFO None 5285657: status FINISHED
2026-07-27 13:48:08 INFO None 5285658: status FINISHED
2026-07-27 13:48:08 INFO None 5285660: status FINISHED
2026-07-27 13:48:11 INFO None 5285661: status FINISHED
2026-07-27 13:48:11 INFO None 5285662: status FINISHED
2026-07-27 13:48:11 INFO None 5285664: status FINISHED
2026-07-27 13:48:11 INFO None 5285665: status FINISHED
2026-07-27 13:48:11 INFO None 5285668: status FINISHED
2026-07-27 13:48:11 INFO None 5285671: status FINISHED
2026-07-27 13:48:11 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:48:11 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:48:11 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:48:11 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:48:11 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:48:11 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:48:11 INFO Jobs still running: ['5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:48:26 INFO None 5285657: status FINISHED
2026-07-27 13:48:26 INFO None 5285658: status FINISHED
2026-07-27 13:48:26 INFO None 5285660: status FINISHED
2026-07-27 13:48:26 INFO None 5285661: status FINISHED
2026-07-27 13:48:26 INFO None 5285662: status FINISHED
2026-07-27 13:48:26 INFO None 5285664: status FINISHED
2026-07-27 13:48:26 INFO None 5285665: status FINISHED
2026-07-27 13:48:26 INFO None 5285668: status FINISHED
2026-07-27 13:48:26 INFO None 5285671: status FINISHED
2026-07-27 13:48:26 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:48:26 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:48:26 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:48:26 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:48:26 INFO None 5285682: status RUNNING/PENDING
2026-07-27 13:48:26 INFO None 5285683: status RUNNING/PENDING
2026-07-27 13:48:26 INFO Jobs still running: ['5285672', '5285673', '5285675', '5285678', '5285682', '5285683']. Waiting...
2026-07-27 13:48:41 INFO None 5285657: status FINISHED
2026-07-27 13:48:41 INFO None 5285658: status FINISHED
2026-07-27 13:48:41 INFO None 5285660: status FINISHED
2026-07-27 13:48:41 INFO None 5285661: status FINISHED
2026-07-27 13:48:41 INFO None 5285662: status FINISHED
2026-07-27 13:48:41 INFO None 5285664: status FINISHED
2026-07-27 13:48:41 INFO None 5285665: status FINISHED
2026-07-27 13:48:41 INFO None 5285668: status FINISHED
2026-07-27 13:48:41 INFO None 5285671: status FINISHED
2026-07-27 13:48:41 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:48:41 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:48:41 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:48:41 INFO None 5285678: status RUNNING/PENDING
2026-07-27 13:48:41 INFO None 5285682: status FINISHED
2026-07-27 13:48:41 INFO None 5285683: status FINISHED
2026-07-27 13:48:41 INFO Jobs still running: ['5285672', '5285673', '5285675', '5285678']. Waiting...
2026-07-27 13:48:56 INFO None 5285657: status FINISHED
2026-07-27 13:48:56 INFO None 5285658: status FINISHED
2026-07-27 13:48:56 INFO None 5285660: status FINISHED
2026-07-27 13:48:56 INFO None 5285661: status FINISHED
2026-07-27 13:48:56 INFO None 5285662: status FINISHED
2026-07-27 13:48:56 INFO None 5285664: status FINISHED
2026-07-27 13:48:56 INFO None 5285665: status FINISHED
2026-07-27 13:48:57 INFO None 5285668: status FINISHED
2026-07-27 13:48:57 INFO None 5285671: status FINISHED
2026-07-27 13:48:57 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:48:57 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:48:57 INFO None 5285675: status RUNNING/PENDING
2026-07-27 13:48:57 INFO None 5285678: status FINISHED
2026-07-27 13:48:57 INFO None 5285682: status FINISHED
2026-07-27 13:48:57 INFO None 5285683: status FINISHED
2026-07-27 13:48:57 INFO Jobs still running: ['5285672', '5285673', '5285675']. Waiting...
2026-07-27 13:49:12 INFO None 5285657: status FINISHED
2026-07-27 13:49:12 INFO None 5285658: status FINISHED
2026-07-27 13:49:12 INFO None 5285660: status FINISHED
2026-07-27 13:49:12 INFO None 5285661: status FINISHED
2026-07-27 13:49:12 INFO None 5285662: status FINISHED
2026-07-27 13:49:12 INFO None 5285664: status FINISHED
2026-07-27 13:49:12 INFO None 5285665: status FINISHED
2026-07-27 13:49:12 INFO None 5285668: status FINISHED
2026-07-27 13:49:12 INFO None 5285671: status FINISHED
2026-07-27 13:49:12 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:49:12 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:49:12 INFO None 5285675: status FINISHED
2026-07-27 13:49:12 INFO None 5285678: status FINISHED
2026-07-27 13:49:12 INFO None 5285682: status FINISHED
2026-07-27 13:49:12 INFO None 5285683: status FINISHED
2026-07-27 13:49:12 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:49:27 INFO None 5285657: status FINISHED
2026-07-27 13:49:27 INFO None 5285658: status FINISHED
2026-07-27 13:49:27 INFO None 5285660: status FINISHED
2026-07-27 13:49:27 INFO None 5285661: status FINISHED
2026-07-27 13:49:27 INFO None 5285662: status FINISHED
2026-07-27 13:49:27 INFO None 5285664: status FINISHED
2026-07-27 13:49:27 INFO None 5285665: status FINISHED
2026-07-27 13:49:27 INFO None 5285668: status FINISHED
2026-07-27 13:49:27 INFO None 5285671: status FINISHED
2026-07-27 13:49:27 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:49:27 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:49:27 INFO None 5285675: status FINISHED
2026-07-27 13:49:27 INFO None 5285678: status FINISHED
2026-07-27 13:49:27 INFO None 5285682: status FINISHED
2026-07-27 13:49:27 INFO None 5285683: status FINISHED
2026-07-27 13:49:27 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:49:42 INFO None 5285657: status FINISHED
2026-07-27 13:49:42 INFO None 5285658: status FINISHED
2026-07-27 13:49:42 INFO None 5285660: status FINISHED
2026-07-27 13:49:42 INFO None 5285661: status FINISHED
2026-07-27 13:49:42 INFO None 5285662: status FINISHED
2026-07-27 13:49:44 INFO None 5285664: status FINISHED
2026-07-27 13:49:44 INFO None 5285665: status FINISHED
2026-07-27 13:49:44 INFO None 5285668: status FINISHED
2026-07-27 13:49:44 INFO None 5285671: status FINISHED
2026-07-27 13:49:44 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:49:44 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:49:44 INFO None 5285675: status FINISHED
2026-07-27 13:49:44 INFO None 5285678: status FINISHED
2026-07-27 13:49:44 INFO None 5285682: status FINISHED
2026-07-27 13:49:45 INFO None 5285683: status FINISHED
2026-07-27 13:49:45 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:50:00 INFO None 5285657: status FINISHED
2026-07-27 13:50:00 INFO None 5285658: status FINISHED
2026-07-27 13:50:00 INFO None 5285660: status FINISHED
2026-07-27 13:50:00 INFO None 5285661: status FINISHED
2026-07-27 13:50:00 INFO None 5285662: status FINISHED
2026-07-27 13:50:00 INFO None 5285664: status FINISHED
2026-07-27 13:50:00 INFO None 5285665: status FINISHED
2026-07-27 13:50:00 INFO None 5285668: status FINISHED
2026-07-27 13:50:00 INFO None 5285671: status FINISHED
2026-07-27 13:50:00 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:50:00 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:50:00 INFO None 5285675: status FINISHED
2026-07-27 13:50:00 INFO None 5285678: status FINISHED
2026-07-27 13:50:00 INFO None 5285682: status FINISHED
2026-07-27 13:50:00 INFO None 5285683: status FINISHED
2026-07-27 13:50:00 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:50:15 INFO None 5285657: status FINISHED
2026-07-27 13:50:15 INFO None 5285658: status FINISHED
2026-07-27 13:50:15 INFO None 5285660: status FINISHED
2026-07-27 13:50:15 INFO None 5285661: status FINISHED
2026-07-27 13:50:15 INFO None 5285662: status FINISHED
2026-07-27 13:50:15 INFO None 5285664: status FINISHED
2026-07-27 13:50:15 INFO None 5285665: status FINISHED
2026-07-27 13:50:15 INFO None 5285668: status FINISHED
2026-07-27 13:50:15 INFO None 5285671: status FINISHED
2026-07-27 13:50:15 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:50:15 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:50:15 INFO None 5285675: status FINISHED
2026-07-27 13:50:15 INFO None 5285678: status FINISHED
2026-07-27 13:50:15 INFO None 5285682: status FINISHED
2026-07-27 13:50:15 INFO None 5285683: status FINISHED
2026-07-27 13:50:15 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:50:30 INFO None 5285657: status FINISHED
2026-07-27 13:50:30 INFO None 5285658: status FINISHED
2026-07-27 13:50:30 INFO None 5285660: status FINISHED
2026-07-27 13:50:30 INFO None 5285661: status FINISHED
2026-07-27 13:50:30 INFO None 5285662: status FINISHED
2026-07-27 13:50:30 INFO None 5285664: status FINISHED
2026-07-27 13:50:30 INFO None 5285665: status FINISHED
2026-07-27 13:50:30 INFO None 5285668: status FINISHED
2026-07-27 13:50:30 INFO None 5285671: status FINISHED
2026-07-27 13:50:30 INFO None 5285672: status RUNNING/PENDING
2026-07-27 13:50:30 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:50:30 INFO None 5285675: status FINISHED
2026-07-27 13:50:30 INFO None 5285678: status FINISHED
2026-07-27 13:50:30 INFO None 5285682: status FINISHED
2026-07-27 13:50:30 INFO None 5285683: status FINISHED
2026-07-27 13:50:30 INFO Jobs still running: ['5285672', '5285673']. Waiting...
2026-07-27 13:50:45 INFO None 5285657: status FINISHED
2026-07-27 13:50:45 INFO None 5285658: status FINISHED
2026-07-27 13:50:45 INFO None 5285660: status FINISHED
2026-07-27 13:50:45 INFO None 5285661: status FINISHED
2026-07-27 13:50:45 INFO None 5285662: status FINISHED
2026-07-27 13:50:45 INFO None 5285664: status FINISHED
2026-07-27 13:50:45 INFO None 5285665: status FINISHED
2026-07-27 13:50:45 INFO None 5285668: status FINISHED
2026-07-27 13:50:46 INFO None 5285671: status FINISHED
2026-07-27 13:50:46 INFO None 5285672: status FINISHED
2026-07-27 13:50:46 INFO None 5285673: status RUNNING/PENDING
2026-07-27 13:50:46 INFO None 5285675: status FINISHED
2026-07-27 13:50:46 INFO None 5285678: status FINISHED
2026-07-27 13:50:46 INFO None 5285682: status FINISHED
2026-07-27 13:50:46 INFO None 5285683: status FINISHED
2026-07-27 13:50:46 INFO Jobs still running: ['5285673']. Waiting...
2026-07-27 13:51:02 INFO None 5285657: status FINISHED
2026-07-27 13:51:02 INFO None 5285658: status FINISHED
2026-07-27 13:51:02 INFO None 5285660: status FINISHED
2026-07-27 13:51:02 INFO None 5285661: status FINISHED
2026-07-27 13:51:02 INFO None 5285662: status FINISHED
2026-07-27 13:51:02 INFO None 5285664: status FINISHED
2026-07-27 13:51:02 INFO None 5285665: status FINISHED
2026-07-27 13:51:02 INFO None 5285668: status FINISHED
2026-07-27 13:51:02 INFO None 5285671: status FINISHED
2026-07-27 13:51:02 INFO None 5285672: status FINISHED
2026-07-27 13:51:02 INFO None 5285673: status FINISHED
2026-07-27 13:51:02 INFO None 5285675: status FINISHED
2026-07-27 13:51:02 INFO None 5285678: status FINISHED
2026-07-27 13:51:03 INFO None 5285682: status FINISHED
2026-07-27 13:51:03 INFO None 5285683: status FINISHED
2026-07-27 13:51:03 INFO Jobs ['5285657', '5285658', '5285660', '5285661', '5285662', '5285664', '5285665', '5285668', '5285671', '5285672', '5285673', '5285675', '5285678', '5285682', '5285683'] have finished
2026-07-27 13:51:03 INFO Checking restart files were created ...
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-27 13:51:03 INFO  Run_model() completed successfully.
2026-07-27 13:51:03 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:51:03 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-27 13:51:03 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 13:51:03 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-27 13:51:03 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:51:03 INFO ---------->>> Running process_satellite_data()
2026-07-27 13:51:03 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-27 13:51:03 INFO ---------->>> Running run_obs_converter()
2026-07-27 13:51:03 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-27 13:51:03 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-27 13:51:03 INFO ---------->>> Running DART
2026-07-27 13:51:03 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO File exists: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-27 13:51:03 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 13:51:03 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 13:51:03 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 13:51:03 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 13:51:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 13:51:03 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 13:51:12 INFO Found: []
2026-07-27 13:51:12 INFO No job id returned by command ./run_filter.bsh
2026-07-27 13:51:12 INFO No monitoring will be performed
2026-07-27 13:51:12 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-27 13:51:12 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020609'
2026-07-27 13:51:12 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 13:51:15 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 13:51:15 INFO run_dart() is DONE.
2026-07-27 13:51:15 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 13:51:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-27 13:51:28 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:51:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-07-27 13:51:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-27 13:51:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:51:42 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-07-27 13:51:42 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-27 13:51:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:51:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-07-27 13:51:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-27 13:52:08 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:52:09 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-07-27 13:52:09 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-27 13:52:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:52:22 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-07-27 13:52:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-27 13:52:35 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:52:36 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-07-27 13:52:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-27 13:52:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:52:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-07-27 13:52:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-27 13:53:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:53:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-07-27 13:53:03 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-27 13:53:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:53:17 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-07-27 13:53:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-27 13:53:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:53:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-07-27 13:53:31 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-27 13:53:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:53:45 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-07-27 13:53:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-27 13:53:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:53:59 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-07-27 13:53:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-27 13:54:12 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:54:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-07-27 13:54:13 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-27 13:54:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:54:26 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-07-27 13:54:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-27 13:54:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 13:54:39 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 13:54:39 INFO [TIME] step_end current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:54:39 INFO [TIME] step_start current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 13:54:39 INFO [TIME] window start=2020-02-06 09:00:00 end=2020-02-06 11:00:00 run_hours=2 has_assimilation=True
2026-07-27 13:54:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:54:40 INFO Hourly dataset computed and listing created
2026-07-27 13:54:47 INFO Hourly dataset computed
2026-07-27 13:54:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:54:48 INFO Hourly dataset computed and listing created
2026-07-27 13:54:51 INFO Hourly dataset computed
2026-07-27 13:54:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:54:52 INFO Hourly dataset computed and listing created
2026-07-27 13:54:54 INFO Hourly dataset computed
2026-07-27 13:54:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:54:55 INFO Hourly dataset computed and listing created
2026-07-27 13:54:57 INFO Hourly dataset computed
2026-07-27 13:54:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:54:58 INFO Hourly dataset computed and listing created
2026-07-27 13:55:01 INFO Hourly dataset computed
2026-07-27 13:55:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:01 INFO Hourly dataset computed and listing created
2026-07-27 13:55:02 INFO Hourly dataset computed
2026-07-27 13:55:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:03 INFO Hourly dataset computed and listing created
2026-07-27 13:55:04 INFO Hourly dataset computed
2026-07-27 13:55:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:05 INFO Hourly dataset computed and listing created
2026-07-27 13:55:06 INFO Hourly dataset computed
2026-07-27 13:55:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:06 INFO Hourly dataset computed and listing created
2026-07-27 13:55:07 INFO Hourly dataset computed
2026-07-27 13:55:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:08 INFO Hourly dataset computed and listing created
2026-07-27 13:55:09 INFO Hourly dataset computed
2026-07-27 13:55:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:10 INFO Hourly dataset computed and listing created
2026-07-27 13:55:11 INFO Hourly dataset computed
2026-07-27 13:55:11 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:12 INFO Hourly dataset computed and listing created
2026-07-27 13:55:12 INFO Hourly dataset computed
2026-07-27 13:55:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:13 INFO Hourly dataset computed and listing created
2026-07-27 13:55:14 INFO Hourly dataset computed
2026-07-27 13:55:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:15 INFO Hourly dataset computed and listing created
2026-07-27 13:55:16 INFO Hourly dataset computed
2026-07-27 13:55:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 13:55:17 INFO Hourly dataset computed and listing created
2026-07-27 13:55:18 INFO Hourly dataset computed
2026-07-27 13:55:18 INFO ---------->>> Running CHIMERE model from 2020-02-06 09:00:00 to 2020-02-06 11:00:00
2026-07-27 13:55:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 13:55:18 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc
2026-07-27 13:55:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 13:55:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:18 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 13:55:18 INFO Queuing job for member 1...
2026-07-27 13:55:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:18 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 13:55:18 INFO Found: ['5285900']
2026-07-27 13:55:23 INFO [TGCC-IRENE] Submitted job with ID:['5285900']
2026-07-27 13:55:23 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:23 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 13:55:23 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc
2026-07-27 13:55:23 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 13:55:23 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:23 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 13:55:23 INFO Queuing job for member 2...
2026-07-27 13:55:23 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:23 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 13:55:24 INFO Found: ['5285901']
2026-07-27 13:55:29 INFO [TGCC-IRENE] Submitted job with ID:['5285901']
2026-07-27 13:55:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 13:55:29 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc
2026-07-27 13:55:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 13:55:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:29 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 13:55:29 INFO Queuing job for member 3...
2026-07-27 13:55:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:29 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 13:55:30 INFO Found: ['5285903']
2026-07-27 13:55:35 INFO [TGCC-IRENE] Submitted job with ID:['5285903']
2026-07-27 13:55:35 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:35 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 13:55:35 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc
2026-07-27 13:55:35 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 13:55:35 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:35 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 13:55:35 INFO Queuing job for member 4...
2026-07-27 13:55:35 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:35 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 13:55:36 INFO Found: ['5285905']
2026-07-27 13:55:41 INFO [TGCC-IRENE] Submitted job with ID:['5285905']
2026-07-27 13:55:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 13:55:41 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc
2026-07-27 13:55:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 13:55:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:41 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 13:55:41 INFO Queuing job for member 5...
2026-07-27 13:55:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:41 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 13:55:41 INFO Found: ['5285906']
2026-07-27 13:55:46 INFO [TGCC-IRENE] Submitted job with ID:['5285906']
2026-07-27 13:55:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 13:55:46 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc
2026-07-27 13:55:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 13:55:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:46 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 13:55:46 INFO Queuing job for member 6...
2026-07-27 13:55:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:46 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 13:55:47 INFO Found: ['5285907']
2026-07-27 13:55:52 INFO [TGCC-IRENE] Submitted job with ID:['5285907']
2026-07-27 13:55:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 13:55:52 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc
2026-07-27 13:55:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 13:55:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:52 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 13:55:52 INFO Queuing job for member 7...
2026-07-27 13:55:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:52 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 13:55:54 INFO Found: ['5285908']
2026-07-27 13:55:59 INFO [TGCC-IRENE] Submitted job with ID:['5285908']
2026-07-27 13:55:59 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:55:59 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 13:55:59 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc
2026-07-27 13:55:59 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 13:55:59 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:55:59 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 13:55:59 INFO Queuing job for member 8...
2026-07-27 13:55:59 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:55:59 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 13:56:01 INFO Found: ['5285914']
2026-07-27 13:56:06 INFO [TGCC-IRENE] Submitted job with ID:['5285914']
2026-07-27 13:56:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 13:56:06 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc
2026-07-27 13:56:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 13:56:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:06 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 13:56:06 INFO Queuing job for member 9...
2026-07-27 13:56:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:06 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 13:56:09 INFO Found: ['5285920']
2026-07-27 13:56:14 INFO [TGCC-IRENE] Submitted job with ID:['5285920']
2026-07-27 13:56:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 13:56:14 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc
2026-07-27 13:56:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 13:56:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:14 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 13:56:14 INFO Queuing job for member 10...
2026-07-27 13:56:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:14 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 13:56:17 INFO Found: ['5285921']
2026-07-27 13:56:22 INFO [TGCC-IRENE] Submitted job with ID:['5285921']
2026-07-27 13:56:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 13:56:22 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc
2026-07-27 13:56:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 13:56:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:22 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 13:56:22 INFO Queuing job for member 11...
2026-07-27 13:56:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:22 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 13:56:24 INFO Found: ['5285922']
2026-07-27 13:56:29 INFO [TGCC-IRENE] Submitted job with ID:['5285922']
2026-07-27 13:56:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 13:56:29 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc
2026-07-27 13:56:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 13:56:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:29 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 13:56:29 INFO Queuing job for member 12...
2026-07-27 13:56:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:29 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 13:56:32 INFO Found: ['5285923']
2026-07-27 13:56:37 INFO [TGCC-IRENE] Submitted job with ID:['5285923']
2026-07-27 13:56:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 13:56:37 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc
2026-07-27 13:56:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 13:56:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:37 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 13:56:37 INFO Queuing job for member 13...
2026-07-27 13:56:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:37 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 13:56:39 INFO Found: ['5285924']
2026-07-27 13:56:44 INFO [TGCC-IRENE] Submitted job with ID:['5285924']
2026-07-27 13:56:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 13:56:44 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc
2026-07-27 13:56:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 13:56:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:44 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 13:56:44 INFO Queuing job for member 14...
2026-07-27 13:56:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:44 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 13:56:47 INFO Found: ['5285925']
2026-07-27 13:56:52 INFO [TGCC-IRENE] Submitted job with ID:['5285925']
2026-07-27 13:56:52 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 13:56:52 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 13:56:52 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc
2026-07-27 13:56:52 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 13:56:52 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 13:56:52 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 13:56:52 INFO Queuing job for member 15...
2026-07-27 13:56:52 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 13:56:52 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 13:56:53 INFO Found: ['5285930']
2026-07-27 13:56:58 INFO [TGCC-IRENE] Submitted job with ID:['5285930']
2026-07-27 13:56:58 INFO Checking job status ...
2026-07-27 13:56:58 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:56:58 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:56:58 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:57:13 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:57:13 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:57:13 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:57:28 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:57:28 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:57:28 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:57:44 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:57:44 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:57:45 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:57:45 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:57:45 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:57:45 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:57:45 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:57:45 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:58:00 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:58:00 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:58:02 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:58:02 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:58:17 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:58:17 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:58:17 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:58:32 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:58:32 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:58:33 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:58:33 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:58:48 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:58:48 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:58:48 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:59:03 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:59:03 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:59:03 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:59:18 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285903: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285905: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285906: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285907: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:59:18 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:59:21 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:59:21 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:59:21 INFO Jobs still running: ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:59:36 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285903: status FINISHED
2026-07-27 13:59:36 INFO None 5285905: status FINISHED
2026-07-27 13:59:36 INFO None 5285906: status FINISHED
2026-07-27 13:59:36 INFO None 5285907: status FINISHED
2026-07-27 13:59:36 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:59:36 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:59:36 INFO Jobs still running: ['5285900', '5285901', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 13:59:51 INFO None 5285900: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285901: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285903: status FINISHED
2026-07-27 13:59:51 INFO None 5285905: status FINISHED
2026-07-27 13:59:51 INFO None 5285906: status FINISHED
2026-07-27 13:59:51 INFO None 5285907: status FINISHED
2026-07-27 13:59:51 INFO None 5285908: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285914: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285920: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285921: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285922: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285923: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285924: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285925: status RUNNING/PENDING
2026-07-27 13:59:51 INFO None 5285930: status RUNNING/PENDING
2026-07-27 13:59:51 INFO Jobs still running: ['5285900', '5285901', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:00:06 INFO None 5285900: status RUNNING/PENDING
2026-07-27 14:00:06 INFO None 5285901: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285903: status FINISHED
2026-07-27 14:00:08 INFO None 5285905: status FINISHED
2026-07-27 14:00:08 INFO None 5285906: status FINISHED
2026-07-27 14:00:08 INFO None 5285907: status FINISHED
2026-07-27 14:00:08 INFO None 5285908: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285920: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285921: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285922: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:00:08 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:00:08 INFO Jobs still running: ['5285900', '5285901', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:00:23 INFO None 5285900: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285901: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285903: status FINISHED
2026-07-27 14:00:24 INFO None 5285905: status FINISHED
2026-07-27 14:00:24 INFO None 5285906: status FINISHED
2026-07-27 14:00:24 INFO None 5285907: status FINISHED
2026-07-27 14:00:24 INFO None 5285908: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285920: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285921: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285922: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:00:24 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:00:24 INFO Jobs still running: ['5285900', '5285901', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:00:39 INFO None 5285900: status FINISHED
2026-07-27 14:00:39 INFO None 5285901: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285903: status FINISHED
2026-07-27 14:00:39 INFO None 5285905: status FINISHED
2026-07-27 14:00:39 INFO None 5285906: status FINISHED
2026-07-27 14:00:39 INFO None 5285907: status FINISHED
2026-07-27 14:00:39 INFO None 5285908: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285920: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285921: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285922: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:00:39 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:00:39 INFO Jobs still running: ['5285901', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:00:54 INFO None 5285900: status FINISHED
2026-07-27 14:00:54 INFO None 5285901: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285903: status FINISHED
2026-07-27 14:00:54 INFO None 5285905: status FINISHED
2026-07-27 14:00:54 INFO None 5285906: status FINISHED
2026-07-27 14:00:54 INFO None 5285907: status FINISHED
2026-07-27 14:00:54 INFO None 5285908: status FINISHED
2026-07-27 14:00:54 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285920: status FINISHED
2026-07-27 14:00:54 INFO None 5285921: status FINISHED
2026-07-27 14:00:54 INFO None 5285922: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:00:54 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:00:54 INFO Jobs still running: ['5285901', '5285914', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:01:10 INFO None 5285900: status FINISHED
2026-07-27 14:01:10 INFO None 5285901: status FINISHED
2026-07-27 14:01:10 INFO None 5285903: status FINISHED
2026-07-27 14:01:10 INFO None 5285905: status FINISHED
2026-07-27 14:01:10 INFO None 5285906: status FINISHED
2026-07-27 14:01:10 INFO None 5285907: status FINISHED
2026-07-27 14:01:10 INFO None 5285908: status FINISHED
2026-07-27 14:01:10 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:01:10 INFO None 5285920: status FINISHED
2026-07-27 14:01:10 INFO None 5285921: status FINISHED
2026-07-27 14:01:10 INFO None 5285922: status RUNNING/PENDING
2026-07-27 14:01:10 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:01:10 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:01:10 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:01:10 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:01:10 INFO Jobs still running: ['5285914', '5285922', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:01:25 INFO None 5285900: status FINISHED
2026-07-27 14:01:25 INFO None 5285901: status FINISHED
2026-07-27 14:01:25 INFO None 5285903: status FINISHED
2026-07-27 14:01:25 INFO None 5285905: status FINISHED
2026-07-27 14:01:25 INFO None 5285906: status FINISHED
2026-07-27 14:01:26 INFO None 5285907: status FINISHED
2026-07-27 14:01:26 INFO None 5285908: status FINISHED
2026-07-27 14:01:26 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:01:26 INFO None 5285920: status FINISHED
2026-07-27 14:01:26 INFO None 5285921: status FINISHED
2026-07-27 14:01:26 INFO None 5285922: status FINISHED
2026-07-27 14:01:26 INFO None 5285923: status RUNNING/PENDING
2026-07-27 14:01:26 INFO None 5285924: status RUNNING/PENDING
2026-07-27 14:01:26 INFO None 5285925: status RUNNING/PENDING
2026-07-27 14:01:28 INFO None 5285930: status RUNNING/PENDING
2026-07-27 14:01:28 INFO Jobs still running: ['5285914', '5285923', '5285924', '5285925', '5285930']. Waiting...
2026-07-27 14:01:43 INFO None 5285900: status FINISHED
2026-07-27 14:01:43 INFO None 5285901: status FINISHED
2026-07-27 14:01:43 INFO None 5285903: status FINISHED
2026-07-27 14:01:43 INFO None 5285905: status FINISHED
2026-07-27 14:01:43 INFO None 5285906: status FINISHED
2026-07-27 14:01:43 INFO None 5285907: status FINISHED
2026-07-27 14:01:43 INFO None 5285908: status FINISHED
2026-07-27 14:01:43 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:01:43 INFO None 5285920: status FINISHED
2026-07-27 14:01:43 INFO None 5285921: status FINISHED
2026-07-27 14:01:43 INFO None 5285922: status FINISHED
2026-07-27 14:01:43 INFO None 5285923: status FINISHED
2026-07-27 14:01:43 INFO None 5285924: status FINISHED
2026-07-27 14:01:43 INFO None 5285925: status FINISHED
2026-07-27 14:01:43 INFO None 5285930: status FINISHED
2026-07-27 14:01:43 INFO Jobs still running: ['5285914']. Waiting...
2026-07-27 14:01:58 INFO None 5285900: status FINISHED
2026-07-27 14:01:58 INFO None 5285901: status FINISHED
2026-07-27 14:01:58 INFO None 5285903: status FINISHED
2026-07-27 14:01:58 INFO None 5285905: status FINISHED
2026-07-27 14:01:58 INFO None 5285906: status FINISHED
2026-07-27 14:01:58 INFO None 5285907: status FINISHED
2026-07-27 14:01:58 INFO None 5285908: status FINISHED
2026-07-27 14:01:58 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:01:58 INFO None 5285920: status FINISHED
2026-07-27 14:01:58 INFO None 5285921: status FINISHED
2026-07-27 14:01:58 INFO None 5285922: status FINISHED
2026-07-27 14:01:58 INFO None 5285923: status FINISHED
2026-07-27 14:01:58 INFO None 5285924: status FINISHED
2026-07-27 14:01:58 INFO None 5285925: status FINISHED
2026-07-27 14:01:58 INFO None 5285930: status FINISHED
2026-07-27 14:01:58 INFO Jobs still running: ['5285914']. Waiting...
2026-07-27 14:02:13 INFO None 5285900: status FINISHED
2026-07-27 14:02:13 INFO None 5285901: status FINISHED
2026-07-27 14:02:13 INFO None 5285903: status FINISHED
2026-07-27 14:02:13 INFO None 5285905: status FINISHED
2026-07-27 14:02:13 INFO None 5285906: status FINISHED
2026-07-27 14:02:13 INFO None 5285907: status FINISHED
2026-07-27 14:02:13 INFO None 5285908: status FINISHED
2026-07-27 14:02:13 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:02:13 INFO None 5285920: status FINISHED
2026-07-27 14:02:13 INFO None 5285921: status FINISHED
2026-07-27 14:02:13 INFO None 5285922: status FINISHED
2026-07-27 14:02:13 INFO None 5285923: status FINISHED
2026-07-27 14:02:13 INFO None 5285924: status FINISHED
2026-07-27 14:02:14 INFO None 5285925: status FINISHED
2026-07-27 14:02:14 INFO None 5285930: status FINISHED
2026-07-27 14:02:14 INFO Jobs still running: ['5285914']. Waiting...
2026-07-27 14:02:29 INFO None 5285900: status FINISHED
2026-07-27 14:02:29 INFO None 5285901: status FINISHED
2026-07-27 14:02:29 INFO None 5285903: status FINISHED
2026-07-27 14:02:29 INFO None 5285905: status FINISHED
2026-07-27 14:02:29 INFO None 5285906: status FINISHED
2026-07-27 14:02:29 INFO None 5285907: status FINISHED
2026-07-27 14:02:29 INFO None 5285908: status FINISHED
2026-07-27 14:02:29 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:02:29 INFO None 5285920: status FINISHED
2026-07-27 14:02:29 INFO None 5285921: status FINISHED
2026-07-27 14:02:29 INFO None 5285922: status FINISHED
2026-07-27 14:02:29 INFO None 5285923: status FINISHED
2026-07-27 14:02:29 INFO None 5285924: status FINISHED
2026-07-27 14:02:29 INFO None 5285925: status FINISHED
2026-07-27 14:02:29 INFO None 5285930: status FINISHED
2026-07-27 14:02:29 INFO Jobs still running: ['5285914']. Waiting...
2026-07-27 14:02:44 INFO None 5285900: status FINISHED
2026-07-27 14:02:44 INFO None 5285901: status FINISHED
2026-07-27 14:02:44 INFO None 5285903: status FINISHED
2026-07-27 14:02:44 INFO None 5285905: status FINISHED
2026-07-27 14:02:44 INFO None 5285906: status FINISHED
2026-07-27 14:02:44 INFO None 5285907: status FINISHED
2026-07-27 14:02:44 INFO None 5285908: status FINISHED
2026-07-27 14:02:44 INFO None 5285914: status RUNNING/PENDING
2026-07-27 14:02:44 INFO None 5285920: status FINISHED
2026-07-27 14:02:44 INFO None 5285921: status FINISHED
2026-07-27 14:02:44 INFO None 5285922: status FINISHED
2026-07-27 14:02:44 INFO None 5285923: status FINISHED
2026-07-27 14:02:44 INFO None 5285924: status FINISHED
2026-07-27 14:02:44 INFO None 5285925: status FINISHED
2026-07-27 14:02:44 INFO None 5285930: status FINISHED
2026-07-27 14:02:44 INFO Jobs still running: ['5285914']. Waiting...
2026-07-27 14:02:59 INFO None 5285900: status FINISHED
2026-07-27 14:02:59 INFO None 5285901: status FINISHED
2026-07-27 14:02:59 INFO None 5285903: status FINISHED
2026-07-27 14:02:59 INFO None 5285905: status FINISHED
2026-07-27 14:02:59 INFO None 5285906: status FINISHED
2026-07-27 14:02:59 INFO None 5285907: status FINISHED
2026-07-27 14:02:59 INFO None 5285908: status FINISHED
2026-07-27 14:02:59 INFO None 5285914: status FINISHED
2026-07-27 14:02:59 INFO None 5285920: status FINISHED
2026-07-27 14:02:59 INFO None 5285921: status FINISHED
2026-07-27 14:02:59 INFO None 5285922: status FINISHED
2026-07-27 14:02:59 INFO None 5285923: status FINISHED
2026-07-27 14:02:59 INFO None 5285924: status FINISHED
2026-07-27 14:02:59 INFO None 5285925: status FINISHED
2026-07-27 14:03:00 INFO None 5285930: status FINISHED
2026-07-27 14:03:00 INFO Jobs ['5285900', '5285901', '5285903', '5285905', '5285906', '5285907', '5285908', '5285914', '5285920', '5285921', '5285922', '5285923', '5285924', '5285925', '5285930'] have finished
2026-07-27 14:03:00 INFO Checking restart files were created ...
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc(1002685915 bytes)
2026-07-27 14:03:00 INFO  Run_model() completed successfully.
2026-07-27 14:03:00 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:03:00 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 11:00:00 days=153072 seconds=39600
2026-07-27 14:03:00 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 14:03:00 INFO [TIME] increment current_time 2020-02-06 09:00:00 -> 2020-02-06 11:00:00
2026-07-27 14:03:00 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:03:00 INFO ---------->>> Running process_satellite_data()
2026-07-27 14:03:00 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12003.nc
2026-07-27 14:03:00 INFO ---------->>> Running run_obs_converter()
2026-07-27 14:03:00 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-27 14:03:00 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_39045_153072.out
2026-07-27 14:03:00 INFO ---------->>> Running DART
2026-07-27 14:03:00 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 14:03:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_1_out_toDART.nc
2026-07-27 14:03:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_1_out_toDART.nc
2026-07-27 14:03:00 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_1_out_toDART.nc
2026-07-27 14:03:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_1_out_toDART.nc
2026-07-27 14:03:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_1_out_toDART.nc
2026-07-27 14:03:01 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_1_out_toDART.nc
2026-07-27 14:03:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_1_out_toDART.nc
2026-07-27 14:03:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_1_out_toDART.nc
2026-07-27 14:03:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_1_out_toDART.nc
2026-07-27 14:03:02 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_1_out_toDART.nc
2026-07-27 14:03:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_1_out_toDART.nc
2026-07-27 14:03:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_1_out_toDART.nc
2026-07-27 14:03:03 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_1_out_toDART.nc
2026-07-27 14:03:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_1_out_toDART.nc
2026-07-27 14:03:04 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_1_out_toDART.nc
2026-07-27 14:03:04 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 14:03:04 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 14:03:04 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 14:03:04 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 14:03:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 14:03:04 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 14:03:21 INFO Found: []
2026-07-27 14:03:21 INFO No job id returned by command ./run_filter.bsh
2026-07-27 14:03:21 INFO No monitoring will be performed
2026-07-27 14:03:21 INFO Moving DART output files to analysis and preassim directories for date 2020020611 if present ...
2026-07-27 14:03:21 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:21 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:21 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020611'
2026-07-27 14:03:22 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:22 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020611'
2026-07-27 14:03:22 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 14:03:22 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 14:03:22 INFO run_dart() is DONE.
2026-07-27 14:03:22 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 14:03:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-27 14:03:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-27 14:03:32 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:32 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-27 14:03:37 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:38 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-27 14:03:43 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:43 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-27 14:03:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:48 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-27 14:03:53 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:53 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-27 14:03:58 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:03:59 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-27 14:04:04 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-27 14:04:09 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:10 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-27 14:04:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-27 14:04:20 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:20 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-27 14:04:25 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:26 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-27 14:04:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:31 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-27 14:04:36 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:36 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-27 14:04:41 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:04:41 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 14:04:41 INFO [TIME] step_end current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:04:41 INFO [TIME] step_start current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 11:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:04:41 INFO [TIME] window start=2020-02-06 11:00:00 end=2020-02-06 13:00:00 run_hours=2 has_assimilation=True
2026-07-27 14:04:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:43 INFO Hourly dataset computed and listing created
2026-07-27 14:04:48 INFO Hourly dataset computed
2026-07-27 14:04:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:48 INFO Hourly dataset computed and listing created
2026-07-27 14:04:49 INFO Hourly dataset computed
2026-07-27 14:04:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:50 INFO Hourly dataset computed and listing created
2026-07-27 14:04:51 INFO Hourly dataset computed
2026-07-27 14:04:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:52 INFO Hourly dataset computed and listing created
2026-07-27 14:04:53 INFO Hourly dataset computed
2026-07-27 14:04:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:54 INFO Hourly dataset computed and listing created
2026-07-27 14:04:54 INFO Hourly dataset computed
2026-07-27 14:04:54 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:55 INFO Hourly dataset computed and listing created
2026-07-27 14:04:56 INFO Hourly dataset computed
2026-07-27 14:04:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:57 INFO Hourly dataset computed and listing created
2026-07-27 14:04:58 INFO Hourly dataset computed
2026-07-27 14:04:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:04:59 INFO Hourly dataset computed and listing created
2026-07-27 14:05:00 INFO Hourly dataset computed
2026-07-27 14:05:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:01 INFO Hourly dataset computed and listing created
2026-07-27 14:05:01 INFO Hourly dataset computed
2026-07-27 14:05:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:02 INFO Hourly dataset computed and listing created
2026-07-27 14:05:03 INFO Hourly dataset computed
2026-07-27 14:05:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:04 INFO Hourly dataset computed and listing created
2026-07-27 14:05:05 INFO Hourly dataset computed
2026-07-27 14:05:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:06 INFO Hourly dataset computed and listing created
2026-07-27 14:05:07 INFO Hourly dataset computed
2026-07-27 14:05:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:07 INFO Hourly dataset computed and listing created
2026-07-27 14:05:08 INFO Hourly dataset computed
2026-07-27 14:05:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:09 INFO Hourly dataset computed and listing created
2026-07-27 14:05:10 INFO Hourly dataset computed
2026-07-27 14:05:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:05:11 INFO Hourly dataset computed and listing created
2026-07-27 14:05:12 INFO Hourly dataset computed
2026-07-27 14:05:12 INFO ---------->>> Running CHIMERE model from 2020-02-06 11:00:00 to 2020-02-06 13:00:00
2026-07-27 14:05:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 14:05:12 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020609_2_ENS1.nc
2026-07-27 14:05:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 14:05:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:12 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 14:05:12 INFO Queuing job for member 1...
2026-07-27 14:05:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:12 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 14:05:13 INFO Found: ['5285991']
2026-07-27 14:05:18 INFO [TGCC-IRENE] Submitted job with ID:['5285991']
2026-07-27 14:05:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 14:05:18 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020609_2_ENS2.nc
2026-07-27 14:05:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 14:05:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:18 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 14:05:18 INFO Queuing job for member 2...
2026-07-27 14:05:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:18 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 14:05:21 INFO Found: ['5285993']
2026-07-27 14:05:26 INFO [TGCC-IRENE] Submitted job with ID:['5285993']
2026-07-27 14:05:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 14:05:26 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020609_2_ENS3.nc
2026-07-27 14:05:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 14:05:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:26 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 14:05:26 INFO Queuing job for member 3...
2026-07-27 14:05:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:26 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 14:05:27 INFO Found: ['5285994']
2026-07-27 14:05:32 INFO [TGCC-IRENE] Submitted job with ID:['5285994']
2026-07-27 14:05:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 14:05:32 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020609_2_ENS4.nc
2026-07-27 14:05:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 14:05:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:32 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 14:05:32 INFO Queuing job for member 4...
2026-07-27 14:05:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:32 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 14:05:32 INFO Found: ['5285996']
2026-07-27 14:05:37 INFO [TGCC-IRENE] Submitted job with ID:['5285996']
2026-07-27 14:05:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 14:05:37 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020609_2_ENS5.nc
2026-07-27 14:05:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 14:05:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:37 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 14:05:37 INFO Queuing job for member 5...
2026-07-27 14:05:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:37 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 14:05:38 INFO Found: ['5285997']
2026-07-27 14:05:43 INFO [TGCC-IRENE] Submitted job with ID:['5285997']
2026-07-27 14:05:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 14:05:43 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020609_2_ENS6.nc
2026-07-27 14:05:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 14:05:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:43 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 14:05:43 INFO Queuing job for member 6...
2026-07-27 14:05:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:43 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 14:05:44 INFO Found: ['5285998']
2026-07-27 14:05:49 INFO [TGCC-IRENE] Submitted job with ID:['5285998']
2026-07-27 14:05:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 14:05:49 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020609_2_ENS7.nc
2026-07-27 14:05:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 14:05:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:49 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 14:05:49 INFO Queuing job for member 7...
2026-07-27 14:05:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:49 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 14:05:50 INFO Found: ['5286000']
2026-07-27 14:05:55 INFO [TGCC-IRENE] Submitted job with ID:['5286000']
2026-07-27 14:05:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:05:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 14:05:55 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020609_2_ENS8.nc
2026-07-27 14:05:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 14:05:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:05:55 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 14:05:55 INFO Queuing job for member 8...
2026-07-27 14:05:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:05:55 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 14:05:55 INFO Found: ['5286001']
2026-07-27 14:06:00 INFO [TGCC-IRENE] Submitted job with ID:['5286001']
2026-07-27 14:06:00 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:00 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 14:06:00 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020609_2_ENS9.nc
2026-07-27 14:06:00 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 14:06:00 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:00 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 14:06:00 INFO Queuing job for member 9...
2026-07-27 14:06:00 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:00 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 14:06:01 INFO Found: ['5286003']
2026-07-27 14:06:06 INFO [TGCC-IRENE] Submitted job with ID:['5286003']
2026-07-27 14:06:06 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:06 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 14:06:06 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020609_2_ENS10.nc
2026-07-27 14:06:06 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 14:06:06 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:06 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 14:06:06 INFO Queuing job for member 10...
2026-07-27 14:06:06 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:06 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 14:06:07 INFO Found: ['5286005']
2026-07-27 14:06:12 INFO [TGCC-IRENE] Submitted job with ID:['5286005']
2026-07-27 14:06:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 14:06:12 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020609_2_ENS11.nc
2026-07-27 14:06:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 14:06:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:12 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 14:06:12 INFO Queuing job for member 11...
2026-07-27 14:06:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:12 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 14:06:14 INFO Found: ['5286006']
2026-07-27 14:06:19 INFO [TGCC-IRENE] Submitted job with ID:['5286006']
2026-07-27 14:06:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 14:06:19 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020609_2_ENS12.nc
2026-07-27 14:06:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 14:06:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:19 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 14:06:19 INFO Queuing job for member 12...
2026-07-27 14:06:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:19 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 14:06:21 INFO Found: ['5286008']
2026-07-27 14:06:26 INFO [TGCC-IRENE] Submitted job with ID:['5286008']
2026-07-27 14:06:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 14:06:26 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020609_2_ENS13.nc
2026-07-27 14:06:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 14:06:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:26 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 14:06:26 INFO Queuing job for member 13...
2026-07-27 14:06:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:26 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 14:06:29 INFO Found: ['5286009']
2026-07-27 14:06:34 INFO [TGCC-IRENE] Submitted job with ID:['5286009']
2026-07-27 14:06:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 14:06:34 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020609_2_ENS14.nc
2026-07-27 14:06:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 14:06:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:34 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 14:06:34 INFO Queuing job for member 14...
2026-07-27 14:06:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:34 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 14:06:36 INFO Found: ['5286011']
2026-07-27 14:06:41 INFO [TGCC-IRENE] Submitted job with ID:['5286011']
2026-07-27 14:06:41 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:06:41 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 14:06:41 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020609_2_ENS15.nc
2026-07-27 14:06:41 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 14:06:41 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:06:41 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 14:06:41 INFO Queuing job for member 15...
2026-07-27 14:06:41 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:06:41 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 14:06:44 INFO Found: ['5286012']
2026-07-27 14:06:49 INFO [TGCC-IRENE] Submitted job with ID:['5286012']
2026-07-27 14:06:49 INFO Checking job status ...
2026-07-27 14:06:49 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:06:49 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:06:49 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:07:04 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:07:04 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:07:06 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:07:06 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:07:21 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:07:21 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:07:21 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:07:22 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:07:22 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:07:37 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:07:37 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:07:37 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:07:52 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:07:52 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:07:52 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:08:09 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:08:09 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:08:09 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:08:24 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:08:24 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:08:25 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:08:25 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:08:40 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:08:40 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:08:40 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:08:40 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:08:40 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:08:40 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:08:42 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:08:42 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:08:57 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:08:57 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:08:57 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:09:12 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:09:12 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:09:12 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:09:27 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:09:27 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:09:27 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286003: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286005: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:09:28 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:09:28 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:09:43 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286003: status FINISHED
2026-07-27 14:09:43 INFO None 5286005: status FINISHED
2026-07-27 14:09:43 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:09:43 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:09:43 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:09:58 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5285998: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286000: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286003: status FINISHED
2026-07-27 14:09:58 INFO None 5286005: status FINISHED
2026-07-27 14:09:58 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:09:58 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:10:00 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:10:00 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:10:15 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:10:15 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:10:15 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:10:15 INFO None 5285996: status RUNNING/PENDING
2026-07-27 14:10:15 INFO None 5285997: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5285998: status FINISHED
2026-07-27 14:10:16 INFO None 5286000: status FINISHED
2026-07-27 14:10:16 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5286003: status FINISHED
2026-07-27 14:10:16 INFO None 5286005: status FINISHED
2026-07-27 14:10:16 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:10:16 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:10:16 INFO Jobs still running: ['5285991', '5285993', '5285994', '5285996', '5285997', '5286001', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:10:31 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:10:31 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:10:31 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:10:31 INFO None 5285996: status FINISHED
2026-07-27 14:10:31 INFO None 5285997: status FINISHED
2026-07-27 14:10:31 INFO None 5285998: status FINISHED
2026-07-27 14:10:31 INFO None 5286000: status FINISHED
2026-07-27 14:10:31 INFO None 5286001: status RUNNING/PENDING
2026-07-27 14:10:31 INFO None 5286003: status FINISHED
2026-07-27 14:10:32 INFO None 5286005: status FINISHED
2026-07-27 14:10:32 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:10:32 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:10:32 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:10:32 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:10:32 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:10:32 INFO Jobs still running: ['5285991', '5285993', '5285994', '5286001', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:10:47 INFO None 5285991: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5285996: status FINISHED
2026-07-27 14:10:47 INFO None 5285997: status FINISHED
2026-07-27 14:10:47 INFO None 5285998: status FINISHED
2026-07-27 14:10:47 INFO None 5286000: status FINISHED
2026-07-27 14:10:47 INFO None 5286001: status FINISHED
2026-07-27 14:10:47 INFO None 5286003: status FINISHED
2026-07-27 14:10:47 INFO None 5286005: status FINISHED
2026-07-27 14:10:47 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:10:47 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:10:47 INFO Jobs still running: ['5285991', '5285993', '5285994', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:11:02 INFO None 5285991: status FINISHED
2026-07-27 14:11:02 INFO None 5285993: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5285996: status FINISHED
2026-07-27 14:11:02 INFO None 5285997: status FINISHED
2026-07-27 14:11:02 INFO None 5285998: status FINISHED
2026-07-27 14:11:02 INFO None 5286000: status FINISHED
2026-07-27 14:11:02 INFO None 5286001: status FINISHED
2026-07-27 14:11:02 INFO None 5286003: status FINISHED
2026-07-27 14:11:02 INFO None 5286005: status FINISHED
2026-07-27 14:11:02 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5286009: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5286011: status RUNNING/PENDING
2026-07-27 14:11:02 INFO None 5286012: status RUNNING/PENDING
2026-07-27 14:11:02 INFO Jobs still running: ['5285993', '5285994', '5286006', '5286008', '5286009', '5286011', '5286012']. Waiting...
2026-07-27 14:11:19 INFO None 5285991: status FINISHED
2026-07-27 14:11:19 INFO None 5285993: status FINISHED
2026-07-27 14:11:19 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:11:19 INFO None 5285996: status FINISHED
2026-07-27 14:11:19 INFO None 5285997: status FINISHED
2026-07-27 14:11:19 INFO None 5285998: status FINISHED
2026-07-27 14:11:19 INFO None 5286000: status FINISHED
2026-07-27 14:11:19 INFO None 5286001: status FINISHED
2026-07-27 14:11:19 INFO None 5286003: status FINISHED
2026-07-27 14:11:19 INFO None 5286005: status FINISHED
2026-07-27 14:11:19 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:11:19 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:11:19 INFO None 5286009: status FINISHED
2026-07-27 14:11:19 INFO None 5286011: status FINISHED
2026-07-27 14:11:19 INFO None 5286012: status FINISHED
2026-07-27 14:11:19 INFO Jobs still running: ['5285994', '5286006', '5286008']. Waiting...
2026-07-27 14:11:34 INFO None 5285991: status FINISHED
2026-07-27 14:11:34 INFO None 5285993: status FINISHED
2026-07-27 14:11:34 INFO None 5285994: status RUNNING/PENDING
2026-07-27 14:11:34 INFO None 5285996: status FINISHED
2026-07-27 14:11:34 INFO None 5285997: status FINISHED
2026-07-27 14:11:34 INFO None 5285998: status FINISHED
2026-07-27 14:11:34 INFO None 5286000: status FINISHED
2026-07-27 14:11:34 INFO None 5286001: status FINISHED
2026-07-27 14:11:34 INFO None 5286003: status FINISHED
2026-07-27 14:11:34 INFO None 5286005: status FINISHED
2026-07-27 14:11:34 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:11:36 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:11:36 INFO None 5286009: status FINISHED
2026-07-27 14:11:36 INFO None 5286011: status FINISHED
2026-07-27 14:11:36 INFO None 5286012: status FINISHED
2026-07-27 14:11:36 INFO Jobs still running: ['5285994', '5286006', '5286008']. Waiting...
2026-07-27 14:11:51 INFO None 5285991: status FINISHED
2026-07-27 14:11:51 INFO None 5285993: status FINISHED
2026-07-27 14:11:51 INFO None 5285994: status FINISHED
2026-07-27 14:11:51 INFO None 5285996: status FINISHED
2026-07-27 14:11:51 INFO None 5285997: status FINISHED
2026-07-27 14:11:51 INFO None 5285998: status FINISHED
2026-07-27 14:11:51 INFO None 5286000: status FINISHED
2026-07-27 14:11:51 INFO None 5286001: status FINISHED
2026-07-27 14:11:51 INFO None 5286003: status FINISHED
2026-07-27 14:11:52 INFO None 5286005: status FINISHED
2026-07-27 14:11:52 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:11:52 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:11:52 INFO None 5286009: status FINISHED
2026-07-27 14:11:52 INFO None 5286011: status FINISHED
2026-07-27 14:11:52 INFO None 5286012: status FINISHED
2026-07-27 14:11:52 INFO Jobs still running: ['5286006', '5286008']. Waiting...
2026-07-27 14:12:07 INFO None 5285991: status FINISHED
2026-07-27 14:12:07 INFO None 5285993: status FINISHED
2026-07-27 14:12:07 INFO None 5285994: status FINISHED
2026-07-27 14:12:07 INFO None 5285996: status FINISHED
2026-07-27 14:12:07 INFO None 5285997: status FINISHED
2026-07-27 14:12:07 INFO None 5285998: status FINISHED
2026-07-27 14:12:07 INFO None 5286000: status FINISHED
2026-07-27 14:12:07 INFO None 5286001: status FINISHED
2026-07-27 14:12:07 INFO None 5286003: status FINISHED
2026-07-27 14:12:07 INFO None 5286005: status FINISHED
2026-07-27 14:12:07 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:12:07 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:12:07 INFO None 5286009: status FINISHED
2026-07-27 14:12:07 INFO None 5286011: status FINISHED
2026-07-27 14:12:07 INFO None 5286012: status FINISHED
2026-07-27 14:12:07 INFO Jobs still running: ['5286006', '5286008']. Waiting...
2026-07-27 14:12:22 INFO None 5285991: status FINISHED
2026-07-27 14:12:22 INFO None 5285993: status FINISHED
2026-07-27 14:12:22 INFO None 5285994: status FINISHED
2026-07-27 14:12:22 INFO None 5285996: status FINISHED
2026-07-27 14:12:22 INFO None 5285997: status FINISHED
2026-07-27 14:12:22 INFO None 5285998: status FINISHED
2026-07-27 14:12:22 INFO None 5286000: status FINISHED
2026-07-27 14:12:22 INFO None 5286001: status FINISHED
2026-07-27 14:12:22 INFO None 5286003: status FINISHED
2026-07-27 14:12:22 INFO None 5286005: status FINISHED
2026-07-27 14:12:22 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:12:22 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:12:22 INFO None 5286009: status FINISHED
2026-07-27 14:12:22 INFO None 5286011: status FINISHED
2026-07-27 14:12:22 INFO None 5286012: status FINISHED
2026-07-27 14:12:22 INFO Jobs still running: ['5286006', '5286008']. Waiting...
2026-07-27 14:12:37 INFO None 5285991: status FINISHED
2026-07-27 14:12:37 INFO None 5285993: status FINISHED
2026-07-27 14:12:37 INFO None 5285994: status FINISHED
2026-07-27 14:12:37 INFO None 5285996: status FINISHED
2026-07-27 14:12:37 INFO None 5285997: status FINISHED
2026-07-27 14:12:37 INFO None 5285998: status FINISHED
2026-07-27 14:12:37 INFO None 5286000: status FINISHED
2026-07-27 14:12:37 INFO None 5286001: status FINISHED
2026-07-27 14:12:37 INFO None 5286003: status FINISHED
2026-07-27 14:12:37 INFO None 5286005: status FINISHED
2026-07-27 14:12:37 INFO None 5286006: status RUNNING/PENDING
2026-07-27 14:12:37 INFO None 5286008: status RUNNING/PENDING
2026-07-27 14:12:37 INFO None 5286009: status FINISHED
2026-07-27 14:12:37 INFO None 5286011: status FINISHED
2026-07-27 14:12:37 INFO None 5286012: status FINISHED
2026-07-27 14:12:37 INFO Jobs still running: ['5286006', '5286008']. Waiting...
2026-07-27 14:12:52 INFO None 5285991: status FINISHED
2026-07-27 14:12:52 INFO None 5285993: status FINISHED
2026-07-27 14:12:52 INFO None 5285994: status FINISHED
2026-07-27 14:12:53 INFO None 5285996: status FINISHED
2026-07-27 14:12:53 INFO None 5285997: status FINISHED
2026-07-27 14:12:53 INFO None 5285998: status FINISHED
2026-07-27 14:12:53 INFO None 5286000: status FINISHED
2026-07-27 14:12:53 INFO None 5286001: status FINISHED
2026-07-27 14:12:53 INFO None 5286003: status FINISHED
2026-07-27 14:12:53 INFO None 5286005: status FINISHED
2026-07-27 14:12:53 INFO None 5286006: status FINISHED
2026-07-27 14:12:53 INFO None 5286008: status FINISHED
2026-07-27 14:12:53 INFO None 5286009: status FINISHED
2026-07-27 14:12:53 INFO None 5286011: status FINISHED
2026-07-27 14:12:53 INFO None 5286012: status FINISHED
2026-07-27 14:12:53 INFO Jobs ['5285991', '5285993', '5285994', '5285996', '5285997', '5285998', '5286000', '5286001', '5286003', '5286005', '5286006', '5286008', '5286009', '5286011', '5286012'] have finished
2026-07-27 14:12:53 INFO Checking restart files were created ...
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc(1002685915 bytes)
2026-07-27 14:12:53 INFO  Run_model() completed successfully.
2026-07-27 14:12:53 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 11:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:12:53 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 13:00:00 days=153072 seconds=46800
2026-07-27 14:12:53 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 14:12:53 INFO [TIME] increment current_time 2020-02-06 11:00:00 -> 2020-02-06 13:00:00
2026-07-27 14:12:53 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:12:53 INFO ---------->>> Running process_satellite_data()
2026-07-27 14:12:53 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12004.nc
2026-07-27 14:12:53 INFO ---------->>> Running run_obs_converter()
2026-07-27 14:12:53 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-27 14:12:53 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_45135_153072.out
2026-07-27 14:12:53 INFO ---------->>> Running DART
2026-07-27 14:12:53 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 14:12:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out_toDART.nc
2026-07-27 14:12:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out_toDART.nc
2026-07-27 14:12:53 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out_toDART.nc
2026-07-27 14:12:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out_toDART.nc
2026-07-27 14:12:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out_toDART.nc
2026-07-27 14:12:54 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out_toDART.nc
2026-07-27 14:12:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out_toDART.nc
2026-07-27 14:12:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out_toDART.nc
2026-07-27 14:12:55 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out_toDART.nc
2026-07-27 14:12:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out_toDART.nc
2026-07-27 14:12:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out_toDART.nc
2026-07-27 14:12:56 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out_toDART.nc
2026-07-27 14:12:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out_toDART.nc
2026-07-27 14:12:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out_toDART.nc
2026-07-27 14:12:57 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020611_2_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out_toDART.nc
2026-07-27 14:12:57 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 14:12:57 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 14:12:57 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 14:12:57 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 14:12:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 14:12:57 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 14:13:16 INFO Found: []
2026-07-27 14:13:16 INFO No job id returned by command ./run_filter.bsh
2026-07-27 14:13:16 INFO No monitoring will be performed
2026-07-27 14:13:16 INFO Moving DART output files to analysis and preassim directories for date 2020020613 if present ...
2026-07-27 14:13:16 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:16 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:16 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:16 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020613'
2026-07-27 14:13:17 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 14:13:17 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 14:13:17 INFO run_dart() is DONE.
2026-07-27 14:13:17 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 14:13:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-27 14:13:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-27 14:13:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-27 14:13:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:33 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-27 14:13:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:39 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-27 14:13:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:44 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-27 14:13:49 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:50 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-27 14:13:55 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:13:55 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-27 14:14:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:01 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-27 14:14:06 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:06 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-27 14:14:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-27 14:14:16 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:17 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-27 14:14:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:22 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-27 14:14:27 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:28 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-27 14:14:33 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-27 14:14:39 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:14:39 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 14:14:39 INFO [TIME] step_end current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:14:39 INFO [TIME] step_start current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 13:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:14:39 INFO [TIME] window start=2020-02-06 13:00:00 end=2020-02-06 14:00:00 run_hours=1 has_assimilation=True
2026-07-27 14:14:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:40 INFO Hourly dataset computed and listing created
2026-07-27 14:14:44 INFO Hourly dataset computed
2026-07-27 14:14:44 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:45 INFO Hourly dataset computed and listing created
2026-07-27 14:14:48 INFO Hourly dataset computed
2026-07-27 14:14:48 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:49 INFO Hourly dataset computed and listing created
2026-07-27 14:14:50 INFO Hourly dataset computed
2026-07-27 14:14:50 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:51 INFO Hourly dataset computed and listing created
2026-07-27 14:14:53 INFO Hourly dataset computed
2026-07-27 14:14:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:54 INFO Hourly dataset computed and listing created
2026-07-27 14:14:56 INFO Hourly dataset computed
2026-07-27 14:14:56 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:57 INFO Hourly dataset computed and listing created
2026-07-27 14:14:58 INFO Hourly dataset computed
2026-07-27 14:14:58 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:14:58 INFO Hourly dataset computed and listing created
2026-07-27 14:14:59 INFO Hourly dataset computed
2026-07-27 14:14:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:00 INFO Hourly dataset computed and listing created
2026-07-27 14:15:00 INFO Hourly dataset computed
2026-07-27 14:15:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:01 INFO Hourly dataset computed and listing created
2026-07-27 14:15:02 INFO Hourly dataset computed
2026-07-27 14:15:02 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:03 INFO Hourly dataset computed and listing created
2026-07-27 14:15:03 INFO Hourly dataset computed
2026-07-27 14:15:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:04 INFO Hourly dataset computed and listing created
2026-07-27 14:15:05 INFO Hourly dataset computed
2026-07-27 14:15:05 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:06 INFO Hourly dataset computed and listing created
2026-07-27 14:15:06 INFO Hourly dataset computed
2026-07-27 14:15:06 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:07 INFO Hourly dataset computed and listing created
2026-07-27 14:15:08 INFO Hourly dataset computed
2026-07-27 14:15:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:09 INFO Hourly dataset computed and listing created
2026-07-27 14:15:09 INFO Hourly dataset computed
2026-07-27 14:15:09 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:15:10 INFO Hourly dataset computed and listing created
2026-07-27 14:15:11 INFO Hourly dataset computed
2026-07-27 14:15:11 INFO ---------->>> Running CHIMERE model from 2020-02-06 13:00:00 to 2020-02-06 14:00:00
2026-07-27 14:15:11 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:11 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 14:15:11 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020611_2_ENS1.nc
2026-07-27 14:15:11 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 14:15:11 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:11 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 14:15:11 INFO Queuing job for member 1...
2026-07-27 14:15:11 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:11 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 14:15:14 INFO Found: ['5286144']
2026-07-27 14:15:19 INFO [TGCC-IRENE] Submitted job with ID:['5286144']
2026-07-27 14:15:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 14:15:19 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020611_2_ENS2.nc
2026-07-27 14:15:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 14:15:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:19 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 14:15:19 INFO Queuing job for member 2...
2026-07-27 14:15:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:19 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 14:15:21 INFO Found: ['5286146']
2026-07-27 14:15:26 INFO [TGCC-IRENE] Submitted job with ID:['5286146']
2026-07-27 14:15:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 14:15:26 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020611_2_ENS3.nc
2026-07-27 14:15:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 14:15:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:26 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 14:15:26 INFO Queuing job for member 3...
2026-07-27 14:15:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:26 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 14:15:29 INFO Found: ['5286147']
2026-07-27 14:15:34 INFO [TGCC-IRENE] Submitted job with ID:['5286147']
2026-07-27 14:15:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 14:15:34 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020611_2_ENS4.nc
2026-07-27 14:15:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 14:15:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:34 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 14:15:34 INFO Queuing job for member 4...
2026-07-27 14:15:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:34 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 14:15:34 INFO Found: ['5286148']
2026-07-27 14:15:39 INFO [TGCC-IRENE] Submitted job with ID:['5286148']
2026-07-27 14:15:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 14:15:39 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020611_2_ENS5.nc
2026-07-27 14:15:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 14:15:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:39 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 14:15:39 INFO Queuing job for member 5...
2026-07-27 14:15:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:39 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 14:15:40 INFO Found: ['5286150']
2026-07-27 14:15:45 INFO [TGCC-IRENE] Submitted job with ID:['5286150']
2026-07-27 14:15:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 14:15:45 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020611_2_ENS6.nc
2026-07-27 14:15:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 14:15:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:45 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 14:15:45 INFO Queuing job for member 6...
2026-07-27 14:15:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:45 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 14:15:46 INFO Found: ['5286151']
2026-07-27 14:15:51 INFO [TGCC-IRENE] Submitted job with ID:['5286151']
2026-07-27 14:15:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 14:15:51 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020611_2_ENS7.nc
2026-07-27 14:15:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 14:15:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:51 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 14:15:51 INFO Queuing job for member 7...
2026-07-27 14:15:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:51 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 14:15:52 INFO Found: ['5286152']
2026-07-27 14:15:57 INFO [TGCC-IRENE] Submitted job with ID:['5286152']
2026-07-27 14:15:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:15:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 14:15:57 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020611_2_ENS8.nc
2026-07-27 14:15:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 14:15:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:15:57 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 14:15:57 INFO Queuing job for member 8...
2026-07-27 14:15:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:15:57 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 14:15:57 INFO Found: ['5286154']
2026-07-27 14:16:02 INFO [TGCC-IRENE] Submitted job with ID:['5286154']
2026-07-27 14:16:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 14:16:02 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020611_2_ENS9.nc
2026-07-27 14:16:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 14:16:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:02 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 14:16:02 INFO Queuing job for member 9...
2026-07-27 14:16:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:02 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 14:16:03 INFO Found: ['5286156']
2026-07-27 14:16:08 INFO [TGCC-IRENE] Submitted job with ID:['5286156']
2026-07-27 14:16:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 14:16:08 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020611_2_ENS10.nc
2026-07-27 14:16:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 14:16:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:08 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 14:16:08 INFO Queuing job for member 10...
2026-07-27 14:16:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:08 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 14:16:09 INFO Found: ['5286157']
2026-07-27 14:16:14 INFO [TGCC-IRENE] Submitted job with ID:['5286157']
2026-07-27 14:16:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 14:16:14 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020611_2_ENS11.nc
2026-07-27 14:16:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 14:16:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:14 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 14:16:14 INFO Queuing job for member 11...
2026-07-27 14:16:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:14 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 14:16:14 INFO Found: ['5286158']
2026-07-27 14:16:19 INFO [TGCC-IRENE] Submitted job with ID:['5286158']
2026-07-27 14:16:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 14:16:19 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020611_2_ENS12.nc
2026-07-27 14:16:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 14:16:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:19 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 14:16:20 INFO Queuing job for member 12...
2026-07-27 14:16:20 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:20 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 14:16:20 INFO Found: ['5286160']
2026-07-27 14:16:25 INFO [TGCC-IRENE] Submitted job with ID:['5286160']
2026-07-27 14:16:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 14:16:25 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020611_2_ENS13.nc
2026-07-27 14:16:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 14:16:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:25 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 14:16:25 INFO Queuing job for member 13...
2026-07-27 14:16:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:25 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 14:16:27 INFO Found: ['5286161']
2026-07-27 14:16:32 INFO [TGCC-IRENE] Submitted job with ID:['5286161']
2026-07-27 14:16:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 14:16:32 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020611_2_ENS14.nc
2026-07-27 14:16:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 14:16:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:32 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 14:16:32 INFO Queuing job for member 14...
2026-07-27 14:16:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:32 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 14:16:35 INFO Found: ['5286164']
2026-07-27 14:16:40 INFO [TGCC-IRENE] Submitted job with ID:['5286164']
2026-07-27 14:16:40 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:16:40 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 14:16:40 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020611_2_ENS15.nc
2026-07-27 14:16:40 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 14:16:40 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:16:40 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 14:16:40 INFO Queuing job for member 15...
2026-07-27 14:16:40 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:16:40 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 14:16:42 INFO Found: ['5286170']
2026-07-27 14:16:47 INFO [TGCC-IRENE] Submitted job with ID:['5286170']
2026-07-27 14:16:47 INFO Checking job status ...
2026-07-27 14:16:47 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:16:48 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:16:48 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:17:03 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:17:03 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:17:03 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:17:05 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:17:05 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:17:20 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:17:20 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:17:21 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:17:21 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:17:36 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:17:36 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:17:36 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:17:51 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:17:51 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:17:51 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:18:08 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:18:08 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:18:08 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:18:23 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:18:23 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:18:24 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:18:26 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:18:26 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:18:26 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:18:26 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:18:26 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:18:41 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:18:41 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:18:41 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:18:56 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:18:56 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:18:58 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:18:58 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:19:13 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286148: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286150: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286154: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:19:13 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:19:13 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:19:28 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:19:28 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:19:28 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:19:28 INFO None 5286148: status FINISHED
2026-07-27 14:19:28 INFO None 5286150: status FINISHED
2026-07-27 14:19:28 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:19:28 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:19:28 INFO None 5286154: status FINISHED
2026-07-27 14:19:29 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:19:29 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:19:29 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286151', '5286152', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:19:45 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286148: status FINISHED
2026-07-27 14:19:45 INFO None 5286150: status FINISHED
2026-07-27 14:19:45 INFO None 5286151: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286152: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286154: status FINISHED
2026-07-27 14:19:45 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:19:45 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:19:45 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286151', '5286152', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:20:00 INFO None 5286144: status RUNNING/PENDING
2026-07-27 14:20:00 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:20:00 INFO None 5286147: status RUNNING/PENDING
2026-07-27 14:20:00 INFO None 5286148: status FINISHED
2026-07-27 14:20:00 INFO None 5286150: status FINISHED
2026-07-27 14:20:00 INFO None 5286151: status FINISHED
2026-07-27 14:20:00 INFO None 5286152: status FINISHED
2026-07-27 14:20:03 INFO None 5286154: status FINISHED
2026-07-27 14:20:03 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:20:03 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:20:03 INFO Jobs still running: ['5286144', '5286146', '5286147', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:20:18 INFO None 5286144: status FINISHED
2026-07-27 14:20:18 INFO None 5286146: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286147: status FINISHED
2026-07-27 14:20:18 INFO None 5286148: status FINISHED
2026-07-27 14:20:18 INFO None 5286150: status FINISHED
2026-07-27 14:20:18 INFO None 5286151: status FINISHED
2026-07-27 14:20:18 INFO None 5286152: status FINISHED
2026-07-27 14:20:18 INFO None 5286154: status FINISHED
2026-07-27 14:20:18 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286164: status RUNNING/PENDING
2026-07-27 14:20:18 INFO None 5286170: status RUNNING/PENDING
2026-07-27 14:20:18 INFO Jobs still running: ['5286146', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170']. Waiting...
2026-07-27 14:20:33 INFO None 5286144: status FINISHED
2026-07-27 14:20:33 INFO None 5286146: status FINISHED
2026-07-27 14:20:33 INFO None 5286147: status FINISHED
2026-07-27 14:20:33 INFO None 5286148: status FINISHED
2026-07-27 14:20:33 INFO None 5286150: status FINISHED
2026-07-27 14:20:33 INFO None 5286151: status FINISHED
2026-07-27 14:20:33 INFO None 5286152: status FINISHED
2026-07-27 14:20:33 INFO None 5286154: status FINISHED
2026-07-27 14:20:35 INFO None 5286156: status RUNNING/PENDING
2026-07-27 14:20:35 INFO None 5286157: status RUNNING/PENDING
2026-07-27 14:20:35 INFO None 5286158: status RUNNING/PENDING
2026-07-27 14:20:35 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:20:35 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:20:35 INFO None 5286164: status FINISHED
2026-07-27 14:20:35 INFO None 5286170: status FINISHED
2026-07-27 14:20:35 INFO Jobs still running: ['5286156', '5286157', '5286158', '5286160', '5286161']. Waiting...
2026-07-27 14:20:50 INFO None 5286144: status FINISHED
2026-07-27 14:20:50 INFO None 5286146: status FINISHED
2026-07-27 14:20:50 INFO None 5286147: status FINISHED
2026-07-27 14:20:50 INFO None 5286148: status FINISHED
2026-07-27 14:20:50 INFO None 5286150: status FINISHED
2026-07-27 14:20:50 INFO None 5286151: status FINISHED
2026-07-27 14:20:50 INFO None 5286152: status FINISHED
2026-07-27 14:20:50 INFO None 5286154: status FINISHED
2026-07-27 14:20:50 INFO None 5286156: status FINISHED
2026-07-27 14:20:50 INFO None 5286157: status FINISHED
2026-07-27 14:20:50 INFO None 5286158: status FINISHED
2026-07-27 14:20:51 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:20:51 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:20:51 INFO None 5286164: status FINISHED
2026-07-27 14:20:51 INFO None 5286170: status FINISHED
2026-07-27 14:20:51 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:21:06 INFO None 5286144: status FINISHED
2026-07-27 14:21:06 INFO None 5286146: status FINISHED
2026-07-27 14:21:06 INFO None 5286147: status FINISHED
2026-07-27 14:21:06 INFO None 5286148: status FINISHED
2026-07-27 14:21:06 INFO None 5286150: status FINISHED
2026-07-27 14:21:06 INFO None 5286151: status FINISHED
2026-07-27 14:21:06 INFO None 5286152: status FINISHED
2026-07-27 14:21:06 INFO None 5286154: status FINISHED
2026-07-27 14:21:06 INFO None 5286156: status FINISHED
2026-07-27 14:21:06 INFO None 5286157: status FINISHED
2026-07-27 14:21:06 INFO None 5286158: status FINISHED
2026-07-27 14:21:06 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:21:06 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:21:06 INFO None 5286164: status FINISHED
2026-07-27 14:21:06 INFO None 5286170: status FINISHED
2026-07-27 14:21:06 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:21:21 INFO None 5286144: status FINISHED
2026-07-27 14:21:21 INFO None 5286146: status FINISHED
2026-07-27 14:21:21 INFO None 5286147: status FINISHED
2026-07-27 14:21:21 INFO None 5286148: status FINISHED
2026-07-27 14:21:21 INFO None 5286150: status FINISHED
2026-07-27 14:21:21 INFO None 5286151: status FINISHED
2026-07-27 14:21:21 INFO None 5286152: status FINISHED
2026-07-27 14:21:21 INFO None 5286154: status FINISHED
2026-07-27 14:21:21 INFO None 5286156: status FINISHED
2026-07-27 14:21:21 INFO None 5286157: status FINISHED
2026-07-27 14:21:21 INFO None 5286158: status FINISHED
2026-07-27 14:21:21 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:21:21 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:21:21 INFO None 5286164: status FINISHED
2026-07-27 14:21:21 INFO None 5286170: status FINISHED
2026-07-27 14:21:21 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:21:38 INFO None 5286144: status FINISHED
2026-07-27 14:21:38 INFO None 5286146: status FINISHED
2026-07-27 14:21:38 INFO None 5286147: status FINISHED
2026-07-27 14:21:38 INFO None 5286148: status FINISHED
2026-07-27 14:21:38 INFO None 5286150: status FINISHED
2026-07-27 14:21:38 INFO None 5286151: status FINISHED
2026-07-27 14:21:38 INFO None 5286152: status FINISHED
2026-07-27 14:21:38 INFO None 5286154: status FINISHED
2026-07-27 14:21:38 INFO None 5286156: status FINISHED
2026-07-27 14:21:38 INFO None 5286157: status FINISHED
2026-07-27 14:21:38 INFO None 5286158: status FINISHED
2026-07-27 14:21:38 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:21:38 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:21:38 INFO None 5286164: status FINISHED
2026-07-27 14:21:38 INFO None 5286170: status FINISHED
2026-07-27 14:21:38 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:21:53 INFO None 5286144: status FINISHED
2026-07-27 14:21:53 INFO None 5286146: status FINISHED
2026-07-27 14:21:53 INFO None 5286147: status FINISHED
2026-07-27 14:21:53 INFO None 5286148: status FINISHED
2026-07-27 14:21:53 INFO None 5286150: status FINISHED
2026-07-27 14:21:54 INFO None 5286151: status FINISHED
2026-07-27 14:21:54 INFO None 5286152: status FINISHED
2026-07-27 14:21:54 INFO None 5286154: status FINISHED
2026-07-27 14:21:54 INFO None 5286156: status FINISHED
2026-07-27 14:21:54 INFO None 5286157: status FINISHED
2026-07-27 14:21:54 INFO None 5286158: status FINISHED
2026-07-27 14:21:54 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:21:54 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:21:54 INFO None 5286164: status FINISHED
2026-07-27 14:21:54 INFO None 5286170: status FINISHED
2026-07-27 14:21:54 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:22:09 INFO None 5286144: status FINISHED
2026-07-27 14:22:09 INFO None 5286146: status FINISHED
2026-07-27 14:22:09 INFO None 5286147: status FINISHED
2026-07-27 14:22:09 INFO None 5286148: status FINISHED
2026-07-27 14:22:09 INFO None 5286150: status FINISHED
2026-07-27 14:22:11 INFO None 5286151: status FINISHED
2026-07-27 14:22:11 INFO None 5286152: status FINISHED
2026-07-27 14:22:11 INFO None 5286154: status FINISHED
2026-07-27 14:22:11 INFO None 5286156: status FINISHED
2026-07-27 14:22:11 INFO None 5286157: status FINISHED
2026-07-27 14:22:11 INFO None 5286158: status FINISHED
2026-07-27 14:22:11 INFO None 5286160: status RUNNING/PENDING
2026-07-27 14:22:11 INFO None 5286161: status RUNNING/PENDING
2026-07-27 14:22:11 INFO None 5286164: status FINISHED
2026-07-27 14:22:11 INFO None 5286170: status FINISHED
2026-07-27 14:22:11 INFO Jobs still running: ['5286160', '5286161']. Waiting...
2026-07-27 14:22:26 INFO None 5286144: status FINISHED
2026-07-27 14:22:26 INFO None 5286146: status FINISHED
2026-07-27 14:22:26 INFO None 5286147: status FINISHED
2026-07-27 14:22:26 INFO None 5286148: status FINISHED
2026-07-27 14:22:26 INFO None 5286150: status FINISHED
2026-07-27 14:22:26 INFO None 5286151: status FINISHED
2026-07-27 14:22:26 INFO None 5286152: status FINISHED
2026-07-27 14:22:26 INFO None 5286154: status FINISHED
2026-07-27 14:22:26 INFO None 5286156: status FINISHED
2026-07-27 14:22:26 INFO None 5286157: status FINISHED
2026-07-27 14:22:26 INFO None 5286158: status FINISHED
2026-07-27 14:22:26 INFO None 5286160: status FINISHED
2026-07-27 14:22:26 INFO None 5286161: status FINISHED
2026-07-27 14:22:26 INFO None 5286164: status FINISHED
2026-07-27 14:22:26 INFO None 5286170: status FINISHED
2026-07-27 14:22:26 INFO Jobs ['5286144', '5286146', '5286147', '5286148', '5286150', '5286151', '5286152', '5286154', '5286156', '5286157', '5286158', '5286160', '5286161', '5286164', '5286170'] have finished
2026-07-27 14:22:26 INFO Checking restart files were created ...
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc(668832435 bytes)
2026-07-27 14:22:26 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc(668832435 bytes)
2026-07-27 14:22:26 INFO  Run_model() completed successfully.
2026-07-27 14:22:26 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 13:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:22:26 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 14:00:00 days=153072 seconds=50400
2026-07-27 14:22:26 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 14:22:26 INFO [TIME] increment current_time 2020-02-06 13:00:00 -> 2020-02-06 14:00:00
2026-07-27 14:22:26 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:22:26 INFO ---------->>> Running process_satellite_data()
2026-07-27 14:22:26 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12005.nc
2026-07-27 14:22:26 INFO ---------->>> Running run_obs_converter()
2026-07-27 14:22:26 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-27 14:22:26 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_51225_153072.out
2026-07-27 14:22:26 INFO ---------->>> Running DART
2026-07-27 14:22:26 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-27 14:22:26 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/chim_ENS1_2020020614_1_out_toDART.nc
2026-07-27 14:22:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/chim_ENS2_2020020614_1_out_toDART.nc
2026-07-27 14:22:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/chim_ENS3_2020020614_1_out_toDART.nc
2026-07-27 14:22:27 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/chim_ENS4_2020020614_1_out_toDART.nc
2026-07-27 14:22:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/chim_ENS5_2020020614_1_out_toDART.nc
2026-07-27 14:22:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/chim_ENS6_2020020614_1_out_toDART.nc
2026-07-27 14:22:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/chim_ENS7_2020020614_1_out_toDART.nc
2026-07-27 14:22:28 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/chim_ENS8_2020020614_1_out_toDART.nc
2026-07-27 14:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/chim_ENS9_2020020614_1_out_toDART.nc
2026-07-27 14:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/chim_ENS10_2020020614_1_out_toDART.nc
2026-07-27 14:22:29 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/chim_ENS11_2020020614_1_out_toDART.nc
2026-07-27 14:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/chim_ENS12_2020020614_1_out_toDART.nc
2026-07-27 14:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/chim_ENS13_2020020614_1_out_toDART.nc
2026-07-27 14:22:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/chim_ENS14_2020020614_1_out_toDART.nc
2026-07-27 14:22:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020613_1_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/chim_ENS15_2020020614_1_out_toDART.nc
2026-07-27 14:22:31 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-27 14:22:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-27 14:22:31 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-27 14:22:31 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-27 14:22:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-27 14:22:31 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-27 14:22:40 INFO Found: []
2026-07-27 14:22:40 INFO No job id returned by command ./run_filter.bsh
2026-07-27 14:22:40 INFO No monitoring will be performed
2026-07-27 14:22:40 INFO Moving DART output files to analysis and preassim directories for date 2020020614 if present ...
2026-07-27 14:22:40 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:40 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_priorinf_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_priorinf_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/analysis/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyICinfl2_0607_15m_low_v2/preassim/2020020614'
2026-07-27 14:22:41 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-27 14:22:41 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-27 14:22:41 INFO run_dart() is DONE.
2026-07-27 14:22:41 INFO ---------->>> Running update_pollutant_in_end()
2026-07-27 14:22:41 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-27 14:22:44 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:22:45 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-27 14:22:48 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:22:49 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-27 14:22:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:22:52 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-27 14:22:56 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:22:56 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-27 14:23:00 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:00 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-27 14:23:03 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:04 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-27 14:23:07 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:07 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-27 14:23:11 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:11 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-27 14:23:15 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:15 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-27 14:23:18 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:19 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-27 14:23:22 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:23 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-27 14:23:26 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:27 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-27 14:23:30 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:30 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-27 14:23:34 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:34 INFO DART's updated NO2 successfully replaced into /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-27 14:23:38 INFO Computing differences between posterior vs. original CHIMERE outputs ...
2026-07-27 14:23:38 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 14:23:38 INFO [TIME] step_end current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:23:38 INFO [TIME] step_start current_time=2020-02-06 14:00:00 simulated_time=2020-02-06 14:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:23:38 INFO [TIME] window start=2020-02-06 14:00:00 end=2020-02-07 00:00:00 run_hours=10 has_assimilation=False
2026-07-27 14:23:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:23:40 INFO Hourly dataset computed and listing created
2026-07-27 14:24:01 INFO Hourly dataset computed
2026-07-27 14:24:01 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:24:02 INFO Hourly dataset computed and listing created
2026-07-27 14:24:21 INFO Hourly dataset computed
2026-07-27 14:24:21 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:24:22 INFO Hourly dataset computed and listing created
2026-07-27 14:24:40 INFO Hourly dataset computed
2026-07-27 14:24:40 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:24:42 INFO Hourly dataset computed and listing created
2026-07-27 14:24:59 INFO Hourly dataset computed
2026-07-27 14:24:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:25:00 INFO Hourly dataset computed and listing created
2026-07-27 14:25:17 INFO Hourly dataset computed
2026-07-27 14:25:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:25:19 INFO Hourly dataset computed and listing created
2026-07-27 14:25:34 INFO Hourly dataset computed
2026-07-27 14:25:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:25:36 INFO Hourly dataset computed and listing created
2026-07-27 14:25:53 INFO Hourly dataset computed
2026-07-27 14:25:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:25:55 INFO Hourly dataset computed and listing created
2026-07-27 14:26:12 INFO Hourly dataset computed
2026-07-27 14:26:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:26:14 INFO Hourly dataset computed and listing created
2026-07-27 14:26:31 INFO Hourly dataset computed
2026-07-27 14:26:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:26:33 INFO Hourly dataset computed and listing created
2026-07-27 14:26:49 INFO Hourly dataset computed
2026-07-27 14:26:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:26:51 INFO Hourly dataset computed and listing created
2026-07-27 14:27:07 INFO Hourly dataset computed
2026-07-27 14:27:07 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:27:08 INFO Hourly dataset computed and listing created
2026-07-27 14:27:25 INFO Hourly dataset computed
2026-07-27 14:27:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:27:27 INFO Hourly dataset computed and listing created
2026-07-27 14:27:45 INFO Hourly dataset computed
2026-07-27 14:27:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:27:46 INFO Hourly dataset computed and listing created
2026-07-27 14:28:04 INFO Hourly dataset computed
2026-07-27 14:28:04 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:28:05 INFO Hourly dataset computed and listing created
2026-07-27 14:28:22 INFO Hourly dataset computed
2026-07-27 14:28:22 INFO ---------->>> Running CHIMERE model from 2020-02-06 14:00:00 to 2020-02-07 00:00:00
2026-07-27 14:28:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 14:28:22 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020613_1_ENS1.nc
2026-07-27 14:28:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 14:28:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:22 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 14:28:22 INFO Queuing job for member 1...
2026-07-27 14:28:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:22 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 14:28:23 INFO Found: ['5286456']
2026-07-27 14:28:28 INFO [TGCC-IRENE] Submitted job with ID:['5286456']
2026-07-27 14:28:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 14:28:28 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020613_1_ENS2.nc
2026-07-27 14:28:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 14:28:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:28 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 14:28:28 INFO Queuing job for member 2...
2026-07-27 14:28:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:28 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 14:28:28 INFO Found: ['5286457']
2026-07-27 14:28:33 INFO [TGCC-IRENE] Submitted job with ID:['5286457']
2026-07-27 14:28:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 14:28:33 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020613_1_ENS3.nc
2026-07-27 14:28:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 14:28:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:33 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 14:28:33 INFO Queuing job for member 3...
2026-07-27 14:28:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:33 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 14:28:34 INFO Found: ['5286458']
2026-07-27 14:28:39 INFO [TGCC-IRENE] Submitted job with ID:['5286458']
2026-07-27 14:28:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 14:28:39 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020613_1_ENS4.nc
2026-07-27 14:28:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 14:28:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:39 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 14:28:39 INFO Queuing job for member 4...
2026-07-27 14:28:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:39 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 14:28:40 INFO Found: ['5286459']
2026-07-27 14:28:45 INFO [TGCC-IRENE] Submitted job with ID:['5286459']
2026-07-27 14:28:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 14:28:45 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020613_1_ENS5.nc
2026-07-27 14:28:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 14:28:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:45 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 14:28:45 INFO Queuing job for member 5...
2026-07-27 14:28:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:45 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 14:28:45 INFO Found: ['5286460']
2026-07-27 14:28:50 INFO [TGCC-IRENE] Submitted job with ID:['5286460']
2026-07-27 14:28:50 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:50 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 14:28:50 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020613_1_ENS6.nc
2026-07-27 14:28:50 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 14:28:50 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:50 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 14:28:50 INFO Queuing job for member 6...
2026-07-27 14:28:50 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:50 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 14:28:51 INFO Found: ['5286461']
2026-07-27 14:28:56 INFO [TGCC-IRENE] Submitted job with ID:['5286461']
2026-07-27 14:28:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:28:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 14:28:56 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020613_1_ENS7.nc
2026-07-27 14:28:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 14:28:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:28:56 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 14:28:56 INFO Queuing job for member 7...
2026-07-27 14:28:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:28:56 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 14:28:57 INFO Found: ['5286462']
2026-07-27 14:29:02 INFO [TGCC-IRENE] Submitted job with ID:['5286462']
2026-07-27 14:29:02 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:02 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 14:29:02 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020613_1_ENS8.nc
2026-07-27 14:29:02 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 14:29:02 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:02 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 14:29:02 INFO Queuing job for member 8...
2026-07-27 14:29:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:02 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 14:29:03 INFO Found: ['5286464']
2026-07-27 14:29:08 INFO [TGCC-IRENE] Submitted job with ID:['5286464']
2026-07-27 14:29:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 14:29:08 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020613_1_ENS9.nc
2026-07-27 14:29:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 14:29:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:08 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 14:29:08 INFO Queuing job for member 9...
2026-07-27 14:29:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:08 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 14:29:08 INFO Found: ['5286465']
2026-07-27 14:29:13 INFO [TGCC-IRENE] Submitted job with ID:['5286465']
2026-07-27 14:29:13 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:13 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 14:29:13 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020613_1_ENS10.nc
2026-07-27 14:29:13 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 14:29:13 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:13 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 14:29:13 INFO Queuing job for member 10...
2026-07-27 14:29:13 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:13 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 14:29:16 INFO Found: ['5286466']
2026-07-27 14:29:21 INFO [TGCC-IRENE] Submitted job with ID:['5286466']
2026-07-27 14:29:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 14:29:21 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020613_1_ENS11.nc
2026-07-27 14:29:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 14:29:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:21 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 14:29:21 INFO Queuing job for member 11...
2026-07-27 14:29:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:21 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 14:29:24 INFO Found: ['5286467']
2026-07-27 14:29:29 INFO [TGCC-IRENE] Submitted job with ID:['5286467']
2026-07-27 14:29:29 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:29 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 14:29:29 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020613_1_ENS12.nc
2026-07-27 14:29:29 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 14:29:29 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:29 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 14:29:29 INFO Queuing job for member 12...
2026-07-27 14:29:29 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:29 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 14:29:31 INFO Found: ['5286469']
2026-07-27 14:29:36 INFO [TGCC-IRENE] Submitted job with ID:['5286469']
2026-07-27 14:29:36 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:36 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 14:29:36 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020613_1_ENS13.nc
2026-07-27 14:29:36 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 14:29:36 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:36 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 14:29:36 INFO Queuing job for member 13...
2026-07-27 14:29:36 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:36 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 14:29:39 INFO Found: ['5286470']
2026-07-27 14:29:44 INFO [TGCC-IRENE] Submitted job with ID:['5286470']
2026-07-27 14:29:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 14:29:44 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020613_1_ENS14.nc
2026-07-27 14:29:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 14:29:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:44 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 14:29:44 INFO Queuing job for member 14...
2026-07-27 14:29:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:44 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 14:29:46 INFO Found: ['5286471']
2026-07-27 14:29:51 INFO [TGCC-IRENE] Submitted job with ID:['5286471']
2026-07-27 14:29:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:29:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 14:29:51 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020613_1_ENS15.nc
2026-07-27 14:29:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 14:29:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:29:51 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 14:29:51 INFO Queuing job for member 15...
2026-07-27 14:29:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:29:51 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 14:29:54 INFO Found: ['5286472']
2026-07-27 14:29:59 INFO [TGCC-IRENE] Submitted job with ID:['5286472']
2026-07-27 14:29:59 INFO Checking job status ...
2026-07-27 14:29:59 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:29:59 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:29:59 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:30:14 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:30:14 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:30:14 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:30:29 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:30:29 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:30:29 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:30:45 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:30:45 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:30:45 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:31:00 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:31:00 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:31:00 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:31:15 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:31:15 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:31:16 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:31:16 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:31:16 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:31:16 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:31:31 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:31:31 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:31:31 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:31:31 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:31:33 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:31:33 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:31:48 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:31:48 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:31:48 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:32:03 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:32:03 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:32:05 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:32:05 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:32:20 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:32:20 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:32:20 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:32:20 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:32:21 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:32:21 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:32:36 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:32:36 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:32:36 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:32:51 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:32:51 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:32:51 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:32:53 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:32:53 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:33:08 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:33:08 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:33:09 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:33:09 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:33:24 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:33:24 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:33:24 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:33:24 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:33:26 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:33:26 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:33:41 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:33:41 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:33:41 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:33:56 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:33:56 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:33:56 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:33:56 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:33:56 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:33:56 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:33:58 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:33:58 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:34:13 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:34:13 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:34:13 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:34:28 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:34:28 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:34:28 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:34:45 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:34:45 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:34:45 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:35:00 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:35:00 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:35:02 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:35:02 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:35:17 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:35:17 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:35:18 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:35:18 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:35:18 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:35:18 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:35:33 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:35:33 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:35:35 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:35:35 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:35:50 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:35:50 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:35:50 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:36:05 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:36:05 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:36:05 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:36:20 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:36:20 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:36:20 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:36:21 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:36:21 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:36:36 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:36:36 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:36:36 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:36:53 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:36:53 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:36:53 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:37:08 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:37:08 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:37:10 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:37:10 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:37:10 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:37:25 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:37:25 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:37:25 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:37:41 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:37:41 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:37:43 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:37:43 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:37:43 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:37:58 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:37:58 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:37:58 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:38:13 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:38:13 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:38:13 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:38:28 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:38:28 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:38:28 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:38:28 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:38:29 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:38:29 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:38:44 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:38:44 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:38:45 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:38:45 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:38:45 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:38:45 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:38:45 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:39:00 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286458: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286459: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286464: status RUNNING/PENDING
2026-07-27 14:39:00 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:39:02 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:39:02 INFO Jobs still running: ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:39:17 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286458: status FINISHED
2026-07-27 14:39:17 INFO None 5286459: status FINISHED
2026-07-27 14:39:17 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286464: status FINISHED
2026-07-27 14:39:17 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:39:17 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:39:17 INFO Jobs still running: ['5286456', '5286457', '5286460', '5286461', '5286462', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:39:32 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286457: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286458: status FINISHED
2026-07-27 14:39:32 INFO None 5286459: status FINISHED
2026-07-27 14:39:32 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286464: status FINISHED
2026-07-27 14:39:32 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:39:32 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:39:34 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:39:34 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:39:34 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:39:34 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:39:34 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:39:34 INFO Jobs still running: ['5286456', '5286457', '5286460', '5286461', '5286462', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:39:50 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286457: status FINISHED
2026-07-27 14:39:50 INFO None 5286458: status FINISHED
2026-07-27 14:39:50 INFO None 5286459: status FINISHED
2026-07-27 14:39:50 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286464: status FINISHED
2026-07-27 14:39:50 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:39:50 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:39:50 INFO Jobs still running: ['5286456', '5286460', '5286461', '5286462', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:40:05 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286457: status FINISHED
2026-07-27 14:40:05 INFO None 5286458: status FINISHED
2026-07-27 14:40:05 INFO None 5286459: status FINISHED
2026-07-27 14:40:05 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286464: status FINISHED
2026-07-27 14:40:05 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:40:05 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:40:05 INFO Jobs still running: ['5286456', '5286460', '5286461', '5286462', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:40:20 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286457: status FINISHED
2026-07-27 14:40:20 INFO None 5286458: status FINISHED
2026-07-27 14:40:20 INFO None 5286459: status FINISHED
2026-07-27 14:40:20 INFO None 5286460: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286462: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286464: status FINISHED
2026-07-27 14:40:20 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286467: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:40:20 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:40:20 INFO Jobs still running: ['5286456', '5286460', '5286461', '5286462', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:40:36 INFO None 5286456: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286457: status FINISHED
2026-07-27 14:40:36 INFO None 5286458: status FINISHED
2026-07-27 14:40:36 INFO None 5286459: status FINISHED
2026-07-27 14:40:36 INFO None 5286460: status FINISHED
2026-07-27 14:40:36 INFO None 5286461: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286462: status FINISHED
2026-07-27 14:40:36 INFO None 5286464: status FINISHED
2026-07-27 14:40:36 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286467: status FINISHED
2026-07-27 14:40:36 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:40:36 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:40:36 INFO Jobs still running: ['5286456', '5286461', '5286465', '5286466', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:40:51 INFO None 5286456: status FINISHED
2026-07-27 14:40:51 INFO None 5286457: status FINISHED
2026-07-27 14:40:51 INFO None 5286458: status FINISHED
2026-07-27 14:40:51 INFO None 5286459: status FINISHED
2026-07-27 14:40:51 INFO None 5286460: status FINISHED
2026-07-27 14:40:51 INFO None 5286461: status FINISHED
2026-07-27 14:40:52 INFO None 5286462: status FINISHED
2026-07-27 14:40:52 INFO None 5286464: status FINISHED
2026-07-27 14:40:52 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:40:52 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:40:52 INFO None 5286467: status FINISHED
2026-07-27 14:40:52 INFO None 5286469: status RUNNING/PENDING
2026-07-27 14:40:54 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:40:54 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:40:54 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:40:54 INFO Jobs still running: ['5286465', '5286466', '5286469', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:41:09 INFO None 5286456: status FINISHED
2026-07-27 14:41:09 INFO None 5286457: status FINISHED
2026-07-27 14:41:09 INFO None 5286458: status FINISHED
2026-07-27 14:41:09 INFO None 5286459: status FINISHED
2026-07-27 14:41:09 INFO None 5286460: status FINISHED
2026-07-27 14:41:09 INFO None 5286461: status FINISHED
2026-07-27 14:41:09 INFO None 5286462: status FINISHED
2026-07-27 14:41:09 INFO None 5286464: status FINISHED
2026-07-27 14:41:09 INFO None 5286465: status RUNNING/PENDING
2026-07-27 14:41:09 INFO None 5286466: status RUNNING/PENDING
2026-07-27 14:41:09 INFO None 5286467: status FINISHED
2026-07-27 14:41:09 INFO None 5286469: status FINISHED
2026-07-27 14:41:09 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:41:09 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:41:09 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:41:09 INFO Jobs still running: ['5286465', '5286466', '5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:41:24 INFO None 5286456: status FINISHED
2026-07-27 14:41:24 INFO None 5286457: status FINISHED
2026-07-27 14:41:24 INFO None 5286458: status FINISHED
2026-07-27 14:41:24 INFO None 5286459: status FINISHED
2026-07-27 14:41:24 INFO None 5286460: status FINISHED
2026-07-27 14:41:24 INFO None 5286461: status FINISHED
2026-07-27 14:41:24 INFO None 5286462: status FINISHED
2026-07-27 14:41:24 INFO None 5286464: status FINISHED
2026-07-27 14:41:24 INFO None 5286465: status FINISHED
2026-07-27 14:41:24 INFO None 5286466: status FINISHED
2026-07-27 14:41:24 INFO None 5286467: status FINISHED
2026-07-27 14:41:26 INFO None 5286469: status FINISHED
2026-07-27 14:41:26 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:41:26 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:41:26 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:41:26 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:41:41 INFO None 5286456: status FINISHED
2026-07-27 14:41:41 INFO None 5286457: status FINISHED
2026-07-27 14:41:41 INFO None 5286458: status FINISHED
2026-07-27 14:41:41 INFO None 5286459: status FINISHED
2026-07-27 14:41:41 INFO None 5286460: status FINISHED
2026-07-27 14:41:41 INFO None 5286461: status FINISHED
2026-07-27 14:41:41 INFO None 5286462: status FINISHED
2026-07-27 14:41:41 INFO None 5286464: status FINISHED
2026-07-27 14:41:41 INFO None 5286465: status FINISHED
2026-07-27 14:41:42 INFO None 5286466: status FINISHED
2026-07-27 14:41:42 INFO None 5286467: status FINISHED
2026-07-27 14:41:42 INFO None 5286469: status FINISHED
2026-07-27 14:41:42 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:41:42 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:41:42 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:41:42 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:41:57 INFO None 5286456: status FINISHED
2026-07-27 14:41:57 INFO None 5286457: status FINISHED
2026-07-27 14:41:57 INFO None 5286458: status FINISHED
2026-07-27 14:41:57 INFO None 5286459: status FINISHED
2026-07-27 14:41:57 INFO None 5286460: status FINISHED
2026-07-27 14:41:57 INFO None 5286461: status FINISHED
2026-07-27 14:41:57 INFO None 5286462: status FINISHED
2026-07-27 14:41:57 INFO None 5286464: status FINISHED
2026-07-27 14:41:57 INFO None 5286465: status FINISHED
2026-07-27 14:41:57 INFO None 5286466: status FINISHED
2026-07-27 14:41:57 INFO None 5286467: status FINISHED
2026-07-27 14:41:57 INFO None 5286469: status FINISHED
2026-07-27 14:41:57 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:41:57 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:41:57 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:41:57 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:42:12 INFO None 5286456: status FINISHED
2026-07-27 14:42:12 INFO None 5286457: status FINISHED
2026-07-27 14:42:12 INFO None 5286458: status FINISHED
2026-07-27 14:42:12 INFO None 5286459: status FINISHED
2026-07-27 14:42:12 INFO None 5286460: status FINISHED
2026-07-27 14:42:12 INFO None 5286461: status FINISHED
2026-07-27 14:42:12 INFO None 5286462: status FINISHED
2026-07-27 14:42:12 INFO None 5286464: status FINISHED
2026-07-27 14:42:12 INFO None 5286465: status FINISHED
2026-07-27 14:42:12 INFO None 5286466: status FINISHED
2026-07-27 14:42:12 INFO None 5286467: status FINISHED
2026-07-27 14:42:12 INFO None 5286469: status FINISHED
2026-07-27 14:42:12 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:42:12 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:42:12 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:42:12 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:42:28 INFO None 5286456: status FINISHED
2026-07-27 14:42:28 INFO None 5286457: status FINISHED
2026-07-27 14:42:28 INFO None 5286458: status FINISHED
2026-07-27 14:42:28 INFO None 5286459: status FINISHED
2026-07-27 14:42:28 INFO None 5286460: status FINISHED
2026-07-27 14:42:28 INFO None 5286461: status FINISHED
2026-07-27 14:42:28 INFO None 5286462: status FINISHED
2026-07-27 14:42:28 INFO None 5286464: status FINISHED
2026-07-27 14:42:28 INFO None 5286465: status FINISHED
2026-07-27 14:42:28 INFO None 5286466: status FINISHED
2026-07-27 14:42:28 INFO None 5286467: status FINISHED
2026-07-27 14:42:28 INFO None 5286469: status FINISHED
2026-07-27 14:42:28 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:42:28 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:42:28 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:42:28 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:42:43 INFO None 5286456: status FINISHED
2026-07-27 14:42:43 INFO None 5286457: status FINISHED
2026-07-27 14:42:44 INFO None 5286458: status FINISHED
2026-07-27 14:42:44 INFO None 5286459: status FINISHED
2026-07-27 14:42:44 INFO None 5286460: status FINISHED
2026-07-27 14:42:44 INFO None 5286461: status FINISHED
2026-07-27 14:42:44 INFO None 5286462: status FINISHED
2026-07-27 14:42:44 INFO None 5286464: status FINISHED
2026-07-27 14:42:44 INFO None 5286465: status FINISHED
2026-07-27 14:42:44 INFO None 5286466: status FINISHED
2026-07-27 14:42:44 INFO None 5286467: status FINISHED
2026-07-27 14:42:44 INFO None 5286469: status FINISHED
2026-07-27 14:42:44 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:42:44 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:42:44 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:42:44 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:42:59 INFO None 5286456: status FINISHED
2026-07-27 14:43:01 INFO None 5286457: status FINISHED
2026-07-27 14:43:01 INFO None 5286458: status FINISHED
2026-07-27 14:43:01 INFO None 5286459: status FINISHED
2026-07-27 14:43:01 INFO None 5286460: status FINISHED
2026-07-27 14:43:01 INFO None 5286461: status FINISHED
2026-07-27 14:43:01 INFO None 5286462: status FINISHED
2026-07-27 14:43:01 INFO None 5286464: status FINISHED
2026-07-27 14:43:01 INFO None 5286465: status FINISHED
2026-07-27 14:43:01 INFO None 5286466: status FINISHED
2026-07-27 14:43:01 INFO None 5286467: status FINISHED
2026-07-27 14:43:01 INFO None 5286469: status FINISHED
2026-07-27 14:43:01 INFO None 5286470: status RUNNING/PENDING
2026-07-27 14:43:01 INFO None 5286471: status RUNNING/PENDING
2026-07-27 14:43:01 INFO None 5286472: status RUNNING/PENDING
2026-07-27 14:43:01 INFO Jobs still running: ['5286470', '5286471', '5286472']. Waiting...
2026-07-27 14:43:16 INFO None 5286456: status FINISHED
2026-07-27 14:43:16 INFO None 5286457: status FINISHED
2026-07-27 14:43:16 INFO None 5286458: status FINISHED
2026-07-27 14:43:16 INFO None 5286459: status FINISHED
2026-07-27 14:43:16 INFO None 5286460: status FINISHED
2026-07-27 14:43:16 INFO None 5286461: status FINISHED
2026-07-27 14:43:16 INFO None 5286462: status FINISHED
2026-07-27 14:43:16 INFO None 5286464: status FINISHED
2026-07-27 14:43:16 INFO None 5286465: status FINISHED
2026-07-27 14:43:16 INFO None 5286466: status FINISHED
2026-07-27 14:43:16 INFO None 5286467: status FINISHED
2026-07-27 14:43:16 INFO None 5286469: status FINISHED
2026-07-27 14:43:16 INFO None 5286470: status FINISHED
2026-07-27 14:43:16 INFO None 5286471: status FINISHED
2026-07-27 14:43:16 INFO None 5286472: status FINISHED
2026-07-27 14:43:16 INFO Jobs ['5286456', '5286457', '5286458', '5286459', '5286460', '5286461', '5286462', '5286464', '5286465', '5286466', '5286467', '5286469', '5286470', '5286471', '5286472'] have finished
2026-07-27 14:43:16 INFO Checking restart files were created ...
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc(3673513755 bytes)
2026-07-27 14:43:16 INFO  Run_model() completed successfully.
2026-07-27 14:43:16 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 14:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:43:16 INFO [TIME] gregorian_conversion simulated_time=2020-02-07 00:00:00 days=153073 seconds=0
2026-07-27 14:43:16 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-27 14:43:16 INFO [TIME] increment current_time 2020-02-06 14:00:00 -> 2020-02-07 00:00:00
2026-07-27 14:43:16 INFO [TIME] after_increment_before_assimilation current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:43:16 INFO ---------->>> Running process_satellite_data()
2026-07-27 14:43:16 INFO [DART] No satellite data found, skipping assimilation
2026-07-27 14:43:16 INFO after_assimilation() skipped
2026-07-27 14:43:16 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-27 14:43:16 INFO [TIME] step_end current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:43:16 INFO [TIME] step_start current_time=2020-02-07 00:00:00 simulated_time=2020-02-07 00:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-27 14:43:16 INFO [TIME] window start=2020-02-07 00:00:00 end=2020-02-07 01:00:00 run_hours=1 has_assimilation=False
2026-07-27 14:43:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:18 INFO Hourly dataset computed and listing created
2026-07-27 14:43:22 INFO Hourly dataset computed
2026-07-27 14:43:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:23 INFO Hourly dataset computed and listing created
2026-07-27 14:43:24 INFO Hourly dataset computed
2026-07-27 14:43:24 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:25 INFO Hourly dataset computed and listing created
2026-07-27 14:43:26 INFO Hourly dataset computed
2026-07-27 14:43:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:27 INFO Hourly dataset computed and listing created
2026-07-27 14:43:28 INFO Hourly dataset computed
2026-07-27 14:43:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:28 INFO Hourly dataset computed and listing created
2026-07-27 14:43:29 INFO Hourly dataset computed
2026-07-27 14:43:29 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:30 INFO Hourly dataset computed and listing created
2026-07-27 14:43:30 INFO Hourly dataset computed
2026-07-27 14:43:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:31 INFO Hourly dataset computed and listing created
2026-07-27 14:43:32 INFO Hourly dataset computed
2026-07-27 14:43:32 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:33 INFO Hourly dataset computed and listing created
2026-07-27 14:43:33 INFO Hourly dataset computed
2026-07-27 14:43:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:34 INFO Hourly dataset computed and listing created
2026-07-27 14:43:35 INFO Hourly dataset computed
2026-07-27 14:43:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:36 INFO Hourly dataset computed and listing created
2026-07-27 14:43:36 INFO Hourly dataset computed
2026-07-27 14:43:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:37 INFO Hourly dataset computed and listing created
2026-07-27 14:43:38 INFO Hourly dataset computed
2026-07-27 14:43:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:39 INFO Hourly dataset computed and listing created
2026-07-27 14:43:39 INFO Hourly dataset computed
2026-07-27 14:43:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:40 INFO Hourly dataset computed and listing created
2026-07-27 14:43:41 INFO Hourly dataset computed
2026-07-27 14:43:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:42 INFO Hourly dataset computed and listing created
2026-07-27 14:43:42 INFO Hourly dataset computed
2026-07-27 14:43:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-27 14:43:43 INFO Hourly dataset computed and listing created
2026-07-27 14:43:44 INFO Hourly dataset computed
2026-07-27 14:43:44 INFO ---------->>> Running CHIMERE model from 2020-02-07 00:00:00 to 2020-02-07 01:00:00
2026-07-27 14:43:44 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:43:44 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1
2026-07-27 14:43:44 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS1/end.2020020614_10_ENS1.nc
2026-07-27 14:43:44 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-27 14:43:44 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:43:44 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-27 14:43:44 INFO Queuing job for member 1...
2026-07-27 14:43:44 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:43:44 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-27 14:43:44 INFO Found: ['5286604']
2026-07-27 14:43:49 INFO [TGCC-IRENE] Submitted job with ID:['5286604']
2026-07-27 14:43:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:43:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2
2026-07-27 14:43:49 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS2/end.2020020614_10_ENS2.nc
2026-07-27 14:43:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-27 14:43:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:43:49 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-27 14:43:49 INFO Queuing job for member 2...
2026-07-27 14:43:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:43:49 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-27 14:43:50 INFO Found: ['5286605']
2026-07-27 14:43:55 INFO [TGCC-IRENE] Submitted job with ID:['5286605']
2026-07-27 14:43:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:43:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3
2026-07-27 14:43:55 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS3/end.2020020614_10_ENS3.nc
2026-07-27 14:43:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-27 14:43:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:43:55 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-27 14:43:55 INFO Queuing job for member 3...
2026-07-27 14:43:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:43:55 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-27 14:43:56 INFO Found: ['5286612']
2026-07-27 14:44:01 INFO [TGCC-IRENE] Submitted job with ID:['5286612']
2026-07-27 14:44:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4
2026-07-27 14:44:01 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS4/end.2020020614_10_ENS4.nc
2026-07-27 14:44:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-27 14:44:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:01 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-27 14:44:01 INFO Queuing job for member 4...
2026-07-27 14:44:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:01 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-27 14:44:02 INFO Found: ['5286614']
2026-07-27 14:44:07 INFO [TGCC-IRENE] Submitted job with ID:['5286614']
2026-07-27 14:44:07 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:07 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5
2026-07-27 14:44:07 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS5/end.2020020614_10_ENS5.nc
2026-07-27 14:44:07 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-27 14:44:07 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:07 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-27 14:44:07 INFO Queuing job for member 5...
2026-07-27 14:44:07 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:07 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-27 14:44:07 INFO Found: ['5286615']
2026-07-27 14:44:12 INFO [TGCC-IRENE] Submitted job with ID:['5286615']
2026-07-27 14:44:12 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:12 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6
2026-07-27 14:44:12 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS6/end.2020020614_10_ENS6.nc
2026-07-27 14:44:12 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-27 14:44:12 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:12 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-27 14:44:12 INFO Queuing job for member 6...
2026-07-27 14:44:12 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:12 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-27 14:44:13 INFO Found: ['5286616']
2026-07-27 14:44:18 INFO [TGCC-IRENE] Submitted job with ID:['5286616']
2026-07-27 14:44:18 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:18 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7
2026-07-27 14:44:18 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS7/end.2020020614_10_ENS7.nc
2026-07-27 14:44:18 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-27 14:44:18 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:18 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-27 14:44:18 INFO Queuing job for member 7...
2026-07-27 14:44:18 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:18 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-27 14:44:19 INFO Found: ['5286617']
2026-07-27 14:44:24 INFO [TGCC-IRENE] Submitted job with ID:['5286617']
2026-07-27 14:44:24 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:24 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8
2026-07-27 14:44:24 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS8/end.2020020614_10_ENS8.nc
2026-07-27 14:44:24 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-27 14:44:24 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:24 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-27 14:44:24 INFO Queuing job for member 8...
2026-07-27 14:44:24 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:24 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-27 14:44:26 INFO Found: ['5286618']
2026-07-27 14:44:31 INFO [TGCC-IRENE] Submitted job with ID:['5286618']
2026-07-27 14:44:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9
2026-07-27 14:44:31 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS9/end.2020020614_10_ENS9.nc
2026-07-27 14:44:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-27 14:44:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:31 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-27 14:44:31 INFO Queuing job for member 9...
2026-07-27 14:44:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:31 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-27 14:44:34 INFO Found: ['5286621']
2026-07-27 14:44:39 INFO [TGCC-IRENE] Submitted job with ID:['5286621']
2026-07-27 14:44:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10
2026-07-27 14:44:39 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS10/end.2020020614_10_ENS10.nc
2026-07-27 14:44:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-27 14:44:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:39 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-27 14:44:39 INFO Queuing job for member 10...
2026-07-27 14:44:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:39 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-27 14:44:41 INFO Found: ['5286622']
2026-07-27 14:44:46 INFO [TGCC-IRENE] Submitted job with ID:['5286622']
2026-07-27 14:44:46 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:46 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11
2026-07-27 14:44:46 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS11/end.2020020614_10_ENS11.nc
2026-07-27 14:44:46 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-27 14:44:46 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:46 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-27 14:44:46 INFO Queuing job for member 11...
2026-07-27 14:44:46 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:46 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-27 14:44:49 INFO Found: ['5286624']
2026-07-27 14:44:54 INFO [TGCC-IRENE] Submitted job with ID:['5286624']
2026-07-27 14:44:54 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:44:54 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12
2026-07-27 14:44:54 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS12/end.2020020614_10_ENS12.nc
2026-07-27 14:44:54 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-27 14:44:54 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:44:54 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-27 14:44:54 INFO Queuing job for member 12...
2026-07-27 14:44:54 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:44:54 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-27 14:44:56 INFO Found: ['5286625']
2026-07-27 14:45:01 INFO [TGCC-IRENE] Submitted job with ID:['5286625']
2026-07-27 14:45:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:45:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13
2026-07-27 14:45:01 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS13/end.2020020614_10_ENS13.nc
2026-07-27 14:45:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-27 14:45:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:45:01 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-27 14:45:01 INFO Queuing job for member 13...
2026-07-27 14:45:01 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:45:01 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-27 14:45:04 INFO Found: ['5286632']
2026-07-27 14:45:09 INFO [TGCC-IRENE] Submitted job with ID:['5286632']
2026-07-27 14:45:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:45:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14
2026-07-27 14:45:09 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS14/end.2020020614_10_ENS14.nc
2026-07-27 14:45:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-27 14:45:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:45:09 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-27 14:45:09 INFO Queuing job for member 14...
2026-07-27 14:45:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:45:09 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-27 14:45:11 INFO Found: ['5286633']
2026-07-27 14:45:16 INFO [TGCC-IRENE] Submitted job with ID:['5286633']
2026-07-27 14:45:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-27 14:45:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15
2026-07-27 14:45:16 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICinfl2_0607_15m_low_v2/ENS15/end.2020020614_10_ENS15.nc
2026-07-27 14:45:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-27 14:45:17 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-27 14:45:17 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-27 14:45:17 INFO Queuing job for member 15...
2026-07-27 14:45:17 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-27 14:45:17 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-27 14:45:19 INFO Found: ['5286634']
2026-07-27 14:45:24 INFO [TGCC-IRENE] Submitted job with ID:['5286634']
2026-07-27 14:45:24 INFO Checking job status ...
2026-07-27 14:45:24 INFO None 5286604: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286605: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286612: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286614: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286615: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286616: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286617: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286618: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286621: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286622: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286624: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286625: status RUNNING/PENDING
2026-07-27 14:45:24 INFO None 5286632: status RUNNING/PENDING
2026-07-27 14:45:27 INFO None 5286633: status RUNNING/PENDING
2026-07-27 14:45:27 INFO None 5286634: status RUNNING/PENDING
2026-07-27 14:45:27 INFO Jobs still running: ['5286604', '5286605', '5286612', '5286614', '5286615', '5286616', '5286617', '5286618', '5286621', '5286622', '5286624', '5286625', '5286632', '5286633', '5286634']. Waiting...
[2026-07-27T14:45:27.222] error: *** JOB 5285432 ON irene4341 CANCELLED AT 2026-07-27T14:45:27 DUE to SIGNAL Terminated ***
