+ SCRIPT_PID=326695
+ /bin/bash -x /tmp/tmp.A8ODhjyksF
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
+ python -u main.py -c config/config_irene_IM.yaml
2026-07-13 20:19:13 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-07-13 20:19:13 INFO [PIPELINE] =======================================
2026-07-13 20:19:13 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-07-13 20:19:13 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-07-13 20:19:13 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-07-13 20:19:13 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260713_201913.log
2026-07-13 20:19:13 INFO [PIPELINE] =======================================
2026-07-13 20:19:13 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-07-13 20:19:13 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-07-13 20:19:13 INFO [STEP] ---- TIME LOOP START ----
2026-07-13 20:19:13 INFO [TIME] step_start current_time=2020-02-06 00:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:19:13 INFO [TIME] window start=2020-02-06 00:00:00 end=2020-02-06 01:00:00 run_hours=1 has_assimilation=False
2026-07-13 20:19:13 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-07-13 20:19:13 INFO Copying EMIS ...
2026-07-13 20:19:14 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens01.nc
2026-07-13 20:19:14 INFO Linking first END ...
2026-07-13 20:19:14 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:14 INFO >> Checking links...
2026-07-13 20:19:14 INFO >> All links are good for ENS1  ...
2026-07-13 20:19:14 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:24 INFO Hourly dataset computed and listing created
2026-07-13 20:19:28 INFO Hourly dataset computed
2026-07-13 20:19:28 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-07-13 20:19:28 INFO Copying EMIS ...
2026-07-13 20:19:28 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens02.nc
2026-07-13 20:19:28 INFO Linking first END ...
2026-07-13 20:19:28 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:28 INFO >> Checking links...
2026-07-13 20:19:28 INFO >> All links are good for ENS2  ...
2026-07-13 20:19:28 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:29 INFO Hourly dataset computed and listing created
2026-07-13 20:19:30 INFO Hourly dataset computed
2026-07-13 20:19:30 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-07-13 20:19:30 INFO Copying EMIS ...
2026-07-13 20:19:31 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens03.nc
2026-07-13 20:19:31 INFO Linking first END ...
2026-07-13 20:19:31 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:31 INFO >> Checking links...
2026-07-13 20:19:31 INFO >> All links are good for ENS3  ...
2026-07-13 20:19:31 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:31 INFO Hourly dataset computed and listing created
2026-07-13 20:19:32 INFO Hourly dataset computed
2026-07-13 20:19:32 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-07-13 20:19:32 INFO Copying EMIS ...
2026-07-13 20:19:33 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens04.nc
2026-07-13 20:19:33 INFO Linking first END ...
2026-07-13 20:19:33 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:33 INFO >> Checking links...
2026-07-13 20:19:33 INFO >> All links are good for ENS4  ...
2026-07-13 20:19:33 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:34 INFO Hourly dataset computed and listing created
2026-07-13 20:19:34 INFO Hourly dataset computed
2026-07-13 20:19:34 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-07-13 20:19:34 INFO Copying EMIS ...
2026-07-13 20:19:35 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens05.nc
2026-07-13 20:19:35 INFO Linking first END ...
2026-07-13 20:19:35 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:35 INFO >> Checking links...
2026-07-13 20:19:35 INFO >> All links are good for ENS5  ...
2026-07-13 20:19:35 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:36 INFO Hourly dataset computed and listing created
2026-07-13 20:19:36 INFO Hourly dataset computed
2026-07-13 20:19:36 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-07-13 20:19:36 INFO Copying EMIS ...
2026-07-13 20:19:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens06.nc
2026-07-13 20:19:37 INFO Linking first END ...
2026-07-13 20:19:37 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:37 INFO >> Checking links...
2026-07-13 20:19:37 INFO >> All links are good for ENS6  ...
2026-07-13 20:19:37 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:38 INFO Hourly dataset computed and listing created
2026-07-13 20:19:38 INFO Hourly dataset computed
2026-07-13 20:19:38 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-07-13 20:19:38 INFO Copying EMIS ...
2026-07-13 20:19:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens07.nc
2026-07-13 20:19:39 INFO Linking first END ...
2026-07-13 20:19:39 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:39 INFO >> Checking links...
2026-07-13 20:19:39 INFO >> All links are good for ENS7  ...
2026-07-13 20:19:39 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:40 INFO Hourly dataset computed and listing created
2026-07-13 20:19:41 INFO Hourly dataset computed
2026-07-13 20:19:41 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-07-13 20:19:41 INFO Copying EMIS ...
2026-07-13 20:19:41 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens08.nc
2026-07-13 20:19:41 INFO Linking first END ...
2026-07-13 20:19:41 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:41 INFO >> Checking links...
2026-07-13 20:19:41 INFO >> All links are good for ENS8  ...
2026-07-13 20:19:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:42 INFO Hourly dataset computed and listing created
2026-07-13 20:19:46 INFO Hourly dataset computed
2026-07-13 20:19:46 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-07-13 20:19:46 INFO Copying EMIS ...
2026-07-13 20:19:47 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens09.nc
2026-07-13 20:19:47 INFO Linking first END ...
2026-07-13 20:19:47 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:47 INFO >> Checking links...
2026-07-13 20:19:47 INFO >> All links are good for ENS9  ...
2026-07-13 20:19:47 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:48 INFO Hourly dataset computed and listing created
2026-07-13 20:19:48 INFO Hourly dataset computed
2026-07-13 20:19:48 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-07-13 20:19:48 INFO Copying EMIS ...
2026-07-13 20:19:49 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens10.nc
2026-07-13 20:19:49 INFO Linking first END ...
2026-07-13 20:19:49 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:49 INFO >> Checking links...
2026-07-13 20:19:49 INFO >> All links are good for ENS10  ...
2026-07-13 20:19:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:50 INFO Hourly dataset computed and listing created
2026-07-13 20:19:50 INFO Hourly dataset computed
2026-07-13 20:19:51 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-07-13 20:19:51 INFO Copying EMIS ...
2026-07-13 20:19:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens11.nc
2026-07-13 20:19:51 INFO Linking first END ...
2026-07-13 20:19:51 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:51 INFO >> Checking links...
2026-07-13 20:19:51 INFO >> All links are good for ENS11  ...
2026-07-13 20:19:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:52 INFO Hourly dataset computed and listing created
2026-07-13 20:19:53 INFO Hourly dataset computed
2026-07-13 20:19:53 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-07-13 20:19:53 INFO Copying EMIS ...
2026-07-13 20:19:53 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens12.nc
2026-07-13 20:19:53 INFO Linking first END ...
2026-07-13 20:19:53 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:53 INFO >> Checking links...
2026-07-13 20:19:53 INFO >> All links are good for ENS12  ...
2026-07-13 20:19:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:54 INFO Hourly dataset computed and listing created
2026-07-13 20:19:55 INFO Hourly dataset computed
2026-07-13 20:19:55 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-07-13 20:19:55 INFO Copying EMIS ...
2026-07-13 20:19:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens13.nc
2026-07-13 20:19:55 INFO Linking first END ...
2026-07-13 20:19:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:55 INFO >> Checking links...
2026-07-13 20:19:55 INFO >> All links are good for ENS13  ...
2026-07-13 20:19:55 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:56 INFO Hourly dataset computed and listing created
2026-07-13 20:19:57 INFO Hourly dataset computed
2026-07-13 20:19:57 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-07-13 20:19:57 INFO Copying EMIS ...
2026-07-13 20:19:57 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens14.nc
2026-07-13 20:19:57 INFO Linking first END ...
2026-07-13 20:19:57 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:19:57 INFO >> Checking links...
2026-07-13 20:19:57 INFO >> All links are good for ENS14  ...
2026-07-13 20:19:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:19:58 INFO Hourly dataset computed and listing created
2026-07-13 20:19:59 INFO Hourly dataset computed
2026-07-13 20:19:59 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-07-13 20:19:59 INFO Copying EMIS ...
2026-07-13 20:20:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Thursday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Thursday.s.ens15.nc
2026-07-13 20:20:00 INFO Linking first END ...
2026-07-13 20:20:00 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020020500_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-07-13 20:20:00 INFO >> Checking links...
2026-07-13 20:20:00 INFO >> All links are good for ENS15  ...
2026-07-13 20:20:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:20:01 INFO Hourly dataset computed and listing created
2026-07-13 20:20:01 INFO Hourly dataset computed
2026-07-13 20:20:01 INFO ---------->>> Running CHIMERE model from 2020-02-06 00:00:00 to 2020-02-06 01:00:00
2026-07-13 20:20:01 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:01 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-13 20:20:01 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020500_24_ENS1.nc
2026-07-13 20:20:01 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-13 20:20:01 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:01 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-13 20:20:02 INFO Queuing job for member 1...
2026-07-13 20:20:02 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:02 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-13 20:20:03 INFO Found: ['5152193']
2026-07-13 20:20:08 INFO [TGCC-IRENE] Submitted job with ID:['5152193']
2026-07-13 20:20:08 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:08 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-13 20:20:08 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020500_24_ENS2.nc
2026-07-13 20:20:08 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-13 20:20:08 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:08 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-13 20:20:08 INFO Queuing job for member 2...
2026-07-13 20:20:08 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:08 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-13 20:20:09 INFO Found: ['5152195']
2026-07-13 20:20:14 INFO [TGCC-IRENE] Submitted job with ID:['5152195']
2026-07-13 20:20:14 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:14 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-13 20:20:14 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020500_24_ENS3.nc
2026-07-13 20:20:14 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-13 20:20:14 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:14 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-13 20:20:14 INFO Queuing job for member 3...
2026-07-13 20:20:14 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:14 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-13 20:20:14 INFO Found: ['5152196']
2026-07-13 20:20:19 INFO [TGCC-IRENE] Submitted job with ID:['5152196']
2026-07-13 20:20:19 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:19 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-13 20:20:19 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020500_24_ENS4.nc
2026-07-13 20:20:19 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-13 20:20:19 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:19 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-13 20:20:19 INFO Queuing job for member 4...
2026-07-13 20:20:19 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:19 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-13 20:20:20 INFO Found: ['5152197']
2026-07-13 20:20:25 INFO [TGCC-IRENE] Submitted job with ID:['5152197']
2026-07-13 20:20:25 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:25 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-13 20:20:25 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020500_24_ENS5.nc
2026-07-13 20:20:25 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-13 20:20:25 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:25 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-13 20:20:25 INFO Queuing job for member 5...
2026-07-13 20:20:25 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:25 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-13 20:20:26 INFO Found: ['5152198']
2026-07-13 20:20:31 INFO [TGCC-IRENE] Submitted job with ID:['5152198']
2026-07-13 20:20:31 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:31 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-13 20:20:31 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020500_24_ENS6.nc
2026-07-13 20:20:31 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-13 20:20:31 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:31 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-13 20:20:31 INFO Queuing job for member 6...
2026-07-13 20:20:31 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:31 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-13 20:20:32 INFO Found: ['5152199']
2026-07-13 20:20:37 INFO [TGCC-IRENE] Submitted job with ID:['5152199']
2026-07-13 20:20:37 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:37 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-13 20:20:37 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020500_24_ENS7.nc
2026-07-13 20:20:37 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-13 20:20:37 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:37 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-13 20:20:37 INFO Queuing job for member 7...
2026-07-13 20:20:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:37 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-13 20:20:38 INFO Found: ['5152200']
2026-07-13 20:20:43 INFO [TGCC-IRENE] Submitted job with ID:['5152200']
2026-07-13 20:20:43 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:43 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-13 20:20:43 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020500_24_ENS8.nc
2026-07-13 20:20:43 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-13 20:20:43 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:43 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-13 20:20:43 INFO Queuing job for member 8...
2026-07-13 20:20:43 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:43 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-13 20:20:44 INFO Found: ['5152201']
2026-07-13 20:20:49 INFO [TGCC-IRENE] Submitted job with ID:['5152201']
2026-07-13 20:20:49 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:49 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-13 20:20:49 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020500_24_ENS9.nc
2026-07-13 20:20:49 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-13 20:20:49 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:49 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-13 20:20:49 INFO Queuing job for member 9...
2026-07-13 20:20:49 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:49 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-13 20:20:50 INFO Found: ['5152203']
2026-07-13 20:20:55 INFO [TGCC-IRENE] Submitted job with ID:['5152203']
2026-07-13 20:20:55 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:20:55 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-13 20:20:55 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020500_24_ENS10.nc
2026-07-13 20:20:55 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-13 20:20:55 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:20:55 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-13 20:20:55 INFO Queuing job for member 10...
2026-07-13 20:20:55 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:20:55 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-13 20:20:56 INFO Found: ['5152204']
2026-07-13 20:21:01 INFO [TGCC-IRENE] Submitted job with ID:['5152204']
2026-07-13 20:23:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:23:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-13 20:23:09 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020500_24_ENS11.nc
2026-07-13 20:23:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-13 20:23:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:23:09 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-13 20:23:09 INFO Queuing job for member 11...
2026-07-13 20:23:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:23:09 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-13 20:23:10 INFO Found: ['5152212']
2026-07-13 20:23:15 INFO [TGCC-IRENE] Submitted job with ID:['5152212']
2026-07-13 20:23:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:23:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-13 20:23:15 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020500_24_ENS12.nc
2026-07-13 20:23:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-13 20:23:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:23:15 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-13 20:23:15 INFO Queuing job for member 12...
2026-07-13 20:23:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:23:15 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-13 20:23:16 INFO Found: ['5152213']
2026-07-13 20:23:21 INFO [TGCC-IRENE] Submitted job with ID:['5152213']
2026-07-13 20:23:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:23:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-13 20:23:21 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020500_24_ENS13.nc
2026-07-13 20:23:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-13 20:23:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:23:21 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-13 20:23:21 INFO Queuing job for member 13...
2026-07-13 20:23:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:23:21 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-13 20:23:21 INFO Found: ['5152214']
2026-07-13 20:23:26 INFO [TGCC-IRENE] Submitted job with ID:['5152214']
2026-07-13 20:23:26 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:23:26 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-13 20:23:26 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020500_24_ENS14.nc
2026-07-13 20:23:26 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-13 20:23:26 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:23:26 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-13 20:23:26 INFO Queuing job for member 14...
2026-07-13 20:23:26 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:23:26 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-13 20:23:27 INFO Found: ['5152215']
2026-07-13 20:23:32 INFO [TGCC-IRENE] Submitted job with ID:['5152215']
2026-07-13 20:23:32 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:23:32 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-13 20:23:32 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020500_24_ENS15.nc
2026-07-13 20:23:32 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-13 20:23:32 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:23:32 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-13 20:23:32 INFO Queuing job for member 15...
2026-07-13 20:23:32 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:23:32 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-13 20:23:33 INFO Found: ['5152216']
2026-07-13 20:23:38 INFO [TGCC-IRENE] Submitted job with ID:['5152216']
2026-07-13 20:23:38 INFO Checking job status ...
2026-07-13 20:23:38 INFO None 5152193: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152197: status FINISHED
2026-07-13 20:23:38 INFO None 5152198: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152199: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152200: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152201: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152203: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152204: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:23:38 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:23:39 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:23:39 INFO Jobs still running: ['5152193', '5152195', '5152196', '5152198', '5152199', '5152200', '5152201', '5152203', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:23:54 INFO None 5152193: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152197: status FINISHED
2026-07-13 20:23:54 INFO None 5152198: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152199: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152200: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152201: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152203: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152204: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:23:54 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:23:54 INFO Jobs still running: ['5152193', '5152195', '5152196', '5152198', '5152199', '5152200', '5152201', '5152203', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:24:10 INFO None 5152193: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152197: status FINISHED
2026-07-13 20:24:10 INFO None 5152198: status FINISHED
2026-07-13 20:24:10 INFO None 5152199: status FINISHED
2026-07-13 20:24:10 INFO None 5152200: status FINISHED
2026-07-13 20:24:10 INFO None 5152201: status FINISHED
2026-07-13 20:24:10 INFO None 5152203: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152204: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:24:10 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:24:10 INFO Jobs still running: ['5152193', '5152195', '5152196', '5152203', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:24:25 INFO None 5152193: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152197: status FINISHED
2026-07-13 20:24:25 INFO None 5152198: status FINISHED
2026-07-13 20:24:25 INFO None 5152199: status FINISHED
2026-07-13 20:24:25 INFO None 5152200: status FINISHED
2026-07-13 20:24:25 INFO None 5152201: status FINISHED
2026-07-13 20:24:25 INFO None 5152203: status FINISHED
2026-07-13 20:24:25 INFO None 5152204: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:24:25 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:24:25 INFO Jobs still running: ['5152193', '5152195', '5152196', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:24:40 INFO None 5152193: status RUNNING/PENDING
2026-07-13 20:24:40 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152197: status FINISHED
2026-07-13 20:24:41 INFO None 5152198: status FINISHED
2026-07-13 20:24:41 INFO None 5152199: status FINISHED
2026-07-13 20:24:41 INFO None 5152200: status FINISHED
2026-07-13 20:24:41 INFO None 5152201: status FINISHED
2026-07-13 20:24:41 INFO None 5152203: status FINISHED
2026-07-13 20:24:41 INFO None 5152204: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:24:41 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:24:41 INFO Jobs still running: ['5152193', '5152195', '5152196', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:24:56 INFO None 5152193: status FINISHED
2026-07-13 20:24:56 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152197: status FINISHED
2026-07-13 20:24:56 INFO None 5152198: status FINISHED
2026-07-13 20:24:56 INFO None 5152199: status FINISHED
2026-07-13 20:24:56 INFO None 5152200: status FINISHED
2026-07-13 20:24:56 INFO None 5152201: status FINISHED
2026-07-13 20:24:56 INFO None 5152203: status FINISHED
2026-07-13 20:24:56 INFO None 5152204: status FINISHED
2026-07-13 20:24:56 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152213: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152214: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152215: status RUNNING/PENDING
2026-07-13 20:24:56 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:24:56 INFO Jobs still running: ['5152195', '5152196', '5152212', '5152213', '5152214', '5152215', '5152216']. Waiting...
2026-07-13 20:25:11 INFO None 5152193: status FINISHED
2026-07-13 20:25:11 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:25:11 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:25:11 INFO None 5152197: status FINISHED
2026-07-13 20:25:11 INFO None 5152198: status FINISHED
2026-07-13 20:25:11 INFO None 5152199: status FINISHED
2026-07-13 20:25:11 INFO None 5152200: status FINISHED
2026-07-13 20:25:11 INFO None 5152201: status FINISHED
2026-07-13 20:25:11 INFO None 5152203: status FINISHED
2026-07-13 20:25:11 INFO None 5152204: status FINISHED
2026-07-13 20:25:11 INFO None 5152212: status RUNNING/PENDING
2026-07-13 20:25:12 INFO None 5152213: status FINISHED
2026-07-13 20:25:12 INFO None 5152214: status FINISHED
2026-07-13 20:25:12 INFO None 5152215: status FINISHED
2026-07-13 20:25:12 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:25:12 INFO Jobs still running: ['5152195', '5152196', '5152212', '5152216']. Waiting...
2026-07-13 20:25:27 INFO None 5152193: status FINISHED
2026-07-13 20:25:27 INFO None 5152195: status RUNNING/PENDING
2026-07-13 20:25:27 INFO None 5152196: status RUNNING/PENDING
2026-07-13 20:25:27 INFO None 5152197: status FINISHED
2026-07-13 20:25:27 INFO None 5152198: status FINISHED
2026-07-13 20:25:27 INFO None 5152199: status FINISHED
2026-07-13 20:25:27 INFO None 5152200: status FINISHED
2026-07-13 20:25:27 INFO None 5152201: status FINISHED
2026-07-13 20:25:27 INFO None 5152203: status FINISHED
2026-07-13 20:25:27 INFO None 5152204: status FINISHED
2026-07-13 20:25:27 INFO None 5152212: status FINISHED
2026-07-13 20:25:27 INFO None 5152213: status FINISHED
2026-07-13 20:25:27 INFO None 5152214: status FINISHED
2026-07-13 20:25:27 INFO None 5152215: status FINISHED
2026-07-13 20:25:27 INFO None 5152216: status RUNNING/PENDING
2026-07-13 20:25:27 INFO Jobs still running: ['5152195', '5152196', '5152216']. Waiting...
2026-07-13 20:25:42 INFO None 5152193: status FINISHED
2026-07-13 20:25:42 INFO None 5152195: status FINISHED
2026-07-13 20:25:42 INFO None 5152196: status FINISHED
2026-07-13 20:25:42 INFO None 5152197: status FINISHED
2026-07-13 20:25:42 INFO None 5152198: status FINISHED
2026-07-13 20:25:42 INFO None 5152199: status FINISHED
2026-07-13 20:25:42 INFO None 5152200: status FINISHED
2026-07-13 20:25:42 INFO None 5152201: status FINISHED
2026-07-13 20:25:42 INFO None 5152203: status FINISHED
2026-07-13 20:25:42 INFO None 5152204: status FINISHED
2026-07-13 20:25:42 INFO None 5152212: status FINISHED
2026-07-13 20:25:42 INFO None 5152213: status FINISHED
2026-07-13 20:25:42 INFO None 5152214: status FINISHED
2026-07-13 20:25:42 INFO None 5152215: status FINISHED
2026-07-13 20:25:42 INFO None 5152216: status FINISHED
2026-07-13 20:25:42 INFO Jobs ['5152193', '5152195', '5152196', '5152197', '5152198', '5152199', '5152200', '5152201', '5152203', '5152204', '5152212', '5152213', '5152214', '5152215', '5152216'] have finished
2026-07-13 20:25:42 INFO Checking restart files were created ...
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc(668832435 bytes)
2026-07-13 20:25:42 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc(668832435 bytes)
2026-07-13 20:25:42 INFO  Run_model() completed successfully.
2026-07-13 20:25:42 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 00:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:25:42 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 01:00:00 days=153072 seconds=3600
2026-07-13 20:25:42 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-13 20:25:42 INFO [TIME] increment current_time 2020-02-06 00:00:00 -> 2020-02-06 01:00:00
2026-07-13 20:25:42 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:25:42 INFO ---------->>> Running process_satellite_data()
2026-07-13 20:25:42 INFO [DART] No satellite data found, skipping assimilation
2026-07-13 20:25:42 INFO after_assimilation() skipped
2026-07-13 20:25:42 INFO /////////////////// Cycle is DONE; starting a new loop!
2026-07-13 20:25:42 INFO [TIME] step_end current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:25:42 INFO [TIME] step_start current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 01:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:25:42 INFO [TIME] window start=2020-02-06 01:00:00 end=2020-02-06 09:00:00 run_hours=8 has_assimilation=True
2026-07-13 20:25:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:25:44 INFO Hourly dataset computed and listing created
2026-07-13 20:25:59 INFO Hourly dataset computed
2026-07-13 20:26:00 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:26:01 INFO Hourly dataset computed and listing created
2026-07-13 20:28:10 INFO Hourly dataset computed
2026-07-13 20:28:10 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:11 INFO Hourly dataset computed and listing created
2026-07-13 20:28:13 INFO Hourly dataset computed
2026-07-13 20:28:13 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:14 INFO Hourly dataset computed and listing created
2026-07-13 20:28:16 INFO Hourly dataset computed
2026-07-13 20:28:16 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:17 INFO Hourly dataset computed and listing created
2026-07-13 20:28:19 INFO Hourly dataset computed
2026-07-13 20:28:19 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:20 INFO Hourly dataset computed and listing created
2026-07-13 20:28:22 INFO Hourly dataset computed
2026-07-13 20:28:22 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:24 INFO Hourly dataset computed and listing created
2026-07-13 20:28:26 INFO Hourly dataset computed
2026-07-13 20:28:26 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:28 INFO Hourly dataset computed and listing created
2026-07-13 20:28:30 INFO Hourly dataset computed
2026-07-13 20:28:30 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:31 INFO Hourly dataset computed and listing created
2026-07-13 20:28:34 INFO Hourly dataset computed
2026-07-13 20:28:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:35 INFO Hourly dataset computed and listing created
2026-07-13 20:28:38 INFO Hourly dataset computed
2026-07-13 20:28:38 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:39 INFO Hourly dataset computed and listing created
2026-07-13 20:28:41 INFO Hourly dataset computed
2026-07-13 20:28:41 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:43 INFO Hourly dataset computed and listing created
2026-07-13 20:28:45 INFO Hourly dataset computed
2026-07-13 20:28:45 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:46 INFO Hourly dataset computed and listing created
2026-07-13 20:28:49 INFO Hourly dataset computed
2026-07-13 20:28:49 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:50 INFO Hourly dataset computed and listing created
2026-07-13 20:28:53 INFO Hourly dataset computed
2026-07-13 20:28:53 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-07-13 20:28:54 INFO Hourly dataset computed and listing created
2026-07-13 20:28:57 INFO Hourly dataset computed
2026-07-13 20:28:57 INFO ---------->>> Running CHIMERE model from 2020-02-06 01:00:00 to 2020-02-06 09:00:00
2026-07-13 20:28:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:28:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1
2026-07-13 20:28:57 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020600_1_ENS1.nc
2026-07-13 20:28:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-07-13 20:28:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:28:57 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-07-13 20:28:57 INFO Queuing job for member 1...
2026-07-13 20:28:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:28:57 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-07-13 20:28:58 INFO Found: ['5152244']
2026-07-13 20:29:03 INFO [TGCC-IRENE] Submitted job with ID:['5152244']
2026-07-13 20:29:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2
2026-07-13 20:29:03 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020600_1_ENS2.nc
2026-07-13 20:29:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-07-13 20:29:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:03 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-07-13 20:29:03 INFO Queuing job for member 2...
2026-07-13 20:29:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:03 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-07-13 20:29:04 INFO Found: ['5152246']
2026-07-13 20:29:09 INFO [TGCC-IRENE] Submitted job with ID:['5152246']
2026-07-13 20:29:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3
2026-07-13 20:29:09 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020600_1_ENS3.nc
2026-07-13 20:29:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-07-13 20:29:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:09 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-07-13 20:29:09 INFO Queuing job for member 3...
2026-07-13 20:29:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:09 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-07-13 20:29:10 INFO Found: ['5152247']
2026-07-13 20:29:15 INFO [TGCC-IRENE] Submitted job with ID:['5152247']
2026-07-13 20:29:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4
2026-07-13 20:29:15 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020600_1_ENS4.nc
2026-07-13 20:29:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-07-13 20:29:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:15 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-07-13 20:29:15 INFO Queuing job for member 4...
2026-07-13 20:29:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:15 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-07-13 20:29:16 INFO Found: ['5152248']
2026-07-13 20:29:21 INFO [TGCC-IRENE] Submitted job with ID:['5152248']
2026-07-13 20:29:21 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:21 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5
2026-07-13 20:29:21 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020600_1_ENS5.nc
2026-07-13 20:29:21 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-07-13 20:29:21 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:21 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-07-13 20:29:21 INFO Queuing job for member 5...
2026-07-13 20:29:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:21 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-07-13 20:29:22 INFO Found: ['5152249']
2026-07-13 20:29:27 INFO [TGCC-IRENE] Submitted job with ID:['5152249']
2026-07-13 20:29:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6
2026-07-13 20:29:27 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020600_1_ENS6.nc
2026-07-13 20:29:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-07-13 20:29:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:27 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-07-13 20:29:27 INFO Queuing job for member 6...
2026-07-13 20:29:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:27 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-07-13 20:29:28 INFO Found: ['5152250']
2026-07-13 20:29:33 INFO [TGCC-IRENE] Submitted job with ID:['5152250']
2026-07-13 20:29:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7
2026-07-13 20:29:33 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020600_1_ENS7.nc
2026-07-13 20:29:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-07-13 20:29:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:33 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-07-13 20:29:33 INFO Queuing job for member 7...
2026-07-13 20:29:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:33 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-07-13 20:29:34 INFO Found: ['5152251']
2026-07-13 20:29:39 INFO [TGCC-IRENE] Submitted job with ID:['5152251']
2026-07-13 20:29:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8
2026-07-13 20:29:39 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020600_1_ENS8.nc
2026-07-13 20:29:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-07-13 20:29:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:39 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-07-13 20:29:39 INFO Queuing job for member 8...
2026-07-13 20:29:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:39 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-07-13 20:29:40 INFO Found: ['5152252']
2026-07-13 20:29:45 INFO [TGCC-IRENE] Submitted job with ID:['5152252']
2026-07-13 20:29:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9
2026-07-13 20:29:45 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020600_1_ENS9.nc
2026-07-13 20:29:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-07-13 20:29:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:45 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-07-13 20:29:45 INFO Queuing job for member 9...
2026-07-13 20:29:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:45 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-07-13 20:29:46 INFO Found: ['5152253']
2026-07-13 20:29:51 INFO [TGCC-IRENE] Submitted job with ID:['5152253']
2026-07-13 20:29:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10
2026-07-13 20:29:51 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020600_1_ENS10.nc
2026-07-13 20:29:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-07-13 20:29:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:51 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-07-13 20:29:51 INFO Queuing job for member 10...
2026-07-13 20:29:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:51 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-07-13 20:29:52 INFO Found: ['5152254']
2026-07-13 20:29:57 INFO [TGCC-IRENE] Submitted job with ID:['5152254']
2026-07-13 20:29:57 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:29:57 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11
2026-07-13 20:29:57 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020600_1_ENS11.nc
2026-07-13 20:29:57 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-07-13 20:29:57 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:29:57 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-07-13 20:29:57 INFO Queuing job for member 11...
2026-07-13 20:29:57 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:29:57 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-07-13 20:29:58 INFO Found: ['5152255']
2026-07-13 20:30:03 INFO [TGCC-IRENE] Submitted job with ID:['5152255']
2026-07-13 20:30:03 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:30:03 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12
2026-07-13 20:30:03 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020600_1_ENS12.nc
2026-07-13 20:30:03 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-07-13 20:30:03 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:30:03 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-07-13 20:30:03 INFO Queuing job for member 12...
2026-07-13 20:30:03 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:30:03 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-07-13 20:30:04 INFO Found: ['5152257']
2026-07-13 20:30:09 INFO [TGCC-IRENE] Submitted job with ID:['5152257']
2026-07-13 20:30:09 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:30:09 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13
2026-07-13 20:30:09 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020600_1_ENS13.nc
2026-07-13 20:30:09 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-07-13 20:30:09 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:30:09 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-07-13 20:30:09 INFO Queuing job for member 13...
2026-07-13 20:30:09 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:30:09 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-07-13 20:30:10 INFO Found: ['5152258']
2026-07-13 20:30:15 INFO [TGCC-IRENE] Submitted job with ID:['5152258']
2026-07-13 20:30:15 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:30:15 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14
2026-07-13 20:30:15 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020600_1_ENS14.nc
2026-07-13 20:30:15 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-07-13 20:30:15 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:30:15 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-07-13 20:30:15 INFO Queuing job for member 14...
2026-07-13 20:30:15 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:30:15 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-07-13 20:30:15 INFO Found: ['5152260']
2026-07-13 20:30:20 INFO [TGCC-IRENE] Submitted job with ID:['5152260']
2026-07-13 20:30:20 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-07-13 20:30:20 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15
2026-07-13 20:30:20 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020600_1_ENS15.nc
2026-07-13 20:30:20 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-07-13 20:30:20 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-07-13 20:30:20 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-07-13 20:30:21 INFO Queuing job for member 15...
2026-07-13 20:30:21 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-07-13 20:30:21 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-07-13 20:30:21 INFO Found: ['5152261']
2026-07-13 20:30:26 INFO [TGCC-IRENE] Submitted job with ID:['5152261']
2026-07-13 20:30:26 INFO Checking job status ...
2026-07-13 20:30:26 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:30:26 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:30:26 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:30:26 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:30:27 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:30:27 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:30:42 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:30:42 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:30:42 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:30:57 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:30:57 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:30:58 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:30:58 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:31:13 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:31:13 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:31:13 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:31:28 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:33:22 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:33:30 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:33:30 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:33:45 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:33:45 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:33:45 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:34:00 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:34:00 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:34:00 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:34:00 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:34:00 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:34:00 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:34:01 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:34:01 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:34:16 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:34:16 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:34:16 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:34:34 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:34:34 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:34:35 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:34:35 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:34:50 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:34:50 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:34:50 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:35:05 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:35:05 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:35:05 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:35:20 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:35:20 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:35:20 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:35:20 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:35:21 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:35:21 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:35:36 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:35:36 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:35:36 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:35:51 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:35:51 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:35:52 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:35:52 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:35:52 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:36:07 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:36:07 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:36:07 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:36:22 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152252: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:38:20 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:38:20 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:38:35 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:38:35 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:38:35 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:38:35 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152252: status FINISHED
2026-07-13 20:38:36 INFO None 5152253: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:38:36 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:38:36 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:38:51 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152247: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152248: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152252: status FINISHED
2026-07-13 20:38:51 INFO None 5152253: status FINISHED
2026-07-13 20:38:51 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152258: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:38:51 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:38:51 INFO Jobs still running: ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261']. Waiting...
2026-07-13 20:39:06 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152247: status FINISHED
2026-07-13 20:39:06 INFO None 5152248: status FINISHED
2026-07-13 20:39:06 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152252: status FINISHED
2026-07-13 20:39:06 INFO None 5152253: status FINISHED
2026-07-13 20:39:06 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:39:06 INFO None 5152257: status RUNNING/PENDING
2026-07-13 20:39:07 INFO None 5152258: status FINISHED
2026-07-13 20:39:07 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:39:07 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:39:07 INFO Jobs still running: ['5152244', '5152246', '5152249', '5152250', '5152251', '5152254', '5152255', '5152257', '5152260', '5152261']. Waiting...
2026-07-13 20:39:22 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152247: status FINISHED
2026-07-13 20:39:22 INFO None 5152248: status FINISHED
2026-07-13 20:39:22 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152252: status FINISHED
2026-07-13 20:39:22 INFO None 5152253: status FINISHED
2026-07-13 20:39:22 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152257: status FINISHED
2026-07-13 20:39:22 INFO None 5152258: status FINISHED
2026-07-13 20:39:22 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:39:22 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:39:22 INFO Jobs still running: ['5152244', '5152246', '5152249', '5152250', '5152251', '5152254', '5152255', '5152260', '5152261']. Waiting...
2026-07-13 20:39:37 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152247: status FINISHED
2026-07-13 20:39:37 INFO None 5152248: status FINISHED
2026-07-13 20:39:37 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152252: status FINISHED
2026-07-13 20:39:37 INFO None 5152253: status FINISHED
2026-07-13 20:39:37 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152257: status FINISHED
2026-07-13 20:39:37 INFO None 5152258: status FINISHED
2026-07-13 20:39:37 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:39:37 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:39:37 INFO Jobs still running: ['5152244', '5152246', '5152249', '5152250', '5152251', '5152254', '5152255', '5152260', '5152261']. Waiting...
2026-07-13 20:39:52 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:39:52 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:39:52 INFO None 5152247: status FINISHED
2026-07-13 20:39:52 INFO None 5152248: status FINISHED
2026-07-13 20:39:52 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:39:52 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:39:52 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:39:52 INFO None 5152252: status FINISHED
2026-07-13 20:39:53 INFO None 5152253: status FINISHED
2026-07-13 20:39:53 INFO None 5152254: status RUNNING/PENDING
2026-07-13 20:39:53 INFO None 5152255: status RUNNING/PENDING
2026-07-13 20:39:53 INFO None 5152257: status FINISHED
2026-07-13 20:39:53 INFO None 5152258: status FINISHED
2026-07-13 20:39:53 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:39:53 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:39:53 INFO Jobs still running: ['5152244', '5152246', '5152249', '5152250', '5152251', '5152254', '5152255', '5152260', '5152261']. Waiting...
2026-07-13 20:40:08 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152247: status FINISHED
2026-07-13 20:40:08 INFO None 5152248: status FINISHED
2026-07-13 20:40:08 INFO None 5152249: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152252: status FINISHED
2026-07-13 20:40:08 INFO None 5152253: status FINISHED
2026-07-13 20:40:08 INFO None 5152254: status FINISHED
2026-07-13 20:40:08 INFO None 5152255: status FINISHED
2026-07-13 20:40:08 INFO None 5152257: status FINISHED
2026-07-13 20:40:08 INFO None 5152258: status FINISHED
2026-07-13 20:40:08 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:40:08 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:40:08 INFO Jobs still running: ['5152244', '5152246', '5152249', '5152250', '5152251', '5152260', '5152261']. Waiting...
2026-07-13 20:40:23 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:40:23 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:40:23 INFO None 5152247: status FINISHED
2026-07-13 20:40:23 INFO None 5152248: status FINISHED
2026-07-13 20:40:23 INFO None 5152249: status FINISHED
2026-07-13 20:40:23 INFO None 5152250: status RUNNING/PENDING
2026-07-13 20:40:23 INFO None 5152251: status RUNNING/PENDING
2026-07-13 20:40:23 INFO None 5152252: status FINISHED
2026-07-13 20:40:23 INFO None 5152253: status FINISHED
2026-07-13 20:40:23 INFO None 5152254: status FINISHED
2026-07-13 20:40:23 INFO None 5152255: status FINISHED
2026-07-13 20:40:23 INFO None 5152257: status FINISHED
2026-07-13 20:40:23 INFO None 5152258: status FINISHED
2026-07-13 20:40:23 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:40:23 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:40:23 INFO Jobs still running: ['5152244', '5152246', '5152250', '5152251', '5152260', '5152261']. Waiting...
2026-07-13 20:40:38 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:40:38 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:40:38 INFO None 5152247: status FINISHED
2026-07-13 20:40:39 INFO None 5152248: status FINISHED
2026-07-13 20:40:39 INFO None 5152249: status FINISHED
2026-07-13 20:40:39 INFO None 5152250: status FINISHED
2026-07-13 20:40:39 INFO None 5152251: status FINISHED
2026-07-13 20:40:39 INFO None 5152252: status FINISHED
2026-07-13 20:40:39 INFO None 5152253: status FINISHED
2026-07-13 20:40:39 INFO None 5152254: status FINISHED
2026-07-13 20:40:39 INFO None 5152255: status FINISHED
2026-07-13 20:40:39 INFO None 5152257: status FINISHED
2026-07-13 20:40:39 INFO None 5152258: status FINISHED
2026-07-13 20:40:39 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:40:39 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:40:39 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:40:54 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:40:54 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:40:54 INFO None 5152247: status FINISHED
2026-07-13 20:40:54 INFO None 5152248: status FINISHED
2026-07-13 20:40:54 INFO None 5152249: status FINISHED
2026-07-13 20:40:54 INFO None 5152250: status FINISHED
2026-07-13 20:40:54 INFO None 5152251: status FINISHED
2026-07-13 20:40:54 INFO None 5152252: status FINISHED
2026-07-13 20:40:54 INFO None 5152253: status FINISHED
2026-07-13 20:40:54 INFO None 5152254: status FINISHED
2026-07-13 20:40:54 INFO None 5152255: status FINISHED
2026-07-13 20:40:54 INFO None 5152257: status FINISHED
2026-07-13 20:40:54 INFO None 5152258: status FINISHED
2026-07-13 20:40:54 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:40:54 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:40:54 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:41:09 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:41:09 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:41:09 INFO None 5152247: status FINISHED
2026-07-13 20:41:09 INFO None 5152248: status FINISHED
2026-07-13 20:41:09 INFO None 5152249: status FINISHED
2026-07-13 20:41:10 INFO None 5152250: status FINISHED
2026-07-13 20:41:10 INFO None 5152251: status FINISHED
2026-07-13 20:41:10 INFO None 5152252: status FINISHED
2026-07-13 20:41:10 INFO None 5152253: status FINISHED
2026-07-13 20:41:10 INFO None 5152254: status FINISHED
2026-07-13 20:41:10 INFO None 5152255: status FINISHED
2026-07-13 20:41:10 INFO None 5152257: status FINISHED
2026-07-13 20:41:10 INFO None 5152258: status FINISHED
2026-07-13 20:41:10 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:41:10 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:41:10 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:41:25 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:41:25 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:41:25 INFO None 5152247: status FINISHED
2026-07-13 20:41:25 INFO None 5152248: status FINISHED
2026-07-13 20:41:25 INFO None 5152249: status FINISHED
2026-07-13 20:41:25 INFO None 5152250: status FINISHED
2026-07-13 20:41:25 INFO None 5152251: status FINISHED
2026-07-13 20:41:25 INFO None 5152252: status FINISHED
2026-07-13 20:41:25 INFO None 5152253: status FINISHED
2026-07-13 20:41:25 INFO None 5152254: status FINISHED
2026-07-13 20:41:25 INFO None 5152255: status FINISHED
2026-07-13 20:41:25 INFO None 5152257: status FINISHED
2026-07-13 20:41:25 INFO None 5152258: status FINISHED
2026-07-13 20:41:25 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:41:25 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:41:25 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:41:40 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:41:40 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:41:40 INFO None 5152247: status FINISHED
2026-07-13 20:41:40 INFO None 5152248: status FINISHED
2026-07-13 20:41:40 INFO None 5152249: status FINISHED
2026-07-13 20:41:40 INFO None 5152250: status FINISHED
2026-07-13 20:41:40 INFO None 5152251: status FINISHED
2026-07-13 20:41:40 INFO None 5152252: status FINISHED
2026-07-13 20:41:40 INFO None 5152253: status FINISHED
2026-07-13 20:41:40 INFO None 5152254: status FINISHED
2026-07-13 20:41:40 INFO None 5152255: status FINISHED
2026-07-13 20:41:40 INFO None 5152257: status FINISHED
2026-07-13 20:41:40 INFO None 5152258: status FINISHED
2026-07-13 20:41:40 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:41:40 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:41:40 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:41:55 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:41:55 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:41:55 INFO None 5152247: status FINISHED
2026-07-13 20:41:55 INFO None 5152248: status FINISHED
2026-07-13 20:41:55 INFO None 5152249: status FINISHED
2026-07-13 20:41:55 INFO None 5152250: status FINISHED
2026-07-13 20:41:55 INFO None 5152251: status FINISHED
2026-07-13 20:41:56 INFO None 5152252: status FINISHED
2026-07-13 20:41:56 INFO None 5152253: status FINISHED
2026-07-13 20:41:56 INFO None 5152254: status FINISHED
2026-07-13 20:41:56 INFO None 5152255: status FINISHED
2026-07-13 20:41:56 INFO None 5152257: status FINISHED
2026-07-13 20:41:56 INFO None 5152258: status FINISHED
2026-07-13 20:41:56 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:41:56 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:41:56 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:42:11 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:42:11 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:42:11 INFO None 5152247: status FINISHED
2026-07-13 20:42:11 INFO None 5152248: status FINISHED
2026-07-13 20:42:11 INFO None 5152249: status FINISHED
2026-07-13 20:42:11 INFO None 5152250: status FINISHED
2026-07-13 20:42:11 INFO None 5152251: status FINISHED
2026-07-13 20:42:11 INFO None 5152252: status FINISHED
2026-07-13 20:42:11 INFO None 5152253: status FINISHED
2026-07-13 20:42:11 INFO None 5152254: status FINISHED
2026-07-13 20:42:11 INFO None 5152255: status FINISHED
2026-07-13 20:42:11 INFO None 5152257: status FINISHED
2026-07-13 20:42:11 INFO None 5152258: status FINISHED
2026-07-13 20:42:11 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:42:11 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:42:11 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:42:26 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:42:26 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:42:26 INFO None 5152247: status FINISHED
2026-07-13 20:42:26 INFO None 5152248: status FINISHED
2026-07-13 20:42:26 INFO None 5152249: status FINISHED
2026-07-13 20:42:26 INFO None 5152250: status FINISHED
2026-07-13 20:42:26 INFO None 5152251: status FINISHED
2026-07-13 20:42:26 INFO None 5152252: status FINISHED
2026-07-13 20:42:26 INFO None 5152253: status FINISHED
2026-07-13 20:42:26 INFO None 5152254: status FINISHED
2026-07-13 20:42:26 INFO None 5152255: status FINISHED
2026-07-13 20:42:26 INFO None 5152257: status FINISHED
2026-07-13 20:42:26 INFO None 5152258: status FINISHED
2026-07-13 20:42:26 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:42:26 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:42:26 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:42:41 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:43:51 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:43:51 INFO None 5152247: status FINISHED
2026-07-13 20:43:51 INFO None 5152248: status FINISHED
2026-07-13 20:43:51 INFO None 5152249: status FINISHED
2026-07-13 20:43:51 INFO None 5152250: status FINISHED
2026-07-13 20:43:52 INFO None 5152251: status FINISHED
2026-07-13 20:43:52 INFO None 5152252: status FINISHED
2026-07-13 20:43:52 INFO None 5152253: status FINISHED
2026-07-13 20:43:52 INFO None 5152254: status FINISHED
2026-07-13 20:43:52 INFO None 5152255: status FINISHED
2026-07-13 20:43:52 INFO None 5152257: status FINISHED
2026-07-13 20:43:52 INFO None 5152258: status FINISHED
2026-07-13 20:43:52 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:43:52 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:43:52 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:44:07 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:44:07 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:44:07 INFO None 5152247: status FINISHED
2026-07-13 20:44:07 INFO None 5152248: status FINISHED
2026-07-13 20:44:07 INFO None 5152249: status FINISHED
2026-07-13 20:44:07 INFO None 5152250: status FINISHED
2026-07-13 20:44:07 INFO None 5152251: status FINISHED
2026-07-13 20:44:07 INFO None 5152252: status FINISHED
2026-07-13 20:44:07 INFO None 5152253: status FINISHED
2026-07-13 20:44:07 INFO None 5152254: status FINISHED
2026-07-13 20:44:07 INFO None 5152255: status FINISHED
2026-07-13 20:44:07 INFO None 5152257: status FINISHED
2026-07-13 20:44:07 INFO None 5152258: status FINISHED
2026-07-13 20:44:07 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:44:07 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:44:07 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:44:22 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:44:22 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:44:22 INFO None 5152247: status FINISHED
2026-07-13 20:44:22 INFO None 5152248: status FINISHED
2026-07-13 20:44:22 INFO None 5152249: status FINISHED
2026-07-13 20:44:22 INFO None 5152250: status FINISHED
2026-07-13 20:44:22 INFO None 5152251: status FINISHED
2026-07-13 20:44:22 INFO None 5152252: status FINISHED
2026-07-13 20:44:22 INFO None 5152253: status FINISHED
2026-07-13 20:44:22 INFO None 5152254: status FINISHED
2026-07-13 20:44:22 INFO None 5152255: status FINISHED
2026-07-13 20:44:22 INFO None 5152257: status FINISHED
2026-07-13 20:44:22 INFO None 5152258: status FINISHED
2026-07-13 20:44:22 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:44:22 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:44:22 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:44:37 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:44:38 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:44:38 INFO None 5152247: status FINISHED
2026-07-13 20:44:38 INFO None 5152248: status FINISHED
2026-07-13 20:44:38 INFO None 5152249: status FINISHED
2026-07-13 20:44:38 INFO None 5152250: status FINISHED
2026-07-13 20:44:38 INFO None 5152251: status FINISHED
2026-07-13 20:44:38 INFO None 5152252: status FINISHED
2026-07-13 20:44:38 INFO None 5152253: status FINISHED
2026-07-13 20:44:38 INFO None 5152254: status FINISHED
2026-07-13 20:44:38 INFO None 5152255: status FINISHED
2026-07-13 20:44:38 INFO None 5152257: status FINISHED
2026-07-13 20:44:38 INFO None 5152258: status FINISHED
2026-07-13 20:44:38 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:44:38 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:44:38 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:44:53 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:44:53 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:44:53 INFO None 5152247: status FINISHED
2026-07-13 20:44:53 INFO None 5152248: status FINISHED
2026-07-13 20:44:53 INFO None 5152249: status FINISHED
2026-07-13 20:44:53 INFO None 5152250: status FINISHED
2026-07-13 20:44:53 INFO None 5152251: status FINISHED
2026-07-13 20:44:53 INFO None 5152252: status FINISHED
2026-07-13 20:44:53 INFO None 5152253: status FINISHED
2026-07-13 20:44:53 INFO None 5152254: status FINISHED
2026-07-13 20:44:53 INFO None 5152255: status FINISHED
2026-07-13 20:44:53 INFO None 5152257: status FINISHED
2026-07-13 20:44:53 INFO None 5152258: status FINISHED
2026-07-13 20:44:53 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:44:53 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:44:53 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:45:08 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:45:08 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:45:08 INFO None 5152247: status FINISHED
2026-07-13 20:45:08 INFO None 5152248: status FINISHED
2026-07-13 20:45:08 INFO None 5152249: status FINISHED
2026-07-13 20:45:08 INFO None 5152250: status FINISHED
2026-07-13 20:45:08 INFO None 5152251: status FINISHED
2026-07-13 20:45:08 INFO None 5152252: status FINISHED
2026-07-13 20:45:08 INFO None 5152253: status FINISHED
2026-07-13 20:45:08 INFO None 5152254: status FINISHED
2026-07-13 20:45:09 INFO None 5152255: status FINISHED
2026-07-13 20:45:09 INFO None 5152257: status FINISHED
2026-07-13 20:45:09 INFO None 5152258: status FINISHED
2026-07-13 20:45:09 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:45:09 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:45:09 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:45:24 INFO None 5152244: status RUNNING/PENDING
2026-07-13 20:45:24 INFO None 5152246: status RUNNING/PENDING
2026-07-13 20:45:24 INFO None 5152247: status FINISHED
2026-07-13 20:45:24 INFO None 5152248: status FINISHED
2026-07-13 20:45:24 INFO None 5152249: status FINISHED
2026-07-13 20:45:24 INFO None 5152250: status FINISHED
2026-07-13 20:45:24 INFO None 5152251: status FINISHED
2026-07-13 20:45:24 INFO None 5152252: status FINISHED
2026-07-13 20:45:24 INFO None 5152253: status FINISHED
2026-07-13 20:45:24 INFO None 5152254: status FINISHED
2026-07-13 20:45:24 INFO None 5152255: status FINISHED
2026-07-13 20:45:24 INFO None 5152257: status FINISHED
2026-07-13 20:45:24 INFO None 5152258: status FINISHED
2026-07-13 20:45:24 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:45:24 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:45:24 INFO Jobs still running: ['5152244', '5152246', '5152260', '5152261']. Waiting...
2026-07-13 20:45:39 INFO None 5152244: status FINISHED
2026-07-13 20:45:39 INFO None 5152246: status FINISHED
2026-07-13 20:45:39 INFO None 5152247: status FINISHED
2026-07-13 20:45:39 INFO None 5152248: status FINISHED
2026-07-13 20:45:39 INFO None 5152249: status FINISHED
2026-07-13 20:45:39 INFO None 5152250: status FINISHED
2026-07-13 20:45:39 INFO None 5152251: status FINISHED
2026-07-13 20:45:39 INFO None 5152252: status FINISHED
2026-07-13 20:45:39 INFO None 5152253: status FINISHED
2026-07-13 20:45:39 INFO None 5152254: status FINISHED
2026-07-13 20:45:39 INFO None 5152255: status FINISHED
2026-07-13 20:45:39 INFO None 5152257: status FINISHED
2026-07-13 20:45:39 INFO None 5152258: status FINISHED
2026-07-13 20:45:39 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:45:39 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:45:39 INFO Jobs still running: ['5152260', '5152261']. Waiting...
2026-07-13 20:45:54 INFO None 5152244: status FINISHED
2026-07-13 20:45:54 INFO None 5152246: status FINISHED
2026-07-13 20:45:54 INFO None 5152247: status FINISHED
2026-07-13 20:45:54 INFO None 5152248: status FINISHED
2026-07-13 20:45:54 INFO None 5152249: status FINISHED
2026-07-13 20:45:54 INFO None 5152250: status FINISHED
2026-07-13 20:45:55 INFO None 5152251: status FINISHED
2026-07-13 20:45:55 INFO None 5152252: status FINISHED
2026-07-13 20:45:55 INFO None 5152253: status FINISHED
2026-07-13 20:45:55 INFO None 5152254: status FINISHED
2026-07-13 20:45:55 INFO None 5152255: status FINISHED
2026-07-13 20:45:55 INFO None 5152257: status FINISHED
2026-07-13 20:45:55 INFO None 5152258: status FINISHED
2026-07-13 20:45:55 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:45:55 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:45:55 INFO Jobs still running: ['5152260', '5152261']. Waiting...
2026-07-13 20:46:10 INFO None 5152244: status FINISHED
2026-07-13 20:46:10 INFO None 5152246: status FINISHED
2026-07-13 20:46:10 INFO None 5152247: status FINISHED
2026-07-13 20:46:10 INFO None 5152248: status FINISHED
2026-07-13 20:46:10 INFO None 5152249: status FINISHED
2026-07-13 20:46:10 INFO None 5152250: status FINISHED
2026-07-13 20:46:10 INFO None 5152251: status FINISHED
2026-07-13 20:46:10 INFO None 5152252: status FINISHED
2026-07-13 20:46:10 INFO None 5152253: status FINISHED
2026-07-13 20:46:10 INFO None 5152254: status FINISHED
2026-07-13 20:46:10 INFO None 5152255: status FINISHED
2026-07-13 20:46:10 INFO None 5152257: status FINISHED
2026-07-13 20:46:10 INFO None 5152258: status FINISHED
2026-07-13 20:46:10 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:46:10 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:46:10 INFO Jobs still running: ['5152260', '5152261']. Waiting...
2026-07-13 20:46:25 INFO None 5152244: status FINISHED
2026-07-13 20:46:25 INFO None 5152246: status FINISHED
2026-07-13 20:46:25 INFO None 5152247: status FINISHED
2026-07-13 20:46:25 INFO None 5152248: status FINISHED
2026-07-13 20:46:25 INFO None 5152249: status FINISHED
2026-07-13 20:46:25 INFO None 5152250: status FINISHED
2026-07-13 20:46:25 INFO None 5152251: status FINISHED
2026-07-13 20:46:25 INFO None 5152252: status FINISHED
2026-07-13 20:46:25 INFO None 5152253: status FINISHED
2026-07-13 20:46:25 INFO None 5152254: status FINISHED
2026-07-13 20:46:25 INFO None 5152255: status FINISHED
2026-07-13 20:46:25 INFO None 5152257: status FINISHED
2026-07-13 20:46:25 INFO None 5152258: status FINISHED
2026-07-13 20:46:25 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:46:25 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:46:25 INFO Jobs still running: ['5152260', '5152261']. Waiting...
2026-07-13 20:46:40 INFO None 5152244: status FINISHED
2026-07-13 20:46:41 INFO None 5152246: status FINISHED
2026-07-13 20:46:41 INFO None 5152247: status FINISHED
2026-07-13 20:46:41 INFO None 5152248: status FINISHED
2026-07-13 20:46:41 INFO None 5152249: status FINISHED
2026-07-13 20:46:41 INFO None 5152250: status FINISHED
2026-07-13 20:46:41 INFO None 5152251: status FINISHED
2026-07-13 20:46:41 INFO None 5152252: status FINISHED
2026-07-13 20:46:41 INFO None 5152253: status FINISHED
2026-07-13 20:46:41 INFO None 5152254: status FINISHED
2026-07-13 20:46:41 INFO None 5152255: status FINISHED
2026-07-13 20:46:41 INFO None 5152257: status FINISHED
2026-07-13 20:46:41 INFO None 5152258: status FINISHED
2026-07-13 20:46:41 INFO None 5152260: status RUNNING/PENDING
2026-07-13 20:46:41 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:46:41 INFO Jobs still running: ['5152260', '5152261']. Waiting...
2026-07-13 20:46:56 INFO None 5152244: status FINISHED
2026-07-13 20:46:56 INFO None 5152246: status FINISHED
2026-07-13 20:46:56 INFO None 5152247: status FINISHED
2026-07-13 20:46:56 INFO None 5152248: status FINISHED
2026-07-13 20:46:56 INFO None 5152249: status FINISHED
2026-07-13 20:46:56 INFO None 5152250: status FINISHED
2026-07-13 20:46:56 INFO None 5152251: status FINISHED
2026-07-13 20:46:56 INFO None 5152252: status FINISHED
2026-07-13 20:46:56 INFO None 5152253: status FINISHED
2026-07-13 20:46:56 INFO None 5152254: status FINISHED
2026-07-13 20:46:56 INFO None 5152255: status FINISHED
2026-07-13 20:46:56 INFO None 5152257: status FINISHED
2026-07-13 20:46:56 INFO None 5152258: status FINISHED
2026-07-13 20:46:56 INFO None 5152260: status FINISHED
2026-07-13 20:46:56 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:46:56 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:47:11 INFO None 5152244: status FINISHED
2026-07-13 20:47:11 INFO None 5152246: status FINISHED
2026-07-13 20:47:11 INFO None 5152247: status FINISHED
2026-07-13 20:47:11 INFO None 5152248: status FINISHED
2026-07-13 20:47:11 INFO None 5152249: status FINISHED
2026-07-13 20:47:11 INFO None 5152250: status FINISHED
2026-07-13 20:47:11 INFO None 5152251: status FINISHED
2026-07-13 20:47:11 INFO None 5152252: status FINISHED
2026-07-13 20:47:11 INFO None 5152253: status FINISHED
2026-07-13 20:47:11 INFO None 5152254: status FINISHED
2026-07-13 20:47:11 INFO None 5152255: status FINISHED
2026-07-13 20:47:11 INFO None 5152257: status FINISHED
2026-07-13 20:47:11 INFO None 5152258: status FINISHED
2026-07-13 20:47:11 INFO None 5152260: status FINISHED
2026-07-13 20:47:11 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:47:11 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:47:27 INFO None 5152244: status FINISHED
2026-07-13 20:47:27 INFO None 5152246: status FINISHED
2026-07-13 20:47:27 INFO None 5152247: status FINISHED
2026-07-13 20:47:27 INFO None 5152248: status FINISHED
2026-07-13 20:47:27 INFO None 5152249: status FINISHED
2026-07-13 20:47:27 INFO None 5152250: status FINISHED
2026-07-13 20:47:27 INFO None 5152251: status FINISHED
2026-07-13 20:47:27 INFO None 5152252: status FINISHED
2026-07-13 20:47:27 INFO None 5152253: status FINISHED
2026-07-13 20:47:27 INFO None 5152254: status FINISHED
2026-07-13 20:47:27 INFO None 5152255: status FINISHED
2026-07-13 20:47:27 INFO None 5152257: status FINISHED
2026-07-13 20:47:27 INFO None 5152258: status FINISHED
2026-07-13 20:47:27 INFO None 5152260: status FINISHED
2026-07-13 20:47:27 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:47:27 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:47:42 INFO None 5152244: status FINISHED
2026-07-13 20:47:42 INFO None 5152246: status FINISHED
2026-07-13 20:47:42 INFO None 5152247: status FINISHED
2026-07-13 20:47:42 INFO None 5152248: status FINISHED
2026-07-13 20:47:42 INFO None 5152249: status FINISHED
2026-07-13 20:47:42 INFO None 5152250: status FINISHED
2026-07-13 20:47:42 INFO None 5152251: status FINISHED
2026-07-13 20:47:42 INFO None 5152252: status FINISHED
2026-07-13 20:47:42 INFO None 5152253: status FINISHED
2026-07-13 20:47:42 INFO None 5152254: status FINISHED
2026-07-13 20:47:42 INFO None 5152255: status FINISHED
2026-07-13 20:47:42 INFO None 5152257: status FINISHED
2026-07-13 20:47:42 INFO None 5152258: status FINISHED
2026-07-13 20:47:42 INFO None 5152260: status FINISHED
2026-07-13 20:47:42 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:47:42 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:47:57 INFO None 5152244: status FINISHED
2026-07-13 20:47:57 INFO None 5152246: status FINISHED
2026-07-13 20:47:57 INFO None 5152247: status FINISHED
2026-07-13 20:47:57 INFO None 5152248: status FINISHED
2026-07-13 20:47:57 INFO None 5152249: status FINISHED
2026-07-13 20:47:57 INFO None 5152250: status FINISHED
2026-07-13 20:47:57 INFO None 5152251: status FINISHED
2026-07-13 20:47:57 INFO None 5152252: status FINISHED
2026-07-13 20:47:57 INFO None 5152253: status FINISHED
2026-07-13 20:47:57 INFO None 5152254: status FINISHED
2026-07-13 20:47:57 INFO None 5152255: status FINISHED
2026-07-13 20:47:57 INFO None 5152257: status FINISHED
2026-07-13 20:47:57 INFO None 5152258: status FINISHED
2026-07-13 20:47:57 INFO None 5152260: status FINISHED
2026-07-13 20:47:57 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:47:57 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:48:12 INFO None 5152244: status FINISHED
2026-07-13 20:48:12 INFO None 5152246: status FINISHED
2026-07-13 20:48:12 INFO None 5152247: status FINISHED
2026-07-13 20:48:12 INFO None 5152248: status FINISHED
2026-07-13 20:48:12 INFO None 5152249: status FINISHED
2026-07-13 20:48:12 INFO None 5152250: status FINISHED
2026-07-13 20:48:13 INFO None 5152251: status FINISHED
2026-07-13 20:48:13 INFO None 5152252: status FINISHED
2026-07-13 20:48:13 INFO None 5152253: status FINISHED
2026-07-13 20:48:13 INFO None 5152254: status FINISHED
2026-07-13 20:48:13 INFO None 5152255: status FINISHED
2026-07-13 20:48:13 INFO None 5152257: status FINISHED
2026-07-13 20:48:13 INFO None 5152258: status FINISHED
2026-07-13 20:48:13 INFO None 5152260: status FINISHED
2026-07-13 20:48:13 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:48:13 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:48:28 INFO None 5152244: status FINISHED
2026-07-13 20:48:28 INFO None 5152246: status FINISHED
2026-07-13 20:48:28 INFO None 5152247: status FINISHED
2026-07-13 20:48:28 INFO None 5152248: status FINISHED
2026-07-13 20:48:28 INFO None 5152249: status FINISHED
2026-07-13 20:48:28 INFO None 5152250: status FINISHED
2026-07-13 20:48:28 INFO None 5152251: status FINISHED
2026-07-13 20:48:28 INFO None 5152252: status FINISHED
2026-07-13 20:48:28 INFO None 5152253: status FINISHED
2026-07-13 20:48:28 INFO None 5152254: status FINISHED
2026-07-13 20:48:28 INFO None 5152255: status FINISHED
2026-07-13 20:48:28 INFO None 5152257: status FINISHED
2026-07-13 20:48:28 INFO None 5152258: status FINISHED
2026-07-13 20:48:28 INFO None 5152260: status FINISHED
2026-07-13 20:48:28 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:48:28 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:48:43 INFO None 5152244: status FINISHED
2026-07-13 20:48:43 INFO None 5152246: status FINISHED
2026-07-13 20:48:43 INFO None 5152247: status FINISHED
2026-07-13 20:48:43 INFO None 5152248: status FINISHED
2026-07-13 20:48:43 INFO None 5152249: status FINISHED
2026-07-13 20:48:43 INFO None 5152250: status FINISHED
2026-07-13 20:48:43 INFO None 5152251: status FINISHED
2026-07-13 20:48:43 INFO None 5152252: status FINISHED (not in squeue)
2026-07-13 20:48:43 INFO None 5152253: status FINISHED (not in squeue)
2026-07-13 20:48:43 INFO None 5152254: status FINISHED
2026-07-13 20:48:43 INFO None 5152255: status FINISHED
2026-07-13 20:48:43 INFO None 5152257: status FINISHED
2026-07-13 20:48:43 INFO None 5152258: status FINISHED
2026-07-13 20:48:43 INFO None 5152260: status FINISHED
2026-07-13 20:48:43 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:48:43 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:48:58 INFO None 5152244: status FINISHED
2026-07-13 20:48:58 INFO None 5152246: status FINISHED
2026-07-13 20:48:59 INFO None 5152247: status FINISHED
2026-07-13 20:48:59 INFO None 5152248: status FINISHED
2026-07-13 20:48:59 INFO None 5152249: status FINISHED
2026-07-13 20:48:59 INFO None 5152250: status FINISHED
2026-07-13 20:48:59 INFO None 5152251: status FINISHED
2026-07-13 20:48:59 INFO None 5152252: status FINISHED (not in squeue)
2026-07-13 20:48:59 INFO None 5152253: status FINISHED (not in squeue)
2026-07-13 20:48:59 INFO None 5152254: status FINISHED
2026-07-13 20:48:59 INFO None 5152255: status FINISHED
2026-07-13 20:48:59 INFO None 5152257: status FINISHED
2026-07-13 20:48:59 INFO None 5152258: status FINISHED
2026-07-13 20:48:59 INFO None 5152260: status FINISHED
2026-07-13 20:48:59 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:48:59 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:49:14 INFO None 5152244: status FINISHED
2026-07-13 20:49:14 INFO None 5152246: status FINISHED
2026-07-13 20:49:14 INFO None 5152247: status FINISHED
2026-07-13 20:49:14 INFO None 5152248: status FINISHED
2026-07-13 20:49:14 INFO None 5152249: status FINISHED
2026-07-13 20:49:14 INFO None 5152250: status FINISHED
2026-07-13 20:49:14 INFO None 5152251: status FINISHED
2026-07-13 20:49:14 INFO None 5152252: status FINISHED (not in squeue)
2026-07-13 20:49:14 INFO None 5152253: status FINISHED (not in squeue)
2026-07-13 20:49:14 INFO None 5152254: status FINISHED
2026-07-13 20:49:14 INFO None 5152255: status FINISHED
2026-07-13 20:49:14 INFO None 5152257: status FINISHED
2026-07-13 20:49:14 INFO None 5152258: status FINISHED
2026-07-13 20:49:14 INFO None 5152260: status FINISHED
2026-07-13 20:49:14 INFO None 5152261: status RUNNING/PENDING
2026-07-13 20:49:14 INFO Jobs still running: ['5152261']. Waiting...
2026-07-13 20:49:29 INFO None 5152244: status FINISHED
2026-07-13 20:49:29 INFO None 5152246: status FINISHED
2026-07-13 20:49:29 INFO None 5152247: status FINISHED
2026-07-13 20:49:29 INFO None 5152248: status FINISHED
2026-07-13 20:49:29 INFO None 5152249: status FINISHED
2026-07-13 20:49:29 INFO None 5152250: status FINISHED
2026-07-13 20:49:29 INFO None 5152251: status FINISHED
2026-07-13 20:49:29 INFO None 5152252: status FINISHED (not in squeue)
2026-07-13 20:49:29 INFO None 5152253: status FINISHED (not in squeue)
2026-07-13 20:49:29 INFO None 5152254: status FINISHED
2026-07-13 20:49:29 INFO None 5152255: status FINISHED
2026-07-13 20:49:29 INFO None 5152257: status FINISHED
2026-07-13 20:49:29 INFO None 5152258: status FINISHED
2026-07-13 20:49:30 INFO None 5152260: status FINISHED
2026-07-13 20:49:30 INFO None 5152261: status FINISHED
2026-07-13 20:49:30 INFO Jobs ['5152244', '5152246', '5152247', '5152248', '5152249', '5152250', '5152251', '5152252', '5152253', '5152254', '5152255', '5152257', '5152258', '5152260', '5152261'] have finished
2026-07-13 20:49:30 INFO Checking restart files were created ...
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/end.2020020601_8_ENS1.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/end.2020020601_8_ENS2.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/end.2020020601_8_ENS3.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/end.2020020601_8_ENS4.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/end.2020020601_8_ENS5.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/end.2020020601_8_ENS6.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/end.2020020601_8_ENS7.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/end.2020020601_8_ENS8.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/end.2020020601_8_ENS9.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/end.2020020601_8_ENS10.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/end.2020020601_8_ENS11.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/end.2020020601_8_ENS12.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/end.2020020601_8_ENS13.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/end.2020020601_8_ENS14.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO ModelType.CHIMERE | restart_file exists for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/end.2020020601_8_ENS15.nc(3005806795 bytes)
2026-07-13 20:49:30 INFO  Run_model() completed successfully.
2026-07-13 20:49:30 INFO [TIME] after_model_set_simulated_time current_time=2020-02-06 01:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:49:30 INFO [TIME] gregorian_conversion simulated_time=2020-02-06 09:00:00 days=153072 seconds=32400
2026-07-13 20:49:30 INFO Saving chimere's output files title timestamp (that is the starting time of the run)
2026-07-13 20:49:30 INFO [TIME] increment current_time 2020-02-06 01:00:00 -> 2020-02-06 09:00:00
2026-07-13 20:49:30 INFO [TIME] after_increment_before_assimilation current_time=2020-02-06 09:00:00 simulated_time=2020-02-06 09:00:00 dt=0 days 01:00:00 end_time=2020-02-07 23:00:00
2026-07-13 20:49:30 INFO ---------->>> Running process_satellite_data()
2026-07-13 20:49:30 INFO Orbit file found: C03/2020/02/E3/S5p_NO2_12002.nc
2026-07-13 20:49:30 INFO ---------->>> Running run_obs_converter()
2026-07-13 20:49:30 INFO Obs sequence file already exists: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-13 20:49:30 INFO [DART] obs_seq created: /ccc/work/cont003/gen7232/demoling/DART/observations/obs_converters/S5P_TROPOMI_L3/data/NO2/C03/2020/02/E3/obs_seq_33067_153072.out
2026-07-13 20:49:30 INFO ---------->>> Running DART
2026-07-13 20:49:30 INFO The timestamp in DART results' titles does not follows chimere's logic: the simulated_time is used
2026-07-13 20:49:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS1/chim_ENS1_2020020609_1_out_toDART.nc
2026-07-13 20:49:30 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS2/chim_ENS2_2020020609_1_out_toDART.nc
2026-07-13 20:49:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS3/chim_ENS3_2020020609_1_out_toDART.nc
2026-07-13 20:49:31 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS4/chim_ENS4_2020020609_1_out_toDART.nc
2026-07-13 20:49:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS5/chim_ENS5_2020020609_1_out_toDART.nc
2026-07-13 20:49:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS6/chim_ENS6_2020020609_1_out_toDART.nc
2026-07-13 20:49:32 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS7/chim_ENS7_2020020609_1_out_toDART.nc
2026-07-13 20:49:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS8/chim_ENS8_2020020609_1_out_toDART.nc
2026-07-13 20:49:33 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS9/chim_ENS9_2020020609_1_out_toDART.nc
2026-07-13 20:49:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS10/chim_ENS10_2020020609_1_out_toDART.nc
2026-07-13 20:49:34 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS11/chim_ENS11_2020020609_1_out_toDART.nc
2026-07-13 20:49:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS12/chim_ENS12_2020020609_1_out_toDART.nc
2026-07-13 20:49:35 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS13/chim_ENS13_2020020609_1_out_toDART.nc
2026-07-13 20:49:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS14/chim_ENS14_2020020609_1_out_toDART.nc
2026-07-13 20:49:36 INFO From /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020601_8_out.nc created /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyEMIdmp_0607_15m_low_v2/ENS15/chim_ENS15_2020020609_1_out_toDART.nc
2026-07-13 20:49:37 INFO Replacement input_template.nml → input.nml completed successfully.
2026-07-13 20:49:37 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_input_list.txt
2026-07-13 20:49:37 INFO Wrote: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/filter_output_list.txt
2026-07-13 20:49:37 INFO Replacement run_filter.template.bsh → run_filter.bsh completed successfully.
2026-07-13 20:49:37 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work
2026-07-13 20:49:37 INFO [CMD] Running: /ccc/work/cont003/gen7232/demoling/DART/models/chimere/work/run_filter.bsh
2026-07-13 20:49:46 INFO Found: []
2026-07-13 20:49:46 INFO No job id returned by command ./run_filter.bsh
2026-07-13 20:49:46 INFO No monitoring will be performed
2026-07-13 20:49:46 INFO Moving DART output files to analysis and preassim directories for date 2020020609 if present ...
2026-07-13 20:49:46 INFO Moved 'analysis_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0013.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0007.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_sd.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0006.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0012.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0010.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'preassim_member_0004.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:46 INFO Moved 'analysis_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:47 INFO Moved 'analysis_member_0011.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0003.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0001.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0014.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_mean.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0015.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0009.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'analysis_member_0008.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/analysis/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0002.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Moved 'preassim_member_0005.nc' to '/ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/preassim/2020020609'
2026-07-13 20:49:47 INFO Computing differences between analysis/preassim means (ana - preassim)...
2026-07-13 20:49:50 ERROR Failed to compute Mean Analysis Increment: "No variable named 'airm'. Variables on the dataset include ['psfc', 'NO2', 'pres', 'NO', 'EMISA', 'EMISB', 'time']"
2026-07-13 20:49:50 INFO run_dart() is DONE.
2026-07-13 20:49:50 INFO ---------->>> Running update_pollutant_in_end()
2026-07-13 20:49:50 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-13 20:49:50 INFO No previous orbit memory found.
2026-07-13 20:49:50 INFO Emission correction applied with pixel-based damping.
2026-07-13 20:49:50 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISA_ENS1_2020020609.nc
2026-07-13 20:49:51 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-13 20:49:51 INFO  > [INFO] Vertical level mismatch detected (DART: 20, EMIS: 7). Scaling LEVEL 0 only.
2026-07-13 20:49:51 INFO No previous orbit memory found.
2026-07-13 20:49:51 INFO Emission correction applied with pixel-based damping.
2026-07-13 20:49:51 INFO Orbit memory saved: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_DART/cyEMIdmp_0607_15m_low_v2/ratio_memory_file_EMISB_ENS1_2020020609.nc
2026-07-13 20:49:52 INFO Emission update completed using pixel-based orbit memory and damping.
2026-07-13 20:49:52 INFO Computing differences between posterior vs. original CHIMERE outputs ...
[2026-07-14T09:18:38.374] error: *** JOB 5152187 ON irene4263 CANCELLED AT 2026-07-14T09:18:38 DUE to SIGNAL Terminated ***
