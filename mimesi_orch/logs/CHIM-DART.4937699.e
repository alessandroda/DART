+ SCRIPT_PID=2285976
+ /bin/bash -x /tmp/tmp.PHRXmX4qbm
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
2026-06-21 19:15:51 INFO 
███    ███ ██ ███    ███ ███████ ███████ ██
████  ████ ██ ████  ████ ██      ██      ██
██ ████ ██ ██ ██ ████ ██ █████   ███████ ██
██  ██  ██ ██ ██  ██  ██ ██           ██ ██
██      ██ ██ ██      ██ ███████ ███████ ██



2026-06-21 19:15:51 INFO [PIPELINE] =======================================
2026-06-21 19:15:51 INFO [PIPELINE] Starting chimere–DART orchestrator
2026-06-21 19:15:51 INFO [PIPELINE] Config file: config/config_irene_IM.yaml
2026-06-21 19:15:51 INFO [PIPELINE] Run dir: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart
2026-06-21 19:15:51 INFO [PIPELINE] Log file: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/mimesi_orchestrator_logs/chimere_DART_20260621_191551.log
2026-06-21 19:15:51 INFO [PIPELINE] =======================================
2026-06-21 19:15:51 INFO Running assimilation with model_type=ModelType.CHIMERE
2026-06-21 19:15:51 INFO Using scheduler=Scheduler.SLURM, queue=rome
2026-06-21 19:15:51 INFO [STEP] ---- TIME LOOP START ----
2026-06-21 19:15:51 INFO [TIME] step_start current_time=2020-02-14 15:00:00 simulated_time=None dt=0 days 01:00:00 end_time=2020-02-15 23:00:00
2026-06-21 19:15:51 INFO [TIME] window start=2020-02-14 15:00:00 end=2020-02-15 00:00:00 run_hours=9 has_assimilation=False
2026-06-21 19:15:51 INFO Creating directories and links for ENS1 to run chimere's parallel part
2026-06-21 19:15:51 INFO Linking EMIS ...
2026-06-21 19:15:51 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens01.nc
2026-06-21 19:15:52 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens01.nc
2026-06-21 19:15:52 INFO Linking END ...
2026-06-21 19:15:52 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021300_24_ENS1.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:15:52 INFO >> Checking links...
2026-06-21 19:16:03 INFO >> All links are good for ENS1  ...
2026-06-21 19:16:03 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:16:12 INFO Hourly dataset computed and listing created
2026-06-21 19:16:23 INFO Hourly dataset computed
2026-06-21 19:16:23 INFO Creating directories and links for ENS2 to run chimere's parallel part
2026-06-21 19:16:23 INFO Linking EMIS ...
2026-06-21 19:16:23 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens02.nc
2026-06-21 19:16:24 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens02.nc
2026-06-21 19:16:24 INFO Linking END ...
2026-06-21 19:16:24 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021300_24_ENS2.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:16:24 INFO >> Checking links...
2026-06-21 19:16:36 INFO >> All links are good for ENS2  ...
2026-06-21 19:16:36 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:16:37 INFO Hourly dataset computed and listing created
2026-06-21 19:16:39 INFO Hourly dataset computed
2026-06-21 19:16:39 INFO Creating directories and links for ENS3 to run chimere's parallel part
2026-06-21 19:16:39 INFO Linking EMIS ...
2026-06-21 19:16:39 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens03.nc
2026-06-21 19:16:40 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens03.nc
2026-06-21 19:16:40 INFO Linking END ...
2026-06-21 19:16:40 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021300_24_ENS3.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:16:40 INFO >> Checking links...
2026-06-21 19:16:51 INFO >> All links are good for ENS3  ...
2026-06-21 19:16:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:16:52 INFO Hourly dataset computed and listing created
2026-06-21 19:16:54 INFO Hourly dataset computed
2026-06-21 19:16:54 INFO Creating directories and links for ENS4 to run chimere's parallel part
2026-06-21 19:16:54 INFO Linking EMIS ...
2026-06-21 19:16:54 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens04.nc
2026-06-21 19:16:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens04.nc
2026-06-21 19:16:55 INFO Linking END ...
2026-06-21 19:16:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021300_24_ENS4.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:16:55 INFO >> Checking links...
2026-06-21 19:17:08 INFO >> All links are good for ENS4  ...
2026-06-21 19:17:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:17:09 INFO Hourly dataset computed and listing created
2026-06-21 19:17:11 INFO Hourly dataset computed
2026-06-21 19:17:11 INFO Creating directories and links for ENS5 to run chimere's parallel part
2026-06-21 19:17:11 INFO Linking EMIS ...
2026-06-21 19:17:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens05.nc
2026-06-21 19:17:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens05.nc
2026-06-21 19:17:12 INFO Linking END ...
2026-06-21 19:17:12 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021300_24_ENS5.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:17:12 INFO >> Checking links...
2026-06-21 19:17:25 INFO >> All links are good for ENS5  ...
2026-06-21 19:17:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:17:26 INFO Hourly dataset computed and listing created
2026-06-21 19:17:28 INFO Hourly dataset computed
2026-06-21 19:17:28 INFO Creating directories and links for ENS6 to run chimere's parallel part
2026-06-21 19:17:28 INFO Linking EMIS ...
2026-06-21 19:17:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens06.nc
2026-06-21 19:17:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens06.nc
2026-06-21 19:17:29 INFO Linking END ...
2026-06-21 19:17:29 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021300_24_ENS6.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:17:29 INFO >> Checking links...
2026-06-21 19:17:42 INFO >> All links are good for ENS6  ...
2026-06-21 19:17:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:17:43 INFO Hourly dataset computed and listing created
2026-06-21 19:17:45 INFO Hourly dataset computed
2026-06-21 19:17:45 INFO Creating directories and links for ENS7 to run chimere's parallel part
2026-06-21 19:17:45 INFO Linking EMIS ...
2026-06-21 19:17:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens07.nc
2026-06-21 19:17:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens07.nc
2026-06-21 19:17:46 INFO Linking END ...
2026-06-21 19:17:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021300_24_ENS7.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:17:46 INFO >> Checking links...
2026-06-21 19:17:59 INFO >> All links are good for ENS7  ...
2026-06-21 19:17:59 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:18:00 INFO Hourly dataset computed and listing created
2026-06-21 19:18:03 INFO Hourly dataset computed
2026-06-21 19:18:03 INFO Creating directories and links for ENS8 to run chimere's parallel part
2026-06-21 19:18:03 INFO Linking EMIS ...
2026-06-21 19:18:03 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens08.nc
2026-06-21 19:18:04 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens08.nc
2026-06-21 19:18:04 INFO Linking END ...
2026-06-21 19:18:04 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021300_24_ENS8.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:18:04 INFO >> Checking links...
2026-06-21 19:18:17 INFO >> All links are good for ENS8  ...
2026-06-21 19:18:17 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:18:18 INFO Hourly dataset computed and listing created
2026-06-21 19:18:20 INFO Hourly dataset computed
2026-06-21 19:18:20 INFO Creating directories and links for ENS9 to run chimere's parallel part
2026-06-21 19:18:20 INFO Linking EMIS ...
2026-06-21 19:18:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens09.nc
2026-06-21 19:18:21 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens09.nc
2026-06-21 19:18:21 INFO Linking END ...
2026-06-21 19:18:21 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021300_24_ENS9.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:18:21 INFO >> Checking links...
2026-06-21 19:18:34 INFO >> All links are good for ENS9  ...
2026-06-21 19:18:34 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:18:35 INFO Hourly dataset computed and listing created
2026-06-21 19:18:37 INFO Hourly dataset computed
2026-06-21 19:18:37 INFO Creating directories and links for ENS10 to run chimere's parallel part
2026-06-21 19:18:37 INFO Linking EMIS ...
2026-06-21 19:18:37 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens10.nc
2026-06-21 19:18:38 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens10.nc
2026-06-21 19:18:38 INFO Linking END ...
2026-06-21 19:18:38 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021300_24_ENS10.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:18:38 INFO >> Checking links...
2026-06-21 19:18:51 INFO >> All links are good for ENS10  ...
2026-06-21 19:18:51 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:18:52 INFO Hourly dataset computed and listing created
2026-06-21 19:18:54 INFO Hourly dataset computed
2026-06-21 19:18:54 INFO Creating directories and links for ENS11 to run chimere's parallel part
2026-06-21 19:18:54 INFO Linking EMIS ...
2026-06-21 19:18:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens11.nc
2026-06-21 19:18:55 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens11.nc
2026-06-21 19:18:55 INFO Linking END ...
2026-06-21 19:18:55 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021300_24_ENS11.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:18:55 INFO >> Checking links...
2026-06-21 19:19:08 INFO >> All links are good for ENS11  ...
2026-06-21 19:19:08 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:19:09 INFO Hourly dataset computed and listing created
2026-06-21 19:19:11 INFO Hourly dataset computed
2026-06-21 19:19:11 INFO Creating directories and links for ENS12 to run chimere's parallel part
2026-06-21 19:19:11 INFO Linking EMIS ...
2026-06-21 19:19:11 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens12.nc
2026-06-21 19:19:12 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens12.nc
2026-06-21 19:19:12 INFO Linking END ...
2026-06-21 19:19:12 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021300_24_ENS12.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:19:12 INFO >> Checking links...
2026-06-21 19:19:25 INFO >> All links are good for ENS12  ...
2026-06-21 19:19:25 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:19:26 INFO Hourly dataset computed and listing created
2026-06-21 19:19:28 INFO Hourly dataset computed
2026-06-21 19:19:28 INFO Creating directories and links for ENS13 to run chimere's parallel part
2026-06-21 19:19:28 INFO Linking EMIS ...
2026-06-21 19:19:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens13.nc
2026-06-21 19:19:29 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens13.nc
2026-06-21 19:19:29 INFO Linking END ...
2026-06-21 19:19:29 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021300_24_ENS13.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:19:29 INFO >> Checking links...
2026-06-21 19:19:42 INFO >> All links are good for ENS13  ...
2026-06-21 19:19:42 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:19:43 INFO Hourly dataset computed and listing created
2026-06-21 19:19:45 INFO Hourly dataset computed
2026-06-21 19:19:45 INFO Creating directories and links for ENS14 to run chimere's parallel part
2026-06-21 19:19:45 INFO Linking EMIS ...
2026-06-21 19:19:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens14.nc
2026-06-21 19:19:46 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens14.nc
2026-06-21 19:19:46 INFO Linking END ...
2026-06-21 19:19:46 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021300_24_ENS14.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:19:46 INFO >> Checking links...
2026-06-21 19:19:57 INFO >> All links are good for ENS14  ...
2026-06-21 19:19:57 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:19:58 INFO Hourly dataset computed and listing created
2026-06-21 19:20:00 INFO Hourly dataset computed
2026-06-21 19:20:00 INFO Creating directories and links for ENS15 to run chimere's parallel part
2026-06-21 19:20:00 INFO Linking EMIS ...
2026-06-21 19:20:00 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Friday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Friday.s.ens15.nc
2026-06-21 19:20:01 INFO Copied ensured: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/EMIS.EUROCOMEX3.02.Saturday.s.nc of /ccc/scratch/cont003/gen7232/demoling/PERTURBED/EMIS_ANT/results_doubled/EMIS.EUROCOMEX3.02.Saturday.s.ens15.nc
2026-06-21 19:20:01 INFO Linking END ...
2026-06-21 19:20:01 INFO Symlink created: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021300_24_ENS15.nc -> /ccc/scratch/cont003/gen7232/demoling/OUT_Chimere/CHIMOUT-EUROCOMEX3_2023_0102_CR_reruned/end.2020021300_24_EUROCOMEX3_2023_0102_CR_reruned.nc
2026-06-21 19:20:01 INFO >> Checking links...
2026-06-21 19:20:12 INFO >> All links are good for ENS15  ...
2026-06-21 19:20:12 INFO Computing BOUNs and exdomouts for the specific hours to run ...
2026-06-21 19:20:13 INFO Hourly dataset computed and listing created
2026-06-21 19:20:16 INFO Hourly dataset computed
2026-06-21 19:20:16 INFO ---------->>> Running CHIMERE model from 2020-02-14 15:00:00 to 2020-02-15 00:00:00
2026-06-21 19:20:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1
2026-06-21 19:20:16 INFO The END file used for ENS1 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021300_24_ENS1.nc
2026-06-21 19:20:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS1.par completed successfully.
2026-06-21 19:20:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:16 INFO Replacement submit_p_template.sh → submit_p_1.sh completed successfully.
2026-06-21 19:20:16 INFO Queuing job for member 1...
2026-06-21 19:20:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:16 INFO [CMD] Running: ccc_msub ./submit_p_1.sh
2026-06-21 19:20:17 INFO Found: ['4937713']
2026-06-21 19:20:22 INFO [TGCC-IRENE] Submitted job with ID:['4937713']
2026-06-21 19:20:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2
2026-06-21 19:20:22 INFO The END file used for ENS2 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021300_24_ENS2.nc
2026-06-21 19:20:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS2.par completed successfully.
2026-06-21 19:20:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:22 INFO Replacement submit_p_template.sh → submit_p_2.sh completed successfully.
2026-06-21 19:20:22 INFO Queuing job for member 2...
2026-06-21 19:20:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:22 INFO [CMD] Running: ccc_msub ./submit_p_2.sh
2026-06-21 19:20:22 INFO Found: ['4937714']
2026-06-21 19:20:27 INFO [TGCC-IRENE] Submitted job with ID:['4937714']
2026-06-21 19:20:27 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:27 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3
2026-06-21 19:20:27 INFO The END file used for ENS3 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021300_24_ENS3.nc
2026-06-21 19:20:27 INFO Replacement chimere.template_ensemble.par → chimere.ENS3.par completed successfully.
2026-06-21 19:20:27 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:27 INFO Replacement submit_p_template.sh → submit_p_3.sh completed successfully.
2026-06-21 19:20:27 INFO Queuing job for member 3...
2026-06-21 19:20:27 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:27 INFO [CMD] Running: ccc_msub ./submit_p_3.sh
2026-06-21 19:20:28 INFO Found: ['4937715']
2026-06-21 19:20:33 INFO [TGCC-IRENE] Submitted job with ID:['4937715']
2026-06-21 19:20:33 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:33 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4
2026-06-21 19:20:33 INFO The END file used for ENS4 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021300_24_ENS4.nc
2026-06-21 19:20:33 INFO Replacement chimere.template_ensemble.par → chimere.ENS4.par completed successfully.
2026-06-21 19:20:33 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:33 INFO Replacement submit_p_template.sh → submit_p_4.sh completed successfully.
2026-06-21 19:20:33 INFO Queuing job for member 4...
2026-06-21 19:20:33 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:33 INFO [CMD] Running: ccc_msub ./submit_p_4.sh
2026-06-21 19:20:34 INFO Found: ['4937716']
2026-06-21 19:20:39 INFO [TGCC-IRENE] Submitted job with ID:['4937716']
2026-06-21 19:20:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5
2026-06-21 19:20:39 INFO The END file used for ENS5 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021300_24_ENS5.nc
2026-06-21 19:20:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS5.par completed successfully.
2026-06-21 19:20:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:39 INFO Replacement submit_p_template.sh → submit_p_5.sh completed successfully.
2026-06-21 19:20:39 INFO Queuing job for member 5...
2026-06-21 19:20:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:39 INFO [CMD] Running: ccc_msub ./submit_p_5.sh
2026-06-21 19:20:40 INFO Found: ['4937717']
2026-06-21 19:20:45 INFO [TGCC-IRENE] Submitted job with ID:['4937717']
2026-06-21 19:20:45 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:45 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6
2026-06-21 19:20:45 INFO The END file used for ENS6 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021300_24_ENS6.nc
2026-06-21 19:20:45 INFO Replacement chimere.template_ensemble.par → chimere.ENS6.par completed successfully.
2026-06-21 19:20:45 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:45 INFO Replacement submit_p_template.sh → submit_p_6.sh completed successfully.
2026-06-21 19:20:45 INFO Queuing job for member 6...
2026-06-21 19:20:45 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:45 INFO [CMD] Running: ccc_msub ./submit_p_6.sh
2026-06-21 19:20:46 INFO Found: ['4937718']
2026-06-21 19:20:51 INFO [TGCC-IRENE] Submitted job with ID:['4937718']
2026-06-21 19:20:51 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:51 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7
2026-06-21 19:20:51 INFO The END file used for ENS7 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021300_24_ENS7.nc
2026-06-21 19:20:51 INFO Replacement chimere.template_ensemble.par → chimere.ENS7.par completed successfully.
2026-06-21 19:20:51 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:51 INFO Replacement submit_p_template.sh → submit_p_7.sh completed successfully.
2026-06-21 19:20:51 INFO Queuing job for member 7...
2026-06-21 19:20:51 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:51 INFO [CMD] Running: ccc_msub ./submit_p_7.sh
2026-06-21 19:20:51 INFO Found: ['4937719']
2026-06-21 19:20:56 INFO [TGCC-IRENE] Submitted job with ID:['4937719']
2026-06-21 19:20:56 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:20:56 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8
2026-06-21 19:20:56 INFO The END file used for ENS8 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021300_24_ENS8.nc
2026-06-21 19:20:56 INFO Replacement chimere.template_ensemble.par → chimere.ENS8.par completed successfully.
2026-06-21 19:20:56 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:20:56 INFO Replacement submit_p_template.sh → submit_p_8.sh completed successfully.
2026-06-21 19:20:56 INFO Queuing job for member 8...
2026-06-21 19:20:56 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:20:56 INFO [CMD] Running: ccc_msub ./submit_p_8.sh
2026-06-21 19:20:59 INFO Found: ['4937720']
2026-06-21 19:21:04 INFO [TGCC-IRENE] Submitted job with ID:['4937720']
2026-06-21 19:21:04 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:04 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9
2026-06-21 19:21:04 INFO The END file used for ENS9 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021300_24_ENS9.nc
2026-06-21 19:21:04 INFO Replacement chimere.template_ensemble.par → chimere.ENS9.par completed successfully.
2026-06-21 19:21:04 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:04 INFO Replacement submit_p_template.sh → submit_p_9.sh completed successfully.
2026-06-21 19:21:04 INFO Queuing job for member 9...
2026-06-21 19:21:04 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:04 INFO [CMD] Running: ccc_msub ./submit_p_9.sh
2026-06-21 19:21:05 INFO Found: ['4937722']
2026-06-21 19:21:10 INFO [TGCC-IRENE] Submitted job with ID:['4937722']
2026-06-21 19:21:10 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:10 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10
2026-06-21 19:21:10 INFO The END file used for ENS10 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021300_24_ENS10.nc
2026-06-21 19:21:10 INFO Replacement chimere.template_ensemble.par → chimere.ENS10.par completed successfully.
2026-06-21 19:21:10 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:10 INFO Replacement submit_p_template.sh → submit_p_10.sh completed successfully.
2026-06-21 19:21:10 INFO Queuing job for member 10...
2026-06-21 19:21:10 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:10 INFO [CMD] Running: ccc_msub ./submit_p_10.sh
2026-06-21 19:21:11 INFO Found: ['4937724']
2026-06-21 19:21:16 INFO [TGCC-IRENE] Submitted job with ID:['4937724']
2026-06-21 19:21:16 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:16 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11
2026-06-21 19:21:16 INFO The END file used for ENS11 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021300_24_ENS11.nc
2026-06-21 19:21:16 INFO Replacement chimere.template_ensemble.par → chimere.ENS11.par completed successfully.
2026-06-21 19:21:16 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:16 INFO Replacement submit_p_template.sh → submit_p_11.sh completed successfully.
2026-06-21 19:21:16 INFO Queuing job for member 11...
2026-06-21 19:21:16 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:16 INFO [CMD] Running: ccc_msub ./submit_p_11.sh
2026-06-21 19:21:17 INFO Found: ['4937725']
2026-06-21 19:21:22 INFO [TGCC-IRENE] Submitted job with ID:['4937725']
2026-06-21 19:21:22 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:22 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12
2026-06-21 19:21:22 INFO The END file used for ENS12 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021300_24_ENS12.nc
2026-06-21 19:21:22 INFO Replacement chimere.template_ensemble.par → chimere.ENS12.par completed successfully.
2026-06-21 19:21:22 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:22 INFO Replacement submit_p_template.sh → submit_p_12.sh completed successfully.
2026-06-21 19:21:22 INFO Queuing job for member 12...
2026-06-21 19:21:22 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:22 INFO [CMD] Running: ccc_msub ./submit_p_12.sh
2026-06-21 19:21:23 INFO Found: ['4937727']
2026-06-21 19:21:28 INFO [TGCC-IRENE] Submitted job with ID:['4937727']
2026-06-21 19:21:28 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:28 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13
2026-06-21 19:21:28 INFO The END file used for ENS13 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021300_24_ENS13.nc
2026-06-21 19:21:28 INFO Replacement chimere.template_ensemble.par → chimere.ENS13.par completed successfully.
2026-06-21 19:21:28 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:28 INFO Replacement submit_p_template.sh → submit_p_13.sh completed successfully.
2026-06-21 19:21:28 INFO Queuing job for member 13...
2026-06-21 19:21:28 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:28 INFO [CMD] Running: ccc_msub ./submit_p_13.sh
2026-06-21 19:21:29 INFO Found: ['4937728']
2026-06-21 19:21:34 INFO [TGCC-IRENE] Submitted job with ID:['4937728']
2026-06-21 19:21:34 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:34 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14
2026-06-21 19:21:34 INFO The END file used for ENS14 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021300_24_ENS14.nc
2026-06-21 19:21:34 INFO Replacement chimere.template_ensemble.par → chimere.ENS14.par completed successfully.
2026-06-21 19:21:34 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:34 INFO Replacement submit_p_template.sh → submit_p_14.sh completed successfully.
2026-06-21 19:21:34 INFO Queuing job for member 14...
2026-06-21 19:21:34 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:34 INFO [CMD] Running: ccc_msub ./submit_p_14.sh
2026-06-21 19:21:34 INFO Found: ['4937729']
2026-06-21 19:21:39 INFO [TGCC-IRENE] Submitted job with ID:['4937729']
2026-06-21 19:21:39 INFO Replacing @TOKENS in CHIMERE .par template file ...
2026-06-21 19:21:39 INFO The output directory (run_dir) is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15
2026-06-21 19:21:39 INFO The END file used for ENS15 is: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021300_24_ENS15.nc
2026-06-21 19:21:39 INFO Replacement chimere.template_ensemble.par → chimere.ENS15.par completed successfully.
2026-06-21 19:21:39 INFO Replacing @TOKENS in CHIMERE template submit script ...
2026-06-21 19:21:39 INFO Replacement submit_p_template.sh → submit_p_15.sh completed successfully.
2026-06-21 19:21:39 INFO Queuing job for member 15...
2026-06-21 19:21:39 INFO [CMD] Entering directory: /ccc/work/cont003/gen7232/demoling/chimere_v2023_dart
2026-06-21 19:21:39 INFO [CMD] Running: ccc_msub ./submit_p_15.sh
2026-06-21 19:21:40 INFO Found: ['4937730']
2026-06-21 19:21:45 INFO [TGCC-IRENE] Submitted job with ID:['4937730']
2026-06-21 19:21:45 INFO Checking job status ...
2026-06-21 19:21:45 INFO None 4937713: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937714: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937715: status FINISHED
2026-06-21 19:21:45 INFO None 4937716: status FINISHED
2026-06-21 19:21:45 INFO None 4937717: status FINISHED
2026-06-21 19:21:45 INFO None 4937718: status FINISHED
2026-06-21 19:21:45 INFO None 4937719: status FINISHED
2026-06-21 19:21:45 INFO None 4937720: status FINISHED
2026-06-21 19:21:45 INFO None 4937722: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937724: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937725: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937727: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937728: status RUNNING/PENDING
2026-06-21 19:21:45 INFO None 4937729: status RUNNING/PENDING
2026-06-21 19:21:46 INFO None 4937730: status RUNNING/PENDING
2026-06-21 19:21:46 INFO Jobs still running: ['4937713', '4937714', '4937722', '4937724', '4937725', '4937727', '4937728', '4937729', '4937730']. Waiting...
2026-06-21 19:22:01 INFO None 4937713: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937714: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937715: status FINISHED
2026-06-21 19:22:01 INFO None 4937716: status FINISHED
2026-06-21 19:22:01 INFO None 4937717: status FINISHED
2026-06-21 19:22:01 INFO None 4937718: status FINISHED
2026-06-21 19:22:01 INFO None 4937719: status FINISHED
2026-06-21 19:22:01 INFO None 4937720: status FINISHED
2026-06-21 19:22:01 INFO None 4937722: status FINISHED
2026-06-21 19:22:01 INFO None 4937724: status FINISHED
2026-06-21 19:22:01 INFO None 4937725: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937727: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937728: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937729: status RUNNING/PENDING
2026-06-21 19:22:01 INFO None 4937730: status RUNNING/PENDING
2026-06-21 19:22:01 INFO Jobs still running: ['4937713', '4937714', '4937725', '4937727', '4937728', '4937729', '4937730']. Waiting...
2026-06-21 19:22:16 INFO None 4937713: status FINISHED
2026-06-21 19:22:16 INFO None 4937714: status FINISHED
2026-06-21 19:22:16 INFO None 4937715: status FINISHED
2026-06-21 19:22:16 INFO None 4937716: status FINISHED
2026-06-21 19:22:16 INFO None 4937717: status FINISHED
2026-06-21 19:22:16 INFO None 4937718: status FINISHED
2026-06-21 19:22:16 INFO None 4937719: status FINISHED
2026-06-21 19:22:16 INFO None 4937720: status FINISHED
2026-06-21 19:22:16 INFO None 4937722: status FINISHED
2026-06-21 19:22:16 INFO None 4937724: status FINISHED
2026-06-21 19:22:16 INFO None 4937725: status FINISHED
2026-06-21 19:22:16 INFO None 4937727: status FINISHED
2026-06-21 19:22:16 INFO None 4937728: status FINISHED
2026-06-21 19:22:16 INFO None 4937729: status FINISHED
2026-06-21 19:22:16 INFO None 4937730: status RUNNING/PENDING
2026-06-21 19:22:16 INFO Jobs still running: ['4937730']. Waiting...
2026-06-21 19:22:31 INFO None 4937713: status FINISHED
2026-06-21 19:22:31 INFO None 4937714: status FINISHED
2026-06-21 19:22:31 INFO None 4937715: status FINISHED
2026-06-21 19:22:31 INFO None 4937716: status FINISHED
2026-06-21 19:22:31 INFO None 4937717: status FINISHED
2026-06-21 19:22:31 INFO None 4937718: status FINISHED
2026-06-21 19:22:31 INFO None 4937719: status FINISHED
2026-06-21 19:22:31 INFO None 4937720: status FINISHED
2026-06-21 19:22:31 INFO None 4937722: status FINISHED
2026-06-21 19:22:31 INFO None 4937724: status FINISHED
2026-06-21 19:22:31 INFO None 4937725: status FINISHED
2026-06-21 19:22:31 INFO None 4937727: status FINISHED
2026-06-21 19:22:31 INFO None 4937728: status FINISHED
2026-06-21 19:22:31 INFO None 4937729: status FINISHED
2026-06-21 19:22:31 INFO None 4937730: status FINISHED
2026-06-21 19:22:31 INFO Jobs ['4937713', '4937714', '4937715', '4937716', '4937717', '4937718', '4937719', '4937720', '4937722', '4937724', '4937725', '4937727', '4937728', '4937729', '4937730'] have finished
2026-06-21 19:22:31 INFO Checking restart files were created ...
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 1: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS1/end.2020021415_9_ENS1.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 2: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS2/end.2020021415_9_ENS2.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 3: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS3/end.2020021415_9_ENS3.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 4: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS4/end.2020021415_9_ENS4.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 5: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS5/end.2020021415_9_ENS5.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 6: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS6/end.2020021415_9_ENS6.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 7: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS7/end.2020021415_9_ENS7.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 8: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS8/end.2020021415_9_ENS8.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 9: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS9/end.2020021415_9_ENS9.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 10: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS10/end.2020021415_9_ENS10.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 11: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS11/end.2020021415_9_ENS11.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 12: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS12/end.2020021415_9_ENS12.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 13: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS13/end.2020021415_9_ENS13.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 14: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS14/end.2020021415_9_ENS14.nc
2026-06-21 19:22:31 WARNING ModelType.CHIMERE | resatrt_file is missing for mem 15: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/end.2020021415_9_ENS15.nc
2026-06-21 19:22:31 INFO Check chimere log file at: /ccc/scratch/cont003/gen7232/demoling/OUT_orch_chimdart/OUT_Chimere/cyICEMIdmp_0615_15m_low_v2/ENS15/ENS15_2020021415.out
2026-06-21 19:22:31 ERROR [PIPELINE] Error: The following chimere ENS run(s) failed (exit code 1): [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
Traceback (most recent call last):
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/main.py", line 103, in <module>
    pipeline.run_pipeline()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/base_pipeline.py", line 138, in run_pipeline
    self.run_model()
  File "/ccc/work/cont003/gen7232/demoling/mimesi_orch/pipelines/chimere2023/pipeline.py", line 329, in run_model
    raise ModelRunError(f"The following chimere ENS run(s) failed (exit code 1): {mems_to_rerun}")
pipeline_errors.ModelRunError: The following chimere ENS run(s) failed (exit code 1): [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
+ exit 0
